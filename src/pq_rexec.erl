%% Where reconstruction tasks run (decision record R1, phase 3): one task interface, several executors. A task is
%% plain data, its result a pure function of it (run_task/2), so the results are the same whichever executor runs
%% them and however the tasks are spread (D1).
%%
%%   local (the default)     this node: a pool with one slot per dirty CPU scheduler (the NIF work runs there).
%%   {nodes, Nodes}          one owner's Erlang nodes, by erpc from this node (a star: the workers need no
%%                           connections to each other, so it scales to 100+ nodes without a full mesh). Each node
%%                           gets as many slots as it has dirty CPU schedulers. A node that fails (no connection,
%%                           a timeout) is dropped and its tasks run elsewhere; with every node gone, here.
%%
%% The federation executor (FVLLMP peers over HTTPS) is to come (R1, phase 5): the same tasks, another transport.
%% Tasks fetch their programs and supports from the coordinator's store by key (pq_rstore).
-module(pq_rexec).
-export([map/3, run_task/2]).

-define(TIMEOUT, 600000).

%% The results of Tasks, in their order.
-spec map([tuple()], local | {nodes, [node()]}, map()) -> [term()].
map([], _Exec, _Opts) -> [];
map(Tasks, Exec, Opts) ->
    %% slots => [Node]: the slots given outright, one entry per slot (for tests, and callers that place work)
    Slots = case maps:get(slots, Opts, none) of none -> slots(Exec, Opts); Given -> Given end,
    loop(lists:enumerate(Tasks), Slots, #{}, #{}, length(Tasks), Opts).

slots(local, Opts) ->
    lists:duplicate(maps:get(workers, Opts, erlang:system_info(dirty_cpu_schedulers_online)), node());
slots({nodes, Nodes}, Opts) ->
    Slots = lists:append([case per_node(N, Opts) of
                              {ok, K} -> lists:duplicate(K, N);
                              error -> []
                          end || N <- lists:usort(Nodes)]),
    case Slots of [] -> slots(local, Opts); _ -> interleave(Slots) end.

per_node(N, Opts) when N =:= node() -> {ok, maps:get(workers, Opts, erlang:system_info(dirty_cpu_schedulers_online))};
per_node(N, Opts) ->
    try erpc:call(N, erlang, system_info, [dirty_cpu_schedulers_online], 10000) of
        K -> {ok, maps:get(per_node, Opts, K)}
    catch
        _:_ -> error
    end.

%% the slots dealt across the nodes in turn, so a batch spreads over all of them
interleave(Slots) ->
    ByNode = lists:foldr(fun(N, M) -> maps:update_with(N, fun(L) -> [N | L] end, [N], M) end, #{}, Slots),
    deal([L || {_, L} <- lists:sort(maps:to_list(ByNode))]).

deal([]) -> [];
deal(Lists) -> [hd(L) || L <- Lists] ++ deal([tl(L) || L <- Lists, tl(L) =/= []]).

%% Queue: [{Index, Task}]; Free: slots (nodes); Running: Pid => {Ref, Index, Task, Node}; Done: Index => Result.
%% Each attempt runs in a monitored process that sends its outcome; a crash of that process (never expected) is
%% raised.
loop(_Queue, _Free, Running, Done, Total, _Opts) when map_size(Done) =:= Total, map_size(Running) =:= 0 ->
    [maps:get(I, Done) || I <- lists:seq(1, Total)];
loop([{I, Task} | Queue], [Node | Free], Running, Done, Total, Opts) ->
    Home = node(),
    Self = self(),
    Timeout = maps:get(task_timeout, Opts, ?TIMEOUT),
    {Pid, Ref} = spawn_monitor(fun() -> Self ! {pq_rexec, self(), attempt(Node, Task, Home, Timeout)} end),
    loop(Queue, Free, Running#{Pid => {Ref, I, Task, Node}}, Done, Total, Opts);
loop(Queue, Free, Running, Done, Total, Opts) ->
    receive
        {pq_rexec, Pid, Outcome} when is_map_key(Pid, Running) ->
            {{Ref, I, Task, Node}, Running2} = maps:take(Pid, Running),
            erlang:demonitor(Ref, [flush]),
            case Outcome of
                {ok, Result} ->
                    loop(Queue, [Node | Free], Running2, Done#{I => Result}, Total, Opts);
                {node_failed, _} ->
                    %% the node is gone: its slots are dropped, the task goes back to the front of the queue
                    Free2 = [N || N <- Free, N =/= Node],
                    Free3 = case Free2 =:= [] andalso [x || {_, {_, _, _, N}} <- maps:to_list(Running2), N =/= Node] =:= [] of
                                true -> slots(local, Opts);            % nothing left: here
                                false -> Free2
                            end,
                    loop([{I, Task} | Queue], Free3, Running2, Done, Total, Opts);
                {task_failed, Class, Why, Stack} ->
                    erlang:raise(Class, Why, Stack)
            end;
        {'DOWN', _Ref, process, Pid, Reason} when is_map_key(Pid, Running) ->
            error({pq_rexec, Reason})
    end.

attempt(Node, Task, Home, _Timeout) when Node =:= Home ->
    try {ok, run_task(Task, Home)}
    catch Class:Why:Stack -> {task_failed, Class, Why, Stack}
    end;
attempt(Node, Task, Home, Timeout) ->
    try erpc:call(Node, ?MODULE, run_task, [Task, Home], Timeout) of
        Result -> {ok, Result}
    catch
        error:{erpc, Why} -> {node_failed, Why};
        error:{exception, Why, Stack} -> {task_failed, error, Why, Stack};
        exit:{exception, Why} -> {task_failed, exit, Why, []};
        throw:Why -> {task_failed, throw, Why, []};
        exit:Why -> {node_failed, Why}
    end.

%%% The tasks

%% {column, ProgramKey, P, R, N, SupportKey, Suffix}: one column of Zippel's interpolation (pq_nif:ff_column/6):
%%   {ok, Coefficients} | bad | singular | bad_prime
%% {line_column, ProgramKey, P, Shift, R, N, SupportKey, Suffix, Ts, Dn, Dd}: one column of the coefficients of a
%%   rational function along lines (pq_nif:ff_line_column/10): {ok, [Coefficients per output]} | bad | singular
%% {linsys_all, BlackBox, P, Points}: a linear system's wanted component for every right-hand side at the points
%%   (pq_nif:ff_linsys_all/3): {ok, [Values | bad]}
%% {values, ProgramKey, P, Points}: the program's outputs at the points (pq_ff:eval/3)
-spec run_task(tuple(), node()) -> term().
run_task({column, ProgKey, P, R, N, SupportKey, Suffix}, Home) ->
    {ok, {ff, Ref, _, _}} = pq_rstore:prepared(ProgKey, Home),
    {ok, Support} = pq_rstore:fetch(SupportKey, Home),
    pq_nif:ff_column(Ref, P, R, N, Support, Suffix);
run_task({line_column, ProgKey, P, Shift, R, N, SupportKey, Suffix, Ts, Dn, Dd}, Home) ->
    {ok, {ff, Ref, _, _}} = pq_rstore:prepared(ProgKey, Home),
    {ok, Support} = pq_rstore:fetch(SupportKey, Home),
    pq_nif:ff_line_column(Ref, P, Shift, R, N, Support, Suffix, Ts, Dn, Dd);
run_task({linsys_all, BB, P, Points}, Home) ->
    {ok, {ff, Ref, _, _}} = pq_rstore:prepared(BB, Home),
    pq_nif:ff_linsys_all(Ref, P, lists:append(Points));
run_task({values, ProgKey, P, Points}, Home) ->
    {ok, FF} = pq_rstore:prepared(ProgKey, Home),
    pq_ff:eval(FF, P, Points).
