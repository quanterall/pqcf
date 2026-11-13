%% Where tasks run (decision records R1, phase 3, and P1): one task interface, several executors. A task is plain
%% data, {pq_task, Module, Args} for a module implementing the pq_task behaviour (or one of the reconstruction's own
%% tuples, pq_ff_task's), its result a pure function of it, so the results are the same whichever executor runs
%% them and however the tasks are spread (D1). map_traced/3 also returns each result's provenance (pq_task).
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
-export([map/3, map_traced/3, run_task/2, run_task_traced/2]).

-define(TIMEOUT, 600000).

%% The results of Tasks, in their order.
-spec map([pq_task:task()], local | {nodes, [node()]}, map()) -> [term()].
map(Tasks, Exec, Opts) ->
    run(Tasks, Exec, Opts#{traced => false}).

%% The results of Tasks, in their order, each with its provenance record: {Result, pq_task:provenance()}. The
%% results are map/3's.
-spec map_traced([pq_task:task()], local | {nodes, [node()]}, map()) -> [{term(), pq_task:provenance()}].
map_traced(Tasks, Exec, Opts) ->
    run(Tasks, Exec, Opts#{traced => true}).

run([], _Exec, _Opts) -> [];
run(Tasks, Exec, Opts) ->
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
    Run = case maps:get(traced, Opts) of true -> run_task_traced; false -> run_task end,
    {Pid, Ref} = spawn_monitor(fun() -> Self ! {pq_rexec, self(), attempt(Node, Run, Task, Home, Timeout)} end),
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

attempt(Node, Run, Task, Home, _Timeout) when Node =:= Home ->
    try {ok, ?MODULE:Run(Task, Home)}
    catch Class:Why:Stack -> {task_failed, Class, Why, Stack}
    end;
attempt(Node, Run, Task, Home, Timeout) ->
    try erpc:call(Node, ?MODULE, Run, [Task, Home], Timeout) of
        Result -> {ok, Result}
    catch
        error:{erpc, Why} -> {node_failed, Why};
        error:{exception, Why, Stack} -> {task_failed, error, Why, Stack};
        exit:{exception, Why} -> {task_failed, exit, Why, []};
        throw:Why -> {task_failed, throw, Why, []};
        exit:Why -> {node_failed, Why}
    end.

%%% The tasks

%% A task's result on this node; Home is the coordinating node (its store holds the tasks' data). See pq_task.
-spec run_task(pq_task:task(), node()) -> term().
run_task(Task, Home) -> pq_task:run(Task, Home).

-spec run_task_traced(pq_task:task(), node()) -> {term(), pq_task:provenance()}.
run_task_traced(Task, Home) -> pq_task:run_traced(Task, Home).
