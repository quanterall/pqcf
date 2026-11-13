%% The task behaviour (pq_task, decision record P1): an application's task module run by the executors (pq_rexec),
%% refused until allowed, with provenance records; the reconstruction's own tuples unchanged as pq_ff_task's tasks.
-module(pq_task_tests).
-include_lib("eunit/include/eunit.hrl").

squares(Ns) -> [{pq_task, pq_task_example, {square, N}} || N <- Ns].

%% a module runs only once allowed; then in task order, whatever the pool size
local_test() ->
    %% start from PQCF's own modules only (the allow-list is the VM's, shared with other tests)
    persistent_term:erase({pq_task, allowed}),
    ?assertError({task_not_allowed, pq_task_example}, pq_task:run(hd(squares([2])), node())),
    ?assertError(_, pq_rexec:map(squares([1, 2, 3]), local, #{})),
    ok = pq_task:allow(pq_task_example),
    ?assert(lists:member(pq_task_example, pq_task:allowed())),
    Want = [N * N || N <- lists:seq(1, 40)],
    ?assertEqual(Want, pq_rexec:map(squares(lists:seq(1, 40)), local, #{workers => 1})),
    ?assertEqual(Want, pq_rexec:map(squares(lists:seq(1, 40)), local, #{workers => 7})),
    %% a failing task raises
    ?assertError(boom, pq_rexec:map([{pq_task, pq_task_example, {fail, boom}}], local, #{})),
    %% resources: declared, or none
    ?assertEqual(#{cpu => 1}, pq_task:resources(hd(squares([5])))),
    ?assertEqual(#{}, pq_task:resources({values, <<"k">>, 101, [[1]]})).

%% the reconstruction's tuples: the same results bare and as {pq_task, pq_ff_task, Tuple}, traced or not
ff_tuples_test() ->
    {ok, Prog} = pq_symbolic:compile([{add, [{sym, x}, {mul, [{sym, y}, {sym, y}]}]}], [x, y], #{exact_constants => true}),
    Key = pq_rstore:put(Prog),
    Bare = [{values, Key, 101, [[X, 3]]} || X <- lists:seq(1, 10)],
    Wrapped = [{pq_task, pq_ff_task, T} || T <- Bare],
    Want = pq_rexec:map(Bare, local, #{}),
    ?assertEqual(Want, pq_rexec:map(Wrapped, local, #{})),
    Traced = pq_rexec:map_traced(Bare, local, #{workers => 3}),
    ?assertEqual(Want, [R || {R, _} <- Traced]).

%% the provenance record: what ran, where, when
provenance_test() ->
    ok = pq_task:allow(pq_task_example),
    [{9, P1}, {16, P2}] = pq_rexec:map_traced(squares([3, 4]), local, #{}),
    ?assertMatch(#{module := pq_task_example, node := N, describe := #{kind := square}} when N =:= node(), P1),
    ?assertEqual(binary:encode_hex(pq_task_example:module_info(md5), lowercase), maps:get(code, P1)),
    ?assert(maps:get(seconds, P1) >= 0.0),
    ?assert(is_integer(maps:get(started, P1))),
    ?assert(is_binary(maps:get(pq_version, P1))),
    %% the task's hash: the same task, the same hash; another task, another
    {_, P1b} = pq_task:run_traced(hd(squares([3])), node()),
    ?assertEqual(maps:get(task, P1), maps:get(task, P1b)),
    ?assertNotEqual(maps:get(task, P1), maps:get(task, P2)),
    ?assertEqual(64, byte_size(maps:get(task, P1))).

nodes_test_() ->
    {setup,
     fun() ->
             os:cmd("epmd -daemon"),
             timer:sleep(300),
             {ok, _} = net_kernel:start(pqtask, #{name_domain => shortnames}),
             Paths = lists:append([["-pa", P] || P <- code:get_path(), string:find(P, "_build") =/= nomatch]),
             [begin {ok, Pid, Node} = peer:start_link(#{name => list_to_atom("pqtk" ++ integer_to_list(K)), args => Paths}),
                    {Pid, Node} end || K <- [1, 2]]
     end,
     fun(Peers) -> [catch peer:stop(P) || {P, _} <- Peers], net_kernel:stop() end,
     fun([{_, N1}, {_, N2}]) ->
             {timeout, 120,
              ?_test(begin
                         ok = pq_task:allow(pq_task_example),
                         Tasks = squares(lists:seq(1, 24)),
                         Want = pq_rexec:map(Tasks, local, #{}),
                         %% the peers have not allowed the module: they refuse it
                         ?assertError({task_not_allowed, pq_task_example},
                                      pq_rexec:map(Tasks, {nodes, [N1]}, #{per_node => 1})),
                         %% once they do, their results are this node's
                         [ok = erpc:call(N, pq_task, allow, [pq_task_example]) || N <- [N1, N2]],
                         ?assertEqual(Want, pq_rexec:map(Tasks, {nodes, [N1, N2]}, #{per_node => 2})),
                         Traced = pq_rexec:map_traced(Tasks, {nodes, [N1, N2]}, #{per_node => 2}),
                         ?assertEqual(Want, [R || {R, _} <- Traced]),
                         Ran = lists:usort([maps:get(node, P) || {_, P} <- Traced]),
                         ?assertEqual(lists:sort([N1, N2]), Ran)
                     end)}
     end}.
