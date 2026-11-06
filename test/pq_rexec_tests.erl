%% Reconstruction across nodes (pq_rexec, pq_rstore; decision record R1, phase 3): the same polynomial on this
%% node, on peer nodes, with a node that can't be reached, and with a node lost in the middle of the run.
-module(pq_rexec_tests).
-include_lib("eunit/include/eunit.hrl").

poly() ->
    rand:seed(exsss, {5, 1, 5}),
    Vars = [w, x, y, z],
    E = {add, [{mul, [{rat, rand:uniform(41) - 21, rand:uniform(9)} |
                      [{pow, {sym, lists:nth(rand:uniform(4), Vars)}, rand:uniform(6)} || _ <- [1, 2, 3]]]}
               || _ <- lists:seq(1, 60)]},
    {ok, Prog} = pq_symbolic:compile([E], Vars, #{exact_constants => true, optimize => horner_cse}),
    {Prog, Vars}.

%% the local executor: results in task order, whatever the pool size
local_test() ->
    {Prog, Vars} = poly(),
    {ok, A} = pq_recon:polynomial(Prog, Vars, #{workers => 1}),
    ?assertEqual({ok, A}, pq_recon:polynomial(Prog, Vars, #{workers => 7})),
    ?assertEqual({ok, A}, pq_recon:polynomial(Prog, Vars, #{workers => 3, batch => 2})),
    %% a task that fails raises, it isn't retried forever
    ?assertError(_, pq_rexec:map([{values, <<"no such key">>, 101, [[1]]}], local, #{})).

nodes_test_() ->
    {setup,
     fun() ->
             os:cmd("epmd -daemon"),
             timer:sleep(300),
             {ok, _} = net_kernel:start(pqrexec, #{name_domain => shortnames}),
             Paths = lists:append([["-pa", P] || P <- code:get_path(), string:find(P, "_build") =/= nomatch]),
             [begin {ok, Pid, Node} = peer:start_link(#{name => list_to_atom("pqrx" ++ integer_to_list(K)), args => Paths}),
                    {Pid, Node} end || K <- [1, 2, 3]]
     end,
     fun(Peers) -> [catch peer:stop(P) || {P, _} <- Peers], net_kernel:stop() end,
     fun([{_, N1}, {_, N2}, {P3, N3}]) ->
             {timeout, 300,
              ?_test(begin
                         {Prog, Vars} = poly(),
                         {ok, Local} = pq_recon:polynomial(Prog, Vars, #{}),
                         %% two peers; this node and a peer
                         ?assertEqual({ok, Local}, pq_recon:polynomial(Prog, Vars, #{executor => {nodes, [N1, N2]}})),
                         ?assertEqual({ok, Local}, pq_recon:polynomial(Prog, Vars, #{executor => {nodes, [node(), N2]}, per_node => 2})),
                         %% a node that can't be reached is left out
                         ?assertEqual({ok, Local}, pq_recon:polynomial(Prog, Vars, #{executor => {nodes, [N1, 'ghost@nohost']}})),
                         %% the store: a peer fetches a blob from this node by its key, once
                         Key = pq_rstore:put({blob, lists:seq(1, 100)}),
                         ?assertEqual({ok, {blob, lists:seq(1, 100)}}, erpc:call(N2, pq_rstore, fetch, [Key, node()])),
                         ?assertEqual({ok, {blob, lists:seq(1, 100)}}, erpc:call(N2, pq_rstore, get, [Key])),
                         ?assertEqual(not_found, erpc:call(N2, pq_rstore, fetch, [pq_rstore:key(nothing), node()])),
                         %% a node lost in the middle of the run: its tasks run elsewhere, the result is the same
                         Self = self(),
                         spawn(fun() -> timer:sleep(30), catch peer:stop(P3), Self ! stopped end),
                         ?assertEqual({ok, Local}, pq_recon:polynomial(Prog, Vars, #{executor => {nodes, [N3, N1]}, per_node => 1})),
                         receive stopped -> ok after 5000 -> ok end,
                         %% the requeue, deterministically: an unreachable node holds slots in the pool, so its
                         %% tasks fail and run on the others
                         Tasks = [{values, pq_rstore:put(Prog), 101, [[X, 2, 3, 4]]} || X <- lists:seq(1, 12)],
                         Want = pq_rexec:map(Tasks, local, #{}),
                         ?assertEqual(Want, pq_rexec:map(Tasks, {nodes, []}, #{slots => ['ghost@nohost', N1, 'ghost@nohost', N2]})),
                         ?assertEqual(Want, pq_rexec:map(Tasks, {nodes, []}, #{slots => ['ghost@nohost', 'ghost@nohost']})),
                         %% every node gone: the work comes back here
                         ?assertEqual({ok, Local}, pq_recon:polynomial(Prog, Vars, #{executor => {nodes, [N3]}}))
                     end)}
     end}.
