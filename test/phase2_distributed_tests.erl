%% Phase 2: distribution (FS 9.5, FR-P02, FR-P05, FR-P08) with peer nodes on this host: identical samples to a local
%% run, node loss mid-job retried elsewhere, the handshake, the memory pre-check.
-module(phase2_distributed_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

setup() ->
    os:cmd("epmd -daemon"),
    timer:sleep(300),
    {ok, _} = net_kernel:start(pqtest, #{name_domain => shortnames}),
    {ok, _} = application:ensure_all_started(pq),
    ok = pq_node_registry:register_node(),
    Paths = lists:append([["-pa", P] || P <- code:get_path(), string:find(P, "_build") =/= nomatch]),
    Peers = [begin
                 {ok, Pid, Node} = peer:start_link(#{name => list_to_atom("pqw" ++ integer_to_list(K)), args => Paths}),
                 {ok, _} = erpc:call(Node, application, ensure_all_started, [pq]),
                 ok = erpc:call(Node, pq_node_registry, register_node, []),
                 {Pid, Node}
             end || K <- [1, 2]],
    wait_members(3, 50),
    Peers.

wait_members(N, 0) -> error({pg_members, length(pg:get_members(pq, workers)), N});
wait_members(N, T) ->
    case length(pg:get_members(pq, workers)) of
        N -> ok;
        _ -> timer:sleep(100), wait_members(N, T - 1)
    end.

cleanup(Peers) ->
    [catch peer:stop(P) || {P, _} <- Peers],
    net_kernel:stop().

c() ->
    {ok, C} = pq_circuit:new(5, {fock, [1, 1, 1, 0, 0]}),
    C#circuit{instrs = [{beamsplitter, J, J + 1, 0.4 + 0.1 * J, 0.2} || J <- lists:seq(0, 3)] ++
                       [{beamsplitter, 0, 1, 0.9, 0.0}, {measure_pnr, [0, 1, 2, 3, 4]}]}.

distributed_test_() ->
    {setup, fun setup/0, fun cleanup/1,
     fun(Peers) ->
             [{"handshake: three compatible nodes",
               ?_test(begin
                          {Nodes, Hetero} = pq_node_registry:eligible_nodes(#{}),
                          ?assertEqual(3, length(Nodes)),
                          ?assertEqual(node(), hd(Nodes)),
                          ?assertNot(Hetero)
                      end)},
              {"distributed samples equal local ones (D1)",
               ?_test(begin
                          Local = samples(#{shots => 4000, seed => 11, max_inflight => 2}),
                          Dist = samples(#{shots => 4000, seed => 11, distributed => true, chunk => 100}),
                          ?assertEqual(Local, Dist)
                      end)},
              {"a node lost mid-job: its chunks run elsewhere, the samples are the same",
               {timeout, 120,
                ?_test(begin
                           Big = big(),
                           Opts = #{shots => 100000, seed => 5, chunk => 200, max_retries => 3},
                           {Tl, Local} = timer:tc(fun() -> samples(Big, Opts#{max_inflight => 4}) end, millisecond),
                           Self = self(),
                           T0 = erlang:monotonic_time(millisecond),
                           spawn(fun() -> Self ! {res, samples(Big, Opts#{distributed => true})} end),
                           Stop = max(20, Tl div 10),
                           timer:sleep(Stop),
                           %% the job must still be running when the node goes, or the test proves nothing
                           receive {res, _} -> error(job_finished_before_the_node_was_stopped) after 0 -> ok end,
                           [{P1, _} | _] = Peers,
                           peer:stop(P1),
                           receive
                               {res, Dist} ->
                                   io:format(user, "~nlocal run ~p ms; node stopped at ~p ms, distributed job done after ~p ms~n",
                                             [Tl, Stop, erlang:monotonic_time(millisecond) - T0]),
                                   ?assertEqual(Local, Dist)
                           after 100000 -> error(timeout)
                           end
                       end)}},
              {"kernel isolation: a node that dies in a kernel takes only itself down (FR-P07)",
               {timeout, 120,
                ?_test(begin
                           Local = samples(#{shots => 3000, seed => 21, max_inflight => 2}),
                           ?assertEqual(Local, samples(#{shots => 3000, seed => 21, isolate_kernels => true, chunk => 100})),
                           ?assertEqual(Local, samples(#{shots => 3000, seed => 21, isolate_kernels => true, chunk => 100,
                                                         isolated_nodes => 2, halt_first_isolated => true, max_retries => 3})),
                           %% the isolated nodes are gone afterwards
                           ?assertEqual([], [N || N <- nodes(), string:prefix(atom_to_list(N), "pq_isolated") =/= nomatch])
                       end)}},
              {"memory pre-check",
               ?_assertMatch({error, {dimension_limit, _, 10}},
                             pq_sampling:sample(boson, c(), #{shots => 1000, seed => 1, max_job_memory => 10}))}]
     end}.

samples(Opts) -> samples(c(), Opts).

samples(C, Opts) ->
    {ok, #result{samples = S}} = pq_sampling:sample(boson, C, Opts),
    S.

%% 8 photons in 12 modes
big() ->
    {ok, C} = pq_circuit:new(12, {fock, lists:duplicate(8, 1) ++ lists:duplicate(4, 0)}),
    C#circuit{instrs = lists:append([[{beamsplitter, J, J + 1, 0.3 + 0.05 * J, 0.1 * L} || J <- lists:seq(L rem 2, 10, 2)]
                                     || L <- lists:seq(0, 11)]) ++ [{measure_pnr, lists:seq(0, 11)}]}.
