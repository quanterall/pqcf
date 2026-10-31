%% FR-F05 trajectory mode: pure-state sampling with drawn lost-photon counts has the density-matrix distribution.
-module(phase3_trajectory_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

setup() -> {ok, _} = application:ensure_all_started(pq), ok.

c(N, In, Is) -> {ok, C} = pq_circuit:new(N, In), C#circuit{instrs = Is}.

trajectory_test_() ->
    {setup, fun setup/0, fun(_) -> ok end,
     [?_test(begin
                 C = c(2, {fock, [1, 1]}, [{loss, 0, 0.6}, {beamsplitter, 0, 1, 0.5, 0.3}, {loss, 1, 0.8},
                                           {measure_pnr, [0, 1]}]),
                 {ok, Exact, _} = pq_api:probabilities(C, #{backend => fock_dm}),
                 {ok, #result{samples = S, metadata = #{backend := fock}}} =
                     pq_simulator:run(pq_fock, C, #{shots => 20000, seed => 6, noise => trajectory}),
                 N = length(S),
                 Hist = lists:foldl(fun(R, A) -> maps:update_with(R, fun(X) -> X + 1 / N end, 1 / N, A) end, #{}, S),
                 {ok, Tv} = pq_stats:tv_distance(Hist, Exact),
                 ?assert(Tv < 0.015),
                 %% D1: workers and chunks do not change the samples
                 {ok, #result{samples = S1}} = pq_simulator:run(pq_fock, C, #{shots => 500, seed => 6, noise => trajectory, max_inflight => 1}),
                 {ok, #result{samples = S4}} = pq_simulator:run(pq_fock, C, #{shots => 500, seed => 6, noise => trajectory, max_inflight => 4}),
                 ?assertEqual(S1, S4),
                 %% without the mode, a pure Fock state refuses loss
                 ?assertMatch({error, {unsupported_instruction, loss, fock}}, pq_simulator:run(pq_fock, C, #{shots => 1}))
             end)]}.
