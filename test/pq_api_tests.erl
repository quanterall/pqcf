%% The facade (FS 7.1): backend choice, run/sample/evolve, exact probabilities on every backend, run_batch.
-module(pq_api_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

c(N, In, Is) -> {ok, C} = pq_api:circuit(N, In), C#circuit{instrs = Is}.

backend_test() ->
    ?assertEqual(gaussian, pq_api:backend_for(c(1, vacuum, [{squeezing, 0, 1, 0}, {measure_pnr, [0]}]), #{})),
    ?assertEqual(fock, pq_api:backend_for(c(1, vacuum, [{kerr, 0, 1}, {measure_pnr, [0]}]), #{})),
    ?assertEqual(fock, pq_api:backend_for(c(2, {fock, [1, 1]}, [{beamsplitter, 0, 1, 1, 0}]), #{})),
    ?assertEqual(fock_dm, pq_api:backend_for(c(1, {fock, [1]}, [{loss, 0, 0.5}]), #{})),
    ?assertEqual(boson, pq_api:backend_for(c(1, {fock, [1]}, []), #{backend => boson})).

probabilities_test() ->
    %% HOM three ways: Fock, boson (permanents), and the Gaussian path for a squeezed state
    H = math:pi() / 4,
    Hom = c(2, {fock, [1, 1]}, [{beamsplitter, 0, 1, H, 0}, {measure_pnr, [0, 1]}]),
    {ok, Pf, _} = pq_api:probabilities(Hom, #{}),
    {ok, Pb, _} = pq_api:probabilities(Hom, #{backend => boson}),
    [?assert(abs(maps:get(K, Pf, 0.0) - maps:get(K, Pb, 0.0)) < 1.0e-14) || K <- [[2, 0], [1, 1], [0, 2]]],
    ?assert(abs(maps:get([2, 0], Pf) - 0.5) < 1.0e-14),
    R = 0.5,
    {ok, Pg, #{truncated_mass := T}} = pq_api:probabilities(c(1, vacuum, [{squeezing, 0, R, 0}, {measure_pnr, [0]}]), #{}),
    ?assert(abs(maps:get([0], Pg) - 1 / math:cosh(R)) < 1.0e-14 andalso T < 1.0e-6),
    {ok, Pt, _} = pq_api:probabilities(c(1, vacuum, [{squeezing, 0, R, 0}, {measure_threshold, [0]}]), #{}),
    ?assert(abs(maps:get([1], Pt) - (1 - 1 / math:cosh(R))) < 1.0e-14),
    ?assertMatch({error, {dimension_limit, _, 2}}, pq_api:probabilities(Hom, #{backend => boson, max_enumeration => 2})).

run_and_batch_test() ->
    ?assertEqual({error, {unsupported_input, shots_required}}, pq_api:sample(c(1, vacuum, [{measure_pnr, [0]}]), #{})),
    C = c(2, {fock, [1, 0]}, [{beamsplitter, 0, 1, {sym, t}, 0}, {measure_pnr, [0, 1]}]),
    {ok, [R0, R1]} = pq_api:run_batch(C, [#{t => 0}, #{t => math:pi() / 2}], #{shots => 50, seed => 1}),
    ?assert(lists:all(fun(S) -> S =:= [1, 0] end, R0#result.samples)),
    ?assert(lists:all(fun(S) -> S =:= [0, 1] end, R1#result.samples)),
    ?assertMatch({error, {unbound_symbols, [t]}}, pq_api:run_batch(C, [#{}], #{shots => 1})),
    {ok, S, []} = pq_api:evolve(c(1, vacuum, [{displacement, 0, 1, 0}]), #{}),
    ?assertEqual({ok, 1.0}, pq_stats:mean_photon_number(S, [0])).

%% threshold detection of a displaced state (the loop torontonian): a coherent state, |alpha| = 0.3, stays dark
%% with probability exp(-0.09)
displaced_threshold_test() ->
    C = c(1, vacuum, [{displacement, 0, 0.3, 0}, {measure_threshold, [0]}]),
    {ok, P, _} = pq_api:probabilities(C, #{}),
    ?assert(abs(maps:get([0], P) - math:exp(-0.09)) < 1.0e-14),
    ?assert(abs(maps:get([1], P) - (1 - math:exp(-0.09))) < 1.0e-14).
