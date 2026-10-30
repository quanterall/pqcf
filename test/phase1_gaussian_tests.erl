%% Phase 1, Gaussian milestone: pq_state (GT-C03 and FR-G07/G09), lowering (FR-C07), capabilities (FS 5.3,
%% FR-M06), and the Gaussian backend end to end (A.3-A.6, FR-G05).
-module(phase1_gaussian_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

%%% GT-C03: basis ordering and dimension (FS 3.4)

basis_order_test() ->
    %% m = 2, C = 2: (0,0) (0,1) (1,0) (0,2) (1,1) (2,0)
    ?assertEqual(6, pq_state:basis_size(2, 2)),
    ?assertEqual([[0, 0], [0, 1], [1, 0], [0, 2], [1, 1], [2, 0]],
                 [element(2, pq_state:pattern(I, 2, 2)) || I <- lists:seq(0, 5)]),
    ?assertEqual(126, pq_state:basis_size(4, 5)),
    ?assertEqual({error, {index_out_of_range, 6}}, pq_state:pattern(6, 2, 2)),
    ?assertMatch({error, {index_out_of_range, _}}, pq_state:index([2, 1], 2, 2)).

basis_round_trip_test_() ->
    [?_test(begin
                D = pq_state:basis_size(M, C),
                Ps = [element(2, pq_state:pattern(I, M, C)) || I <- lists:seq(0, D - 1)],
                ?assertEqual(lists:seq(0, D - 1), [element(2, pq_state:index(P, M, C)) || P <- Ps]),
                %% total ascending, then lexicographic ascending
                Keys = [{lists:sum(P), P} || P <- Ps],
                ?assertEqual(lists:sort(Keys), Keys),
                ?assertEqual(D, length(lists:usort(Ps)))
            end) || {M, C} <- [{1, 7}, {2, 6}, {3, 5}, {4, 4}, {6, 3}]].

%%% States (FR-G09, FR-G07, A.12)

states_test() ->
    {M0, C0} = pq_state:mean_cov(pq_state:vacuum(2, 2)),
    ?assertEqual([0.0, 0.0, 0.0, 0.0], M0),
    ?assertEqual([[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]], C0),
    {Mc, _} = pq_state:mean_cov(pq_state:coherent(1, [{complex, 0.7, 0.2}], 2)),
    ?assertEqual([1.4, 0.4], [round(X * 1.0e12) / 1.0e12 || X <- Mc]),
    R = 0.5,
    {_, [[Vx, Cxp], [_, Vp]]} = pq_state:mean_cov(pq_state:squeezed_vacuum(1, [R], 2)),
    ?assert(abs(Vx - math:exp(-2 * R)) < 1.0e-15 andalso abs(Vp - math:exp(2 * R)) < 1.0e-14 andalso abs(Cxp) < 1.0e-15),
    {_, [[T | _] | _]} = pq_state:mean_cov(pq_state:thermal(2, [1.5, 0], 2)),
    ?assertEqual(4.0, T),
    ?assertEqual(ok, pq_state:check_physical(pq_state:thermal(2, [1.5, 0], 2))),
    ?assertMatch({error, {non_physical_state, _}}, pq_state:check_physical(pq_state:from_mean_cov([0, 0], [[0.5, 0], [0, 0.5]], 2))),
    S = pq_state:coherent(3, [1, {complex, 0, 2}, 3], 2),
    ?assertEqual(pq_state:mean_cov(S), pq_state:mean_cov(pq_state:from_xpxp(pq_state:to_xpxp(S), 2))),
    {Xp, _} = pq_state:to_xpxp(S),
    ?assertEqual([2.0, 0.0, 0.0, 4.0, 6.0, 0.0], Xp),
    {ok, Red} = pq_state:reduce(S, [2, 0]),
    ?assertEqual([6.0, 2.0, 0.0, 0.0], element(1, pq_state:mean_cov(Red))).

%%% Lowering (FR-C07) and capabilities (FS 5.3)

c(N, Instrs) -> {ok, C} = pq_circuit:new(N), C#circuit{instrs = Instrs}.

lower_test() ->
    C = (c(2, [{beamsplitter, 0, 1, {sym, t}, {mul, [2, {sym, p}]}}, {unitary, [0, 1], {matrix, [[{cos, {sym, a}}, 0], [0, 1]]}},
               {measure_homodyne, 0, {sym, p}}]))#circuit{params = #{t => {mul, [{rat, 1, 4}, {const, pi}]}}},
    {ok, L} = pq_circuit:lower(C, #{p => 0.5, a => 0, zz => 1}),
    [{beamsplitter, 0, 1, T, P}, {unitary, [0, 1], {cmat, 2, 2, _}}, {measure_homodyne, 0, 0.5}] = L#circuit.instrs,
    ?assert(abs(T - math:pi() / 4) < 1.0e-15 andalso P == 1.0),
    ?assertMatch(#{lowered := true, warnings := [{unused_bindings, [zz]}]}, L#circuit.meta),
    ?assertEqual(#{}, L#circuit.params),
    ?assertEqual({error, {unbound_symbols, [a, p]}}, pq_circuit:lower(C, #{})),
    ?assertMatch({error, {unsupported_input, #{complex_parameter := _}}},
                 pq_circuit:lower(c(1, [{phaseshift, 0, {complex, 0, 1}}]), #{})),
    %% a symbolic unitary that is not unitary under the binding is caught after lowering
    ?assertMatch({error, {not_unitary, _}},
                 pq_circuit:lower(c(2, [{unitary, [0, 1], {matrix, [[{sym, a}, 0], [0, 1]]}}]), #{a => 2})).

custom_test() ->
    application:set_env(pq, custom_instructions, [pq_custom_example]),
    try
        {ok, L} = pq_circuit:lower(c(2, [{custom, half_mirror, #{modes => [0, 1]}}]), #{}),
        ?assertMatch([{beamsplitter, 0, 1, _, +0.0}], L#circuit.instrs)
    after
        application:unset_env(pq, custom_instructions)
    end.

capabilities_test() ->
    ?assertEqual({error, {unsupported_instruction, kerr, gaussian}},
                 pq_circuit:check_capabilities(c(1, [{kerr, 0, 1}]), pq_gaussian)),
    ?assertMatch({error, {unsupported_instruction, measure_pnr, #{placement := midcircuit}}},
                 pq_circuit:check_capabilities(c(2, [{measure_pnr, [0]}, {phaseshift, 1, 0}]), pq_gaussian)),
    ?assertEqual(ok, pq_circuit:check_capabilities(c(2, [{measure_homodyne, 0, 0}, {phaseshift, 1, 0}]), pq_gaussian)).

%%% The Gaussian backend (A.3-A.6, FR-G05)

means(C, Opts) ->
    {ok, S, []} = pq_simulator:evolve(pq_gaussian, C, Opts),
    pq_state:mean_cov(S).

passive_checks_test() ->
    In = {gaussian, pq_state:coherent(2, [1, {complex, 0, 0.5}], 2)},
    {ok, C0} = pq_circuit:new(2, In),
    %% A.3: phase shifter x' = cos x - sin p
    Ph = 0.3,
    {[X, _, P, _], _} = means(C0#circuit{instrs = [{phaseshift, 0, Ph}]}, #{}),
    ?assert(abs(X - (math:cos(Ph) * 2)) < 1.0e-15 andalso abs(P - math:sin(Ph) * 2) < 1.0e-15),
    %% beamsplitter at Phi = 0: x1' = cos x1 - sin x2, x2' = sin x1 + cos x2 (and p the same)
    T = 0.4,
    {[X1, X2, P1, P2], _} = means(C0#circuit{instrs = [{beamsplitter, 0, 1, T, 0}]}, #{}),
    ?assert(abs(X1 - math:cos(T) * 2) < 1.0e-15 andalso abs(X2 - math:sin(T) * 2) < 1.0e-15),
    ?assert(abs(P1 + math:sin(T) * 1) < 1.0e-15 andalso abs(P2 - math:cos(T) * 1) < 1.0e-15),
    %% the beamsplitter's phase: alpha_b' = e^(i phi) sin t alpha_a for alpha = (1, 0)
    {ok, C1} = pq_circuit:new(2, {gaussian, pq_state:coherent(2, [1, 0], 2)}),
    {[_, Xb, _, Pb], _} = means(C1#circuit{instrs = [{beamsplitter, 0, 1, T, 1.1}]}, #{}),
    ?assert(abs(Xb - 2 * math:cos(1.1) * math:sin(T)) < 1.0e-15 andalso abs(Pb - 2 * math:sin(1.1) * math:sin(T)) < 1.0e-15),
    ?assertEqual({error, {unsupported_input, fock}},
                 pq_simulator:evolve(pq_gaussian, (c(1, []))#circuit{input = {fock, [1]}}, #{})).

homodyne_test() ->
    R = 0.5,
    C = c(2, [{squeezing, 0, R, 0}, {displacement, 1, 0.5, math:pi() / 2}, {loss, 1, 0.5},
              {measure_homodyne, 0, 0}, {measure_homodyne, 1, math:pi() / 2}]),
    {ok, #result{kind = homodyne, measured = [0, 1], samples = S, metadata = Meta}} =
        pq_simulator:run(pq_gaussian, C, #{shots => 100000, seed => 42}),
    ?assertMatch(#{seed := 42, shots_requested := 100000, shots_accepted := 100000, backend := gaussian}, Meta),
    N = length(S),
    Mean = fun(F) -> lists:sum([F(Row) || Row <- S]) / N end,
    VarX = Mean(fun([A, _]) -> A * A end),
    %% mode 1: displaced by 0.5 i (p mean 2 * 0.5 = 1), loss 0.5 (mean sqrt(0.5)), measured along p
    MeanP = Mean(fun([_, B]) -> B end),
    ?assert(abs(VarX / math:exp(-2 * R) - 1) < 0.02),
    ?assert(abs(MeanP - math:sqrt(0.5)) < 0.01),
    %% the same seed gives the same samples; a missing seed is drawn and recorded
    {ok, #result{samples = S2}} = pq_simulator:run(pq_gaussian, C, #{shots => 5, seed => 42}),
    ?assertEqual(lists:sublist(S, 5), S2),
    {ok, #result{metadata = #{seed := Drawn}}} = pq_simulator:run(pq_gaussian, C, #{shots => 1}),
    ?assert(is_integer(Drawn)),
    ?assertMatch({error, {unsupported_instruction, measure_pnr, _}},
                 pq_simulator:run(pq_gaussian, c(1, [{measure_pnr, [0]}]), #{shots => 1})).
