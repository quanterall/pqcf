%% Requirements found unimplemented in the Phase 3 requirement audit: FR-B06 (collision-free boson sampling as
%% post-selection; boson sampling also ignored postselect), FR-G06 (physicality after every instruction with
%% check_physical), FR-F01 (Fock basis states, superpositions within the cutoff, density matrices from pure states).
-module(phase3_requirements_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

c(N, In, Is) -> {ok, C} = pq_circuit:new(N, In), C#circuit{instrs = Is}.

boson() ->
    c(4, {fock, [1, 1, 1, 0]}, [{beamsplitter, J, J + 1, 0.4 + 0.15 * J + 0.1 * L, 0.3 * L - 0.2 * J}
                                || L <- lists:seq(0, 3), J <- lists:seq(L rem 2, 2, 2)] ++ [{measure_pnr, [0, 1, 2, 3]}]).

requirements_test_() ->
    {setup, fun() -> {ok, _} = application:ensure_all_started(pq) end, fun(_) -> ok end,
     [{timeout, 120, {"FR-B06: collision-free samples follow the conditional distribution", fun collision_free/0}},
      {timeout, 60, {"FR-B06: accepted mode, and combined with postselect", fun collision_free_modes/0}},
      {"boson sampling honours postselect (it was ignored)", fun boson_postselect/0},
      {"FR-G06: check_physical after every instruction", fun physical/0},
      {"FR-F01: basis states and superpositions", fun superposition/0},
      {"FR-F01: density matrices from pure states and mixtures", fun mixture/0}]}.

collision_free() ->
    C = boson(),
    {ok, Exact, _} = pq_api:probabilities(C, #{backend => boson}),
    Free = maps:filter(fun(K, _) -> lists:max(K) =< 1 end, Exact),
    Z = lists:sum(maps:values(Free)),
    Cond = maps:map(fun(_, V) -> V / Z end, Free),
    {ok, #result{samples = S, metadata = M}} = pq_api:sample(boson, C, #{shots => 1000000, seed => 21, collision_free => true}),
    ?assert(lists:all(fun(X) -> lists:max(X) =< 1 end, S)),
    #{attempts := 1000000, shots_accepted := N, acceptance_rate := A, collision_free := true} = M,
    ?assertEqual(length(S), N),
    %% the acceptance rate estimates the collision-free probability (binomial, 6 sigma)
    ?assert(abs(A - Z) < 6 * math:sqrt(Z * (1 - Z) / 1000000)),
    {Pv, Info} = pq_chi2:p_value(Cond, pq_chi2:counts(S)),
    io:format(user, "~n  collision-free: acceptance ~.4f (exact ~.4f), p = ~.3g ~p", [A, Z, Pv, Info]),
    ?assert(Pv > 1.0e-4).

collision_free_modes() ->
    C = boson(),
    {ok, #result{samples = S, metadata = #{attempts := At}}} =
        pq_api:sample(boson, C, #{shots => 5000, seed => 3, collision_free => true, shots_mode => accepted}),
    ?assertEqual(5000, length(S)),
    ?assert(At > 5000),
    %% the accepted samples are the collision-free ones among the plain run's, in order (the attempt index is the
    %% stream index)
    {ok, #result{samples = Plain}} = pq_api:sample(boson, C, #{shots => At, seed => 3}),
    ?assertEqual(lists:sublist([X || X <- Plain, lists:max(X) =< 1], 5000), S),
    {ok, #result{samples = S2}} = pq_api:sample(boson, C, #{shots => 20000, seed => 4, collision_free => true,
                                                            postselect => {eq, 3, 0}}),
    ?assert(lists:all(fun([_, _, _, D] = X) -> D =:= 0 andalso lists:max(X) =< 1 end, S2)),
    ?assert(length(S2) > 0).

boson_postselect() ->
    {ok, #result{samples = S, metadata = M}} = pq_api:sample(boson, boson(), #{shots => 10000, seed => 5, postselect => {eq, 0, 1}}),
    ?assert(lists:all(fun([A | _]) -> A =:= 1 end, S)),
    ?assert(length(S) < 10000),
    ?assertEqual(10000, maps:get(attempts, M)),
    ?assertMatch({error, {postselection_exhausted, _}},
                 pq_api:sample(boson, boson(), #{shots => 10, seed => 5, postselect => {eq, 0, 4}, shots_mode => accepted,
                                                 max_attempts => 200})).

physical() ->
    %% a state that is not physical (covariance below the uncertainty bound), passed without the init check
    {ok, C0} = pq_circuit:new(1),
    Bad = pq_state:from_mean_cov([0.0, 0.0], [[0.5, 0.0], [0.0, 0.5]], 2),
    ?assertMatch({error, {non_physical_state, #{after_instruction := {phaseshift, 0, 0.3}}}},
                 pq_gaussian:apply({phaseshift, 0, 0.3}, Bad, #{check_physical => true})),
    {ok, _, []} = pq_gaussian:apply({phaseshift, 0, 0.3}, Bad, #{}),
    %% a physical run with the check on is unchanged
    C = C0#circuit{instrs = [{squeezing, 0, 0.5, 0}, {displacement, 0, 0.3, 0.1}, {loss, 0, 0.6}]},
    {ok, S1, _} = pq_api:evolve(C, #{check_physical => true}),
    {ok, S2, _} = pq_api:evolve(C, #{}),
    ?assertEqual(S1, S2).

superposition() ->
    {ok, B} = pq_state:fock_state([1, 0], 2),
    ?assertEqual(6, element(2, B#fock_state.amps)),
    %% (|1,0> + i |0,1>) / sqrt 2 through a 50:50 beamsplitter
    {ok, S} = pq_state:fock_superposition([{1, [1, 0]}, {{complex, 0, 1}, [0, 1]}], 1),
    {ok, {cvec, A}} = pq_core:unpack(S#fock_state.amps),
    ?assert(abs(lists:sum([X * X + Y * Y || {complex, X, Y} <- A]) - 1) < 1.0e-15),
    C = c(2, {fock_state, S}, [{beamsplitter, 0, 1, math:pi() / 4, 0}, {measure_pnr, [0, 1]}]),
    {ok, P, _} = pq_api:probabilities(C, #{}),
    %% B(pi/4, 0) takes a0 -> (a0 + a1)/sqrt2, a1 -> (-a0 + a1)/sqrt2: amplitudes (1 - i)/2 and (1 + i)/2
    ?assert(abs(maps:get([1, 0], P) - 0.5) < 1.0e-12),
    ?assert(abs(maps:get([0, 1], P) - 0.5) < 1.0e-12),
    {ok, Sm} = pq_state:fock_superposition([{1, [1, 0]}, {1, [1, 0]}, {2, [0, 1]}], 1),
    {ok, {cvec, [_, {complex, A1, _}, {complex, A2, _}]}} = pq_core:unpack(Sm#fock_state.amps),
    ?assert(abs(A1 - A2) < 1.0e-15),
    ?assertMatch({error, {unsupported_input, _}}, pq_state:fock_state([2, 1], 2)),
    ?assertMatch({error, {bad_shape, _}}, pq_state:fock_superposition([{1, [1, 0]}, {1, [1]}], 2)),
    ?assertMatch({error, {bad_shape, zero_norm}}, pq_state:fock_superposition([{1, [1, 0]}, {-1, [1, 0]}], 2)).

mixture() ->
    {ok, A} = pq_state:fock_state([1, 0], 1),
    {ok, B} = pq_state:fock_state([0, 1], 1),
    {ok, #fock_dm{rho = Ra}} = pq_state:fock_dm(A),
    {ok, {cmat, RowsA}} = pq_core:unpack(Ra),
    {ok, Ia} = pq_state:index([1, 0], 2, 1),
    {ok, Ib} = pq_state:index([0, 1], 2, 1),
    ?assertEqual({complex, 1.0, 0.0}, lists:nth(Ia + 1, lists:nth(Ia + 1, RowsA))),
    {ok, #fock_dm{rho = R}} = pq_state:fock_dm([{3, A}, {1, B}]),
    {ok, {cmat, Rows}} = pq_core:unpack(R),
    Diag = [Re || {K, Row} <- lists:enumerate(1, Rows), {complex, Re, _} <- [lists:nth(K, Row)]],
    ?assertEqual(0.75, lists:nth(Ia + 1, Diag)),
    ?assertEqual(0.25, lists:nth(Ib + 1, Diag)),
    ?assertEqual(1.0, lists:sum(Diag)),
    {ok, Mx} = pq_state:fock_dm([{3, A}, {1, B}]),
    {ok, Pp} = pq_nif:fock_pnr_probabilities(Mx#fock_dm.rho, [0, 1], 2, 1),
    [{[0, 1], P01}, {[1, 0], P10}] = lists:sort([X || {_, V} = X <- Pp, V > 0]),
    ?assert(abs(P01 - 0.25) < 1.0e-15 andalso abs(P10 - 0.75) < 1.0e-15),
    ?assertMatch({error, {bad_shape, mixture}}, pq_state:fock_dm([{-1, A}])).
