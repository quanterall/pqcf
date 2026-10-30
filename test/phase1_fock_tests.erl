%% Phase 1, Fock milestone: the pure-state Fock backend end to end, GT-B01 (bunched probabilities against the
%% permanent kernel), the Clements mesh replayed (FR-C10), truncation (FS 3.5), statistics (FR-X01).
-module(phase1_fock_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

c(N, In, Instrs) -> {ok, C} = pq_circuit:new(N, In), C#circuit{instrs = Instrs}.

%% a random k x k unitary by Gram-Schmidt (test data only)
haar(K, Seed) ->
    rand:seed(exsss, {Seed, 2, 3}),
    Cols = [[{rand:normal(), rand:normal()} || _ <- lists:seq(1, K)] || _ <- lists:seq(1, K)],
    Orth = lists:foldl(fun(V, Acc) ->
                               W = lists:foldl(fun(Q, X) -> sub(X, scale(Q, dot(Q, X))) end, V, Acc),
                               Acc ++ [scale(W, {1 / math:sqrt(re(dot(W, W))), 0.0})]
                       end, [], Cols),
    %% columns -> rows
    Rows = [[lists:nth(I, Col) || Col <- Orth] || I <- lists:seq(1, K)],
    {ok, U} = pq_core:pack({cmat, [[{complex, R, Im} || {R, Im} <- Row] || Row <- Rows]}),
    U.

mul({A, B}, {C, D}) -> {A * C - B * D, A * D + B * C}.
dot(X, Y) -> lists:foldl(fun({{A, B}, Z}, {S, T}) -> {P, Q} = mul({A, -B}, Z), {S + P, T + Q} end, {0.0, 0.0}, lists:zip(X, Y)).
scale(V, S) -> [mul(Z, S) || Z <- V].
sub(X, Y) -> [{A - C, B - D} || {{A, B}, {C, D}} <- lists:zip(X, Y)].
re({A, _}) -> A.

hom_test() ->
    C = c(2, {fock, [1, 1]}, [{beamsplitter, 0, 1, math:pi() / 4, 0}, {measure_pnr, [0, 1]}]),
    {ok, #result{kind = pnr, measured = [0, 1], samples = S, metadata = M}} =
        pq_simulator:run(pq_fock, C, #{shots => 4000, seed => 1}),
    ?assertNot(lists:member([1, 1], S)),
    ?assert(abs(length([x || [2, 0] <- S]) / 4000 - 0.5) < 0.03),
    ?assertMatch(#{cutoff := 2, backend := fock}, M).

%% GT-B01: Fock-backend probabilities equal |perm|^2 / (S! T!) from the permanent kernel, bunched input
gt_b01_test() ->
    U = haar(4, 7),
    In = [2, 1, 1, 0],
    {ok, S, []} = pq_simulator:evolve(pq_fock, c(4, {fock, In}, [{unitary, [0, 1, 2, 3], U}]), #{}),
    {ok, PP} = pq_nif:fock_pnr_probabilities(S#fock_state.amps, [0, 1, 2, 3], 4, 4),
    Four = [{P, Pr} || {P, Pr} <- PP, lists:sum(P) =:= 4],
    ?assertEqual(35, length(Four)),
    [begin
         {ok, {Want, _}} = pq_nif:bs_probability(U, In, P, #{}),
         ?assert(abs(Pr - Want) < 1.0e-12)
     end || {P, Pr} <- Four],
    ?assert(abs(lists:sum([Pr || {_, Pr} <- PP]) - 1) < 1.0e-12).

%% FR-C10: the mesh's instructions on the Gaussian backend map coherent amplitudes alpha -> U alpha
mesh_replay_test() ->
    K = 5,
    U = haar(K, 11),
    {ok, Instrs} = pq_circuit:mesh(U),
    ?assertEqual(K * (K - 1) div 2 + K, length(Instrs)),
    Alpha = [{complex, 0.1 * J, 0.05 * (K - J)} || J <- lists:seq(1, K)],
    {ok, S, []} = pq_simulator:evolve(pq_gaussian, c(K, {gaussian, pq_state:coherent(K, Alpha, 2)}, Instrs), #{}),
    {Mu, _} = pq_state:mean_cov(S),
    {ok, {cmat, Rows}} = pq_core:unpack(U),
    Want = [lists:foldl(fun({{complex, A, B}, {complex, C, D}}, {X, Y}) -> {X + A * C - B * D, Y + A * D + B * C} end,
                        {0.0, 0.0}, lists:zip(Row, Alpha)) || Row <- Rows],
    Got = lists:zip(lists:sublist(Mu, K), lists:nthtail(K, Mu)),
    [?assert(abs(Gx - 2 * Wx) + abs(Gp - 2 * Wp) < 1.0e-12) || {{Gx, Gp}, {Wx, Wp}} <- lists:zip(Got, Want)].

truncation_test() ->
    C = c(1, vacuum, [{displacement, 0, 1.5, 0}, {measure_pnr, [0]}]),
    ?assertEqual({error, {unsupported_input, cutoff_required}}, pq_simulator:run(pq_fock, C, #{shots => 1})),
    {ok, #result{metadata = #{warnings := [{truncation, #{deficit := D, instr_index := 0}}]}}} =
        pq_simulator:run(pq_fock, C, #{shots => 10, cutoff => 3}),
    ?assert(D > 0.1 andalso D < 0.5),
    ?assertMatch({error, {truncation, _}}, pq_simulator:run(pq_fock, C, #{shots => 1, cutoff => 3, strict => true})),
    {ok, #result{metadata = #{warnings := []}}} = pq_simulator:run(pq_fock, C, #{shots => 1, cutoff => 40}),
    ?assertMatch({error, {dimension_limit, 1001, 10}},
                 pq_simulator:run(pq_fock, c(4, vacuum, [{measure_pnr, [0]}]), #{cutoff => 10, max_fock_dim => 10})),
    ?assertEqual({error, {unsupported_instruction, loss, fock}},
                 pq_simulator:run(pq_fock, c(1, {fock, [1]}, [{loss, 0, 0.5}, {measure_pnr, [0]}]), #{})),
    ?assertMatch({error, {unsupported_instruction, midcircuit_measurement, _}},
                 pq_simulator:run(pq_fock, c(2, {fock, [1, 0]}, [{measure_pnr, [0]}, {phaseshift, 1, 0}]), #{})).

detectors_test() ->
    C = c(3, {fock, [1, 1, 0]}, [{beamsplitter, 0, 1, math:pi() / 4, 0}, {measure_threshold, [0, 1]}, {measure_pnr, [2]}]),
    {ok, #result{kind = mixed, measured = [0, 1, 2], samples = S, metadata = M}} =
        pq_simulator:run(pq_fock, C, #{shots => 500, seed => 3}),
    ?assertEqual([threshold, threshold, pnr], maps:get(column_kinds, M)),
    ?assertEqual([[0, 1, 0], [1, 0, 0]], lists:usort(S)),
    {ok, #result{samples = S2}} = pq_simulator:run(pq_fock, C, #{shots => 500, seed => 3}),
    ?assertEqual(S, S2).

mean_photon_number_test() ->
    {ok, N1} = pq_stats:mean_photon_number(pq_state:coherent(2, [{complex, 0.6, 0.8}, 0], 2), [0]),
    ?assert(abs(N1 - 1.0) < 1.0e-14),
    R = 0.7,
    {ok, N2} = pq_stats:mean_photon_number(pq_state:squeezed_vacuum(1, [R], 2), [0]),
    ?assert(abs(N2 - math:pow(math:sinh(R), 2)) < 1.0e-14),
    {ok, S, []} = pq_simulator:evolve(pq_fock, c(1, vacuum, [{squeezing, 0, R, 0}]), #{cutoff => 60}),
    {ok, N3} = pq_stats:mean_photon_number(S, [0]),
    ?assert(abs(N3 - math:pow(math:sinh(R), 2)) < 1.0e-12),
    {ok, Res} = pq_simulator:run(pq_fock, c(2, {fock, [2, 0]}, [{measure_pnr, [0, 1]}]), #{shots => 3}),
    ?assertEqual({ok, [2.0, 0.0]}, pq_stats:mean_photon_number(Res, [0, 1])).
