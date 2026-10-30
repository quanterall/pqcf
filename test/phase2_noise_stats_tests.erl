%% Phase 2: noise insertion (FS 5.8) and statistics (FS 5.12, GT-X01).
-module(phase2_noise_stats_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

c(N, In, Instrs) -> {ok, C} = pq_circuit:new(N, In), C#circuit{instrs = Instrs}.

noise_test() ->
    C = c(2, vacuum, [{squeezing, 0, 0.5, 0}, {beamsplitter, 0, 1, 0.3, 0}, {measure_pnr, [0, 1]}]),
    {ok, L} = pq_noise:insert_loss(C, #{0 => 0.9, 1 => {sym, eta}}, #{}),
    ?assertEqual([{squeezing, 0, 0.5, 0}, {loss, 0, 0.9}, {beamsplitter, 0, 1, 0.3, 0}, {loss, 0, 0.9},
                  {loss, 1, {sym, eta}}, {measure_pnr, [0, 1]}], L#circuit.instrs),
    {ok, D} = pq_noise:detector_efficiency(C, 0.8),
    ?assertEqual([{loss, 0, 0.8}, {loss, 1, 0.8}, {measure_pnr, [0, 1]}], lists:nthtail(2, D#circuit.instrs)),
    %% detector efficiency on a Gaussian state: mean photon number scales by eta
    {ok, S, []} = pq_simulator:evolve(pq_gaussian, D, #{}),
    {ok, N0} = pq_stats:mean_photon_number(pq_state:squeezed_vacuum(1, [0.5], 2), [0]),
    {ok, N1} = pq_stats:mean_photon_number(S, [0, 1]),
    ?assert(abs(N1 - 0.8 * N0) < 1.0e-14).

%% GT-X01
gt_x01_test() ->
    {ok, Gc} = pq_stats:g2(pq_state:coherent(1, [{complex, 0.8, 0.3}], 2), 0),
    ?assert(abs(Gc - 1) < 1.0e-9),
    {ok, Gt} = pq_stats:g2(pq_state:thermal(1, [0.7], 2), 0),
    ?assert(abs(Gt - 2) < 1.0e-12),
    R = 0.6,
    {ok, Gs} = pq_stats:g2(pq_state:squeezed_vacuum(1, [{complex, R * math:cos(0.4), R * math:sin(0.4)}], 2), 0),
    ?assert(abs(Gs - (3 + 1 / math:pow(math:sinh(R), 2))) < 1.0e-12),
    %% the same squeezed state in Fock space
    {ok, Sq, []} = pq_simulator:evolve(pq_fock, c(1, vacuum, [{squeezing, 0, R, 0.4}]), #{cutoff => 60}),
    {ok, Gf} = pq_stats:g2(Sq, 0),
    ?assert(abs(Gf - Gs) < 1.0e-9),
    {ok, F1, []} = pq_simulator:evolve(pq_fock, c(1, {fock, [1]}, []), #{}),
    ?assertEqual({ok, 0.0}, pq_stats:g2(F1, 0)),
    ?assertEqual({error, {undefined_statistic, zero_mean}}, pq_stats:g2(pq_state:vacuum(1, 2), 0)),
    ?assertEqual({ok, 0.0}, pq_stats:tv_distance(#{[0] => 0.5, [1] => 0.5}, #{[1] => 0.5, [0] => 0.5})),
    ?assertEqual({ok, 0.5}, pq_stats:tv_distance(#{[0] => 1.0}, #{[0] => 0.5, [1] => 0.5})).

result_stats_test() ->
    C = c(2, {fock, [1, 1]}, [{beamsplitter, 0, 1, math:pi() / 4, 0}, {measure_pnr, [0, 1]}]),
    {ok, R} = pq_simulator:run(pq_fock, C, #{shots => 4000, seed => 2}),
    {ok, H} = pq_stats:histogram(R, [0, 1]),
    ?assertEqual([[0, 2], [2, 0]], lists:sort(maps:keys(H))),
    {ok, Tv} = pq_stats:tv_distance(R, #{[2, 0] => 0.5, [0, 2] => 0.5}),
    ?assert(Tv < 0.03),
    {ok, G} = pq_stats:g2(R, 0),
    ?assert(abs(G - 1) < 0.1),                         % half the shots 2 photons: <n(n-1)> = 1, <n> = 1
    {ok, A, []} = pq_simulator:evolve(pq_fock, c(2, {fock, [1, 0]}, [{beamsplitter, 0, 1, 0.3, 0}]), #{}),
    {ok, B, []} = pq_simulator:evolve(pq_fock, c(2, {fock, [1, 0]}, [{beamsplitter, 0, 1, 0.3, 0}]), #{}),
    ?assert(abs(element(2, pq_stats:fidelity(A, B)) - 1) < 1.0e-14),
    Rho = pq_fock_dm:from_pure(B),
    ?assert(abs(element(2, pq_stats:fidelity(A, Rho)) - 1) < 1.0e-14),
    {ok, Z, []} = pq_simulator:evolve(pq_fock, c(2, {fock, [0, 1]}, []), #{cutoff => 1}),
    ?assert(abs(element(2, pq_stats:fidelity(A, Z)) - math:pow(math:sin(0.3), 2)) < 1.0e-14).
