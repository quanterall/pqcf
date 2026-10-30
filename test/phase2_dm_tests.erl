%% Phase 2: density matrices and loss (FR-F05, A.11). GT-N01: loss on Gaussian states against Fock density
%% matrices; promote_dm; max_dm_dim; sampling from rho.
-module(phase2_dm_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

c(N, In, Instrs) -> {ok, C} = pq_circuit:new(N, In), C#circuit{instrs = Instrs}.

%% GT-N01: the photon statistics after loss agree between the Gaussian backend (exact) and the density matrix
gt_n01_test_() ->
    Cases = [{1, [{squeezing, 0, 0.4, 0.3}, {displacement, 0, 0.5, 1.0}, {loss, 0, 0.6}], 24},
             {2, [{squeezing, 0, 0.3, 0}, {beamsplitter, 0, 1, 0.7, 0.2}, {loss, 1, 0.5}, {loss, 0, 0.8}], 10}],
    [?_test(begin
                {ok, G, []} = pq_simulator:evolve(pq_gaussian, c(N, vacuum, Instrs), #{}),
                {ok, Rho, []} = pq_simulator:evolve(pq_fock_dm, c(N, vacuum, Instrs), #{cutoff => Cut}),
                Modes = lists:seq(0, N - 1),
                {ok, PP} = pq_nif:fock_pnr_probabilities(Rho#fock_dm.rho, Modes, N, Cut),
                Low = [{P, Pr} || {P, Pr} <- PP, lists:sum(P) =< 4],
                ?assert(length(Low) >= 5),
                [begin
                     {ok, Want} = pq_nif:gaussian_probability(G#gaussian_state.mean, G#gaussian_state.cov, 2.0, P, pnr),
                     ?assert(abs(Pr - Want) < 1.0e-9)
                 end || {P, Pr} <- Low]
            end) || {N, Instrs, Cut} <- Cases].

promote_test() ->
    C = c(1, {fock, [2]}, [{loss, 0, 0.7}, {measure_pnr, [0]}]),
    ?assertEqual({error, {unsupported_instruction, loss, fock}}, pq_simulator:run(pq_fock, C, #{shots => 1})),
    {ok, #fock_dm{rho = R}, []} = pq_simulator:evolve(pq_fock, C, #{promote_dm => true}),
    {ok, [{[0], P0}, {[1], P1}, {[2], P2}]} = pq_nif:fock_pnr_probabilities(R, [0], 1, 2),
    ?assert(abs(P0 - 0.09) + abs(P1 - 0.42) + abs(P2 - 0.49) < 1.0e-14),
    {ok, #result{samples = S}} = pq_simulator:run(pq_fock_dm, C, #{shots => 20000, seed => 8}),
    ?assert(abs(length([x || [1] <- S]) / 20000 - 0.42) < 0.015).

limits_test() ->
    ?assertMatch({error, {dimension_limit, 5005, 4096}},
                 pq_simulator:run(pq_fock_dm, c(6, vacuum, [{measure_pnr, [0]}]), #{cutoff => 9, shots => 1})),
    %% the trace deficit of a displacement under a small cutoff is a truncation warning
    {ok, #result{metadata = #{warnings := [{truncation, _}]}}} =
        pq_simulator:run(pq_fock_dm, c(1, vacuum, [{displacement, 0, 1.5, 0}, {measure_pnr, [0]}]), #{cutoff => 3, shots => 1}).
