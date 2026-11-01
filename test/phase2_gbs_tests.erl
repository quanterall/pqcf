%% Phase 2, Gaussian boson sampling through the Gaussian backend (FS 5.7, A.9, A.10; AN 6). The kernels are
%% checked against thewalrus in the Rust tests (gt_g01).
-module(phase2_gbs_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

c(N, Instrs) -> {ok, C} = pq_circuit:new(N), C#circuit{instrs = Instrs}.

pnr_test() ->
    R = 0.5,
    C = c(2, [{squeezing, 0, R, 0}, {measure_pnr, [0, 1]}]),
    {ok, #result{kind = pnr, measured = [0, 1], samples = S, metadata = M}} =
        pq_simulator:run(pq_gaussian, C, #{shots => 20000, seed => 4}),
    ?assert(lists:all(fun([A, 0]) -> A rem 2 =:= 0; (_) -> false end, S)),
    P0 = length([x || [0, 0] <- S]) / 20000,
    ?assert(abs(P0 - 1 / math:cosh(R)) < 0.01),
    ?assertMatch(#{truncated_mass := T, warnings := []} when T < 1.0e-6, M).

threshold_and_caps_test() ->
    R = 1.0,
    C = c(1, [{squeezing, 0, R, 0}, {measure_threshold, [0]}]),
    {ok, #result{kind = threshold, samples = S}} = pq_simulator:run(pq_gaussian, C, #{shots => 20000, seed => 2}),
    Click = length([x || [1] <- S]) / 20000,
    ?assert(abs(Click - (1 - 1 / math:cosh(R))) < 0.01),
    %% tight caps: the lost mass is reported, and is an error under strict
    Cp = c(1, [{squeezing, 0, R, 0}, {measure_pnr, [0]}]),
    {ok, #result{metadata = #{truncated_mass := T, warnings := [{gbs_truncation, _}]}}} =
        pq_simulator:run(pq_gaussian, Cp, #{shots => 2000, seed => 2, max_photons => 2}),
    ?assert(T > 0.05),
    ?assertMatch({error, {gbs_truncation, _}},
                 pq_simulator:run(pq_gaussian, Cp, #{shots => 200, seed => 2, max_photons => 2, strict => true})),
    %% displaced threshold detection (the loop torontonian, FS 15 item 2): a coherent state clicks with
    %% probability 1 - exp(-|alpha|^2), |alpha| = 0.5 here
    {ok, #result{kind = threshold, samples = Clicks}} =
        pq_simulator:run(pq_gaussian, c(1, [{displacement, 0, 0.5, 0}, {measure_threshold, [0]}]), #{shots => 20000, seed => 3}),
    Pc = 1 - math:exp(-0.25),
    F = length([x || [1] <- Clicks]) / 20000,
    ?assert(abs(F - Pc) < 5 * math:sqrt(Pc * (1 - Pc) / 20000)).

amat_test() ->
    %% coherent state: A = 0, gamma = conj(beta), Norm = exp(|alpha|^2 / ... ) such that P(0) = exp(-|alpha|^2)
    S = pq_state:coherent(1, [{complex, 0.7, 0.2}], 2),
    {ok, {{cmat, 2, 2, _} = A, {cvec, 2, _}, Norm}} = pq_nif:gaussian_to_amat(S#gaussian_state.mean, S#gaussian_state.cov, 2.0),
    {ok, {cmat, Rows}} = pq_core:unpack(A),
    ?assert(lists:all(fun({complex, Re, Im}) -> abs(Re) + abs(Im) < 1.0e-15 end, lists:append(Rows))),
    ?assert(abs(1 / Norm - math:exp(-0.53)) < 1.0e-15),
    {ok, P1} = pq_nif:gaussian_probability(S#gaussian_state.mean, S#gaussian_state.cov, 2.0, [1], pnr),
    ?assert(abs(P1 - 0.53 * math:exp(-0.53)) < 1.0e-15).
