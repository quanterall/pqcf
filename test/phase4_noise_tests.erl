%% Phase 4, item 4 (FS 15 item 4): noise models (docs/decisions/N1-noise-models.md). Partial distinguishability of
%% boson-sampling photons and dark counts of the detectors. Exit: the Hong-Ou-Mandel dip against the visibility,
%% closed forms for dark counts, chi-square of the samplers.
-module(phase4_noise_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

c(N, In, Is) -> {ok, C} = pq_circuit:new(N, In), C#circuit{instrs = Is}.

setup_test_() ->
    {setup, fun() -> {ok, _} = application:ensure_all_started(pq) end, fun(_) -> ok end,
     [{timeout, 300, fun hom_dip/0}, {timeout, 300, fun mixture_against_gram/0}, {timeout, 300, fun general_gram/0},
      {timeout, 120, fun chunking/0}, {timeout, 300, fun dark_threshold/0}, {timeout, 300, fun dark_pnr/0},
      {timeout, 120, fun dark_postselect/0}, fun refusals/0, fun thermal_gaussian/0, {timeout, 120, fun thermal_dm/0},
      fun thermal_refusals/0, {timeout, 300, fun multi_photon/0}]}.

hom() -> c(2, {fock, [1, 1]}, [{beamsplitter, 0, 1, math:pi() / 4, 0}, {measure_pnr, [0, 1]}]).

%% P(1,1) = (1 - V)/2 exactly; 400,000 shots at V = 0.6 fit the exact distribution
hom_dip() ->
    [begin
         {ok, P, _} = pq_api:probabilities(hom(), #{visibility => V}),
         ?assert(abs(maps:get([1, 1], P) - (1 - V) / 2) < 1.0e-14),
         ?assert(abs(maps:get([2, 0], P) - (1 + V) / 4) < 1.0e-14)
     end || V <- [0.0, 0.3, 0.6, 1.0]],
    {ok, P, _} = pq_api:probabilities(hom(), #{visibility => 0.6}),
    {ok, #result{samples = S, metadata = #{backend := boson}}} = pq_api:sample(hom(), #{shots => 400000, seed => 11, visibility => 0.6}),
    {Pv, _} = pq_chi2:p_value(P, pq_chi2:counts(S)),
    ?debugFmt("HOM at visibility 0.6: (1,1) fraction ~.4f (exact 0.2), chi-square p = ~.3f",
              [length([x || [1, 1] <- S]) / 400000, Pv]),
    ?assert(Pv > 1.0e-4).

%% three photons with different visibilities on a four-mode mesh: the mixture sampler against the Gram formula
mesh() ->
    [{beamsplitter, 0, 1, 0.7, 0.2}, {beamsplitter, 1, 2, 0.5, -0.4}, {beamsplitter, 2, 3, 0.9, 0.1},
     {beamsplitter, 0, 1, 0.4, 0.8}, {beamsplitter, 1, 2, 1.1, 0.0}, {phaseshift, 1, 0.3}, {beamsplitter, 2, 3, 0.3, 0.6}].

mixture_against_gram() ->
    C = c(4, {fock, [1, 1, 1, 0]}, mesh() ++ [{measure_pnr, [0, 1, 2, 3]}]),
    Vis = [0.9, 0.5, 0.7, 1.0],
    {ok, P, _} = pq_api:probabilities(C, #{visibility => Vis}),
    ?assert(abs(lists:sum(maps:values(P)) - 1) < 1.0e-12),
    %% the limits: V = 1 is the boson distribution, V = 0 the distinguishable one
    {ok, Pb, _} = pq_api:probabilities(C, #{backend => boson}),
    {ok, P1, _} = pq_api:probabilities(C, #{visibility => 1.0}),
    [?assert(abs(maps:get(K, P1, 0.0) - V) < 1.0e-13) || K := V <- Pb],
    {ok, #result{samples = S}} = pq_api:sample(C, #{shots => 1000000, seed => 5, visibility => Vis}),
    {Pv, _} = pq_chi2:p_value(P, pq_chi2:counts(S)),
    ?debugFmt("three photons, visibilities ~p: chi-square p = ~.3f over 1,000,000 shots", [Vis, Pv]),
    ?assert(Pv > 1.0e-4),
    %% power: the same samples against fully indistinguishable photons
    {Pv0, _} = pq_chi2:p_value(Pb, pq_chi2:counts(S)),
    ?assert(Pv0 < 1.0e-10).

%% a general Gram matrix with complex overlaps (two photons in one mode included): sampling against probabilities
general_gram() ->
    C = c(3, {fock, [2, 1, 0]}, [{beamsplitter, 0, 1, 0.6, 0.3}, {beamsplitter, 1, 2, 0.8, -0.2}, {beamsplitter, 0, 1, 0.3, 0.5},
                                 {measure_pnr, [0, 1, 2]}]),
    %% the Gram matrix of three normalised vectors, so positive semidefinite
    Vs = [[{1.0, 0.0}, {0.0, 0.0}], [{0.6, 0.0}, {0.0, 0.8}], [{0.8, 0.0}, {0.36, 0.48}]],
    Dot = fun(A, B) -> lists:foldl(fun({{Ar, Ai}, {Br, Bi}}, {R, I}) -> {R + Ar * Br + Ai * Bi, I + Ar * Bi - Ai * Br} end,
                                   {0.0, 0.0}, lists:zip(A, B)) end,
    G = [[case Dot(A, B) of {R, I} -> {complex, R, I} end || B <- Vs] || A <- Vs],
    {ok, P, _} = pq_api:probabilities(C, #{gram => G}),
    ?assert(abs(lists:sum(maps:values(P)) - 1) < 1.0e-12),
    {ok, #result{samples = S}} = pq_api:sample(C, #{shots => 300000, seed => 9, gram => G}),
    {Pv, _} = pq_chi2:p_value(P, pq_chi2:counts(S)),
    ?debugFmt("general Gram matrix: chi-square p = ~.3f", [Pv]),
    ?assert(Pv > 1.0e-4),
    %% a matrix that is not a Gram matrix (overlap above 1) gives negative probabilities: refused
    Bad = [[1, 0.99, 0.99], [0.99, 1, -0.99], [0.99, -0.99, 1]],
    ?assertMatch({error, {bad_shape, #{gram := not_positive_semidefinite}}}, pq_api:probabilities(C, #{gram => Bad})),
    ?assertMatch({error, {bad_shape, #{gram := not_hermitian_unit_diagonal}}},
                 pq_api:probabilities(C, #{gram => [[1, 0.5, 0], [0.4, 1, 0], [0, 0, 1]]})).

%% D1: the same samples for any chunking, with partial distinguishability and dark counts
chunking() ->
    C = c(4, {fock, [1, 1, 1, 0]}, mesh() ++ [{measure_pnr, [0, 1, 2, 3]}]),
    O = #{shots => 20000, seed => 77, visibility => [0.9, 0.5, 0.7, 1.0], dark_counts => 0.05},
    {ok, #result{samples = A}} = pq_api:sample(C, O#{parallel => false}),
    {ok, #result{samples = B}} = pq_api:sample(C, O#{max_inflight => 8, chunk => 37}),
    ?assertEqual(A, B).

%% threshold dark counts: an empty detector clicks with probability 1 - exp(-Nu); a lit one always clicks
dark_threshold() ->
    Nu = 0.2,
    C = c(2, {fock, [1, 0]}, [{beamsplitter, 0, 1, 0.5, 0}, {measure_threshold, [0, 1]}]),
    {ok, P0, _} = pq_api:probabilities(C, #{}),
    {ok, P, _} = pq_api:probabilities(C, #{dark_counts => Nu}),
    Q = 1 - math:exp(-Nu),
    ?assert(abs(maps:get([1, 1], P) - (maps:get([1, 0], P0) * Q + maps:get([0, 1], P0) * Q)) < 1.0e-14),
    ?assert(abs(maps:get([1, 0], P) - maps:get([1, 0], P0) * (1 - Q)) < 1.0e-14),
    ?assert(maps:get([0, 0], P, 0.0) < 1.0e-30),
    {ok, #result{samples = S}} = pq_api:sample(C, #{shots => 400000, seed => 3, dark_counts => Nu}),
    {Pv, _} = pq_chi2:p_value(P, pq_chi2:counts(S)),
    ?debugFmt("threshold dark counts: chi-square p = ~.3f", [Pv]),
    ?assert(Pv > 1.0e-4),
    %% per-mode means
    {ok, Pm, _} = pq_api:probabilities(C, #{dark_counts => #{1 => Nu}}),
    ?assert(abs(maps:get([1, 1], Pm) - maps:get([1, 0], P0) * Q) < 1.0e-14),
    ?assert(abs(maps:get([0, 1], Pm) - maps:get([0, 1], P0)) < 1.0e-14).

%% PNR dark counts: n + Poisson(Nu); the mean count rises by Nu per detector, the distribution is the convolution
dark_pnr() ->
    Nu = 0.3,
    C = c(2, {fock, [2, 1]}, [{beamsplitter, 0, 1, 0.7, 0.2}, {measure_pnr, [0, 1]}]),
    {ok, P0, _} = pq_api:probabilities(C, #{}),
    {ok, P, Meta} = pq_api:probabilities(C, #{dark_counts => Nu}),
    ?assert(maps:get(dark_count_truncated_mass, Meta) < 1.0e-12),
    %% P(0, 3) only from (0, 3) with no dark counts
    ?assert(abs(maps:get([0, 3], P) - maps:get([0, 3], P0) * math:exp(-2 * Nu)) < 1.0e-15),
    Mean0 = lists:sum([V * hd(K) || K := V <- P0]),
    Mean = lists:sum([V * hd(K) || K := V <- P]),
    ?assert(abs(Mean - Mean0 - Nu) < 1.0e-10),
    {ok, #result{samples = S}} = pq_api:sample(C, #{shots => 400000, seed => 8, dark_counts => Nu}),
    {Pv, _} = pq_chi2:p_value(P, pq_chi2:counts(S)),
    ?debugFmt("PNR dark counts: chi-square p = ~.3f", [Pv]),
    ?assert(Pv > 1.0e-4),
    %% the same on boson sampling
    {ok, #result{samples = Sb}} = pq_api:sample(boson, C, #{shots => 400000, seed => 8, dark_counts => Nu}),
    {Pvb, _} = pq_chi2:p_value(P, pq_chi2:counts(Sb)),
    ?assert(Pvb > 1.0e-4).

%% the predicate sees the noisy outcomes: collision-free acceptance falls with dark counts
dark_postselect() ->
    C = c(2, {fock, [1, 1]}, [{beamsplitter, 0, 1, 0.3, 0}, {measure_pnr, [0, 1]}]),
    {ok, #result{metadata = #{acceptance_rate := A0}}} = pq_api:sample(boson, C, #{shots => 50000, seed => 1, collision_free => true}),
    {ok, #result{metadata = #{acceptance_rate := A1}, samples = S}} =
        pq_api:sample(boson, C, #{shots => 50000, seed => 1, collision_free => true, dark_counts => 0.5}),
    ?assert(A1 < A0 - 0.1),
    ?assertEqual([], [X || X <- S, lists:max(X) > 1]).

refusals() ->
    F = c(2, {fock, [1, 1]}, [{kerr, 0, 0.1}, {beamsplitter, 0, 1, 0.3, 0}, {measure_pnr, [0, 1]}]),
    ?assertMatch({error, {unsupported_input, #{distinguishability := _}}}, pq_api:sample(F, #{shots => 10, visibility => 0.5, backend => fock})),
    ?assertMatch({error, {unsupported_input, #{distinguishability := _}}}, pq_api:probabilities(F, #{visibility => 0.5, backend => fock})),
    ?assertMatch({error, {bad_shape, #{visibility := _}}}, pq_api:sample(hom(), #{shots => 10, visibility => 1.5})),
    ?assertMatch({error, {bad_shape, _}}, pq_api:sample(hom(), #{shots => 10, visibility => 0.5, algorithm => metropolis})),
    ?assertMatch({error, {bad_shape, #{dark_counts := _}}}, pq_api:sample(hom(), #{shots => 10, dark_counts => -1})),
    Mid = c(2, {fock, [1, 0]}, [{beamsplitter, 0, 1, 0.3, 0}, {measure_pnr, [0]}, {phaseshift, 1, 0.1}, {measure_pnr, [1]}]),
    ?assertMatch({error, {unsupported_input, #{dark_counts := midcircuit}}}, pq_api:sample(Mid, #{shots => 10, dark_counts => 0.1})).

%%% 4b: thermal loss and multi-photon source terms

%% Gaussian: <n> -> eta <n> + (1 - eta) Nth, the covariance eta Sigma + (1 - eta)(2 Nth + 1) hbar/2 on the mode
thermal_gaussian() ->
    {Eta, Nth} = {0.7, 0.5},
    C = c(1, vacuum, [{squeezing, 0, 0.4, 0}, {displacement, 0, 0.6, 0.3}, {thermal_loss, 0, Eta, Nth}]),
    C0 = C#circuit{instrs = lists:droplast(C#circuit.instrs)},
    {ok, S0, _} = pq_api:evolve(C0, #{}),
    {ok, S, _} = pq_api:evolve(C, #{}),
    {ok, N0} = pq_stats:mean_photon_number(S0, [0]),
    {ok, N1} = pq_stats:mean_photon_number(S, [0]),
    ?assert(abs(N1 - (Eta * N0 + (1 - Eta) * Nth)) < 1.0e-13),
    {ok, {rmat, [[Sxx0, _], _]}} = pq_core:unpack(S0#gaussian_state.cov),
    {ok, {rmat, [[Sxx, _], _]}} = pq_core:unpack(S#gaussian_state.cov),
    H = S#gaussian_state.hbar,
    ?assert(abs(Sxx - (Eta * Sxx0 + (1 - Eta) * (2 * Nth + 1) * H / 2)) < 1.0e-13),
    %% Nth = 0 is pure loss
    {ok, Sl, _} = pq_api:evolve(C0#circuit{instrs = C0#circuit.instrs ++ [{loss, 0, Eta}]}, #{}),
    {ok, St, _} = pq_api:evolve(C0#circuit{instrs = C0#circuit.instrs ++ [{thermal_loss, 0, Eta, 0.0}]}, #{}),
    ?assertEqual(Sl, St).

%% the density-matrix backend (amplifier after loss) against the Gaussian backend's photon statistics
thermal_dm() ->
    Is = [{squeezing, 0, 0.3, 0}, {displacement, 0, 0.5, 0.2}, {beamsplitter, 0, 1, 0.6, 0.1},
          {thermal_loss, 0, 0.8, 0.3}, {thermal_loss, 1, 0.6, 0.2}, {measure_pnr, [0, 1]}],
    C = c(2, vacuum, Is),
    gaussian = pq_api:backend_for(C, #{}),
    {ok, Pg, _} = pq_api:probabilities(C, #{max_photons => 6}),
    {ok, Pd, Meta} = pq_api:probabilities(C, #{backend => fock_dm, cutoff => 18}),
    Worst = lists:max([abs(V - maps:get(K, Pd, 0.0)) || K := V <- Pg]),
    ?debugFmt("thermal loss, density matrix (cutoff 18) against Gaussian: max difference ~.2e, warnings ~p", [Worst, maps:get(warnings, Meta)]),
    ?assert(Worst < 1.0e-9),
    %% a Fock input goes to the density-matrix backend by default
    F = c(1, {fock, [1]}, [{thermal_loss, 0, 0.6, 0.8}, {measure_pnr, [0]}]),
    fock_dm = pq_api:backend_for(F, #{}),
    {ok, Pf, _} = pq_api:probabilities(F, #{cutoff => 40}),
    Mean = lists:sum([V * N || [N] := V <- Pf]),
    ?assert(abs(Mean - (0.6 + 0.4 * 0.8)) < 1.0e-10),
    %% JSON round trip
    {ok, J} = pq_circuit:to_json(F),
    ?assertEqual({ok, F}, pq_circuit:from_json(J)).

thermal_refusals() ->
    ?assertMatch({error, {domain_violation, #{param := nth}}}, pq_circuit:validate(c(1, vacuum, [{thermal_loss, 0, 0.5, -1.0}]))),
    ?assertMatch({error, {domain_violation, #{param := eta}}}, pq_circuit:validate(c(1, vacuum, [{thermal_loss, 0, 1.5, 0.1}]))),
    F = c(1, {fock, [1]}, [{thermal_loss, 0, 0.6, 0.8}, {measure_pnr, [0]}]),
    ?assertMatch({error, _}, pq_api:probabilities(F, #{backend => fock})),
    {ok, _, _} = pq_api:probabilities(F, #{backend => fock, promote_dm => true, cutoff => 10}),
    %% not periodic: no shift rule; the adjoint is for pure states
    G = c(1, {fock, [1]}, [{thermal_loss, 0, {sym, e}, 0.1}]),
    ?assertMatch({error, {differentiate_unsupported, _}}, pq_gradient:gradient(G, {mean_photon_number, [0]}, #{e => 0.5}, #{mode => parameter_shift})),
    ?assertMatch({error, {differentiate_unsupported, _}}, pq_gradient:gradient(G, {mean_photon_number, [0]}, #{e => 0.5}, #{mode => adjoint})).

%% multi-photon terms: g2 = 2q/(1+q)^2 from a Hanbury Brown-Twiss split, and the sampler against the mixture
multi_photon() ->
    Q = 0.08,
    Hbt = c(2, {fock, [1, 0]}, [{beamsplitter, 0, 1, math:pi() / 4, 0}, {measure_pnr, [0, 1]}]),
    {ok, P, _} = pq_api:probabilities(Hbt, #{multi_photon => Q}),
    ?assert(abs(maps:get([1, 1], P) - Q / 2) < 1.0e-14),
    Mean = lists:sum([V * (A + B) || [A, B] := V <- P]),
    G2 = 4 * lists:sum([V * A * B || [A, B] := V <- P]) / (Mean * Mean),
    ?assert(abs(G2 - 2 * Q / math:pow(1 + Q, 2)) < 1.0e-13),
    {ok, #result{samples = S}} = pq_api:sample(Hbt, #{shots => 300000, seed => 2, multi_photon => Q}),
    {Pv, _} = pq_chi2:p_value(P, pq_chi2:counts(S)),
    ?assert(Pv > 1.0e-4),
    %% with partial distinguishability on a three-photon mesh
    C = c(4, {fock, [1, 1, 1, 0]}, mesh() ++ [{measure_pnr, [0, 1, 2, 3]}]),
    O = #{multi_photon => [0.05, 0.1, 0.0, 0.0], visibility => 0.85},
    {ok, Pm, _} = pq_api:probabilities(C, O),
    ?assert(abs(lists:sum(maps:values(Pm)) - 1) < 1.0e-12),
    {ok, #result{samples = Sm}} = pq_api:sample(C, O#{shots => 500000, seed => 21}),
    {Pvm, _} = pq_chi2:p_value(Pm, pq_chi2:counts(Sm)),
    ?debugFmt("multi-photon terms: HBT g2 ~.5f (closed form ~.5f), mesh chi-square p = ~.3f", [G2, 2 * Q / math:pow(1 + Q, 2), Pvm]),
    ?assert(Pvm > 1.0e-4),
    %% D1
    {ok, #result{samples = A}} = pq_api:sample(C, O#{shots => 5000, seed => 3, parallel => false}),
    {ok, #result{samples = B}} = pq_api:sample(C, O#{shots => 5000, seed => 3, max_inflight => 4, chunk => 53}),
    ?assertEqual(A, B),
    ?assertMatch({error, {unsupported_input, _}}, pq_api:sample(C, #{shots => 10, multi_photon => 0.1, backend => fock})).
