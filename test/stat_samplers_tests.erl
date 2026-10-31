%% FS 12.4: every sampler, on instances whose exact distribution is enumerable, 1,000,000 shots with fixed seeds,
%% chi-square goodness of fit p > 1e-4 (low-probability bins pooled to an expected count of at least 5). The exact
%% distributions come from pq_api:probabilities on another path (permanents, the density matrix, the hafnian), and
%% a wrong distribution must fail (the power check).
-module(stat_samplers_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

-define(SHOTS, 1000000).
-define(ALPHA, 1.0e-4).

c(N, In, Is) -> {ok, C} = pq_circuit:new(N, In), C#circuit{instrs = Is}.

mesh(N) -> [{beamsplitter, J, J + 1, 0.4 + 0.15 * J + 0.1 * L, 0.3 * L - 0.2 * J} || L <- lists:seq(0, N - 1), J <- lists:seq(L rem 2, N - 2, 2)].

boson() -> c(4, {fock, [1, 1, 1, 0]}, mesh(4) ++ [{measure_pnr, [0, 1, 2, 3]}]).

%% the samples' fit to the exact probabilities of Exact (a circuit, or a map)
fit(Exact, Run) ->
    P = case Exact of
            #{} -> Exact;
            {C, Opts} -> {ok, M, _} = pq_api:probabilities(C, Opts), M
        end,
    {ok, #result{samples = S}} = Run(),
    ?assertEqual(?SHOTS, length(S)),
    {Pv, Info} = pq_chi2:p_value(P, pq_chi2:counts(S)),
    io:format(user, "~n  p = ~.3g ~p", [Pv, Info]),
    Pv.

pass(Exact, Run) -> ?assert(fit(Exact, Run) > ?ALPHA).

stat_test_() ->
    {setup, fun() -> {ok, _} = application:ensure_all_started(pq) end, fun(_) -> ok end,
     [{Name, {timeout, 600, F}} || {Name, F} <-
      [{"chi-square survival function against scipy", fun sf/0},
       {"boson sampling, Clifford-Clifford", fun() ->
            pass({boson(), #{backend => boson}}, fun() -> pq_api:sample(boson, boson(), #{shots => ?SHOTS, seed => 101}) end) end},
       {"boson sampling, brute force", fun() ->
            pass({boson(), #{backend => boson}},
                 fun() -> pq_api:sample(boson, boson(), #{shots => ?SHOTS, seed => 102, algorithm => brute_force}) end) end},
       {"boson sampling, Metropolis (independence sampler, thinned)", fun() ->
            pass({boson(), #{backend => boson}},
                 fun() -> pq_api:sample(boson, boson(), #{shots => ?SHOTS, seed => 112, algorithm => metropolis}) end) end},
       {"lossy boson sampling (loss at the inputs) against the density matrix", fun lossy/0},
       {"Gaussian boson sampling, photon counting", fun() ->
            pass({gbs(measure_pnr), #{}}, fun() -> pq_api:sample(gbs(measure_pnr), #{shots => ?SHOTS, seed => 103}) end) end},
       {"Gaussian boson sampling, threshold detectors", fun() ->
            pass({gbs(measure_threshold), #{}}, fun() -> pq_api:sample(gbs(measure_threshold), #{shots => ?SHOTS, seed => 104}) end) end},
       {"Fock backend, photon counting", fun() ->
            pass({fock(), #{}}, fun() -> pq_api:sample(fock(), #{shots => ?SHOTS, seed => 105}) end) end},
       {"density matrix, photon counting after loss", fun() ->
            pass({dm(), #{}}, fun() -> pq_api:sample(dm(), #{shots => ?SHOTS, seed => 106}) end) end},
       {"trajectory loss against the density matrix", fun() ->
            pass({dm(), #{backend => fock_dm}},
                 fun() -> pq_api:sample(fock, dm(), #{shots => ?SHOTS, seed => 107, noise => trajectory}) end) end},
       {"mid-circuit photon counting against the terminal joint distribution", fun midcircuit/0},
       {"homodyne against the normal distribution of the quadrature", fun homodyne/0},
       {"power: the samples of one instance fail against another's distribution", fun power/0}]]}.

sf() ->
    Ref = [{3.841458820694124, 1, 0.04999999999999994}, {18.307038053275146, 10, 0.050000000000000024},
           {0.5, 3, 0.9188914116546758}, {100.0, 50, 3.4549313829848465e-05}, {250.0, 200, 0.0093791316688261},
           {40.0, 10, 1.694474393006737e-05}, {1.0e-3, 2, 0.9995001249791693}, {1500.0, 1400, 0.031503720989392366}],
    [?assert(abs(pq_chi2:sf(X, K) - P) =< 1.0e-10 * P) || {X, K, P} <- Ref].

lossy() ->
    Etas = [0.6, 0.9, 0.75, 1.0],
    Exact = c(4, {fock, [1, 1, 1, 0]}, [{loss, M, E} || {M, E} <- lists:enumerate(0, Etas)] ++ mesh(4) ++ [{measure_pnr, [0, 1, 2, 3]}]),
    pass({Exact, #{}}, fun() -> pq_api:sample(boson, boson(), #{shots => ?SHOTS, seed => 108, loss => {input, Etas}}) end).

%% photon counting with a displacement (the loop hafnian); threshold detectors need a zero mean (the torontonian)
gbs(measure_pnr) ->
    c(3, vacuum, [{squeezing, 0, 0.5, 0}, {squeezing, 1, 0.3, 0.4}, {displacement, 2, 0.2, 0.1},
                  {beamsplitter, 0, 1, 0.7, 0.2}, {beamsplitter, 1, 2, 0.5, -0.3}, {measure_pnr, [0, 1, 2]}]);
gbs(measure_threshold) ->
    c(3, vacuum, [{squeezing, 0, 0.5, 0}, {squeezing, 1, 0.3, 0.4}, {squeezing, 2, 0.6, -0.2},
                  {beamsplitter, 0, 1, 0.7, 0.2}, {beamsplitter, 1, 2, 0.5, -0.3}, {measure_threshold, [0, 1, 2]}]).

fock() -> c(3, {fock, [2, 1, 0]}, [{beamsplitter, 0, 1, 0.6, 0.1}, {kerr, 1, 0.7}, {beamsplitter, 1, 2, 0.9, 0.4},
                                    {phaseshift, 2, 0.3}, {beamsplitter, 0, 2, 0.5, 0.0}, {measure_pnr, [0, 1, 2]}]).

dm() -> c(2, {fock, [2, 1]}, [{beamsplitter, 0, 1, 0.6, 0.1}, {loss, 0, 0.7}, {kerr, 1, 0.4}, {loss, 1, 0.5},
                              {beamsplitter, 0, 1, 0.4, 0.2}, {measure_pnr, [0, 1]}]).

midcircuit() ->
    Mid = c(3, {fock, [1, 1, 1]}, [{beamsplitter, 0, 1, 0.7, 0.2}, {measure_pnr, [0]}, {beamsplitter, 1, 2, 0.5, 0.1},
                                   {measure_pnr, [1, 2]}]),
    %% the measured mode is not used again, so the joint distribution is that of measuring all three at the end
    End = c(3, {fock, [1, 1, 1]}, [{beamsplitter, 0, 1, 0.7, 0.2}, {beamsplitter, 1, 2, 0.5, 0.1}, {measure_pnr, [0, 1, 2]}]),
    pass({End, #{}}, fun() -> pq_api:sample(Mid, #{shots => ?SHOTS, seed => 109}) end).

homodyne() ->
    C = c(1, vacuum, [{squeezing, 0, 0.4, 0.3}, {displacement, 0, 0.5, 0.2}, {measure_homodyne, 0, 0.6}]),
    %% x cos phi + p sin phi of the state before the measurement
    {ok, S, _} = pq_api:evolve(C#circuit{instrs = lists:droplast(C#circuit.instrs)}, #{}),
    {[Mx, Mp], [[Sxx, Sxp], [_, Spp]]} = pq_state:mean_cov(S),
    {Co, Si} = {math:cos(0.6), math:sin(0.6)},
    Mu = Mx * Co + Mp * Si,
    Sd = math:sqrt(Sxx * Co * Co + 2 * Sxp * Co * Si + Spp * Si * Si),
    Edges = [Mu + Sd * (K / 10) || K <- lists:seq(-40, 40)],
    Cdf = fun(X) -> 0.5 * (1 + math:erf((X - Mu) / (Sd * math:sqrt(2)))) end,
    Bin = fun(X) -> length([E || E <- Edges, E =< X]) end,
    Probs = maps:from_list([{K, case K of
                                    0 -> Cdf(hd(Edges));
                                    81 -> 1 - Cdf(lists:last(Edges));
                                    _ -> Cdf(lists:nth(K + 1, Edges)) - Cdf(lists:nth(K, Edges))
                                end} || K <- lists:seq(0, 81)]),
    {ok, #result{kind = homodyne, samples = Rows}} = pq_api:sample(C, #{shots => ?SHOTS, seed => 110}),
    {Pv, Info} = pq_chi2:p_value(Probs, pq_chi2:counts([Bin(X) || [X] <- Rows])),
    io:format(user, "~n  p = ~.3g ~p", [Pv, Info]),
    ?assert(Pv > ?ALPHA).

%% the test can fail: Clifford-Clifford samples against the distribution of distinguishable photons on the same
%% circuit (the permanents of |U|^2), and a 1 % shift of one beamsplitter angle
power() ->
    C = boson(),
    U = pq_sampling:passive_unitary(element(2, pq_circuit:lower(C#circuit{instrs = lists:droplast(C#circuit.instrs)}, #{}))),
    {ok, {cmat, Rows}} = pq_core:unpack(U),
    Abs2 = [[{complex, Re * Re + Im * Im, 0.0} || {complex, Re, Im} <- R] || R <- Rows],
    {ok, A} = pq_core:pack({cmat, Abs2}),
    {ok, Exact, _} = pq_api:probabilities(C, #{backend => boson}),
    Perm = fun(P) -> {ok, {Z, _}} = pq_nif:permanent_repeated(A, P, [1, 1, 1, 0], #{}), re(Z) end,
    Dist = maps:from_list([{P, Perm(P) / fact(P)} || P <- maps:keys(Exact)]),
    ?assert(abs(lists:sum(maps:values(Dist)) - 1) < 1.0e-12),
    Run = fun() -> pq_api:sample(boson, C, #{shots => ?SHOTS, seed => 111}) end,
    ?assert(fit(Dist, Run) < 1.0e-10),
    Shifted = C#circuit{instrs = [setelement(4, hd(C#circuit.instrs), 1.01 * element(4, hd(C#circuit.instrs))) | tl(C#circuit.instrs)]},
    ?assert(fit({Shifted, #{backend => boson}}, Run) < ?ALPHA).

re({R, _Im}) -> R.

fact(P) -> lists:foldl(fun(N, A) -> A * lists:foldl(fun erlang:'*'/2, 1, lists:seq(1, N)) end, 1, P).
