%% Holographic renormalisation (decision record SY2, phase 5): the free energy of the AdS_5 black brane from its
%% renormalised on-shell action. With ds^2 = (L^2/z^2)(f dtau^2 + dx^2 + dz^2/f), f = 1 - z^4/z_h^4, Euclidean:
%%   I = -(1/16 pi G) int sqrt(g) (R - 2 Lambda) - (1/8 pi G) int sqrt(gamma) K + (1/8 pi G) int sqrt(gamma) 3/L
%% (de Haro, Skenderis, Solodukhin, hep-th/0002230, https://arxiv.org/abs/hep-th/0002230), the boundary at z = eps,
%% eps -> 0. Per unit volume F = T I = I / beta; with T = 1/(pi z_h) and L^3/G = 2 N^2/pi, F = -pi^2 N^2 T^4 / 8.
-module(pq_holo_tests).
-include_lib("eunit/include/eunit.hrl").

-define(S(X), {sym, X}).

can(E) ->
    {ok, C} = pq_symbolic:canonicalise(E, #{domains => maps:from_list([{V, positive} || V <- [l, z, zh, eps, g5, n, t]])}),
    C.
d(E, X) -> {ok, R} = pq_symbolic:differentiate(E, X), can(R).
at(E, X, V) -> {ok, R} = pq_symbolic:subs(E, #{X => V}), can(R).
inv(E) -> {pow, E, -1}.
pi() -> {const, pi}.

free_energy_test_() ->
    {timeout, 600, fun() ->
        Z = ?S(z), L = ?S(l), Zh = ?S(zh), G5 = ?S(g5), Eps = ?S(eps),
        F = {add, [1, {mul, [-1, {pow, Z, 4}, {pow, Zh, -4}]}]},
        A = {mul, [{pow, L, 2}, {pow, Z, -2}]},
        Ds = [{mul, [A, F]}, A, A, A, {mul, [A, inv(F)]}],
        Metric = [[case I of J -> X; _ -> 0 end || J <- lists:seq(1, 5)] || {I, X} <- lists:enumerate(Ds)],
        {ok, #{scalar := Rs}} = pq_geometry:curvature(Metric, [tau, x, y, w, z]),
        Lambda = {mul, [-6, {pow, L, -2}]},
        % det g reduced to normal form first (f / f cancels), then the root
        {ok, Det} = pq_symbolic:cancel({mul, Ds}),
        SqrtG = can({pow, Det, {rat, 1, 2}}),
        ?assertEqual(can({mul, [{pow, L, 5}, {pow, Z, -5}]}), SqrtG),
        % the bulk term, z from eps to z_h
        Bulk = can({mul, [{rat, -1, 16}, inv(pi()), inv(G5), SqrtG, {add, [Rs, {mul, [-2, Lambda]}]}]}),
        {ok, Prim} = pq_symbolic:integrate(Bulk, z),
        IBulk = can({add, [at(Prim, z, Zh), {mul, [-1, at(Prim, z, Eps)]}]}),
        % the boundary z = eps: unit normal toward the boundary n^z = -z sqrt(f)/L; K = (1/sqrt g) d_z(sqrt(g) n^z)
        Nz = {mul, [-1, Z, {pow, F, {rat, 1, 2}}, inv(L)]},
        % sqrt(gamma) = sqrt(A^4) sqrt(f), both positive: the monomial root exact
        {ok, A4} = pq_symbolic:cancel({pow, A, 4}),
        SqrtGamma = can({mul, [{pow, A4, {rat, 1, 2}}, {pow, F, {rat, 1, 2}}]}),
        K = can({mul, [inv(SqrtG), d({mul, [SqrtG, Nz]}, z)]}),
        GH = {mul, [{rat, -1, 8}, inv(pi()), inv(G5), SqrtGamma, K]},
        CT = {mul, [{rat, 1, 8}, inv(pi()), inv(G5), {mul, [3, inv(L)]}, SqrtGamma]},
        Total = can({add, [IBulk, at(can({add, [GH, CT]}), z, Eps)]}),
        % the divergences cancel: the limit is finite
        {ok, Fe} = pq_symbolic:limit(Total, eps, 0),
        ?assert(pq_symbolic:equiv(Fe, {mul, [{rat, -1, 16}, inv(pi()), inv(G5), {pow, L, 3}, {pow, Zh, -4}]}, #{}) =:= {ok, true}, Fe),
        % T = 1/(pi z_h), L^3/G = 2 N^2/pi
        {ok, Ft} = pq_symbolic:subs(Fe, #{zh => inv({mul, [pi(), ?S(t)]}), g5 => {mul, [pi(), {pow, L, 3}, {rat, 1, 2}, {pow, ?S(n), -2}]}}),
        ?assertEqual({ok, true}, pq_symbolic:equiv(Ft, {mul, [{rat, -1, 8}, {pow, pi(), 2}, {pow, ?S(n), 2}, {pow, ?S(t), 4}]}, #{}))
    end}.

integrate_test() ->
    X = ?S(x),
    % x^2 + 1/x + 1/x^2 -> x^3/3 + log x - 1/x: at x = 2, 8/3 + log 2 - 1/2
    {ok, I} = pq_symbolic:integrate({add, [{pow, X, 2}, {pow, X, -1}, {pow, X, -2}]}, x),
    {ok, V} = pq_mp:evaluate(I, #{x => 2}, 30),
    {ok, W} = pq_mp:evaluate({add, [{rat, 13, 6}, {log, 2}]}, #{}, 30),
    ?assertEqual(W, V),
    ?assertMatch({error, {symbolic_unsupported, #{integrate := _}}}, pq_symbolic:integrate({exp, X}, x)).

%%% The holographic Weyl anomaly (Henningson, Skenderis, hep-th/9806087, https://arxiv.org/abs/hep-th/9806087):
%%% <T> = (c W^2 - a E_4) / (16 pi^2) with a = c = N^2/4 for N = 4 super-Yang-Mills. The bulk metric in
%%% Fefferman-Graham form, L^2 [d rho^2 / (4 rho^2) + g(rho)_ij dx^i dx^j / rho], g(rho) = g0 (1 + rho m + rho^2 tau) for an
%%% Einstein boundary g0; m and tau (which fixes Tr g_(4)) from Einstein's equations order by order; the volume
%%% coefficient a_4 (rho^2 in sqrt(det g(rho)) / sqrt(det g0)) gives <T> = -(L^3 / (2 pi G)) a_4 (the log divergence of the
%%% on-shell action under a constant rescaling of g0). S^4 (W = 0) gives a; S^2 x S^2 (equal radii, W^2 /= 0) then c.

%% the lowest-order coefficient of E's series in rho at exponent K (log power 0)
coef(E, K) ->
    {ok, {Ts, _}} = pq_series:terms(E, rho, K + 1, #{}),
    case [C || {Ex, 0, C} <- Ts, Ex =:= K] of [C] -> C; [] -> 0 end.

%% X with X linear: alpha + beta X = 0
solve_linear(E, X) ->
    {ok, B} = pq_symbolic:differentiate(E, X),
    {ok, A} = pq_symbolic:subs(E, #{X => 0}),
    {ok, R} = pq_symbolic:cancel({mul, [-1, A, {pow, B, -1}]}),
    R.

%% {m, tau, <T>, E_4, W^2} for the boundary metric G0 (diagonal, Einstein) in Coords
fefferman_graham(G0Diag, Coords) ->
    L = ?S(l), Rho = ?S(rho),
    Mult = {add, [1, {mul, [Rho, ?S(m)]}, {mul, [{pow, Rho, 2}, ?S(tau)]}]},
    Ds = [{mul, [{pow, L, 2}, {rat, 1, 4}, {pow, Rho, -2}]} |
          [{mul, [{pow, L, 2}, X, Mult, {pow, Rho, -1}]} || X <- G0Diag]],
    N = length(Ds),
    Metric = [[case I of J -> X; _ -> 0 end || J <- lists:seq(1, N)] || {I, X} <- lists:enumerate(Ds)],
    {ok, #{ricci := Ric}} = pq_geometry:curvature(Metric, [rho | Coords]),
    E = fun(I) -> {ok, V} = pq_symbolic:cancel({add, [lists:nth(I, lists:nth(I, Ric)), {mul, [4, {pow, L, -2}, lists:nth(I, Ds)]}]}), V end,
    % E_xx at rho^0 fixes m; then E_rho rho at rho^0 fixes tau
    M = solve_linear(coef(E(2), 0), m),
    {ok, Err} = pq_symbolic:subs(E(1), #{m => M}),
    Tau = solve_linear(coef(Err, 0), tau),
    % a_4: rho^2 in sqrt(det g(rho)/det g0) = (1 + rho m + rho^2 tau)^(n/2), n = 4
    {ok, A4} = pq_symbolic:cancel({add, [{pow, M, 2}, {mul, [2, Tau]}]}),
    T = {mul, [-1, {pow, L, 3}, {pow, {mul, [2, {const, pi}, ?S(g5)]}, -1}, A4]},
    % the boundary's curvature: E_4 = Riem^2 - 4 Ric^2 + R^2, W^2 = Riem^2 - 2 Ric^2 + R^2/3
    N0 = length(G0Diag),
    G0 = [[case I of J -> X; _ -> 0 end || J <- lists:seq(1, N0)] || {I, X} <- lists:enumerate(G0Diag)],
    {ok, #{ricci := R0, inverse := Gi, scalar := Sc, kretschmann := K}} = pq_geometry:curvature(G0, Coords),
    Ric2 = {add, [{mul, [lists:nth(I, lists:nth(I, Gi)), lists:nth(J, lists:nth(J, Gi)), lists:nth(J, lists:nth(I, R0)), lists:nth(J, lists:nth(I, R0))]}
                  || I <- lists:seq(1, N0), J <- lists:seq(1, N0)]},
    {ok, E4} = pq_symbolic:cancel({add, [K, {mul, [-4, Ric2]}, {pow, Sc, 2}]}),
    {ok, W2} = pq_symbolic:cancel({add, [K, {mul, [-2, Ric2]}, {mul, [{rat, 1, 3}, {pow, Sc, 2}]}]}),
    {M, Tau, T, E4, W2}.

weyl_anomaly_test_() ->
    {timeout, 1200, fun() ->
        R = ?S(r), Pi = {const, pi},
        Units = #{g5 => {mul, [Pi, {pow, ?S(l), 3}, {rat, 1, 2}, {pow, ?S(n), -2}]}},   % L^3/G = 2 N^2 / pi
        Conf = fun(Xs) -> {mul, [4, {pow, R, 2}, {pow, {add, [1 | [{pow, ?S(X), 2} || X <- Xs]]}, -2}]} end,
        % S^4 of radius r: conformally flat, so <T> = -a E_4 / (16 pi^2)
        S4 = Conf([x, y, u, v]),
        {M1, _, T1, E41, W21} = fefferman_graham([S4, S4, S4, S4], [x, y, u, v]),
        ?assert(pq_symbolic:equiv(M1, {mul, [{rat, -1, 2}, {pow, R, -2}]}, #{}) =:= {ok, true}),   % g_(2) = -(R_ij - R g/6)/2
        ?assertEqual({ok, true}, pq_symbolic:equiv(W21, 0, #{})),
        {ok, A} = pq_symbolic:canonicalise({mul, [-16, {pow, Pi, 2}, T1, {pow, E41, -1}]}),
        {ok, Ac} = pq_symbolic:subs(A, Units),
        ?assertEqual({ok, true}, pq_symbolic:equiv(Ac, {mul, [{rat, 1, 4}, {pow, ?S(n), 2}]}, #{})),
        % S^2 x S^2, both of radius r: W^2 /= 0; c = (16 pi^2 <T> + a E_4) / W^2
        S2a = Conf([x, y]), S2b = Conf([u, v]),
        {M2, _, T2, E42, W22} = fefferman_graham([S2a, S2a, S2b, S2b], [x, y, u, v]),
        ?assert(pq_symbolic:equiv(M2, {mul, [{rat, -1, 6}, {pow, R, -2}]}, #{}) =:= {ok, true}),
        ?assertNot(pq_symbolic:equiv(W22, 0, #{}) =:= {ok, true}),
        {ok, C} = pq_symbolic:canonicalise({mul, [{add, [{mul, [16, {pow, Pi, 2}, T2]}, {mul, [Ac, E42]}]}, {pow, W22, -1}]}),
        {ok, Cc} = pq_symbolic:subs(C, Units),
        ?assertEqual({ok, true}, pq_symbolic:equiv(Cc, {mul, [{rat, 1, 4}, {pow, ?S(n), 2}]}, #{}))
    end}.

%%% Gubser-Nellore (arXiv:0804.0434, https://arxiv.org/abs/0804.0434), comparison/holographic/gubser.erl

load_gubser() ->
    Dir = filename:dirname(code:which(?MODULE)),
    Src = hd([F || F <- [filename:join([Dir, "..", "comparison", "holographic", "gubser.erl"]),
                         filename:join([Dir, "..", "..", "..", "..", "comparison", "holographic", "gubser.erl"]),
                         "comparison/holographic/gubser.erl"], filelib:is_regular(F)]),
    {ok, gubser, Bin} = compile:file(Src, [binary]),
    {module, gubser} = code:load_binary(gubser, Src, Bin).

gubser_test_() ->
    {timeout, 900, fun() ->
        load_gubser(),
        % Chamblin-Reall: V = -12 e^(phi/2), G = -1/(3 gamma) = -2/3 kept exactly by the master equation
        Ex = {exp, {mul, [{rat, 1, 2}, ?S(phi)]}},
        {ok, [G1, Gp1, _, _]} = pq_jet:integrate(gubser:rhs({mul, [-12, Ex]}, {mul, [-6, Ex]}, {mul, [-3, Ex]}, none), phi,
                                                [g, gp, i1, i2], 3, [{rat, -2, 3}, 0, 0, 0], 1, #{digits => 30}),
        ?assert(abs(pq_mp:to_float(G1) + 2 / 3) < 1.0e-25),
        ?assert(abs(pq_mp:to_float(Gp1)) < 1.0e-25),
        % the paper's parameters: Delta = 3.93; c_s^2 -> 1/3 in the UV; 0.070171 at phi_H = 4 (an independent SciPy
        % solution agrees to 6 digits there)
        Dm4 = gubser:delta_minus_4({606, 1000}, {206, 100}),
        ?assert(abs(pq_mp:to_float(Dm4) + 4 - 3.93) < 0.005),
        Pot = gubser:pot({rat, 606, 1000}, {rat, 206, 100}),
        {C1, _} = gubser:cs2(Pot, Dm4, 0.2, 0.001, 20),
        {C4, _} = gubser:cs2(Pot, Dm4, 4.0, 0.001, 20),
        ?assert(abs(C1 - 1 / 3) < 2.0e-5),
        ?assert(abs(C4 - 0.070171) < 2.0e-6)
    end}.

%%% Improved holographic QCD (arXiv:0707.1349), comparison/holographic/ihqcd.erl: the model's consistency
%%% d log W / d lambda = -(4/3) X / lambda (2.22), and pq_jet:diff against a central difference

ihqcd_model_test() ->
    Dir = filename:dirname(code:which(?MODULE)),
    Src = hd([F || F <- [filename:join([Dir, "..", "comparison", "holographic", "ihqcd.erl"]),
                         filename:join([Dir, "..", "..", "..", "..", "comparison", "holographic", "ihqcd.erl"]),
                         "comparison/holographic/ihqcd.erl"], filelib:is_regular(F)]),
    {ok, ihqcd, Bin} = compile:file(Src, [binary]),
    {module, ihqcd} = code:load_binary(ihqcd, Src, Bin),
    {X, W} = ihqcd:x_w({rat, 21, 5}),
    Ev = fun(E, L) -> {ok, S} = pq_mp:evaluate(E, #{l => L}, 30), pq_mp:to_float(pq_mp:from_string(S, 120)) end,
    L = 0.3, H = 1.0e-6,
    DlogW = (math:log(Ev(W, L + H)) - math:log(Ev(W, L - H))) / (2 * H),
    ?assert(abs(DlogW + 4 / 3 * Ev(X, L) / L) < 1.0e-8),
    DX = pq_jet:diff(X, l),
    ?assert(abs(Ev(DX, L) - (Ev(X, L + H) - Ev(X, L - H)) / (2 * H)) < 1.0e-8).
