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
