%% pq_geometry (decision record SY2, phase 4): curvature of the backgrounds of holographic QCD and of Schwarzschild.
-module(pq_geometry_tests).
-include_lib("eunit/include/eunit.hrl").

-define(L, {sym, l}).

eq(A, B) -> pq_symbolic:equiv(A, B, #{}) =:= {ok, true}.

diag(Ds) ->
    N = length(Ds),
    [[case I of J -> lists:nth(I, Ds); _ -> 0 end || J <- lists:seq(1, N)] || I <- lists:seq(1, N)].

%% AdS_5 in Poincare coordinates: ds^2 = (L^2/z^2)(-dt^2 + dx^2 + dy^2 + dw^2 + dz^2): R = -20/L^2,
%% R_ab = -(4/L^2) g_ab, Weyl = 0, Kretschmann = 40/L^4
ads5_test_() ->
    {timeout, 300, fun() ->
        Z = {sym, z},
        F = {mul, [{pow, ?L, 2}, {pow, Z, -2}]},
        G = diag([{mul, [-1, F]}, F, F, F, F]),
        {ok, #{scalar := R, ricci := Ric, weyl_zero := W, kretschmann := K}} = pq_geometry:curvature(G, [t, x, y, w, z]),
        ?assert(eq(R, {mul, [-20, {pow, ?L, -2}]})),
        [?assert(eq(lists:nth(J, lists:nth(I, Ric)), {mul, [-4, {pow, ?L, -2}, lists:nth(J, lists:nth(I, G))]}))
         || I <- lists:seq(1, 5), J <- lists:seq(1, 5)],
        ?assert(W),
        ?assert(eq(K, {mul, [40, {pow, ?L, -4}]}))
    end}.

%% the AdS_5 black brane, f = 1 - z^4/z_h^4: G_ab + Lambda g_ab = 0 with Lambda = -6/L^2; R = -20/L^2 still; the
%% Hawking temperature |f'(z_h)| / (4 pi) = 1 / (pi z_h)
black_brane_test_() ->
    {timeout, 300, fun() ->
        Z = {sym, z}, Zh = {sym, zh},
        Fz = {add, [1, {mul, [-1, {pow, Z, 4}, {pow, Zh, -4}]}]},
        A = {mul, [{pow, ?L, 2}, {pow, Z, -2}]},
        G = diag([{mul, [-1, A, Fz]}, A, A, A, {mul, [A, {pow, Fz, -1}]}]),
        {ok, #{einstein := Ein, scalar := R, weyl_zero := W}} = pq_geometry:curvature(G, [t, x, y, w, z]),
        Lambda = {mul, [-6, {pow, ?L, -2}]},
        [?assert(eq({add, [lists:nth(J, lists:nth(I, Ein)), {mul, [Lambda, lists:nth(J, lists:nth(I, G))]}]}, 0))
         || I <- lists:seq(1, 5), J <- lists:seq(1, 5)],
        ?assert(eq(R, {mul, [-20, {pow, ?L, -2}]})),
        ?assertNot(W),
        {ok, Df} = pq_symbolic:differentiate(Fz, z),
        {ok, DfH} = pq_symbolic:subs(Df, #{z => Zh}),
        T = {mul, [-1, DfH, {rat, 1, 4}, {pow, {const, pi}, -1}]},
        ?assert(eq(T, {mul, [{pow, {const, pi}, -1}, {pow, Zh, -1}]}))
    end}.

%% Schwarzschild with x = cos theta: Ricci flat, Kretschmann = 48 M^2 / r^6
schwarzschild_test_() ->
    {timeout, 300, fun() ->
        R = {sym, r}, X = {sym, x}, M = {sym, m},
        F = {add, [1, {mul, [-2, M, {pow, R, -1}]}]},
        G = diag([{mul, [-1, F]}, {pow, F, -1},
                  {mul, [{pow, R, 2}, {pow, {add, [1, {mul, [-1, {pow, X, 2}]}]}, -1}]},
                  {mul, [{pow, R, 2}, {add, [1, {mul, [-1, {pow, X, 2}]}]}]}]),
        {ok, #{ricci := Ric, kretschmann := K}} = pq_geometry:curvature(G, [t, r, x, phi]),
        ?assert(lists:all(fun(E) -> eq(E, 0) end, lists:append(Ric))),
        ?assert(eq(K, {mul, [48, {pow, M, 2}, {pow, R, -6}]}))
    end}.
