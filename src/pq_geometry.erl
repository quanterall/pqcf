%% Curvature of a metric given in components (decision record SY2, phase 4): the inverse metric, Christoffel
%% symbols, Riemann, Ricci, the scalar curvature, Einstein, Weyl and the Kretschmann scalar, in any dimension.
%%
%% The components are rational functions of the coordinates and the parameters, held in normal form (pq_sym_ratfun)
%% so that every cancellation is exact and every zero is recognised; a metric that is not rational in its
%% coordinates (sin theta, say) is written in rational ones (x = cos theta). Results come back as canonical
%% expressions.
%%
%% Conventions: Gamma^a_bc = (1/2) g^ad (d_b g_dc + d_c g_db - d_d g_bc); R^a_bcd = d_c Gamma^a_db - d_d Gamma^a_cb
%% + Gamma^a_ce Gamma^e_db - Gamma^a_de Gamma^e_cb; R_bd = R^a_bad; R = g^bd R_bd (Misner, Thorne, Wheeler).
-module(pq_geometry).
-export([curvature/2]).

%% Metric: a list of rows of expressions; Coords: the coordinate symbols, in the metric's order. A map of canonical
%% expressions:
%%   inverse => rows; christoffel => #{{A, B, C} => Gamma^A_BC} (nonzero ones);
%%   riemann => #{{A, B, C, D} => R^A_BCD} (nonzero); ricci => rows; scalar => R; einstein => rows;
%%   weyl_zero => boolean() (all components of C^a_bcd vanish); kretschmann => R_abcd R^abcd
%% Indices count from 1.
-spec curvature([[term()]], [atom()]) -> {ok, map()} | {error, term()}.
curvature(Metric, Coords) ->
    try
        N = length(Coords),
        length(Metric) =:= N andalso lists:all(fun(R) -> length(R) =:= N end, Metric)
            orelse throw({geometry, bad_shape}),
        Vars = lists:usort(Coords ++ lists:append([pq_symbolic:free_symbols(E) || R <- Metric, E <- R])),
        put({?MODULE, vars}, Vars),
        G = maps:from_list([{{I, J}, rf(E, Vars)} || {I, Row} <- lists:enumerate(Metric), {J, E} <- lists:enumerate(Row)]),
        Gi = inverse(G, N),
        Ix = lists:seq(1, N),
        Pos = maps:from_list([{C, index(C, Vars)} || C <- Coords]),
        Dx = fun(F, K) -> deriv(F, maps:get(lists:nth(K, Coords), Pos)) end,
        % d_c g_ab, once
        Dg = maps:from_list([{{A, B, C}, Dx(maps:get({A, B}, G), C)} || A <- Ix, B <- Ix, C <- Ix]),
        Gam = maps:from_list([{{A, B, C}, sum([mul(maps:get({A, D}, Gi),
                                                   add(add(maps:get({D, C, B}, Dg), maps:get({D, B, C}, Dg)), neg(maps:get({B, C, D}, Dg))))
                                               || D <- Ix], half())}
                              || A <- Ix, B <- Ix, C <- Ix]),
        Dgam = maps:from_list([{{A, B, C, D}, Dx(maps:get({A, B, C}, Gam), D)} || A <- Ix, B <- Ix, C <- Ix, D <- Ix]),
        Riem = maps:from_list([{{A, B, C, D},
                                add(add(maps:get({A, D, B, C}, Dgam), neg(maps:get({A, C, B, D}, Dgam))),
                                    sum([add(mul(maps:get({A, C, E}, Gam), maps:get({E, D, B}, Gam)),
                                             neg(mul(maps:get({A, D, E}, Gam), maps:get({E, C, B}, Gam)))) || E <- Ix], one()))}
                               || A <- Ix, B <- Ix, C <- Ix, D <- Ix, C < D]),
        R = fun(A, B, C, D) when C < D -> maps:get({A, B, C, D}, Riem);
               (A, B, C, D) when C > D -> neg(maps:get({A, B, D, C}, Riem));
               (_, _, _, _) -> zero() end,
        Ric = maps:from_list([{{B, D}, sum([R(A, B, A, D) || A <- Ix], one())} || B <- Ix, D <- Ix]),
        Scal = sum([mul(maps:get({B, D}, Gi), maps:get({B, D}, Ric)) || B <- Ix, D <- Ix], one()),
        Ein = maps:from_list([{{B, D}, add(maps:get({B, D}, Ric), neg(mul(half(), mul(Scal, maps:get({B, D}, G)))))} || B <- Ix, D <- Ix]),
        % R_abcd = g_ae R^e_bcd, and R^abcd = g^bf g^cg g^dh R^a_fgh
        Low = fun(A, B, C, D) -> sum([mul(maps:get({A, E}, G), R(E, B, C, D)) || E <- Ix], one()) end,
        Lows = maps:from_list([{{A, B, C, D}, Low(A, B, C, D)} || A <- Ix, B <- Ix, C <- Ix, D <- Ix]),
        Up = fun(A, B, C, D) ->
                     sum([mul(mul(maps:get({B, F}, Gi), mul(maps:get({C, Gg}, Gi), maps:get({D, H}, Gi))), R(A, F, Gg, H))
                          || F <- Ix, Gg <- Ix, H <- Ix], one())
             end,
        Kre = sum([mul(maps:get({A, B, C, D}, Lows), Up(A, B, C, D)) || A <- Ix, B <- Ix, C <- Ix, D <- Ix], one()),
        WeylZero = weyl_zero(Lows, G, Ric, Scal, N),
        Out = fun(F) -> ast(F, Vars) end,
        Rows = fun(M) -> [[Out(maps:get({I, J}, M)) || J <- Ix] || I <- Ix] end,
        {ok, #{inverse => Rows(Gi),
               christoffel => maps:from_list([{K, Out(V)} || K := V <- Gam, not is_zero(V)]),
               riemann => maps:from_list([{K, Out(V)} || K := V <- Riem, not is_zero(V)]),
               ricci => Rows(Ric), scalar => Out(Scal), einstein => Rows(Ein),
               weyl_zero => WeylZero, kretschmann => Out(Kre)}}
    catch
        throw:{geometry, Why} -> {error, {symbolic_unsupported, #{geometry => Why}}};
        throw:{symbolic_unsupported, rational} -> {error, {symbolic_unsupported, #{geometry => not_rational}}}
    after
        erase({?MODULE, vars})
    end.

%% C_abcd = R_abcd - (g_ac R_bd - g_ad R_bc - g_bc R_ad + g_bd R_ac)/(n-2) + R (g_ac g_bd - g_ad g_bc)/((n-1)(n-2))
weyl_zero(_, _, _, _, N) when N < 4 -> true;
weyl_zero(Lows, G, Ric, Scal, N) ->
    Ix = lists:seq(1, N),
    K1 = q({1, N - 2}), K2 = q({1, (N - 1) * (N - 2)}),
    lists:all(fun({A, B, C, D}) ->
                      T1 = add(add(mul(g(G, {A, C}), g(Ric, {B, D})), neg(mul(g(G, {A, D}), g(Ric, {B, C})))),
                               add(neg(mul(g(G, {B, C}), g(Ric, {A, D}))), mul(g(G, {B, D}), g(Ric, {A, C})))),
                      T2 = add(mul(g(G, {A, C}), g(G, {B, D})), neg(mul(g(G, {A, D}), g(G, {B, C})))),
                      Cw = add(add(maps:get({A, B, C, D}, Lows), neg(mul(K1, T1))), mul(K2, mul(Scal, T2))),
                      is_zero(Cw)
              end, [{A, B, C, D} || A <- Ix, B <- Ix, C <- Ix, D <- Ix, A < B, C < D]).

g(M, K) -> maps:get(K, M).

%% the inverse by Gauss-Jordan elimination on exact rational functions
inverse(G, N) ->
    Ix = lists:seq(1, N),
    Aug = [[maps:get({I, J}, G) || J <- Ix] ++ [case I of J -> one(); _ -> zero() end || J <- Ix] || I <- Ix],
    Rows = gj(Aug, 1, N),
    maps:from_list([{{I, J}, lists:nth(N + J, lists:nth(I, Rows))} || I <- Ix, J <- Ix]).

gj(Rows, K, N) when K > N -> Rows;
gj(Rows, K, N) ->
    {Before, After} = lists:split(K - 1, Rows),
    case lists:splitwith(fun(R) -> is_zero(lists:nth(K, R)) end, After) of
        {_, []} -> throw({geometry, singular_metric});
        {Zs, [Piv | Rest]} ->
            Inv = inv(lists:nth(K, Piv)),
            PivN = [mul(Inv, E) || E <- Piv],
            Others = Before ++ Zs ++ Rest,
            Red = [begin F = lists:nth(K, R), [add(E, neg(mul(F, P))) || {E, P} <- lists:zip(R, PivN)] end || R <- Others],
            {B2, A2} = lists:split(K - 1, Red),
            gj(B2 ++ [PivN] ++ A2, K + 1, N)
    end.

%%% Rational functions (pq_sym_ratfun normal forms) over the fixed variable list

nv() -> length(get({?MODULE, vars})).
zero() -> {#{}, #{lists:duplicate(nv(), 0) => 1}}.
one() -> q({1, 1}).
half() -> q({1, 2}).
q({A, B}) ->
    Z = lists:duplicate(nv(), 0),
    {#{Z => case B of 1 -> A; _ -> {A, B} end}, #{Z => 1}}.   % A /= 0

rf(E, Vars) ->
    {ok, C} = pq_symbolic:canonicalise(E),
    pq_sym_ratfun:from_ast(C, Vars).

add(A, B) -> pq_sym_ratfun:add(A, B).
mul(A, B) -> pq_sym_ratfun:mul(A, B).
neg(A) -> mul(q({-1, 1}), A).
inv(A) -> pq_sym_ratfun:inv(A).
is_zero({N, _}) -> map_size(N) =:= 0.

sum(L, Factor) ->
    S = lists:foldl(fun(X, Acc) -> add(Acc, X) end, zero(), L),
    case Factor =:= one() of true -> S; false -> mul(Factor, S) end.

%% d(N/D)/dx_I = (N' D - N D') / D^2
deriv({N, D}, I) ->
    case map_size(N) of
        0 -> zero();
        _ ->
            Dn = pq_sym_gcd:derivative(N, I),
            Dd = pq_sym_gcd:derivative(D, I),
            Num = padd(pq_sym_gcd:mul(Dn, D), pneg(pq_sym_gcd:mul(N, Dd))),
            case map_size(Num) of
                0 -> zero();
                _ -> mul({Num, #{lists:duplicate(nv(), 0) => 1}}, inv({pq_sym_gcd:mul(D, D), #{lists:duplicate(nv(), 0) => 1}}))
            end
    end.

padd(A, B) ->
    maps:filter(fun(_, C) -> C =/= 0 end,
                maps:fold(fun(K, V, Acc) -> Acc#{K => cadd(maps:get(K, Acc, 0), V)} end, A, B)).
pneg(A) -> maps:map(fun(_, C) -> cneg(C) end, A).
cneg({A, B}) -> {-A, B};
cneg(A) -> -A.
cadd(A, B) when is_integer(A), is_integer(B) -> A + B;
cadd(A, B) ->
    {X, Y} = qq(A), {Z, W} = qq(B),
    N = X * W + Z * Y, M = Y * W,
    G = pq_sym_exact:gcd(abs(N), M),
    case M div G of 1 -> N div G; M2 -> {N div G, M2} end.
qq(A) when is_integer(A) -> {A, 1};
qq(Q) -> Q.

index(X, L) -> length(lists:takewhile(fun(Y) -> Y =/= X end, L)) + 1.

ast(F, Vars) ->
    {ok, C} = pq_symbolic:canonicalise(pq_sym_ratfun:to_ast(F, Vars)),
    C.
