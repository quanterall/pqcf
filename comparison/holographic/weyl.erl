%% The holographic Weyl anomaly of N = 4 super-Yang-Mills in PQCF (decision record SY2, phase 5): a = c = N^2/4 from
%% the Fefferman-Graham expansion on two boundaries, S^4 and S^2 x S^2. Writes WEYL.md here. The same steps as
%% test/pq_holo_tests.erl.
%%
%%   docker run --rm -u $(id -u):$(id -g) -e PQCF_COMMIT=$(git rev-parse --short HEAD) -v "$PWD":/w -w /w pqcf-dev \
%%     sh -c 'rebar3 compile >/dev/null && erlc -o /tmp comparison/holographic/weyl.erl &&
%%     erl -noshell -pa _build/default/lib/*/ebin -pa /tmp -eval "weyl:main(\"comparison/holographic\"), halt()."'
-module(weyl).
-export([main/1]).

-define(S(X), {sym, X}).

pp(E) ->
    case catch pq_symbolic:pretty(E) of
        B when is_binary(B) -> B;
        _ -> io_lib:format("~0p", [E])
    end.

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

main(Dir) ->
    R = ?S(r), Pi = {const, pi},
    Units = #{g5 => {mul, [Pi, {pow, ?S(l), 3}, {rat, 1, 2}, {pow, ?S(n), -2}]}},
    Conf = fun(Xs) -> {mul, [4, {pow, R, 2}, {pow, {add, [1 | [{pow, ?S(X), 2} || X <- Xs]]}, -2}]} end,
    S4 = Conf([x, y, u, v]),
    {Ta, {M1, Tau1, T1, E41, W21}} = timer:tc(fun() -> fefferman_graham([S4, S4, S4, S4], [x, y, u, v]) end),
    {ok, A} = pq_symbolic:canonicalise({mul, [-16, {pow, Pi, 2}, T1, {pow, E41, -1}]}),
    {ok, Ac} = pq_symbolic:subs(A, Units),
    S2a = Conf([x, y]), S2b = Conf([u, v]),
    {Tb, {M2, Tau2, T2, E42, W22}} = timer:tc(fun() -> fefferman_graham([S2a, S2a, S2b, S2b], [x, y, u, v]) end),
    {ok, C} = pq_symbolic:canonicalise({mul, [{add, [{mul, [16, {pow, Pi, 2}, T2]}, {mul, [Ac, E42]}]}, {pow, W22, -1}]}),
    {ok, Cc} = pq_symbolic:subs(C, Units),
    Want = {mul, [{rat, 1, 4}, {pow, ?S(n), 2}]},
    Ok = {pq_symbolic:equiv(Ac, Want, #{}), pq_symbolic:equiv(Cc, Want, #{})},
    Md = ["# The holographic Weyl anomaly: a = c = N^2/4\n\n",
          "Produced by `weyl.erl` (decision record SY2, phase 5), PQCF ", os:getenv("PQCF_COMMIT", "(commit not given)"),
          ". Every value below is PQCF's output.\n\n",
          "## Method\n\n",
          "Henningson and Skenderis, [hep-th/9806087](https://arxiv.org/abs/hep-th/9806087); de Haro, Skenderis and\n",
          "Solodukhin, [hep-th/0002230](https://arxiv.org/abs/hep-th/0002230).\n\n",
          "1. **Bulk metric:** L^2 [drho^2/(4 rho^2) + g(rho)_ij dx^i dx^j / rho], with\n",
          "   g(rho) = g0 (1 + rho m + rho^2 tau) for an Einstein boundary g0.\n",
          "2. **Coefficients:** the 5-dimensional Einstein equations R_ab + (4/L^2) G_ab = 0 are solved order by order.\n",
          "   The boundary components at rho^0 fix m, so g_(2) = m g0; the radial component at rho^0 then fixes tau,\n",
          "   and with it Tr g_(4).\n",
          "3. **Volume coefficient:** a_4 is the rho^2 coefficient of sqrt(det g(rho) / det g0) = (1 + rho m + rho^2 tau)^2.\n",
          "4. **Anomaly:** the log divergence of the on-shell action under a constant rescaling of g0 gives\n",
          "   <T> = -(L^3/(2 pi G)) a_4, compared with (c W^2 - a E_4)/(16 pi^2).\n",
          "5. **Separating a from c:** S^4 has W = 0, so it gives a. S^2 x S^2 with equal radii has W^2 /= 0, so it then\n",
          "   gives c. Both are in stereographic coordinates.\n",
          "6. **Units:** L^3/G = 2 N^2/pi.\n\n",
          "## Results\n\n",
          "| boundary (radius r) | m | tau | <T> | E_4 | W^2 |\n| --- | --- | --- | --- | --- | --- |\n",
          io_lib:format("| S^4 | `~s` | `~s` | `~s` | `~s` | `~s` |~n", [pp(M1), pp(Tau1), pp(T1), pp(E41), pp(W21)]),
          io_lib:format("| S^2 x S^2 | `~s` | `~s` | `~s` | `~s` | `~s` |~n~n", [pp(M2), pp(Tau2), pp(T2), pp(E42), pp(W22)]),
          "- a = `", pp(Ac), "`\n- c = `", pp(Cc), "`\n\n",
          "Against N^2/4, PQCF's `equiv` answers `", io_lib:format("~p", [Ok]), "`.\n",
          io_lib:format("~nTimes: S^4 ~.1f s, S^2 x S^2 ~.1f s.~n", [Ta / 1.0e6, Tb / 1.0e6])],
    ok = file:write_file(filename:join(Dir, "WEYL.md"), unicode:characters_to_binary(Md)),
    io:format("~s", [Md]).
