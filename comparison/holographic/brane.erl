%% The free energy of the AdS_5 black brane by holographic renormalisation, in PQCF (decision record SY2, phase 5).
%% Writes BRANE.md here. The same steps as test/pq_holo_tests.erl, with the intermediate expressions shown.
%%
%%   docker run --rm -u $(id -u):$(id -g) -e PQCF_COMMIT=$(git rev-parse --short HEAD) -v "$PWD":/w -w /w pqcf-dev \
%%     sh -c 'rebar3 compile >/dev/null && erlc -o /tmp comparison/holographic/brane.erl &&
%%     erl -noshell -pa _build/default/lib/*/ebin -pa /tmp -eval "brane:main(\"comparison/holographic\"), halt()."'
-module(brane).
-export([main/1]).

-define(S(X), {sym, X}).

can(E) ->
    {ok, C} = pq_symbolic:canonicalise(E, #{domains => maps:from_list([{V, positive} || V <- [l, z, zh, eps, g5, n, t]])}),
    C.
d(E, X) -> {ok, R} = pq_symbolic:differentiate(E, X), can(R).
at(E, X, V) -> {ok, R} = pq_symbolic:subs(E, #{X => V}), can(R).
inv(E) -> {pow, E, -1}.
pi() -> {const, pi}.
pp(E) ->
    case catch pq_symbolic:pretty(E) of
        B when is_binary(B) -> B;
        {ok, B} when is_binary(B) -> B;
        _ -> io_lib:format("~0p", [E])
    end.

main(Dir) ->
    Z = ?S(z), L = ?S(l), Zh = ?S(zh), G5 = ?S(g5), Eps = ?S(eps),
    F = {add, [1, {mul, [-1, {pow, Z, 4}, {pow, Zh, -4}]}]},
    A = {mul, [{pow, L, 2}, {pow, Z, -2}]},
    Ds = [{mul, [A, F]}, A, A, A, {mul, [A, inv(F)]}],
    Metric = [[case I of J -> X; _ -> 0 end || J <- lists:seq(1, 5)] || {I, X} <- lists:enumerate(Ds)],
    {T0, {ok, #{scalar := Rs, kretschmann := Kr}}} = timer:tc(fun() -> pq_geometry:curvature(Metric, [tau, x, y, w, z]) end),
    Lambda = {mul, [-6, {pow, L, -2}]},
    {ok, Det} = pq_symbolic:cancel({mul, Ds}),
    SqrtG = can({pow, Det, {rat, 1, 2}}),
    Bulk = can({mul, [{rat, -1, 16}, inv(pi()), inv(G5), SqrtG, {add, [Rs, {mul, [-2, Lambda]}]}]}),
    {ok, Prim} = pq_symbolic:integrate(Bulk, z),
    IBulk = can({add, [at(Prim, z, Zh), {mul, [-1, at(Prim, z, Eps)]}]}),
    Nz = {mul, [-1, Z, {pow, F, {rat, 1, 2}}, inv(L)]},
    {ok, A4} = pq_symbolic:cancel({pow, A, 4}),
    SqrtGamma = can({mul, [{pow, A4, {rat, 1, 2}}, {pow, F, {rat, 1, 2}}]}),
    K = can({mul, [inv(SqrtG), d({mul, [SqrtG, Nz]}, z)]}),
    GH = at(can({mul, [{rat, -1, 8}, inv(pi()), inv(G5), SqrtGamma, K]}), z, Eps),
    CT = at(can({mul, [{rat, 3, 8}, inv(pi()), inv(G5), inv(L), SqrtGamma]}), z, Eps),
    Total = can({add, [IBulk, GH, CT]}),
    {ok, SB} = pq_symbolic:series(IBulk, eps, 1),
    {ok, SG} = pq_symbolic:series(GH, eps, 1),
    {ok, SC} = pq_symbolic:series(CT, eps, 1),
    {ok, Fe} = pq_symbolic:limit(Total, eps, 0),
    {ok, Ft} = pq_symbolic:subs(Fe, #{zh => inv({mul, [pi(), ?S(t)]}), g5 => {mul, [pi(), {pow, L, 3}, {rat, 1, 2}, {pow, ?S(n), -2}]}}),
    {ok, Fc} = pq_symbolic:canonicalise(Ft),
    Ok = pq_symbolic:equiv(Fc, {mul, [{rat, -1, 8}, {pow, pi(), 2}, {pow, ?S(n), 2}, {pow, ?S(t), 4}]}, #{}),
    Md = ["# The AdS_5 black brane's free energy by holographic renormalisation\n\n",
          "Produced by `brane.erl` (decision record SY2, phase 5), PQCF ", os:getenv("PQCF_COMMIT", "(commit not given)"),
          ". Each step below is computed by PQCF; the expressions are its output.\n\n",
          "## Set-up\n\n",
          "- **Metric** (Euclidean): ds^2 = (L^2/z^2)(f dtau^2 + dx^2 + dy^2 + dw^2 + dz^2/f), with f = 1 - z^4/z_h^4.\n",
          "- **Action:** I = -(1/16 pi G) int sqrt(g) (R - 2 Lambda) - (1/8 pi G) int sqrt(gamma) K\n",
          "  + (1/8 pi G) int sqrt(gamma) 3/L, with the boundary at z = eps (de Haro, Skenderis, Solodukhin,\n",
          "  [hep-th/0002230](https://arxiv.org/abs/hep-th/0002230)).\n",
          "- **Free energy:** per unit volume, F = T I = I / beta.\n\n",
          "## Steps\n\n",
          "1. **Curvature** (`pq_geometry`, ", io_lib:format("~.2f s", [T0 / 1.0e6]), "): R = `", pp(Rs),
          "`. The Kretschmann scalar is `", pp(Kr), "`.\n",
          "2. **Volume element and bulk term:** sqrt(g) = `", pp(SqrtG), "`. The bulk term integrated from eps to\n",
          "   z_h is `", pp(IBulk), "`.\n",
          "3. **Extrinsic curvature:** K = `", pp(K), "`.\n",
          "4. **Behaviour as eps -> 0, term by term:**\n",
          "   - bulk: `", pp(SB), "`;\n",
          "   - Gibbons-Hawking: `", pp(SG), "`;\n",
          "   - counterterm: `", pp(SC), "`.\n\n",
          "   The eps^-4 divergences cancel.\n",
          "5. **Finite limit:** F = `", pp(Fe), "`.\n",
          "6. **In field-theory units:** with T = 1/(pi z_h) and L^3/G = 2 N^2/pi, F = `", pp(Fc), "`.\n\n",
          "## Result\n\n",
          "The expected value is -pi^2 N^2 T^4 / 8, and PQCF's `equiv` answers: `", io_lib:format("~p", [Ok]), "`.\n"],
    ok = file:write_file(filename:join(Dir, "BRANE.md"), unicode:characters_to_binary(Md)),
    io:format("~s", [Md]).
