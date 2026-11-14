%% The soft-wall and hard-wall models of holographic QCD from the action to the vector-meson spectrum, by PQCF alone
%% (decision record SY2, phase 3), checked against the closed forms of the papers. Writes RESULTS.md here.
%%
%%   docker run --rm -u $(id -u):$(id -g) -e PQCF_COMMIT=$(git rev-parse --short HEAD) -v "$PWD":/w -w /w pqcf-dev sh -c 'rebar3 compile >/dev/null &&
%%     erlc -o /tmp comparison/holographic/walls.erl && erl -noshell -pa _build/default/lib/*/ebin -pa /tmp \
%%     -eval "walls:main(\"comparison/holographic\"), halt()."'
-module(walls).
-export([main/1]).

-define(D, 30).

main(Dir) ->
    Z = {sym, z},
    L = fun(Dil) -> {mul, [{rat, 1, 2}, {exp, {mul, [-1, Dil]}}, {pow, Z, -1},
                           {add, [{pow, {sym, dpsi}, 2}, {mul, [-1, {sym, m2}, {pow, {sym, psi}, 2}]}]}]} end,
    {T1, {ok, #{potential := Vs}}} = timer:tc(fun() -> pq_ode:schrodinger(L({pow, Z, 2}), z, psi, dpsi, m2) end),
    {T2, {ok, Es}} = timer:tc(fun() -> pq_ode:eigenvalues(Vs, z, #{right => {dirichlet, 10}, range => {1, 17}, grid => 64, digits => ?D}) end),
    {T3, {ok, #{potential := Vh}}} = timer:tc(fun() -> pq_ode:schrodinger(L(0), z, psi, dpsi, m2) end),
    {T4, {ok, Eh}} = timer:tc(fun() -> pq_ode:eigenvalues(Vh, z, #{right => {robin, 1, {1, 2}, 1}, range => {1, 80}, grid => 80, digits => ?D}) end),
    Exact = [pq_mp:to_string(pq_mp:from(4 * (N + 1), 120), ?D) || N <- lists:seq(0, length(Es) - 1)],
    Zeros = [pq_mp:to_string(pq_mp:mul(J, J, 130), ?D) || K <- lists:seq(1, length(Eh)), J <- [pq_mp:besselj_zero(0, K, 130)]],
    Row = fun(N, A, B) -> io_lib:format("| ~b | `~s` | `~s` | ~s |~n", [N, A, B, case A =:= B of true -> "equal"; false -> "differ" end]) end,
    Md = ["# Holographic QCD: the soft and hard walls\n\n",
          "Produced by `walls.erl` (decision record SY2, phase 3), PQCF ", pqcf_commit(), ", ", today(), ". ",
          "Each model goes from its action to its spectrum in PQCF alone:\n\n",
          "1. **Equation of motion:** the Euler-Lagrange equation of the radial Lagrangian of a vector mode,\n",
          "   (1/2) e^(-Phi) (1/z) (psi'^2 - m^2 psi^2).\n",
          "2. **Schrödinger form:** the Liouville transformation gives -phi'' + V(z) phi = m^2 phi, symbolically\n",
          "   (`pq_ode:schrodinger/5`).\n",
          "3. **Spectrum:** the eigenvalues to ", integer_to_list(?D), " digits (`pq_ode:eigenvalues/3`): a Frobenius\n",
          "   series at z = 0, Taylor-series integration in arbitrary precision, then the Illinois method.\n\n",
          "## Soft wall\n\n",
          "Karch, Katz, Son and Stephanov, [hep-ph/0602229](https://arxiv.org/abs/hep-ph/0602229). The dilaton is\n",
          "Phi = z^2 (c = 1), and the decaying solution is imposed at z = 10.\n\n",
          "- **Potential derived:** `", io_lib:format("~p", [Vs]), "`, which is z^2 + 3/(4 z^2).\n",
          "- **Prediction:** m_n^2 = 4 (n + 1).\n\n",
          "| n | m_n^2, PQCF | 4 (n + 1) | |\n| --- | --- | --- | --- |\n",
          [Row(N, A, B) || {N, {A, B}} <- lists:enumerate(0, lists:zip(Es, Exact))],
          io_lib:format("~nTimes: potential ~.2f s, spectrum ~.2f s.~n~n", [T1 / 1.0e6, T2 / 1.0e6]),
          "## Hard wall\n\n",
          "Erlich, Katz, Son and Stephanov, [hep-ph/0501128](https://arxiv.org/abs/hep-ph/0501128). There is no\n",
          "dilaton, 0 < z =< z_m with z_m = 1, and psi'(z_m) = 0, which is phi'(1) + phi(1)/2 = 0.\n\n",
          "- **Potential derived:** `", io_lib:format("~p", [Vh]), "`, which is 3/(4 z^2).\n",
          "- **Prediction:** m_n z_m = j_(0,n), the zeros of J_0. Here they come from `pq_mp:besselj_zero/3`, which\n",
          "  applies Newton's method to the Bessel series: a computation independent of the ODE solver.\n\n",
          "| n | (m_n z_m)^2, PQCF | j_(0,n)^2 | |\n| --- | --- | --- | --- |\n",
          [Row(N, A, B) || {N, {A, B}} <- lists:enumerate(1, lists:zip(Eh, Zeros))],
          io_lib:format("~nTimes: potential ~.2f s, spectrum ~.2f s.~n", [T3 / 1.0e6, T4 / 1.0e6])],
    ok = file:write_file(filename:join(Dir, "RESULTS.md"), Md),
    io:format("~s", [Md]).

pqcf_commit() -> os:getenv("PQCF_COMMIT", "(commit not given)").
today() -> {{Y, M, D}, _} = calendar:universal_time(), io_lib:format("~4..0b-~2..0b-~2..0b", [Y, M, D]).
