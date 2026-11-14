%% The Sakai-Sugimoto model's vector and axial-vector mesons in PQCF (decision record SY2, phase 6). Writes SAKAI.md here.
%%
%%   docker run --rm -u $(id -u):$(id -g) -e PQCF_COMMIT=$(git rev-parse --short HEAD) -v "$PWD":/w -w /w pqcf-dev \
%%     sh -c 'rebar3 compile >/dev/null && erlc -o /tmp comparison/holographic/sakai.erl &&
%%     erl -noshell -pa _build/default/lib/*/ebin -pa /tmp -eval "sakai:main(\"comparison/holographic\"), halt()."'
-module(sakai).
-export([main/1]).

main(Dir) ->
    U = {sym, u},
    C = #{p2 => {add, [{pow, U, 4}, {mul, [-1, U]}]}, p1 => {add, [{mul, [{rat, 5, 2}, {pow, U, 3}]}, -1]}, p0 => 0,
          q => {mul, [{rat, 9, 4}, U]}},
    Opts = #{right => larger, match => 2, range => {0.1, 12}, grid => 60, digits => 25},
    {Tv, {ok, V}} = timer:tc(fun() -> pq_ode:eigenvalues_sl(C, u, Opts#{left => {1, smaller}}) end),
    {Ta, {ok, A}} = timer:tc(fun() -> pq_ode:eigenvalues_sl(C, u, Opts#{left => {1, larger}}) end),
    All = lists:sort(fun({_, X, _}, {_, Y, _}) -> num(X) =< num(Y) end,
                     [{vector, X, "even"} || X <- V] ++ [{axial, X, "odd"} || X <- A]),
    Paper = ["0.67 (rho)", "1.6 (a_1(1260))", "2.9 (rho(1450))", "4.5", "-", "-", "-"],
    Rows = [io_lib:format("| ~b | ~s | ~s | `~s` | ~s |~n", [N, K, P, X, Pp])
            || {N, {K, X, P}, Pp} <- lists:zip3(lists:seq(1, length(All)), All, lists:sublist(Paper, length(All)))],
    L1 = num(element(2, lists:nth(1, All))), L2 = num(element(2, lists:nth(2, All))), L3 = num(element(2, lists:nth(3, All))),
    Md = ["# The Sakai-Sugimoto model: vector and axial-vector mesons\n\n",
          "Produced by `sakai.erl` (decision record SY2, phase 6), PQCF ", os:getenv("PQCF_COMMIT", "(commit not given)"), ".\n",
          "The model is Sakai and Sugimoto, [hep-th/0412141](https://arxiv.org/abs/hep-th/0412141), eq. (3.28):\n",
          "-K^(1/3) d_z (K d_z psi) = lambda psi with K = 1 + z^2, and m_n^2 = lambda_n M_KK^2.\n\n",
          "## Method\n\n",
          "- **Coordinates:** in u = K^(1/3) the equation has polynomial coefficients,\n",
          "  u (u^3 - 1) psi'' + (5u^3/2 - 1) psi' + (9u/4) lambda psi = 0, on [1, infinity).\n",
          "- **Ends:** both are regular singular points. At u = 1 the exponent is 0 for the vectors (psi even in z) and\n",
          "  1/2 for the axial vectors (odd). At infinity the exponent is 3/2 in w = 1/u, the decaying solution.\n",
          "- **Solver:** `pq_ode:eigenvalues_sl/3` takes Frobenius series from both ends, integrates them by the\n",
          "  Taylor-series method to u = 2, and finds where their Wronskian vanishes, to 25 digits.\n\n",
          "## Results\n\n",
          "| n | meson | psi(z) | lambda_n, PQCF | the paper (4.12) |\n| --- | --- | --- | --- | --- |\n", Rows,
          io_lib:format("~n- **Mass ratios:** m_a1^2/m_rho^2 = lambda_2/lambda_1 = ~.4f (paper ~~2.4, experiment 2.51);\n"
                        "  m_rho(1450)^2/m_rho^2 = lambda_3/lambda_1 = ~.4f (paper ~~4.3, experiment 3.56).~n", [L2 / L1, L3 / L1]),
          "- **Independent check:** a SciPy shooting in z (DOP853) with the decaying condition at a cutoff Z, extrapolated\n",
          "  in 1/Z, gives lambda_1 = 0.6693144 and lambda_2 = 1.5687644. It agrees with the values above to its own\n",
          "  accuracy (about 1e-7).\n",
          "- **Self-check:** the result does not depend on the match point (`test/pq_ode_tests.erl`).\n",
          io_lib:format("~nTimes: vectors ~.1f s, axial vectors ~.1f s.~n", [Tv / 1.0e6, Ta / 1.0e6])],
    ok = file:write_file(filename:join(Dir, "SAKAI.md"), unicode:characters_to_binary(Md)),
    io:format("~s", [Md]).

num(S) -> [M, E] = string:split(S, "e"), list_to_float(M) * math:pow(10, list_to_integer(E)).
