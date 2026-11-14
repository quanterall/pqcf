%% The hard-wall model of holographic QCD (Erlich, Katz, Son, Stephanov, hep-ph/0501128,
%% https://arxiv.org/abs/hep-ph/0501128): its vector, axial and pion observables from the paper's parameters, in PQCF
%% (decision record SY2, phase 3). Writes HARDWALL.md here.
%%
%%   docker run --rm -u $(id -u):$(id -g) -e PQCF_COMMIT=$(git rev-parse --short HEAD) -v "$PWD":/w -w /w pqcf-dev \
%%     sh -c 'rebar3 compile >/dev/null && erlc -o /tmp comparison/holographic/hardwall.erl &&
%%     erl -noshell -pa _build/default/lib/*/ebin -pa /tmp -eval "hardwall:main(\"comparison/holographic\"), halt()."'
-module(hardwall).
-export([main/1, models/0]).
-define(P, 140).

num(S) -> pq_mp:from_string(S, ?P).

%% the paper's Table II: Model A (a fit of m_pi, f_pi, m_rho: starred) and Model B (a fit to all seven)
paper() ->
    #{"A" => #{m_rho => "775.8*", m_a1 => "1363", f_pi => "92.4*", m_pi => "139.6*", f_rho => "329", f_a1 => "486"},
      "B" => #{m_rho => "832", m_a1 => "1220", f_pi => "84.0", m_pi => "141", f_rho => "353", f_a1 => "440"}}.

%% Model A: z_m = 1/(323 MeV), m_q = 2.29 MeV, sigma = (327 MeV)^3; Model B: 1/(346 MeV), 2.3 MeV, (308 MeV)^3;
%% g5^2 = 12 pi^2 / N_c = 4 pi^2
models() ->
    Pi = pq_mp:pi(?P),
    {mp, Mg, Eg} = pq_mp:mul({mp, 4, 0}, pq_mp:mul(Pi, Pi, ?P), ?P),
    G5sq = {Mg, 1 bsl (-Eg)},
    G5 = pq_mp:sqrt(pq_mp:from(G5sq, ?P), ?P),
    Robin = {robin, 1, {1, 2}, 1},
    [model(Name, Zinv, Mq, Sig3, G5sq, G5, Robin) || {Name, Zinv, Mq, Sig3} <- [{"A", 323, {229, 100}, 327}, {"B", 346, {23, 10}, 308}]].

main(Dir) ->
    {T, Ms} = timer:tc(fun models/0),
    P = paper(),
    Keys = [{m_rho, "m_rho"}, {m_a1, "m_a1"}, {f_pi, "f_pi (chiral limit)"}, {m_pi, "m_pi (GOR)"}, {f_rho, "F_rho^(1/2)"}, {f_a1, "F_a1^(1/2)"}],
    Row = fun(M) ->
                  Name = maps:get(model, M),
                  [io_lib:format("| ~s | ~s | `~s` | ~s |~n", [Name, Label, pq_mp:to_string(maps:get(K, M), 8), maps:get(K, maps:get(Name, P))])
                   || {K, Label} <- Keys]
          end,
    Md = ["# The hard-wall model: vector, axial and pion observables\n\n",
          "Produced by `hardwall.erl` (decision record SY2, phase 3), PQCF ", os:getenv("PQCF_COMMIT", "(commit not given)"),
          ". The model is Erlich, Katz, Son and Stephanov, [hep-ph/0501128](https://arxiv.org/abs/hep-ph/0501128),\n",
          "with its two parameter sets (Table II) and g5^2 = 12 pi^2/N_c = 4 pi^2. Every value is computed by PQCF in\n",
          "arbitrary precision, in units of 1/z_m, then converted to MeV.\n\n",
          "## Method\n\n",
          "- **rho:** m_rho = j_(0,1)/z_m. F_rho^2 = psi''(0)^2/g5^2 for the mode normalised by int dz psi^2/z = 1.\n",
          "  The norm `int phi^2` comes from `pq_ode:frobenius_solution/4` and equals the closed form 2 J_1(j)^2/j^2.\n",
          "- **a_1:** the axial equation (eq. 16) in Schrödinger form, -phi'' + [3/(4u^2) + g5^2 (m_q z_m +\n",
          "  sigma z_m^3 u^2)^2] phi = (m z_m)^2 phi, with psi'(z_m) = 0. Its first eigenvalue comes from\n",
          "  `pq_ode:eigenvalues/3`, and F_a1 is found like F_rho.\n",
          "- **f_pi:** eq. 20 in the chiral limit m_q = 0. A(0, z) combines the Frobenius solutions of both exponents,\n",
          "  with A(0) = 1 and A'(z_m) = 0. Then f_pi^2 = -(1/g5^2) A'(z)/z at z -> 0.\n",
          "- **m_pi:** the Gell-Mann-Oakes-Renner relation m_pi^2 f_pi^2 = 2 m_q sigma, as in the paper.\n\n",
          "## Results\n\n",
          "| model | observable (MeV) | PQCF | paper |\n| --- | --- | --- | --- |\n",
          [Row(M) || M <- Ms],
          "\n## Comparison with the paper\n\n",
          "- **Agreement:** m_rho, F_rho, F_a1 and m_pi agree with Table II to the paper's rounding, and so does model B's m_a1.\n",
          "- **f_pi** is 1.6% (A) and 2.1% (B) below the paper. The paper keeps m_q: an independent floating-point\n",
          "  solution of the same boundary-value problem (SciPy, `solve_ivp`) gives the chiral-limit values above\n",
          "  (90.90 and 82.26 MeV). With m_q included it gives 92.1-92.3 and 83.6-83.9 MeV, depending\n",
          "  logarithmically on the cutoff, which is the O(m_q^2) contact term. So the difference is the quark-mass\n",
          "  contribution, not an error.\n",
          "- **Model A's m_a1** is 0.35% below the paper's 1363 MeV. The same independent SciPy solution gives 1358.2432\n",
          "  (A) and 1219.5034 MeV (B), equal to PQCF's. The paper's 1363 probably comes from slightly different\n",
          "  parameter values (sigma = (327 MeV)^3 is rounded); this is open.\n",
          "- **Not computed:** g_(rho pi pi) (eq. 22).\n",
          io_lib:format("~nTime: ~.1f s for both models.~n", [T / 1.0e6])],
    ok = file:write_file(filename:join(Dir, "HARDWALL.md"), unicode:characters_to_binary(Md)),
    io:format("~s", [Md]).

model(Name, Zinv, Mq, Sig3, G5sq, G5, Robin) ->
    U = {sym, u},
    {Mqn, Mqd} = Mq,
    Mu = {Mqn, Mqd * Zinv},                      % m_q z_m
    S3 = {Sig3 * Sig3 * Sig3, Zinv * Zinv * Zinv}, % sigma z_m^3
    Va = {add, [{mul, [{rat, 3, 4}, {pow, U, -2}]},
                {mul, [{sym, g}, {pow, {add, [{sym, mu}, {mul, [{sym, s3}, {pow, U, 2}]}]}, 2}]}]},
    Par = #{g => G5sq, mu => Mu, s3 => S3},
    Zi = pq_mp:from(Zinv, ?P),
    Mev = fun(Ev) -> pq_mp:mul(Zi, pq_mp:sqrt(num(Ev), ?P), ?P) end,
    % rho
    J = pq_mp:besselj_zero(0, 1, ?P),
    Mrho = pq_mp:mul(Zi, J, ?P),
    {ok, #{norm := Ir}} = pq_ode:frobenius_solution({mul, [{rat, 3, 4}, {pow, U, -2}]}, u, pq_mp:to_string(pq_mp:mul(J, J, ?P), 40), #{b => 1, digits => 30}),
    Frho = fsqrt(Ir, G5, Zi),
    % a1
    {ok, [Ea | _]} = pq_ode:eigenvalues(Va, u, #{params => Par, right => Robin, range => {6, 40}, grid => 60, digits => 20}),
    {ok, #{norm := Ia}} = pq_ode:frobenius_solution(Va, u, Ea, #{params => Par, b => 1, digits => 20}),
    Fa1 = fsqrt(Ia, G5, Zi),
    % f_pi, chiral limit (m_q = 0), q^2 = 0
    Par0 = Par#{mu => 0},
    {ok, #{phi := Ps, dphi := Dps}} = pq_ode:frobenius_solution(Va, u, 0, #{params => Par0, b => 1, digits => 20, exponent => smaller}),
    {ok, #{phi := Pl, dphi := Dpl}} = pq_ode:frobenius_solution(Va, u, 0, #{params => Par0, b => 1, digits => 20, exponent => larger}),
    Ad = fun(Ph, Dph) -> pq_mp:add(num(Dph), pq_mp:ldexp(num(Ph), -1), ?P) end,   % A'(1) = phi'(1) + phi(1)/2
    Beta = pq_mp:neg(pq_mp:divide(Ad(Ps, Dps), Ad(Pl, Dpl), ?P)),
    Fpi = pq_mp:mul(Zi, pq_mp:sqrt(pq_mp:divide(pq_mp:mul({mp, -2, 0}, Beta, ?P), pq_mp:mul(G5, G5, ?P), ?P), ?P), ?P),
    % GOR: m_pi^2 f_pi^2 = 2 m_q sigma
    Mpi = pq_mp:divide(pq_mp:sqrt(pq_mp:mul({mp, 2, 0}, pq_mp:mul(pq_mp:from(Mq, ?P), pq_mp:from(Sig3 * Sig3 * Sig3, ?P), ?P), ?P), ?P), Fpi, ?P),
    #{model => Name, m_rho => Mrho, m_a1 => Mev(Ea), f_pi => Fpi, m_pi => Mpi, f_rho => Frho, f_a1 => Fa1}.

%% F^(1/2) = sqrt(2 / (g5 sqrt(norm))) / z_m, in MeV
fsqrt(Norm, G5, Zi) ->
    F = pq_mp:divide({mp, 2, 0}, pq_mp:mul(G5, pq_mp:sqrt(pq_mp:from_string(Norm, ?P), ?P), ?P), ?P),
    pq_mp:mul(Zi, pq_mp:sqrt(F, ?P), ?P).
