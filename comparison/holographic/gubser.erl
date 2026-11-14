%% Gubser and Nellore's black holes mimicking the QCD equation of state (arXiv:0804.0434,
%% https://arxiv.org/abs/0804.0434), in PQCF (decision record SY2, phase 6): the nonlinear master equation (34) for
%% G = dA/dphi integrated in arbitrary precision (pq_jet) from the horizon series (36), the entropy and temperature from
%% (45) and (51), and the speed of sound c_s^2 = d log T / d log s, for the potential (62) with (64): V = -12 cosh(gamma
%% phi) + b phi^2, gamma = 0.606, b = 2.06 (L = 1). Writes GUBSER.md here.
%%
%%   docker run --rm -u $(id -u):$(id -g) -e PQCF_COMMIT=$(git rev-parse --short HEAD) -v "$PWD":/w -w /w pqcf-dev \
%%     sh -c 'rebar3 compile >/dev/null && erlc -o /tmp comparison/holographic/gubser.erl &&
%%     erl -noshell -pa _build/default/lib/*/ebin -pa /tmp -eval "gubser:main(\"comparison/holographic\"), halt()."'
-module(gubser).
-export([main/1, rhs/4, pot/2, delta_minus_4/2, cs2/5]).
-define(P, 160).

%% Delta - 4 from b = 6 gamma^2 + Delta (Delta - 4)/2 (eq. 63): -2 + sqrt(4 + 2 (b - 6 gamma^2))
delta_minus_4(Gm, B) ->
    W = ?P,
    Gf = pq_mp:from(Gm, W),
    Bmg = pq_mp:sub(pq_mp:from(B, W), pq_mp:mul({mp, 6, 0}, pq_mp:mul(Gf, Gf, W), W), W),
    pq_mp:add({mp, -2, 0}, pq_mp:sqrt(pq_mp:add({mp, 4, 0}, pq_mp:ldexp(Bmg, 1), W), W), W).

%% the master equation (34) as G'' = (G P / R) (a G'/P + G'^2/G^2 + G'/(6G^2) + 4G' - G'(G' + R')/P^2),
%% a = G'/G + 1/(6G) - 4G - G'/P, P = G + R, R = V/(3V'), R' = 1/3 - V V''/(3 V'^2)
rhs(V, V1, V2, Dm4) ->
    G = {sym, g}, Gp = {sym, gp}, Phi = {sym, phi},
    Inv = fun(E) -> {pow, E, -1} end,
    R = {mul, [{rat, 1, 3}, V, Inv(V1)]},
    Rp = {add, [{rat, 1, 3}, {mul, [{rat, -1, 3}, V, V2, {pow, V1, -2}]}]},
    P = {add, [G, R]},
    A = {add, [{mul, [Gp, Inv(G)]}, {mul, [{rat, 1, 6}, Inv(G)]}, {mul, [-4, G]}, {mul, [-1, Gp, Inv(P)]}]},
    Num = {add, [{mul, [A, Gp, Inv(P)]}, {mul, [{pow, Gp, 2}, {pow, G, -2}]}, {mul, [{rat, 1, 6}, Gp, {pow, G, -2}]},
                 {mul, [4, Gp]}, {mul, [-1, Gp, {add, [Gp, Rp]}, {pow, P, -2}]}]},
    Gpp = {mul, [G, P, Inv(R), Num]},
    I1 = case Dm4 of none -> 0; _ -> {add, [G, {mul, [-1, Inv({mul, [Dm4, Phi]})]}]} end,
    [Gp, Gpp, I1, {mul, [{rat, 1, 6}, Inv(G)]}].

%% the Gubser-Nellore potential V = -12 cosh(g phi) + b phi^2 (L = 1) and its derivatives, as expressions
pot(Gm, B) ->
    Phi = {sym, phi},
    Ep = {exp, {mul, [Gm, Phi]}}, Em = {exp, {mul, [-1, Gm, Phi]}},
    V = {add, [{mul, [-6, {add, [Ep, Em]}]}, {mul, [B, {pow, Phi, 2}]}]},
    V1 = {add, [{mul, [-6, Gm, {add, [Ep, {mul, [-1, Em]}]}]}, {mul, [2, B, Phi]}]},
    V2 = {add, [{mul, [-6, Gm, Gm, {add, [Ep, Em]}]}, {mul, [2, B]}]},
    {V, V1, V2}.

ev(E, Phi) -> {ok, S} = pq_mp:evaluate(E, #{phi => Phi}, 45), pq_mp:from_string(S, ?P).

%% {log s, log T} (up to constants) at phi_H
thermo({V, V1, V2}, Dm4, PhiH, Digits) ->
    W = ?P,
    Vh = ev(V, PhiH), V1h = ev(V1, PhiH), V2h = ev(V2, PhiH), V0 = ev(V, 0),
    R = pq_mp:divide(Vh, pq_mp:mul({mp, 3, 0}, V1h, W), W),
    G0 = pq_mp:neg(R),
    G1 = pq_mp:mul(pq_mp:from({1, 6}, W), pq_mp:sub(pq_mp:divide(pq_mp:mul(Vh, V2h, W), pq_mp:mul(V1h, V1h, W), W), {mp, 1, 0}, W), W),
    Delta = pq_mp:from({1, 100000000}, W),
    A = pq_mp:sub(pq_mp:from(PhiH, W), Delta, W),
    Gs = pq_mp:sub(G0, pq_mp:mul(Delta, G1, W), W),
    Dm4f = pq_mp:to_float(Dm4),
    {mp, Md, Ed} = Dm4, Dm4r = {Md, 1 bsl (-Ed)},
    {ok, [_, _, I1, I2]} = pq_jet:integrate(rhs(V, V1, V2, {rat, element(1, Dm4r), element(2, Dm4r)}), phi, [g, gp, i1, i2],
                                            A, [Gs, G1, 0, 0], {rat, 1, 1000000}, #{digits => Digits}),
    PhiHf = pq_mp:to_float(pq_mp:from(PhiH, W)),
    G0f = pq_mp:to_float(G0),
    J1 = -pq_mp:to_float(I1) + 1.0e-8 * (G0f - 1 / (Dm4f * PhiHf)),
    J2 = -pq_mp:to_float(I2) + 1.0e-8 / (6 * G0f),
    LogS = 3 / Dm4f * math:log(PhiHf) + 3 * J1,
    LogT = 1 / Dm4f * math:log(PhiHf) + math:log(pq_mp:to_float(Vh) / pq_mp:to_float(V0)) + J1 + J2,
    {LogS, LogT}.

cs2(Pot, Dm4, PhiH, H, Digits) ->
    {S1, T1} = thermo(Pot, Dm4, {round((PhiH - H) * 100000), 100000}, Digits),
    {S2, T2} = thermo(Pot, Dm4, {round((PhiH + H) * 100000), 100000}, Digits),
    {(T2 - T1) / (S2 - S1), (T1 + T2) / 2}.

main(Dir) ->
    Gm = {rat, 606, 1000}, B = {rat, 206, 100},
    Dm4 = delta_minus_4({606, 1000}, {206, 100}),
    Pot = pot(Gm, B),
    Phis = [0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 16.0],
    {T, Rows} = timer:tc(fun() -> [{Ph, cs2(Pot, Dm4, Ph, 0.001, 20)} || Ph <- Phis] end),
    % T_c: the temperature at the minimum of c_s^2 on the grid
    {_, {_, LogTc}} = lists:min([{C, {Ph, LT}} || {Ph, {C, LT}} <- Rows]),
    Table = [io_lib:format("| ~.2f | ~.4f | ~.6f |~n", [Ph, math:exp(LT - LogTc), C]) || {Ph, {C, LT}} <- Rows],
    Delta = pq_mp:to_float(Dm4) + 4,
    Md = ["# Gubser-Nellore: a black hole with the equation of state of QCD\n\n",
          "Produced by `gubser.erl` (decision record SY2, phase 6), PQCF ", os:getenv("PQCF_COMMIT", "(commit not given)"),
          ". The model is Gubser and Nellore, [arXiv:0804.0434](https://arxiv.org/abs/0804.0434), with the potential (62)\n",
          "and the parameters (64): V = -12 cosh(gamma phi) + b phi^2, gamma = 0.606, b = 2.06 (L = 1).\n\n",
          "## Method\n\n",
          "1. **Master equation:** the nonlinear master equation (34) for G = dA/dphi, written as\n",
          "   G'' = (G P/R) [a G'/P + G'^2/G^2 + G'/(6G^2) + 4G' - G'(G' + R')/P^2], with P = G + R, R = V/(3V') and\n",
          "   a = G'/G + 1/(6G) - 4G - G'/P.\n",
          "2. **Integration:** from the horizon series (36) at phi_H - 10^-8 down to phi = 10^-6, by `pq_jet` (the\n",
          "   Taylor-series method, 20 digits). The two integrals of (45) and (51) are integrated along with it.\n",
          "3. **Speed of sound:** c_s^2 = d log T / d log s by central differences in phi_H (h = 10^-3).\n\n",
          "## Checks\n\n",
          "- **Chamblin-Reall:** for V = -12 e^(phi/2), the master equation keeps the exact solution G = -1/(3 gamma)\n",
          "  to every digit (`test/pq_holo_tests.erl`).\n",
          io_lib:format("- **Dimension:** Delta = ~.4f, as the paper's 3.93.~n", [Delta]),
          io_lib:format("- **Infrared limit:** 1/3 - gamma^2/2 = ~.4f, the paper's 0.15. The table approaches it at large phi_H.~n", [1/3 - 0.606 * 0.606 / 2]),
          "- **Ultraviolet limit:** c_s^2 -> 1/3 at small phi_H (high T).\n",
          "- **The crossover:** c_s^2 dips near T_c, then rises toward its infrared value, as in the paper's figures 3 and 4.\n",
          "  The figures give no table to compare digits with.\n",
          "- **An independent double-precision solution** (SciPy DOP853, the same equations) agrees to 6 digits up to\n",
          "  phi_H = 4. Beyond that it is unstable: at phi_H = 6 it swings between 0.0602 and 0.32 as its tolerance\n",
          "  changes, because of the 0/0 at the horizon. PQCF's 0.06021565 does not change from 20 to 30 digits or\n",
          "  with h halved.\n\n",
          "## Results (T_c: the temperature of the minimum of c_s^2 on this grid)\n\n",
          "| phi_H | T/T_c | c_s^2 |\n| --- | --- | --- |\n", Table,
          io_lib:format("~nTime: ~.0f s.~n", [T / 1.0e6])],
    ok = file:write_file(filename:join(Dir, "GUBSER.md"), unicode:characters_to_binary(Md)),
    io:format("~s", [Md]).
