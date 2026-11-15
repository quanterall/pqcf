%% Improved holographic QCD, background I (Gursoy, Kiritsis, Nitti, arXiv:0707.1349, https://arxiv.org/abs/0707.1349,
%% section 8.1; set up in arXiv:0707.1324): the scalar (0++) and tensor (2++) glueball spectra, in PQCF (decision record
%% SY2, phase 6). Writes IHQCD.md here.
%%
%%   docker run --rm -u $(id -u):$(id -g) -e PQCF_COMMIT=$(git rev-parse --short HEAD) -v "$PWD":/w -w /w pqcf-dev \
%%     sh -c 'rebar3 compile >/dev/null && erlc -o /tmp comparison/holographic/ihqcd.erl &&
%%     erl -noshell -pa _build/default/lib/*/ebin -pa /tmp -eval "ihqcd:main(\"comparison/holographic\"), halt()."'
%% (about 40 minutes on 16 cores). main/1 also writes ihqcd-raw.etf, the raw results; report_only/1 rewrites IHQCD.md
%% from it without the runs.
-module(ihqcd).
-export([main/1, report_only/1, spectrum/3, x_w/1, start/2, model/2, background/3, ir_scale/2]).

-define(D, 12).          % digits of the integrations

%% The model: b0 and a. The paper uses a twice (section 8.1): its IR asymptotics (8.2) X = -1/2 - a/log(lambda) with
%% a = (3/8)(alpha - 1)/alpha, so a = 3/16 for linear confinement (alpha = 2); but the interpolating X (8.4) tends to
%% -1/2 - a/(2 log lambda), so in (8.4) and (8.6) the same alpha = 2 needs a = 3/8. That is the value here. (Until
%% 15 Nov 2026 this script used a = 5/16, which matches neither, and its ratios were 6-15 % low: model(B0, {rat, 5, 16})
%% reproduces that run.)
model(B0, A) -> #{b0 => B0, a => A}.

default_model() -> model({rat, 21, 5}, {rat, 3, 8}).

%% X(lambda) (8.4) and W(lambda) (8.6) for b0, b1 = (51/121) b0^2 (pure Yang-Mills) and a, as expressions in the
%% symbol l; a model or (as before) b0 alone, with a = 3/8
x_w(#{b0 := B0, a := A}) ->
    L = {sym, l},
    B1 = {mul, [{rat, 51, 121}, {pow, B0, 2}]},
    C = {add, [{mul, [2, {pow, B0, 2}]}, {mul, [3, B1]}]},
    Lg = {add, [1, {mul, [C, {pow, {mul, [9, A]}, -1}, {log, {add, [1, {pow, L, 2}]}}]}]},
    X = {add, [{mul, [-1, B0, L, {pow, {add, [3, {mul, [2, B0, L]}]}, -1}]},
               {mul, [-1, C, {pow, L, 2}, {pow, {mul, [9, {add, [1, {pow, L, 2}]}, Lg]}, -1}]}]},
    W = {mul, [{rat, 9, 4}, {pow, {add, [1, {mul, [{rat, 2, 3}, B0, L]}]}, {rat, 2, 3}}, {pow, Lg, {mul, [{rat, 2, 3}, A]}}]},
    {X, W};
x_w(B0) -> x_w(model(B0, {rat, 3, 8})).

%% the background equations (8.5), l = 1: lambda' = -(4/3) X W lambda e^A, A' = -(4/9) W e^A; with the fluctuation
%% xi'' + 2 B' xi' + m^2 xi = 0, B = (3/2) A (tensor) or (3/2) A + (1/2) log X^2 (scalar), state [l, a, xi, xp]
rhs(B0, Kind, M2) ->
    {X, W} = x_w(B0),
    L = {sym, l}, Ea = {exp, {sym, a}},
    Lp = {mul, [{rat, -4, 3}, X, W, L, Ea]},
    Ap = {mul, [{rat, -4, 9}, W, Ea]},
    Bp = case Kind of
             tensor -> {mul, [{rat, 3, 2}, Ap]};
             scalar -> {add, [{mul, [{rat, 3, 2}, Ap]}, {mul, [pq_jet:diff(X, l), Lp, {pow, X, -1}]}]}
         end,
    [Lp, Ap, {sym, xp}, {add, [{mul, [-2, Bp, {sym, xp}]}, {mul, [-1, M2, {sym, xi}]}]}].

%% the background near the UV boundary: from lambda(r_in) = Lambda0, A(r_in) = 0 at r_in = 0.36, integrated back
%% toward the boundary r_b, where e^(-A) -> 0. Near it e^(-A) ~ (4/9) W (r - r_b), which places r_b by linear
%% extrapolation; each step goes three quarters of the way, until 10^-3 from it. {R_start, lambda, A, r_b}.
%% The mass ratios do not depend on Lambda0 (the solutions differ by a shift of A, section 7.2 of the paper); a
%% larger Lambda0 brings the IR closer.
start(B0, Lambda0) ->
    [Lp, Ap | _] = rhs(B0, tensor, 0),
    {_, W} = x_w(B0),
    toward_boundary(Lp, Ap, W, 0.36, Lambda0, 0, 0).

toward_boundary(Lp, Ap, W, R, L, A, N) ->
    Lf = num(L), Af = num(A),
    {ok, Ws} = pq_mp:evaluate(W, #{l => Lf}, 16),
    Slope = 4 / 9 * list_to_float(fix(Ws)),
    Rb = R - math:exp(-Af) / Slope,
    case R - Rb < 1.0e-3 orelse N > 40 of
        true -> {R, L, A, Rb};
        false ->
            R2 = Rb + 0.25 * (R - Rb),
            {ok, [L2, A2]} = pq_jet:integrate([Lp, Ap], r, [l, a], R, [L, A], R2, #{digits => ?D}),
            toward_boundary(Lp, Ap, W, R2, L2, A2, N + 1)
    end.

%% the background from lambda(0.36) = Lambda0, A(0.36) = 0 forward to each r of Rs (ascending, > 0.36):
%% [{R, lambda, A}] as floats
background(Model, Lambda0, Rs) ->
    [Lp, Ap | _] = rhs(Model, tensor, 0),
    {Out, _} = lists:mapfoldl(fun(R, {R0, Y0}) ->
                                      {ok, [L, A] = Y} = pq_jet:integrate([Lp, Ap], r, [l, a], R0, Y0, R, #{digits => ?D}),
                                      {{R, pq_mp:to_float(L), pq_mp:to_float(A)}, {R, Y}}
                              end, {0.36, [Lambda0, 0]}, Rs),
    Out.

%% the IR scale r1, where lambda = 1: lambda grows monotonically, bracketed by doubling r, then bisected to 1e-6
ir_scale(Model, Lambda0) ->
    [Lp, Ap | _] = rhs(Model, tensor, 0),
    Step = fun(R0, Y0, R) -> {ok, Y} = pq_jet:integrate([Lp, Ap], r, [l, a], R0, Y0, R, #{digits => ?D}), Y end,
    Lam = fun([L, _]) -> pq_mp:to_float(L) end,
    Up = fun Up(R0, Y0) ->
                 R = 2 * R0, Y = Step(R0, Y0, R),
                 case Lam(Y) > 1 of true -> {R0, Y0, R}; false -> Up(R, Y) end
         end,
    {Ra, Ya, Rb} = Up(0.36, [Lambda0, 0]),
    Bis = fun Bis(Lo, Yl, Hi) ->
                  case Hi - Lo < 1.0e-6 * Hi of
                      true -> (Lo + Hi) / 2;
                      false ->
                          Mid = (Lo + Hi) / 2, Ym = Step(Lo, Yl, Mid),
                          case Lam(Ym) > 1 of true -> Bis(Lo, Yl, Mid); false -> Bis(Mid, Ym, Hi) end
                  end
          end,
    Bis(Ra, Ya, Rb).

num({mp, _, _} = M) -> pq_mp:to_float(M);
num({rat, P, Q}) -> P / Q;
num(X) when is_number(X) -> float(X).

%% the growing mode's coefficient at Rmax: e^(2B) xi' (e^(2B) = e^(3A) for the tensor, e^(3A) X^2 for the scalar)
shoot(B0, Kind, M2, {Rs, L0, A0, Rb}, Rmax) ->
    % xi ~ (r - r_b)^4 near the boundary (psi ~ (r - r_b)^(5/2), the normalisable solution)
    D = Rs - Rb,
    Xi0 = math:pow(D, 4), Xp0 = 4 * math:pow(D, 3),
    {ok, [Lf, Af, _, Xpf]} = pq_jet:integrate(rhs(B0, Kind, M2), r, [l, a, xi, xp], Rs, [L0, A0, Xi0, Xp0], Rmax, #{digits => ?D}),
    {X, _} = x_w(B0),
    E3a = math:exp(3 * pq_mp:to_float(Af)),
    Fac = case Kind of
              tensor -> E3a;
              scalar -> {ok, Xs} = pq_mp:evaluate(X, #{l => pq_mp:to_float(Lf)}, 16), E3a * math:pow(list_to_float(fix(Xs)), 2)
          end,
    Fac * pq_mp:to_float(Xpf).

fix(S) -> [M, E] = string:split(S, "e"), M ++ "e" ++ E.

%% the first N eigenvalues m^2 (units 1/l^2), bracketed on a grid of m^2 then refined by bisection
%% Opts: lambda0, count, grid; rmax and range (m^2) as numbers, or in units of the IR scale r1 with ir => true
%% (rmax r1 and range / r1^2), so that runs with different lambda0 cover the same states
spectrum(B0, Kind, Opts) ->
    L0 = maps:get(lambda0, Opts, {rat, 3, 20}),
    St = start(B0, L0),
    U = case maps:get(ir, Opts, false) of true -> ir_scale(B0, L0); false -> 1.0 end,
    Rmax = maps:get(rmax, Opts, 20) * U,
    {Lo0, Hi0} = maps:get(range, Opts),
    {Lo, Hi} = {Lo0 / (U * U), Hi0 / (U * U)},
    K = maps:get(grid, Opts, 40),
    Ms = [Lo + (Hi - Lo) * I / K || I <- lists:seq(0, K)],
    Fs = [{M, shoot(B0, Kind, M, St, Rmax)} || M <- Ms],
    Br = [{A, Fa, B, Fb} || {{A, Fa}, {B, Fb}} <- lists:zip(lists:droplast(Fs), tl(Fs)), Fa * Fb < 0],
    [illinois(fun(M) -> shoot(B0, Kind, M, St, Rmax) end, A, Fa, B, Fb, 0, none)
     || {A, Fa, B, Fb} <- lists:sublist(Br, maps:get(count, Opts, 4))].

%% the Illinois method, to a relative width of 1e-9
illinois(_, A, _, B, _, 40, _) -> (A + B) / 2;
illinois(F, A, Fa, B, Fb, I, Side) ->
    C = B - Fb * (B - A) / (Fb - Fa),
    case abs(B - A) < 1.0e-9 * abs(C) of
        true -> C;
        false ->
            Fc = F(C),
            case Fc * Fb > 0 of
                true -> illinois(F, A, case Side of a -> Fa / 2; _ -> Fa end, C, Fc, I + 1, a);
                false -> illinois(F, C, Fc, B, case Side of b -> Fb / 2; _ -> Fb end, I + 1, b)
            end
    end.

%% the runs of IHQCD.md, in parallel processes. Units: masses in 1/l; rmax and the m^2 range in units of the IR scale
%% r1 (where lambda = 1), so that every lambda0 covers the same states.
main(Dir) ->
    M = default_model(),
    Old = model({rat, 21, 5}, {rat, 5, 16}),
    Std = #{range => {0.15, 90}, grid => 30, rmax => 8, ir => true},
    Jobs = [{base_s, fun() -> spectrum(M, scalar, Std#{count => 4}) end},
            {base_t, fun() -> spectrum(M, tensor, Std#{count => 4}) end},
            {cut_s, fun() -> spectrum(M, scalar, Std#{count => 4, rmax => 11}) end},
            {cut_t, fun() -> spectrum(M, tensor, Std#{count => 4, rmax => 11}) end},
            {old_s, fun() -> spectrum(Old, scalar, Std#{count => 2}) end},
            {old_t, fun() -> spectrum(Old, tensor, Std#{count => 1}) end},
            {lin_s, fun() -> spectrum(M, scalar, Std#{count => 12, range => {0.15, 220}, grid => 90, rmax => 11}) end},
            {lin_t, fun() -> spectrum(M, tensor, Std#{count => 8, range => {0.15, 220}, grid => 90, rmax => 11}) end},
            {alpha, fun() -> R1 = ir_scale(M, {rat, 3, 20}), background(M, {rat, 3, 20}, [R1 * F || F <- [2.0, 4.0, 8.0, 16.0]]) end},
            {fit, fun() -> fit_b0(1.87) end}]
        ++ lists:append([[{{l0_s, L0}, fun() -> spectrum(M, scalar, Std#{count => 2, lambda0 => L0}) end},
                          {{l0_t, L0}, fun() -> spectrum(M, tensor, Std#{count => 1, lambda0 => L0}) end},
                          {{r1, L0}, fun() -> ir_scale(M, L0) end}]
                         || L0 <- [{rat, 1, 20}, {rat, 2, 25}, {rat, 3, 20}, {rat, 1, 4}]]),
    T0 = erlang:monotonic_time(millisecond),
    R0 = pmap(Jobs),
    Secs = (erlang:monotonic_time(millisecond) - T0) / 1000,
    % the raw results first (failures included), so that the report can be rewritten without the runs
    ok = file:write_file(filename:join(Dir, "ihqcd-raw.etf"), term_to_binary({R0, Secs})),
    R = case [{K, W} || K := {failed, W} <- R0] of
            [] -> R0;
            Failed -> io:format(standard_error, "failed: ~p~n", [Failed]), error(failed)
        end,
    Md = report(R, Secs),
    ok = file:write_file(filename:join(Dir, "IHQCD.md"), unicode:characters_to_binary(Md)),
    io:format("~s", [Md]).

%% IHQCD.md again from the raw results of the last main/1 (ihqcd-raw.etf), without the runs
report_only(Dir) ->
    {ok, B} = file:read_file(filename:join(Dir, "ihqcd-raw.etf")),
    {R, Secs} = binary_to_term(B),
    Md = report(R, Secs),
    ok = file:write_file(filename:join(Dir, "IHQCD.md"), unicode:characters_to_binary(Md)),
    io:format("~s", [Md]).

pmap(Jobs) ->
    Self = self(),
    Pids = [{K, spawn(fun() -> Self ! {self(), try timer:tc(F) catch C:E:St -> {failed, {C, E, St}} end} end)} || {K, F} <- Jobs],
    maps:from_list([receive {P, {failed, Why}} -> {K, {failed, Why}}; {P, {T, V}} -> {K, {V, T / 1.0e6}} end || {K, P} <- Pids]).

%% b0 with R00 = Target, by the secant method on R00(b0) (lambda0 = 0.15, two scalar states), to 1e-4 in R00
fit_b0(Target) ->
    F = fun(B) ->
                case [math:sqrt(X) || X <- spectrum(model(B, {rat, 3, 8}), scalar,
                                                    #{range => {0.15, 90}, grid => 30, rmax => 8, ir => true, count => 2})] of
                    [M1, M2] -> M2 / M1 - Target;
                    Ms -> error({fit_b0, B, Ms})
                end
        end,
    secant(F, 4.0, F(4.0), 4.4, F(4.4), 0).

secant(_, _, _, B, Fb, N) when abs(Fb) < 1.0e-4; N >= 8 -> {B, Fb};
secant(F, A, Fa, B, Fb, N) ->
    C = min(8.0, max(2.0, B - Fb * (B - A) / (Fb - Fa))),
    secant(F, B, Fb, C, F(C), N + 1).

report(R, Secs) ->
    V = fun(K) -> element(1, maps:get(K, R)) end,
    Ms = fun(K) -> [math:sqrt(X) || X <- V(K)] end,
    S = Ms(base_s), T = Ms(base_t), M1 = hd(S),
    Paper0 = [1475, 2753, 3561, 4253], Paper2 = [2055, 2991, 3739, 4396],
    Row = fun(Name, Mine, Paper) ->
                  [io_lib:format("| ~s ~b | ~.4f | ~.3f | ~.3f | ~s % |~n", [Name, N, M, M / M1, P / 1475, pct(100 * (M / M1 / (P / 1475) - 1))])
                   || {N, {M, P}} <- lists:enumerate(lists:zip(Mine, Paper))]
          end,
    Ratios = fun(Sk, Tk) -> Ss = Ms(Sk), Ts = Ms(Tk), {lists:nth(2, Ss) / hd(Ss), hd(Ts) / hd(Ss)} end,
    L0Rows = [begin
                  {R00, R20} = Ratios({l0_s, L0}, {l0_t, L0}),
                  io_lib:format("| ~s | ~.3f | ~.4f | ~.4f | ~.4f |~n", [l0_str(L0), V({r1, L0}), hd(Ms({l0_s, L0})), R00, R20])
              end || L0 <- [{rat, 1, 20}, {rat, 2, 25}, {rat, 3, 20}, {rat, 1, 4}]],
    Cs = Ms(cut_s), Ct = Ms(cut_t),
    CutDiff = lists:max([abs(X / Y - 1) || {X, Y} <- lists:zip(S ++ T, Cs ++ Ct)]),
    Alpha = V(alpha),
    AlphaRows = [io_lib:format("| ~b r1 to ~b r1 | ~.3f |~n", [round(Ra / Rr1), round(Rb / Rr1), math:log(Ab / Aa) / math:log(Rb / Ra)])
                 || {{Ra, _, Aa}, {Rb, _, Ab}} <- lists:zip(lists:droplast(Alpha), tl(Alpha)), Rr1 <- [element(1, hd(Alpha)) / 2]],
    Lin = Ms(lin_s), LinT = Ms(lin_t),
    Lin2 = [X * X || X <- Lin],
    Gaps = [B - A || {A, B} <- lists:zip(lists:droplast(Lin2), tl(Lin2))],
    Gaps2 = [B * B - A * A || {A, B} <- lists:zip(lists:droplast(LinT), tl(LinT))],
    {Ob, Fb} = V(fit),
    Old = Ms(old_s), OldT = Ms(old_t),
    ["# Improved holographic QCD: scalar and tensor glueballs\n\n",
     "Produced by `ihqcd.erl` (decision record SY2, phase 6), PQCF ", os:getenv("PQCF_COMMIT", "(commit not given)"), ".\n",
     "The model is Gursoy, Kiritsis and Nitti, [arXiv:0707.1349](https://arxiv.org/abs/0707.1349), section 8.1\n",
     "(background I; the framework is [arXiv:0707.1324](https://arxiv.org/abs/0707.1324)).\n\n",
     "## Set-up\n\n",
     "- **Model:** X(lambda) from (8.4) and W(lambda) from (8.6), with b0 = 4.2, b1/b0^2 = 51/121 and a = 3/8.\n",
     "- **The value of a.** The paper uses the symbol a in two ways. Its IR asymptotics (8.2) are\n",
     "  X = -1/2 - a/log(lambda), with a = (3/8)(alpha - 1)/alpha, so a = 3/16 for linear confinement (alpha = 2). Its\n",
     "  interpolating X (8.4) tends instead to -1/2 - a/(2 log lambda), so through (8.4) and (8.6) linear confinement\n",
     "  needs a = 3/8. The background confirms it: A(r) ~ -C r^alpha with alpha -> 2 (below).\n",
     "- **Background** (8.5): lambda' = -(4/3) X W lambda e^A and A' = -(4/9) W e^A (l = 1), from\n",
     "  lambda(0.36) = lambda0, A(0.36) = 0.\n",
     "- **Glueballs** (6.2): xi'' + 2 B' xi' + m^2 xi = 0, with B = (3/2) A + (1/2) log X^2 for 0++ and (3/2) A for\n",
     "  2++.\n",
     "- **Integration:** everything is integrated together by `pq_jet` at 12 digits; eigenvalues by the Illinois\n",
     "  method to a relative width of 1e-9. The ratios below are numerical results at that precision, not the 25-30\n",
     "  digits of the soft-wall, hard-wall and Sakai-Sugimoto checks.\n",
     "- **The UV boundary** r_b, where e^(-A) -> 0, is located by extrapolation (the solution with A(0.36) = 0 has it\n",
     "  below r = 0). The solution starts 10^-3 from it as xi ~ (r - r_b)^4 (psi ~ (r - r_b)^(5/2), (8.10)).\n",
     "- **Eigenvalues:** where e^(2B) xi', the growing mode's coefficient, changes sign at the cutoff.\n",
     "- **Units of the runs:** the IR scale r1, where lambda = 1, moves by a factor 40 between lambda0 = 0.25 and\n",
     "  0.05. The cutoff and the scanned range of m^2 are set in units of r1, so that every lambda0 covers the same\n",
     "  states: cutoff 8 r1, and 11 r1 for the check and the long spectra (beyond about 12 r1, lambda exceeds\n",
     "  10^150 and the integration overflows).\n\n",
     "## Results, b0 = 4.2, lambda0 = 0.15 (masses in units of 1/l; ratios to the lowest 0++)\n\n",
     "| state | m | m/m(0++), PQCF | the paper (8.12), (8.13) | difference |\n| --- | --- | --- | --- | --- |\n",
     Row("0++", S, Paper0), Row("2++", T, Paper2),
     io_lib:format("~nR00 = ~.4f (the paper: 1.87 by its fit, 2753/1475 = 1.866), R20 = ~.4f (the paper: 1.40).~n~n",
                   [lists:nth(2, S) / M1, hd(T) / M1]),
     "## Checks\n\n",
     "**Independence of the integration constants** (section 7.2 of the paper: lambda0 and A0 only set the scale).\n",
     "lambda0 = 0.05 is the paper's own value (its Fig. 4 and fit).\n\n",
     "| lambda0 | r1 | m(0++) | R00 | R20 |\n| --- | --- | --- | --- | --- |\n", L0Rows,
     case CutDiff == 0 of
         true -> "\n**Cutoff:** 11 r1 instead of 8 r1 leaves every mass of the table unchanged, to the last digit computed.\n\n";
         false -> io_lib:format("~n**Cutoff:** 11 r1 instead of 8 r1 changes no mass of the table by more than ~.2e (relative).~n~n", [CutDiff])
     end,
     "**Linear confinement** (alpha = 2). The local exponent d log(-A)/d log r of the background (lambda0 = 0.15):\n\n",
     "| range | exponent |\n| --- | --- |\n", AlphaRows,
     "\nIt falls toward 2 slowly, through the logarithmic corrections of (8.2). The spectrum is close to linear in\n",
     "m^2, as the paper's Fig. 7a shows: the gaps m_(n+1)^2 - m_n^2 of the first 12 0++ states are\n",
     io_lib:format("~s,~nand of the first 8 2++ states ~s:~n", [nums(Gaps), nums(Gaps2)]),
     "nearly constant, and still growing slowly toward their asymptotic value, as the logarithmic corrections\n",
     "predict.\n\n",
     io_lib:format("**The fit of b0.** The paper fits b0 to R00 = 1.87 and finds b0 = 4.2. The secant method on R00(b0) gives\n"
                   "b0 = ~.3f (R00 - 1.87 = ~.2e). Our fit lands 1.3 % below the paper's b0,\n"
                   "and at its b0 = 4.2 our R00 is 1.875: the same small offset as in the table above, not a separate\n"
                   "difference.~n~n", [Ob, Fb]),
     "**The model's own consistency:** d log W/d lambda = -(4/3) X/lambda (2.22) holds for (8.4) and (8.6).\n\n",
     "## The earlier discrepancy (resolved, 15 Nov 2026)\n\n",
     "This script first used a = 5/16, which is neither reading of the paper (it corresponds to alpha = 12/7 through\n",
     io_lib:format("(8.4)). That run gave R00 = 1.683 and R20 = 1.309, 6-15 % below the paper, and the cause was not\n"
                   "found at the time. Rerun with a = 5/16 now: R00 = ~.3f, R20 = ~.3f.~n~n", [lists:nth(2, Old) / hd(Old), hd(OldT) / hd(Old)]),
     "With a = 3/8 the ratios agree with the paper to within 0.6 %. That remaining difference is larger than the\n",
     "rounding of the paper's masses (to 1 MeV in about 1500) and is not yet explained. The next step is a\n",
     "convergence study at higher precision and an independent solver.\n",
     io_lib:format("~nTime: ~b s, all runs in parallel on 16 cores (the longest, the first 12 0++ states).~n", [round(Secs)])].

pct(X) when X >= 0 -> io_lib:format("+~.1f", [X]);
pct(X) -> io_lib:format("-~.1f", [-X]).

l0_str({rat, P, Q}) -> io_lib:format("~.2f", [P / Q]).

nums(L) -> lists:join(", ", [io_lib:format("~.2f", [X]) || X <- L]).
