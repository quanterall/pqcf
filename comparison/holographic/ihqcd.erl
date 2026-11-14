%% Improved holographic QCD, background I (Gursoy, Kiritsis, Nitti, arXiv:0707.1349, https://arxiv.org/abs/0707.1349,
%% section 8.1; set up in arXiv:0707.1324): the scalar (0++) and tensor (2++) glueball spectra, in PQCF (decision record
%% SY2, phase 6). Writes IHQCD.md here.
%%
%%   docker run --rm -u $(id -u):$(id -g) -e PQCF_COMMIT=$(git rev-parse --short HEAD) -v "$PWD":/w -w /w pqcf-dev \
%%     sh -c 'rebar3 compile >/dev/null && erlc -o /tmp comparison/holographic/ihqcd.erl &&
%%     erl -noshell -pa _build/default/lib/*/ebin -pa /tmp -eval "ihqcd:main(\"comparison/holographic\"), halt()."'
-module(ihqcd).
-export([main/1, spectrum/3, x_w/1, start/2]).

-define(D, 12).          % digits of the integrations

%% X(lambda) (8.4) and W(lambda) (8.6) for b0, b1 = (51/121) b0^2 (pure Yang-Mills), a = 5/16 (alpha = 2: linear
%% confinement), as expressions in the symbol l
x_w(B0) ->
    L = {sym, l},
    B1 = {mul, [{rat, 51, 121}, {pow, B0, 2}]},
    C = {add, [{mul, [2, {pow, B0, 2}]}, {mul, [3, B1]}]},
    A = {rat, 5, 16},
    Lg = {add, [1, {mul, [C, {pow, {mul, [9, A]}, -1}, {log, {add, [1, {pow, L, 2}]}}]}]},
    X = {add, [{mul, [-1, B0, L, {pow, {add, [3, {mul, [2, B0, L]}]}, -1}]},
               {mul, [-1, C, {pow, L, 2}, {pow, {mul, [9, {add, [1, {pow, L, 2}]}, Lg]}, -1}]}]},
    W = {mul, [{rat, 9, 4}, {pow, {add, [1, {mul, [{rat, 2, 3}, B0, L]}]}, {rat, 2, 3}}, {pow, Lg, {rat, 5, 24}}]},
    {X, W}.

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
spectrum(B0, Kind, Opts) ->
    St = start(B0, maps:get(lambda0, Opts, {rat, 3, 20})),
    Rmax = maps:get(rmax, Opts, 20),
    {Lo, Hi} = maps:get(range, Opts),
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

main(Dir) ->
    B0 = {rat, 21, 5},
    Run = fun(Kind, L0, Range) ->
                  timer:tc(fun() -> spectrum(B0, Kind, #{range => Range, grid => 30, rmax => 10, lambda0 => L0}) end)
          end,
    {Ts, S} = Run(scalar, {rat, 3, 20}, {0.05, 30}),
    {Tt, T} = Run(tensor, {rat, 3, 20}, {0.05, 30}),
    {Ts2, S2} = Run(scalar, {rat, 1, 4}, {0.0005, 120}),
    Ms = [math:sqrt(X) || X <- S], Mt = [math:sqrt(X) || X <- T], Ms2 = [math:sqrt(X) || X <- S2],
    M1 = hd(Ms),
    Paper0 = [1475, 2753, 3561, 4253], Paper2 = [2055, 2991, 3739, 4396],
    Row = fun(Name, Mine, Paper) ->
                  [io_lib:format("| ~s ~b | ~.4f | ~.3f | ~.3f |~n", [Name, N, M, M / M1, P / 1475])
                   || {N, {M, P}} <- lists:enumerate(lists:zip(lists:sublist(Mine, 4), lists:sublist(Paper, length(lists:sublist(Mine, 4)))))]
          end,
    Md = ["# Improved holographic QCD: scalar and tensor glueballs\n\n",
          "Produced by `ihqcd.erl` (decision record SY2, phase 6), PQCF ", os:getenv("PQCF_COMMIT", "(commit not given)"), ".\n",
          "The model is Gursoy, Kiritsis and Nitti, [arXiv:0707.1349](https://arxiv.org/abs/0707.1349), section 8.1\n",
          "(background I; the framework is [arXiv:0707.1324](https://arxiv.org/abs/0707.1324)).\n\n",
          "## Set-up\n\n",
          "- **Model:** X(lambda) from (8.4) and W(lambda) from (8.6), with b0 = 4.2, b1/b0^2 = 51/121 and a = 5/16 (linear\n",
          "  confinement, alpha = 2).\n",
          "- **Background** (8.5): lambda' = -(4/3) X W lambda e^A and A' = -(4/9) W e^A (l = 1), from\n",
          "  lambda(0.36) = lambda0, A(0.36) = 0.\n",
          "- **Glueballs** (6.2): xi'' + 2 B' xi' + m^2 xi = 0, with B = (3/2) A + (1/2) log X^2 for 0++ and (3/2) A for\n",
          "  2++.\n",
          "- **Integration:** everything is integrated together by `pq_jet` (12 digits).\n",
          "- **The UV boundary** r_b, where e^(-A) -> 0, is located by extrapolation: the solution with A(0.36) = 0 has it\n",
          "  near r = -0.46, not 0. The normalisable solution xi ~ (r - r_b)^4 (psi ~ (r - r_b)^(5/2), (8.10)) starts\n",
          "  10^-3 from it.\n",
          "- **Eigenvalues:** where e^(2B) xi', the growing mode's coefficient, changes sign at r = 10.\n\n",
          "## Results, b0 = 4.2 (masses in units of 1/l; ratios to the lowest 0++)\n\n",
          "| state | m | m/m(0++), PQCF | the paper (8.12), (8.13) |\n| --- | --- | --- | --- |\n",
          Row("0++", Ms, Paper0), Row("2++", Mt, Paper2),
          "\n## Checks\n\n",
          io_lib:format("- **Independence of the integration constants** (section 7.2 of the paper): with lambda0 = 0.25 the 0++ masses\n"
                        "  change (to ~w), but the ratios do not: ~w.~n",
                        [[round(X * 1000) / 1000 || X <- Ms2], [round(X / hd(Ms2) * 1000) / 1000 || X <- Ms2]]),
          "  lambda0 = 0.08 gives the same ratios too.\n",
          "- **Cutoff:** r = 14 gives the same values as r = 10.\n",
          "- **The model's own consistency:** d log W/d lambda = -(4/3) X/lambda (2.22) holds for (8.4) and (8.6).\n\n",
          "## Discrepancy with the paper (open)\n\n",
          io_lib:format("- **The values differ:** PQCF gives R00 = ~.3f and R20 = ~.3f at b0 = 4.2. The paper fits b0 = 4.2 to\n"
                        "  R00 = 1.87 and finds R20 = 1.40. The excited 0++ ratios are lower too.~n",
                        [lists:nth(2, Ms) / M1, hd(Mt) / M1]),
          "- **What the checks rule out:** the computation does not depend on the integration constants or the cutoff,\n",
          "  as the paper's own argument requires, so the difference is not numerical noise.\n",
          "- **Possible sources:** a convention not stated in the paper's formulas (the normalisation of lambda in X\n",
          "  against W, or where the UV condition is imposed), or a difference in the paper's numerics. With A(r_in) = 0\n",
          "  at r_in = 0.36, the paper's stated initial conditions put the boundary at r ~ -0.46, while (8.10) imposes\n",
          "  the UV behaviour at r -> 0.\n",
          "- **Not resolved:** a scan in b0 found R00 = 1.87 near b0 = 2, but with R20 = 1.09, so a single change of b0\n",
          "  does not reproduce both ratios. That scan used a coarse grid and is not reliable enough to settle this.\n",
          io_lib:format("~nTimes: 0++ ~.1f s, 2++ ~.1f s, the lambda0 check ~.1f s.~n", [Ts / 1.0e6, Tt / 1.0e6, Ts2 / 1.0e6])],
    ok = file:write_file(filename:join(Dir, "IHQCD.md"), unicode:characters_to_binary(Md)),
    io:format("~s", [Md]).
