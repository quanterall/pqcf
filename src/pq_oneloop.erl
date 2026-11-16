%% One-loop scalar and tensor integrals (plan "theoretical physics instruments", F3, first part): the one- and
%% two-point functions at arbitrary precision, for real kinematics.
%%
%% Conventions are LoopTools' (Hahn and Perez-Victoria, arXiv:hep-ph/9807565) and Denner's (arXiv:0709.1075): the
%% finite parts, Delta = 2/(4 - d) - gamma_E + log(4 pi) set to 0, the scale mu^2 (default 1),
%%   A0(m^2)           = m^2 (1 - log(m^2 / mu^2)),
%%   B0(p^2, m1^2, m2^2) = - int_0^1 dx log((x m1^2 + (1 - x) m2^2 - x (1 - x) p^2 - i eps) / mu^2),
%% B1, B00, B11 by Passarino-Veltman reduction to A0 and B0, and DB0 = dB0/dp^2. A massless B0 at p^2 = 0 is
%% scaleless and 0 (its UV and IR poles cancel).
%%
%% B0 in closed form (Denner, eq. (4.23)), with s = p^2, a = m1^2, b = m2^2 (all > 0) and
%% lambda = s^2 + a^2 + b^2 - 2 s a - 2 s b - 2 a b:
%%   B0 = 2 - log(m1 m2) + (a - b)/s log(m2/m1) + sqrt(lambda)/s log((a + b - s + sqrt(lambda)) / (2 m1 m2)),
%% the last term real below the pseudo-threshold, -sqrt(-lambda)/s theta between the thresholds
%% (cos theta = (a + b - s) / (2 m1 m2), theta in (0, pi)), and with imaginary part pi sqrt(lambda)/s above
%% (m1 + m2)^2. A dimensionless integral depends on s/mu^2, a/mu^2, b/mu^2 only, so everything is computed at mu^2 = 1.
%%
%% Arguments are exact: integers or fractions {N, D}; mu2 also. Results are {Re, Im} of pq_mp numbers at P bits, or
%% decimal strings (value/3).
-module(pq_oneloop).
-export([a0/2, a0/3, b0/4, b0/5, b1/4, b1/5, b00/4, b00/5, b11/4, b11/5, db0/4, db0/5, c0/7, d0/11, value/3]).

-define(Z, {mp, 0, 0}).

a0(M2, P) -> a0(M2, P, #{}).
a0(M2, P, Opts) ->
    Mu2 = q(maps:get(mu2, Opts, 1)),
    W = P + 32,
    {R, I} = a0_(qdiv(q(M2), Mu2), W),
    round_c({mul(R, mp(Mu2, W), W), mul(I, mp(Mu2, W), W)}, P).

b0(S, A, B, P) -> b0(S, A, B, P, #{}).
b0(S, A, B, P, Opts) -> at_scale(fun b0_/4, S, A, B, P, Opts, 0).

b1(S, A, B, P) -> b1(S, A, B, P, #{}).
b1(S, A, B, P, Opts) -> at_scale(fun b1_/4, S, A, B, P, Opts, 0).

%% B00 has mass dimension 2, DB0 -2: rescaled by mu^2 and mu^-2
b00(S, A, B, P) -> b00(S, A, B, P, #{}).
b00(S, A, B, P, Opts) -> at_scale(fun b00_/4, S, A, B, P, Opts, 1).

b11(S, A, B, P) -> b11(S, A, B, P, #{}).
b11(S, A, B, P, Opts) -> at_scale(fun b11_/4, S, A, B, P, Opts, 0).

db0(S, A, B, P) -> db0(S, A, B, P, #{}).
db0(S, A, B, P, Opts) -> at_scale(fun db0_/4, S, A, B, P, Opts, -1).

%% {Re, Im} as decimal strings to D digits: value(Fun, Args, D), e.g. value(b0, [5, 1, 1], 30)
value(F, Args, D) ->
    P = ceil(D * 3.3219280948873622) + 8,
    {R, I} = apply(?MODULE, F, Args ++ [P]),
    {pq_mp:to_string(R, D), pq_mp:to_string(I, D)}.

at_scale(F, S, A, B, P, Opts, Dim) ->
    Mu2 = q(maps:get(mu2, Opts, 1)),
    W = P + 32,
    {R, I} = F(qdiv(q(S), Mu2), qdiv(q(A), Mu2), qdiv(q(B), Mu2), W),
    Scale = mp(qpow(Mu2, Dim), W),
    round_c({mul(R, Scale, W), mul(I, Scale, W)}, P).

%%% The functions at mu^2 = 1, on exact rationals {N, D}, at W bits

a0_({0, _}, _) -> {?Z, ?Z};
a0_(M2, W) -> {mul(mp(M2, W), sub(one(), log_(M2, W), W), W), ?Z}.

b0_({0, _}, {0, _}, {0, _}, _) -> {?Z, ?Z};
b0_({0, _}, A, B, W) -> {b0_zero(A, B, W), ?Z};
b0_(S, {0, _}, {0, _}, W) ->
    % 2 - log(-s - i eps)
    Re = sub(from(2, W), log_(qabs(S), W), W),
    {Re, case sgn(S) > 0 of true -> pq_mp:pi(W); false -> ?Z end};
b0_(S, {0, _}, B, W) -> b0_(S, B, {0, 1}, W);
b0_(S, A, {0, _}, W) ->
    % 2 - log a + (a - s)/s log((a - s - i eps) / a)
    D = qsub(A, S),
    Base = sub(from(2, W), log_(A, W), W),
    case sgn(D) of
        0 -> {Base, ?Z};
        Sg ->
            C = mp(qdiv(D, S), W),
            L = log_(qdiv(qabs(D), A), W),
            {add(Base, mul(C, L, W), W), case Sg < 0 of true -> neg(mul(C, pq_mp:pi(W), W)); false -> ?Z end}
    end;
b0_(S, A, B, W) ->
    Base = add(sub(from(2, W), half(log_(qmul(A, B), W)), W),
               mul(mp(qdiv(qsub(A, B), S), W), half(log_(qdiv(B, A), W)), W), W),
    {Lr, Li} = lterm(S, A, B, W),                      % sqrt(lambda) L, as {Re, Im}
    Inv = mp(qdiv({1, 1}, S), W),
    {add(Base, mul(Lr, Inv, W), W), mul(Li, Inv, W)}.

%% sqrt(lambda) log((a + b - s + sqrt(lambda)) / (2 m1 m2)), with lambda < 0 continued to -sqrt(-lambda) theta, and the
%% branch above threshold
lterm(S, A, B, W) ->
    Lam = lambda(S, A, B),
    Sig = qsub(qadd(A, B), S),
    case sgn(Lam) of
        0 -> {?Z, ?Z};
        1 ->
            Sq = pq_mp:sqrt(mp(Lam, W), W),
            U = divide(add(mp(qabs(Sig), W), case sgn(Sig) of -1 -> neg(Sq); _ -> Sq end, W),
                       mul(from(2, W), pq_mp:sqrt(mp(qmul(A, B), W), W), W), W),
            % above threshold u < 0: |u| = (s - a - b - sqrt lambda)/(2 m1 m2), log u = log |u| + i pi
            case sgn(Sig) of
                -1 -> {mul(Sq, pq_mp:log(U, W), W), mul(Sq, pq_mp:pi(W), W)};
                _ -> {mul(Sq, pq_mp:log(U, W), W), ?Z}
            end;
        -1 ->
            Sq = pq_mp:sqrt(mp(qneg(Lam), W), W),
            {neg(mul(Sq, theta(Sq, Sig, W), W)), ?Z}
    end.

%% theta in (0, pi) with tan theta = sqrt(-lambda) / (a + b - s)
theta(Sq, Sig, W) ->
    case sgn(Sig) of
        0 -> pq_mp:ldexp(pq_mp:pi(W), -1);
        1 -> pq_mp:atan(divide(Sq, mp(Sig, W), W), W);
        -1 -> sub(pq_mp:pi(W), pq_mp:atan(divide(Sq, mp(qneg(Sig), W), W), W), W)
    end.

%% B0 at p^2 = 0, not both masses 0
b0_zero({0, _}, B, W) -> sub(one(), log_(B, W), W);
b0_zero(A, {0, _}, W) -> sub(one(), log_(A, W), W);
b0_zero(A, A, W) -> neg(log_(A, W));
b0_zero(A, B, W) ->
    % 1 - (a log a - b log b) / (a - b)
    sub(one(), divide(sub(mul(mp(A, W), log_(A, W), W), mul(mp(B, W), log_(B, W), W), W), mp(qsub(A, B), W), W), W).

%% B1 = (A0(a) - A0(b) - (s - b + a) B0) / (2 s)
b1_({0, _}, _, _, _) -> error({pq_oneloop, zero_momentum});
b1_(S, A, B, W) ->
    {A0a, _} = a0_(A, W), {A0b, _} = a0_(B, W),
    F = mp(qadd(qsub(S, B), A), W),
    {R0, I0} = b0_(S, A, B, W),
    K = mp(qdiv({1, 1}, qmul({2, 1}, S)), W),
    {mul(sub(sub(A0a, A0b, W), mul(F, R0, W), W), K, W), neg(mul(mul(F, I0, W), K, W))}.

%% B00 = (A0(b) + 2 a B0 + (s - b + a) B1 + a + b - s/3) / 6
b00_(S, A, B, W) ->
    {A0b, _} = a0_(B, W),
    {R0, I0} = b0_(S, A, B, W),
    {R1, I1} = case S of {0, _} -> b1_zero(A, B, W); _ -> b1_(S, A, B, W) end,
    F = mp(qadd(qsub(S, B), A), W),
    TwoA = mp(qmul({2, 1}, A), W),
    C = mp(qsub(qadd(A, B), qdiv(S, {3, 1})), W),
    Six = mp({1, 6}, W),
    {mul(add(add(add(A0b, mul(TwoA, R0, W), W), mul(F, R1, W), W), C, W), Six, W),
     mul(add(mul(TwoA, I0, W), mul(F, I1, W), W), Six, W)}.

%% B1 at p^2 = 0 (for B00 there), from its integral: B1 = int_0^1 x log(x b + (1 - x) a) dx (Denner's (4.7))
b1_zero(A, B, W) -> {int_x_log(B, A, W), ?Z}.

%% int_0^1 x log(x a + (1 - x) b) dx, a, b >= 0, not both 0
int_x_log({0, _}, B, W) -> sub(half(log_(B, W)), mp({3, 4}, W), W);
int_x_log(A, {0, _}, W) -> sub(half(log_(A, W)), mp({1, 4}, W), W);
int_x_log(A, A, W) -> half(log_(A, W));
int_x_log(A, B, W) ->
    % with c = a - b: int_0^1 x log(b + c x) dx
    %   = [(x^2/2 - b^2/(2 c^2)) log(b + c x) - x^2/4 + b x/(2 c)]_0^1
    C = qsub(A, B),
    T1 = mul(mp(qsub({1, 2}, qdiv(qmul(B, B), qmul({2, 1}, qmul(C, C)))), W), log_(A, W), W),
    T0 = mul(mp(qneg(qdiv(qmul(B, B), qmul({2, 1}, qmul(C, C)))), W), log_(B, W), W),
    sub(add(sub(T1, T0, W), mp(qdiv(B, qmul({2, 1}, C)), W), W), mp({1, 4}, W), W).

%% B11 = (A0(b) - 2 (s - b + a) B1 - a B0 - (a + b - s/3)/2) / (3 s)
b11_({0, _}, _, _, _) -> error({pq_oneloop, zero_momentum});
b11_(S, A, B, W) ->
    {A0b, _} = a0_(B, W),
    {R0, I0} = b0_(S, A, B, W),
    {R1, I1} = b1_(S, A, B, W),
    F2 = mp(qmul({2, 1}, qadd(qsub(S, B), A)), W),
    Am = mp(A, W),
    C = mp(qdiv(qsub(qadd(A, B), qdiv(S, {3, 1})), {2, 1}), W),
    K = mp(qdiv({1, 1}, qmul({3, 1}, S)), W),
    {mul(sub(sub(sub(A0b, mul(F2, R1, W), W), mul(Am, R0, W), W), C, W), K, W),
     mul(neg(add(mul(F2, I1, W), mul(Am, I0, W), W)), K, W)}.

%% dB0/ds = -(a - b)/s^2 log(m2/m1) - 1/s + (s (a + b) - (a - b)^2) / s^2 * L / sqrt(lambda)
db0_({0, _}, _, _, _) -> error({pq_oneloop, zero_momentum});
db0_(S, {0, _}, {0, _}, W) -> {neg(mp(qdiv({1, 1}, S), W)), ?Z};
db0_(S, {0, _}, B, W) -> db0_(S, B, {0, 1}, W);
db0_(S, A, {0, _}, W) ->
    % -a/s^2 log((a - s - i eps) / a) - 1/s; on shell (s = a) infrared divergent
    D = qsub(A, S),
    sgn(D) =:= 0 andalso error({pq_oneloop, infrared_divergent}),
    C = mp(qneg(qdiv(A, qmul(S, S))), W),
    L = log_(qdiv(qabs(D), A), W),
    {sub(mul(C, L, W), mp(qdiv({1, 1}, S), W), W), case sgn(D) < 0 of true -> neg(mul(C, pq_mp:pi(W), W)); false -> ?Z end};
db0_(S, A, B, W) ->
    S2 = qmul(S, S),
    T1 = mul(mp(qneg(qdiv(qsub(A, B), S2)), W), half(log_(qdiv(B, A), W)), W),
    T2 = neg(mp(qdiv({1, 1}, S), W)),
    K = mp(qdiv(qsub(qmul(S, qadd(A, B)), qmul(qsub(A, B), qsub(A, B))), S2), W),
    {Lr, Li} = l_over_sqrt(S, A, B, W),
    {add(add(T1, T2, W), mul(K, Lr, W), W), mul(K, Li, W)}.

%% L / sqrt(lambda): log u / sqrt(lambda), theta / sqrt(-lambda) between the thresholds, 1/(a + b - s) at the
%% pseudo-threshold; at the threshold (m1 + m2)^2 dB0 diverges
l_over_sqrt(S, A, B, W) ->
    Lam = lambda(S, A, B),
    Sig = qsub(qadd(A, B), S),
    case sgn(Lam) of
        0 when element(1, Sig) > 0 -> {mp(qdiv({1, 1}, Sig), W), ?Z};
        0 -> error({pq_oneloop, threshold});
        1 ->
            {Lr, Li} = lterm(S, A, B, W),
            L2 = mp(Lam, W),
            {divide(Lr, L2, W), divide(Li, L2, W)};
        -1 ->
            Sq = pq_mp:sqrt(mp(qneg(Lam), W), W),
            {divide(theta(Sq, Sig, W), Sq, W), ?Z}
    end.

lambda(S, A, B) ->
    qsub(qadd(qadd(qmul(S, S), qmul(A, A)), qmul(B, B)),
         qmul({2, 1}, qadd(qadd(qmul(S, A), qmul(S, B)), qmul(A, B)))).

%%% The three-point function
%%
%% C0(p1^2, p2^2, p3^2; m1^2, m2^2, m3^2) in LoopTools' order (p3 = p1 + p2; propagators q^2 - m1^2,
%% (q + p1)^2 - m2^2, (q + p1 + p2)^2 - m3^2), UV finite and independent of mu:
%%   C0 = int_0^1 h(x) dx,  h(x) = - int_0^x dy 1 / (D(x, y) - i eps),
%%   D = a x^2 + b y^2 + c x y + d x + e y + f,
%%   a = p1^2, b = p2^2, c = p3^2 - p1^2 - p2^2, d = m2^2 - m1^2 - p1^2, e = m3^2 - m2^2 + p1^2 - p3^2, f = m1^2.
%% The y integral in closed form: D = b (y - r1)(y - r2) as a quadratic in y; real roots get the prescription of
%% D - i eps, r_k + i eps sign(D'(r_k)), so log(u - i eps s) = log|u| - i pi s [u < 0]; complex roots need none (the
%% logs stay on one sheet along y in [0, x]); b = 0 is linear. The x integral by pq_quad, split where a root crosses
%% an end of [0, x] (D(x, 0) = 0, D(x, x) = 0) and where the discriminant in y vanishes: there h is singular
%% (logarithms, inverse square roots), elsewhere analytic. No change of variables into the complex plane, so every
%% region (lambda < 0 included, where 't Hooft and Veltman's real substitution does not exist) is one case.
%% All p^2 = 0: the closed form -phi[m1^2, m2^2, m3^2], the second divided difference of phi(z) = z log z
%% (Hermite-Genocchi). Refused: infrared-divergent configurations, soft (a massless propagator between two on-shell
%% lines) and collinear (a massless external line between two massless propagators).
c0(P1, P2, P3, M1, M2, M3, P) ->
    [Q1, Q2, Q3, Ma, Mb, Mc] = [q(X) || X <- [P1, P2, P3, M1, M2, M3]],
    ir_check(Q1, Q2, Q3, Ma, Mb, Mc),
    W = pq_quad:bits(P + 16),
    R = case [Q1, Q2, Q3] of
            [{0, _}, {0, _}, {0, _}] -> {neg(divided2(Ma, Mb, Mc, W)), ?Z};
            _ -> c0_int(Q1, Q2, Q3, Ma, Mb, Mc, P, W)
        end,
    round_c(R, P).

%% soft: a massless propagator whose two neighbouring lines are on shell; collinear: a massless external line
%% (p^2 = 0) between two massless propagators
ir_check(P1, P2, P3, M1, M2, M3) ->
    Soft = [I || {I, M, Pa, Ma, Pb, Mb} <- [{1, M1, P1, M2, P3, M3}, {2, M2, P1, M1, P2, M3}, {3, M3, P3, M1, P2, M2}],
                 sgn(M) =:= 0, Pa =:= Ma, Pb =:= Mb],
    Coll = [Pij || {Pij, Mi, Mj} <- [{P1, M1, M2}, {P2, M2, M3}, {P3, M1, M3}], sgn(Pij) =:= 0, sgn(Mi) =:= 0, sgn(Mj) =:= 0],
    Soft =:= [] andalso Coll =:= [] orelse error({pq_oneloop, infrared_divergent}).

c0_int(P1, P2, P3, M1, M2, M3, P, W) ->
    A = P1, B = P2, C = qsub(qsub(P3, P1), P2),
    D = qsub(qsub(M2, M1), P1), E = qadd(qsub(qsub(M3, M2), P3), P1), F = M1,
    % where h is singular: D(x, 0) = 0, D(x, x) = 0, the discriminant (c x + e)^2 - 4 b (a x^2 + d x + f) = 0,
    % and for b = 0 the zero of the y coefficient c x + e
    Edges = [{A, D, F}, {qadd(qadd(A, B), C), qadd(D, E), F},
             {qsub(qmul(C, C), qmul({4, 1}, qmul(A, B))), qsub(qmul({2, 1}, qmul(C, E)), qmul({4, 1}, qmul(B, D))),
              qsub(qmul(E, E), qmul({4, 1}, qmul(B, F)))}],
    Breaks = lists:append([real_roots(Qd, W) || Qd <- Edges]) ++
             case sgn(B) of 0 -> real_roots({{0, 1}, C, E}, W); _ -> [] end,
    % D on the two edges, y = 0 and y = x, through their zeros: accurate next to them, where the roots in y reach 0 or x
    G0 = {A, D, F}, G1 = {qadd(qadd(A, B), C), qadd(D, E), F},
    Ctx = #{g0 => {G0, real_roots(G0, W)}, g1 => {G1, real_roots(G1, W)}},
    H = fun(X) -> cneg(inner(X, A, B, C, D, E, F, Ctx, W)) end,
    {R, I} = pq_quad:integrate(H, 0, 1, Breaks, P + 16),
    {pq_mp:from(R, W), pq_mp:from(I, W)}.

%% int_0^x dy / (b y^2 + beta y + gamma - i eps), beta = c x + e, gamma = a x^2 + d x + f, x real in (0, 1).
%% The logarithms need -r_k and x - r_k; their products are gamma/b and D(x, x)/b, computed through the edges' zeros,
%% so the smaller of each pair is the product over the larger: no cancellation where a root meets 0 or x.
inner(X, _A, B, C, _D, E, _F, #{g0 := {G0, Z0}, g1 := {G1, Z1}}, W) ->
    Beta = add(mul(mp(C, W), X, W), mp(E, W), W),
    Gamma = qeval_r(G0, Z0, X, W),                       % D(x, 0)
    Dxx = qeval_r(G1, Z1, X, W),                         % D(x, x)
    inner_q(X, B, Beta, Gamma, Dxx, undefined, W).

%% int_0^X dz / (b z^2 + beta z + gamma - i eps): b exact, beta, gamma = D(0) and Dxx = D(X) given (accurate where
%% they vanish), and the discriminant beta^2 - 4 b gamma if the caller has it more accurately than by the subtraction
inner_q(X, B, Beta, Gamma, Dxx, Disc0, W) ->
    case sgn(B) of
        0 ->
            % b = 0: D(x, y) = beta y + gamma, and D(x, x) = beta x + gamma
            case is_zero(Beta) of
                true -> cdiv({X, ?Z}, {Gamma, ?Z}, W);
                false -> cscale(csub(logr(Dxx, W), logr(Gamma, W), W), divide(one(), Beta, W), W)
            end;
        _ ->
            Bm = mp(B, W),
            Disc = case Disc0 of
                       undefined -> sub(mul(Beta, Beta, W), mul(from(4, W), mul(Bm, Gamma, W), W), W);
                       _ -> Disc0
                   end,
            case pq_mp:cmp(Disc, ?Z) of
                1 ->
                    Sq = pq_mp:sqrt(Disc, W),
                    % the stable pair: q = -(beta + sign(beta) sqrt(disc))/2, roots q/b and gamma/q
                    Qv = neg(pq_mp:ldexp(case pq_mp:cmp(Beta, ?Z) < 0 of true -> sub(Beta, Sq, W); false -> add(Beta, Sq, W) end, -1)),
                    {R1, R2} = case is_zero(Qv) of
                                   true -> {?Z, ?Z};
                                   false -> {divide(Qv, Bm, W), divide(Gamma, Qv, W)}
                               end,
                    D1 = mul(Bm, sub(R1, R2, W), W),
                    S1 = pq_mp:cmp(D1, ?Z),
                    {N1, N2} = pair(neg(R1), neg(R2), divide(Gamma, Bm, W), W),          % -r1, -r2
                    {X1, X2} = pair(sub(X, R1, W), sub(X, R2, W), divide(Dxx, Bm, W), W), % x - r1, x - r2
                    L1 = csub(logs(X1, S1, W), logs(N1, S1, W), W),
                    L2 = csub(logs(X2, -S1, W), logs(N2, -S1, W), W),
                    cscale(csub(L1, L2, W), divide(one(), D1, W), W);
                -1 ->
                    % complex roots r, conj(r): 1/D = (1/(y - r) - 1/(y - conj r)) / (b (r - conj r))
                    Inv2b = divide(one(), mul(from(2, W), Bm, W), W),
                    Re = neg(mul(Beta, Inv2b, W)),
                    Im = mul(pq_mp:sqrt(neg(Disc), W), rabs(Inv2b), W),
                    Rr = {Re, Im}, Rc = {Re, neg(Im)},
                    La = csub(clog(csub({X, ?Z}, Rr, W), W), clog(cneg(Rr), W), W),
                    Lb = csub(clog(csub({X, ?Z}, Rc, W), W), clog(cneg(Rc), W), W),
                    cdiv(csub(La, Lb, W), {?Z, mul(Bm, pq_mp:ldexp(Im, 1), W)}, W);
                0 ->
                    % a double root r: (1/(-r) - 1/(x - r)) / b
                    R = neg(divide(Beta, mul(from(2, W), Bm, W), W)),
                    {divide(sub(divide(one(), neg(R), W), divide(one(), sub(X, R, W), W), W), Bm, W), ?Z}
            end
    end.

%% a quadratic {t^2, t, 1} (exact) at a real t, through its real zeros when it has any (accurate next to them)
qeval_r({{0, _}, {0, _}, C}, _, _, W) -> mp(C, W);
qeval_r({{0, _}, B, _}, [R], T, W) -> mul(mp(B, W), sub(T, R, W), W);
qeval_r({A, _, _}, [R1, R2], T, W) -> mul(mp(A, W), mul(sub(T, R1, W), sub(T, R2, W), W), W);
qeval_r({A, _, _}, [R], T, W) -> mul(mp(A, W), mul(sub(T, R, W), sub(T, R, W), W), W);
qeval_r({A, B, C}, [], T, W) -> add(mul(add(mul(mp(A, W), T, W), mp(B, W), W), T, W), mp(C, W), W).

%% two numbers whose product Prod is known: the larger as given, the smaller as Prod over the larger
pair(U1, U2, Prod, W) ->
    case pq_mp:cmp(rabs(U1), rabs(U2)) >= 0 of
        true -> {U1, small(Prod, U1, W)};
        false -> {small(Prod, U2, W), U2}
    end.
small(_, {mp, 0, _}, _) -> ?Z;
small(Prod, Big, W) -> divide(Prod, Big, W).

%% log(u - i eps s) of a real u, s = +-1
logs(U, S, W) ->
    case pq_mp:cmp(U, ?Z) < 0 of
        true -> {pq_mp:log(neg(U), W), case S > 0 of true -> neg(pq_mp:pi(W)); false -> pq_mp:pi(W) end};
        false -> {pq_mp:log(U, W), ?Z}
    end.


is_zero({mp, 0, _}) -> true;
is_zero(_) -> false.

%% the real zeros of a quadratic {t^2, t, 1} with exact coefficients (pq_quad keeps those inside (0, 1))
real_roots({{0, _}, {0, _}, _}, _) -> [];
real_roots({{0, _}, B, C}, W) -> [mp(qneg(qdiv(C, B)), W)];
real_roots({A, B, C}, W) ->
    Disc = qsub(qmul(B, B), qmul({4, 1}, qmul(A, C))),
    case sgn(Disc) of
        -1 -> [];
        0 -> [mp(qneg(qdiv(B, qmul({2, 1}, A))), W)];
        1 ->
            Sq = pq_mp:sqrt(mp(Disc, W), W),
            Inv = mp(qdiv({1, 1}, qmul({2, 1}, A)), W),
            [mul(sub(neg(mp(B, W)), Sq, W), Inv, W), mul(add(neg(mp(B, W)), Sq, W), Inv, W)]
    end.

%% log(v - i eps) of a real v, and of a complex value (real ones get the -i eps)
logr(V, W) ->
    case pq_mp:cmp(V, ?Z) < 0 of
        true -> {pq_mp:log(neg(V), W), neg(pq_mp:pi(W))};
        false -> {pq_mp:log(V, W), ?Z}
    end.

%% the second divided difference of z log z at three points (masses squared, exact)
divided2(A, B, C, W) ->
    case lists:sort([A, B, C]) of
        [X, X, X] -> sgn(X) =:= 0 andalso error({pq_oneloop, infrared_divergent}),
                     mp(qdiv({1, 1}, qmul({2, 1}, X)), W);
        [X, X, Y] -> confluent2(X, Y, W);
        [X, Y, Y] -> confluent2(Y, X, W);
        [X, Y, Z] ->
            Phi = fun(V) -> zlogz(V, W) end,
            add(add(divide(Phi(X), mp(qmul(qsub(X, Y), qsub(X, Z)), W), W),
                    divide(Phi(Y), mp(qmul(qsub(Y, X), qsub(Y, Z)), W), W), W),
                divide(Phi(Z), mp(qmul(qsub(Z, X), qsub(Z, Y)), W), W), W)
    end.

%% phi[x, x, y] = (phi[x, y] - phi'(x)) / (y - x), phi'(x) = log x + 1
confluent2(X, Y, W) ->
    sgn(X) =:= 0 andalso error({pq_oneloop, infrared_divergent}),
    Dxy = divide(sub(zlogz(Y, W), zlogz(X, W), W), mp(qsub(Y, X), W), W),
    divide(sub(Dxy, add(log_(X, W), one(), W), W), mp(qsub(Y, X), W), W).

zlogz({0, _}, _) -> ?Z;
zlogz(X, W) -> mul(mp(X, W), log_(X, W), W).

%%% The four-point function
%%
%% D0(p1^2, p2^2, p3^2, p4^2, s12, s23; m1^2, m2^2, m3^2, m4^2) in LoopTools' order (propagators q^2 - m1^2,
%% (q + p1)^2 - m2^2, (q + p1 + p2)^2 - m3^2, (q - p4)^2 - m4^2; s12 = (p1 + p2)^2, s23 = (p2 + p3)^2):
%%   D0 = int over the 3-simplex of 1 / (x^T S x - i eps)^2,  S_ij = (m_i^2 + m_j^2 - p_ij^2) / 2,
%% with p_12 = p1, p_23 = p2, p_34 = p3, p_14 = p4, p_13 = s12, p_24 = s23. One analytic integration of 1/F^2 leaves
%% singularities that are not integrable, so D0 is reduced by the dimension-shift relation (Bern, Dixon and Kosower,
%% hep-ph/9306240; Binoth, Guillet and Heinrich, hep-ph/9911342):
%%   D0 = -(1/2) sum_i b_i C0^(i) - (1/2) B E,   b = S^-1 (1, 1, 1, 1),  B = sum_i b_i,
%% C0^(i) the triangle without propagator i and E = int over the simplex of 1/(F - i eps) (the six-dimensional box,
%% finite). E as C0: x1 = 1 - x, x2 = x - y, x3 = y - z, x4 = z, the z integral in closed form (inner_q), y and x by
%% pq_quad, split where F vanishes on the faces z = 0 and z = y or the discriminant in z vanishes (quadratics in y),
%% and in x where those curves turn, meet y = 0 or y = x, or meet each other (resultants, pq_rroots).
%% Refused: det S = 0, where the reduction does not apply (the leading Landau singularity, but also degenerate cases
%% where D0 is finite, such as zero momenta with equal masses), and divergent triangles with b_i =/= 0.
d0(P1, P2, P3, P4, S12, S23, M1, M2, M3, M4, P) ->
    Pij = #{{1, 2} => q(P1), {2, 3} => q(P2), {3, 4} => q(P3), {1, 4} => q(P4), {1, 3} => q(S12), {2, 4} => q(S23)},
    Ms = [q(M) || M <- [M1, M2, M3, M4]],
    Sm = [[qdiv(qsub(qadd(lists:nth(I, Ms), lists:nth(J, Ms)), pij(Pij, I, J)), {2, 1}) || J <- lists:seq(1, 4)]
          || I <- lists:seq(1, 4)],
    Bv = case solve(Sm, [{1, 1} || _ <- Sm]) of
             singular -> error({pq_oneloop, degenerate_kinematics});
             V -> V
         end,
    Pw = P + 16,
    W = pq_quad:bits(Pw),
    Tri = [begin
               [J, K, L] = [N || N <- lists:seq(1, 4), N =/= I],
               {Rt, It} = c0(pij(Pij, J, K), pij(Pij, K, L), pij(Pij, J, L),
                             lists:nth(J, Ms), lists:nth(K, Ms), lists:nth(L, Ms), Pw),
               cscale({pq_mp:from(Rt, W), pq_mp:from(It, W)}, mp(Bi, W), W)
           end || {I, Bi} <- lists:enumerate(Bv), sgn(Bi) =/= 0],
    Bsum = lists:foldl(fun qadd/2, {0, 1}, Bv),
    {Er, Ei} = box6(Sm, Pw),
    Ebox = cscale({pq_mp:from(Er, W), pq_mp:from(Ei, W)}, mp(Bsum, W), W),
    Sum = lists:foldl(fun(T, Acc) -> cadd2(Acc, T, W) end, Ebox, Tri),
    round_c(cscale(Sum, mp({-1, 2}, W), W), P).

pij(_, I, I) -> {0, 1};
pij(Pij, I, J) -> maps:get({min(I, J), max(I, J)}, Pij).

%% E = int over the 3-simplex of 1 / (x^T S x - i eps)
box6(Sm, P) ->
    Dot = fun(X, Y) -> lists:foldl(fun qadd/2, {0, 1},
                                   [qmul(qmul(Xi, lists:nth(J, lists:nth(I, Sm))), Yj)
                                    || {I, Xi} <- lists:enumerate(X), {J, Yj} <- lists:enumerate(Y)]) end,
    Z = {0, 1}, O = {1, 1}, Mo = {-1, 1},
    V = [Z, Z, Mo, O], U0 = [O, Z, Z, Z], U1 = [Mo, O, Z, Z], U2 = [Z, Mo, O, Z],
    Two = fun(X) -> qmul({2, 1}, X) end,
    Bz = Dot(V, V),
    {B0, Bx, By} = {Two(Dot(U0, V)), Two(Dot(U1, V)), Two(Dot(U2, V))},
    {C00, Cx, Cy} = {Dot(U0, U0), Two(Dot(U0, U1)), Two(Dot(U0, U2))},
    {Cxx, Cxy, Cyy} = {Dot(U1, U1), Two(Dot(U1, U2)), Dot(U2, U2)},
    F4 = {4, 1},
    % quadratics in y, coefficients polynomials in x: F on z = 0, F on z = y, the discriminant in z
    G0 = {[Cyy], [Cy, Cxy], [C00, Cx, Cxx]},
    G1 = {[qadd(qadd(Cyy, By), Bz)], [qadd(Cy, B0), qadd(Cxy, Bx)], [C00, Cx, Cxx]},
    Dz = {[qsub(qmul(By, By), qmul(F4, qmul(Bz, Cyy)))],
          [qsub(qmul({2, 1}, qmul(By, B0)), qmul(F4, qmul(Bz, Cy))), qsub(qmul({2, 1}, qmul(By, Bx)), qmul(F4, qmul(Bz, Cxy)))],
          [qsub(qmul(B0, B0), qmul(F4, qmul(Bz, C00))), qsub(qmul({2, 1}, qmul(B0, Bx)), qmul(F4, qmul(Bz, Cx))),
           qsub(qmul(Bx, Bx), qmul(F4, qmul(Bz, Cxx)))]},
    Curves = [G0, G1, Dz] ++ case sgn(Bz) of 0 -> [{[Z], [By], [B0, Bx]}]; _ -> [] end,
    Wn = pq_quad:bits(P, 3),
    OuterBreaks = lists:append([pq_rroots:roots(Pl, 0, 1, Wn) || Pl <- outer_polys(Curves)]),
    % the discriminant in z, Dz(x, y) = a (y - re)^2 - disc_y(x) / (4 a), with disc_y(x) (exact, degree 2) through its
    % zeros: accurate where it vanishes, where the subtraction beta^2 - 4 b gamma would lose everything
    {[Ad], Bd, Cd} = Dz,
    DiscY = pq_rroots:sub(pq_rroots:mul(Bd, Bd), pq_rroots:scale(pq_rroots:mul([Ad], Cd), 4)),
    DiscYq = list_to_tuple(lists:reverse(lists:sublist(DiscY ++ [Z, Z, Z], 3))),
    DiscYr = real_roots(DiscYq, Wn),
    G = fun(X) ->
                Ev = fun({A, Bp, Cp}) -> {peval(A, X, Wn), peval(Bp, X, Wn), peval(Cp, X, Wn)} end,
                [Q0, Q1 | _] = Qs = [Ev(Cv) || Cv <- Curves],
                R0 = mp_roots(Q0, Wn), R1 = mp_roots(Q1, Wn),
                % the real roots, and the real parts of complex ones close to the real line (near-singularities)
                Near = pq_mp:ldexp(X, -4),
                Breaks = lists:append([mp_roots(Qv, Wn) ++ near_roots(Qv, Near, Wn) || Qv <- Qs]),
                BetaX = add(mp(B0, Wn), mul(mp(Bx, Wn), X, Wn), Wn),
                DzAt = dz_at(Ad, Bd, DiscYq, DiscYr, X, Wn),
                H = fun(Y) ->
                            Beta = add(BetaX, mul(mp(By, Wn), Y, Wn), Wn),
                            inner_q(Y, Bz, Beta, mp_qeval(Q0, R0, Y, Wn), mp_qeval(Q1, R1, Y, Wn), DzAt(Y), Wn)
                    end,
                pq_quad:integrate(H, 0, X, Breaks, P, #{on_no_convergence => last, max_level => 11, ends => 3})
        end,
    pq_quad:integrate(G, 0, 1, OuterBreaks, P, #{parallel => true, ends => 3}).

%% Dz(x, y) at fixed x as a function of y, without the cancellation of beta^2 - 4 b gamma: a (y - re)^2 - d / (4 a),
%% re = -b_d(x) / (2 a), d = disc_y(x) through its zeros (through its real roots in y when d > 0); undefined (the
%% subtraction) when a = 0
dz_at({0, _}, _, _, _, _, _) -> fun(_) -> undefined end;
dz_at(Ad, Bd, DiscYq, DiscYr, X, W) ->
    A = mp(Ad, W),
    Re = neg(divide(peval(Bd, X, W), pq_mp:ldexp(A, 1), W)),
    Dv = qeval_r(DiscYq, DiscYr, X, W),
    case pq_mp:cmp(Dv, ?Z) > 0 of
        true ->
            % real roots re -+ sqrt(d) / (2 |a|): through them
            H = divide(pq_mp:sqrt(Dv, W), rabs(pq_mp:ldexp(A, 1)), W),
            {R1, R2} = {sub(Re, H, W), add(Re, H, W)},
            fun(Y) -> mul(A, mul(sub(Y, R1, W), sub(Y, R2, W), W), W) end;
        false ->
            K = divide(Dv, pq_mp:ldexp(A, 2), W),
            fun(Y) -> T = sub(Y, Re, W), sub(mul(A, mul(T, T, W), W), K, W) end
    end.

%% the polynomials in x whose roots are the outer breakpoints
outer_polys(Curves) ->
    Each = lists:append([[pq_rroots:sub(pq_rroots:mul(B, B), pq_rroots:scale(pq_rroots:mul(A, C), 4)),   % double root
                          C,                                                                            % root at y = 0
                          pq_rroots:add(pq_rroots:add(pq_rroots:mul([0, 0, 1], A), pq_rroots:mul([0, 1], B)), C)]  % at y = x
                         ++ case lists:all(fun(K) -> sgn(K) =:= 0 end, A) of true -> [B]; false -> [] end
                         || {A, B, C} <- Curves]),
    Pairs = [resultant(Q1, Q2) || {I, Q1} <- lists:enumerate(Curves), {J, Q2} <- lists:enumerate(Curves), I < J],
    [Pl || Pl <- Each ++ Pairs, Pl =/= []].

%% the resultant in y of a1 y^2 + b1 y + c1 and a2 y^2 + b2 y + c2
resultant({A1, B1, C1}, {A2, B2, C2}) ->
    M = fun pq_rroots:mul/2, S = fun pq_rroots:sub/2,
    S(M(S(M(A1, C2), M(A2, C1)), S(M(A1, C2), M(A2, C1))), M(S(M(A1, B2), M(A2, B1)), S(M(B1, C2), M(B2, C1)))).

peval(Poly, X, W) -> lists:foldr(fun(C, Acc) -> add(mul(Acc, X, W), mp(C, W), W) end, ?Z, Poly).

%% the real roots of a y^2 + b y + c (pq_mp coefficients), the stable pair
mp_roots({A, B, C}, W) ->
    case is_zero(A) of
        true -> case is_zero(B) of true -> []; false -> [neg(divide(C, B, W))] end;
        false ->
            Disc = sub(mul(B, B, W), mul(from(4, W), mul(A, C, W), W), W),
            case pq_mp:cmp(Disc, ?Z) of
                -1 -> [];
                _ ->
                    Sq = pq_mp:sqrt(Disc, W),
                    Qv = neg(pq_mp:ldexp(case pq_mp:cmp(B, ?Z) < 0 of true -> sub(B, Sq, W); false -> add(B, Sq, W) end, -1)),
                    case is_zero(Qv) of
                        true -> [?Z];
                        false -> [divide(Qv, A, W), divide(C, Qv, W)]
                    end
            end
    end.

%% the real part of a complex pair of roots of a y^2 + b y + c whose imaginary part is below Near
near_roots({A, B, C}, Near, W) ->
    case is_zero(A) of
        true -> [];
        false ->
            Disc = sub(mul(B, B, W), mul(from(4, W), mul(A, C, W), W), W),
            case pq_mp:cmp(Disc, ?Z) < 0 of
                false -> [];
                true ->
                    Im = divide(pq_mp:sqrt(neg(Disc), W), rabs(pq_mp:ldexp(A, 1)), W),
                    case pq_mp:cmp(Im, Near) < 0 of
                        true -> [neg(divide(B, pq_mp:ldexp(A, 1), W))];
                        false -> []
                    end
            end
    end.

%% a y^2 + b y + c at y through its real roots
mp_qeval({A, B, C}, Roots, Y, W) ->
    case {is_zero(A), Roots} of
        {true, [R]} -> mul(B, sub(Y, R, W), W);
        {true, []} -> C;
        {false, [R1, R2]} -> mul(A, mul(sub(Y, R1, W), sub(Y, R2, W), W), W);
        {false, _} -> add(mul(add(mul(A, Y, W), B, W), Y, W), C, W)
    end.

%% S x = rhs over the rationals (Gauss-Jordan); singular if det S = 0
solve(Sm, Rhs) ->
    Rows = [Row ++ [R] || {Row, R} <- lists:zip(Sm, Rhs)],
    gj(Rows, 1, length(Sm)).

gj(Rows, K, N) when K > N -> [lists:last(R) || R <- Rows];
gj(Rows, K, N) ->
    case [R || R <- lists:nthtail(K - 1, Rows), sgn(lists:nth(K, R)) =/= 0] of
        [] -> singular;
        [Piv | _] ->
            Others = lists:sublist(Rows, K - 1) ++ (lists:nthtail(K - 1, Rows) -- [Piv]),
            Pn = [qdiv(X, lists:nth(K, Piv)) || X <- Piv],
            Elim = [[qsub(X, qmul(lists:nth(K, R), Y)) || {X, Y} <- lists:zip(R, Pn)] || R <- Others],
            {Before, After} = lists:split(K - 1, Elim),
            gj(Before ++ [Pn] ++ After, K + 1, N)
    end.

cadd2({A, B}, {C, D}, W) -> {add(A, C, W), add(B, D, W)}.

%%% Complex numbers {Re, Im} of pq_mp numbers

csub({A, B}, {C, D}, W) -> {sub(A, C, W), sub(B, D, W)}.
cneg({A, B}) -> {neg(A), neg(B)}.
cscale({A, B}, S, W) -> {mul(A, S, W), mul(B, S, W)}.
cdiv({A, B}, {C, D}, W) ->
    N = add(mul(C, C, W), mul(D, D, W), W),
    {divide(add(mul(A, C, W), mul(B, D, W), W), N, W), divide(sub(mul(B, C, W), mul(A, D, W), W), N, W)}.
clog({A, B}, W) ->
    {pq_mp:ldexp(pq_mp:log(add(mul(A, A, W), mul(B, B, W), W), W), -1), atan2(B, A, W)}.
rabs({mp, M, E}) -> {mp, abs(M), E}.

atan2(Y, X, W) ->
    case {pq_mp:cmp(X, ?Z), pq_mp:cmp(Y, ?Z)} of
        {1, _} -> pq_mp:atan(divide(Y, X, W), W);
        {-1, S} when S >= 0 -> add(pq_mp:atan(divide(Y, X, W), W), pq_mp:pi(W), W);
        {-1, _} -> sub(pq_mp:atan(divide(Y, X, W), W), pq_mp:pi(W), W);
        {0, 1} -> pq_mp:ldexp(pq_mp:pi(W), -1);
        {0, -1} -> neg(pq_mp:ldexp(pq_mp:pi(W), -1));
        {0, 0} -> error({pq_oneloop, log_of_zero})
    end.

%%% Exact rationals {N, D}, D > 0, reduced

q(N) when is_integer(N) -> {N, 1};
q({N, D}) when is_integer(N), is_integer(D), D =/= 0 -> norm_q(N, D);
q(X) -> error({pq_oneloop, {not_exact, X}}).

norm_q(N, D) when D < 0 -> norm_q(-N, -D);
norm_q(N, D) -> G = max(1, gcd(abs(N), D)), {N div G, D div G}.

gcd(A, 0) -> A;
gcd(A, B) -> gcd(B, A rem B).

qadd({A, B}, {C, D}) -> norm_q(A * D + C * B, B * D).
qsub(X, {C, D}) -> qadd(X, {-C, D}).
qmul({A, B}, {C, D}) -> norm_q(A * C, B * D).
qdiv({A, B}, {C, D}) when C =/= 0 -> norm_q(A * D, B * C).
qneg({A, B}) -> {-A, B}.
qabs({A, B}) -> {abs(A), B}.
qpow(_, 0) -> {1, 1};
qpow(X, K) when K > 0 -> qmul(X, qpow(X, K - 1));
qpow(X, K) -> qdiv({1, 1}, qpow(X, -K)).
sgn({0, _}) -> 0;
sgn({N, _}) when N > 0 -> 1;
sgn(_) -> -1.

%%% pq_mp shorthands

mp({N, D}, W) -> pq_mp:from({N, D}, W).
from(N, W) -> pq_mp:from(N, W).
one() -> {mp, 1, 0}.
add(X, Y, W) -> pq_mp:add(X, Y, W).
sub(X, Y, W) -> pq_mp:sub(X, Y, W).
mul(X, Y, W) -> pq_mp:mul(X, Y, W).
divide(X, Y, W) -> pq_mp:divide(X, Y, W).
neg(X) -> pq_mp:neg(X).
half(X) -> pq_mp:ldexp(X, -1).
log_(Q, W) -> pq_mp:log(mp(Q, W), W).
round_c({R, I}, P) -> {pq_mp:from(R, P), pq_mp:from(I, P)}.
