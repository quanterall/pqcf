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
-export([a0/2, a0/3, b0/4, b0/5, b1/4, b1/5, b00/4, b00/5, b11/4, b11/5, db0/4, db0/5, value/3]).

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
