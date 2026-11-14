%% Arbitrary-precision real numbers (decision record SY2, phase 2): binary floating point on Erlang's integers, the
%% elementary and special functions the series and the eigenvalue solvers need, and the evaluation of expressions
%% to a requested number of decimal digits.
%%
%% A number is {mp, M, E}: M 2^E, |M| below 2^P for the working precision P (bits), rounded to nearest. Functions
%% take P explicitly; they add their own guard bits and round the result to P.
%%
%% Methods (each checked against mpmath to 50 digits, test/fixtures/gt_mp.json):
%% - pi: Machin's formula; log 2: 2 atanh(1/3); Euler's gamma: Brent and McMillan's algorithm B1; zeta(k):
%%   Borwein's algorithm 2 in exact rationals;
%% - exp: reduction by log 2 and halving, Taylor, squaring; log: Newton on exp; sin, cos: reduction modulo pi/2,
%%   Taylor; atan by halving;
%% - Gamma: shift up, Stirling's series with exact Bernoulli numbers, reflection below 1/2; polygamma likewise;
%% - Bessel J_nu (nu rational) and Y_n (n an integer): their power series with guard bits for the cancellation (x up
%%   to a few hundred); their zeros: McMahon's asymptotic start, Newton;
%% - Laguerre L_n^alpha: the sum; 2F1(a, b; c; z) for |z| < 1: the series.
-module(pq_mp).
-export([from/2, from_string/2, to_float/1, to_string/2, add/3, sub/3, mul/3, divide/3, neg/1, sqrt/2, cmp/2, ldexp/2,
         pi/1, ln2/1, euler_gamma/1, zeta/2, exp/2, log/2, sin/2, cos/2, gamma/2, polygamma/3,
         besselj/3, bessely/3, besselj_zero/3, laguerre/4, hyp2f1/5, evaluate/3, pslq/3]).

-export_type([mp/0]).
-type mp() :: {mp, integer(), integer()}.

%%% Core

norm(0, _, _) -> {mp, 0, 0};
norm(M, E, P) ->
    L = bitlen(abs(M)),
    case L > P of
        true ->
            S = L - P,
            % round to nearest, ties away from zero
            Q = (abs(M) + (1 bsl (S - 1))) bsr S,
            {mp, sign(M) * Q, E + S};
        false -> {mp, M, E}
    end.

sign(M) when M < 0 -> -1;
sign(_) -> 1.

bitlen(0) -> 0;
bitlen(N) -> bitlen(N, 0).
bitlen(N, Acc) when N >= 1 bsl 64 -> bitlen(N bsr 64, Acc + 64);
bitlen(N, Acc) when N > 0 -> bitlen(N bsr 1, Acc + 1);
bitlen(0, Acc) -> Acc.

%% an integer, a fraction {rat, N, D} or {N, D}, a float, or an mp, at precision P
-spec from(term(), pos_integer()) -> mp().
from({mp, M, E}, P) -> norm(M, E, P);
from(N, _) when is_integer(N) -> {mp, N, 0};
from({rat, N, D}, P) -> from({N, D}, P);
from({N, D}, P) when is_integer(N), is_integer(D), D > 0 ->
    S = P + bitlen(D) + 2,
    norm((N bsl S) div D, -S, P);
from(F, P) when is_float(F) ->
    {M, E} = float_parts(F),
    norm(M, E, P).

float_parts(F) when F == 0.0 -> {0, 0};
float_parts(F) ->
    <<S:1, Ex:11, Frac:52>> = <<F/float>>,
    {M, E} = case Ex of 0 -> {Frac, -1074}; _ -> {Frac bor (1 bsl 52), Ex - 1075} end,
    {case S of 1 -> -M; 0 -> M end, E}.

%% a decimal string "d.ddd[e±k]" at precision P
-spec from_string(string(), pos_integer()) -> mp().
from_string(S, P) ->
    {Mant, Exp} = case string:split(S, "e") of
                      [M, E] -> {M, list_to_integer(E)};
                      [M] -> {M, 0}
                  end,
    {Sign, Digits} = case Mant of [$- | R] -> {-1, R}; R -> {1, R} end,
    {Int, Frac} = case string:split(Digits, ".") of [I, F] -> {I, F}; [I] -> {I, ""} end,
    N = Sign * list_to_integer(Int ++ Frac),
    E2 = Exp - length(Frac),
    case E2 >= 0 of
        true -> from(N * pow10(E2), P);
        false -> from({N, pow10(-E2)}, P)
    end.

-spec to_float(mp()) -> float().
to_float({mp, 0, _}) -> 0.0;
to_float({mp, M, E}) ->
    L = bitlen(abs(M)),
    S = max(0, L - 60),
    (M bsr S) * math:pow(2, E + S).

neg({mp, M, E}) -> {mp, -M, E}.

add({mp, 0, _}, B, P) -> from(B, P);
add(A, {mp, 0, _}, P) -> from(A, P);
add({mp, Ma, Ea} = A, {mp, Mb, Eb} = B, P) ->
    % a term far below the other's last bit only affects the rounding: kept to P + 4 bits of it
    Ta = Ea + bitlen(abs(Ma)), Tb = Eb + bitlen(abs(Mb)),
    if
        Ta - Tb > P + 4 -> add(A, trunc_to(B, Ta - P - 4), P);
        Tb - Ta > P + 4 -> add(trunc_to(A, Tb - P - 4), B, P);
        true ->
            E = min(Ea, Eb),
            norm((Ma bsl (Ea - E)) + (Mb bsl (Eb - E)), E, P)
    end.

%% A rounded to exponent at least Emin (a tiny remainder kept as one sticky unit)
trunc_to({mp, M, E}, Emin) when E >= Emin -> {mp, M, E};
trunc_to({mp, M, E}, Emin) ->
    S = Emin - E,
    Q = abs(M) bsr S,
    {mp, sign(M) * (Q bor 1), Emin}.

sub(A, B, P) -> add(A, neg(B), P).
mul({mp, Ma, Ea}, {mp, Mb, Eb}, P) -> norm(Ma * Mb, Ea + Eb, P).

divide(_, {mp, 0, _}, _) -> error(badarith);
divide({mp, Ma, Ea}, {mp, Mb, Eb}, P) ->
    S = max(0, P + bitlen(abs(Mb)) - bitlen(abs(Ma)) + 2),
    norm((Ma bsl S) div Mb, Ea - Eb - S, P).

ldexp({mp, M, E}, K) -> {mp, M, E + K}.

cmp(A, B) ->
    {mp, M, _} = sub(A, B, 64),
    if M > 0 -> 1; M < 0 -> -1; true -> 0 end.

sqrt({mp, 0, _}, _) -> {mp, 0, 0};
sqrt({mp, M, E}, P) when M > 0 ->
    S0 = max(0, 2 * P + 4 - bitlen(M)),
    S = case (E - S0) band 1 of 0 -> S0; 1 -> S0 + 1 end,
    norm(isqrt(M bsl S), (E - S) div 2, P).

isqrt(0) -> 0;
isqrt(N) -> isqrt(N, 1 bsl ((bitlen(N) + 1) div 2)).
isqrt(N, X) ->
    Y = (X + N div X) bsr 1,
    case Y >= X of true -> X; false -> isqrt(N, Y) end.

is_zero({mp, M, _}) -> M =:= 0.
mag({mp, 0, _}) -> -1 bsl 30;
mag({mp, M, E}) -> E + bitlen(abs(M)).   % 2^(mag - 1) =< |x| < 2^mag

%% D decimal digits (significant), as a string "d.ddd...e±k"
-spec to_string(mp(), pos_integer()) -> string().
to_string({mp, 0, _}, _) -> "0";
to_string({mp, M, E} = X, D) ->
    % 10^K x rounded to D digits: K chosen from the magnitude
    Dec = floor((mag(X) - 1) * 0.30102999566398120),
    K = D - 1 - Dec,
    N0 = scaled_int(M, E, K),
    {N, K2} = case length(integer_to_list(abs(N0))) of
                  L when L > D -> {round_div(N0, 10), K - 1};
                  L when L < D -> {scaled_int(M, E, K + 1), K + 1};
                  _ -> {N0, K}
              end,
    Ds = integer_to_list(abs(N)),
    Exp = length(Ds) - 1 - K2,
    [$- || N < 0] ++ [hd(Ds), $.] ++ tl(Ds) ++ "e" ++ integer_to_list(Exp).

%% round(M 2^E 10^K)
scaled_int(M, E, K) ->
    {Num, Den} = case K >= 0 of true -> {M * pow10(K), 1}; false -> {M, pow10(-K)} end,
    case E >= 0 of
        true -> round_div(Num bsl E, Den);
        false -> round_div(Num, Den bsl (-E))
    end.

round_div(A, B) when A >= 0 -> (2 * A + B) div (2 * B);
round_div(A, B) -> -round_div(-A, B).

pow10(K) -> pow(10, K).
pow(_, 0) -> 1;
pow(B, K) when K rem 2 =:= 0 -> X = pow(B, K div 2), X * X;
pow(B, K) -> B * pow(B, K - 1).

%%% Constants (cached per precision in the calling process)

cached(Key, P, F) ->
    K = {?MODULE, Key, P},
    case get(K) of
        undefined -> V = F(), put(K, V), V;
        V -> V
    end.

%% pi = 16 atan(1/5) - 4 atan(1/239), in fixed point at P + 32 bits
pi(P) ->
    cached(pi, P, fun() ->
                          W = P + 32,
                          V = 16 * atan_inv(5, W) - 4 * atan_inv(239, W),
                          norm(V, -W, P)
                  end).

%% atan(1/N) 2^W, by its series
atan_inv(N, W) ->
    X = (1 bsl W) div N,
    N2 = N * N,
    atan_inv(X, N2, 1, X, 1).
atan_inv(0, _, _, Acc, _) -> Acc;
atan_inv(T, N2, K, Acc, S) ->
    T2 = T div N2,
    S2 = -S,
    case T2 div (2 * K + 1) of
        0 -> Acc;
        Q -> atan_inv(T2, N2, K + 1, Acc + S2 * Q, S2)
    end.

%% log 2 = 2 atanh(1/3) = 2 sum 1 / ((2k + 1) 3^(2k + 1))
ln2(P) ->
    cached(ln2, P, fun() ->
                           W = P + 32,
                           norm(2 * atanh_inv(3, W), -W, P)
                   end).

atanh_inv(N, W) ->
    X = (1 bsl W) div N,
    atanh_inv(X, N * N, 1, X).
atanh_inv(T, N2, K, Acc) ->
    T2 = T div N2,
    case T2 div (2 * K + 1) of
        0 -> Acc;
        Q -> atanh_inv(T2, N2, K + 1, Acc + Q)
    end.

%% Euler's gamma by Brent and McMillan (B1): with n large, A = sum_k (n^k / k!)^2 (H_k - log n), B = sum_k (n^k/k!)^2,
%% gamma = A / B to within e^(-4n)
euler_gamma(P) ->
    cached(euler_gamma, P, fun() ->
                                   W = P + 64,
                                   N = ceil(W * 0.6931471805599453 / 4) + 2,
                                   Kmax = ceil(3.6 * N) + 10,
                                   One = 1 bsl W,
                                   {mp, Ml, El} = log(from(N, W), W),
                                   LogN = shift(Ml, El + W),     % log n in fixed point
                                   N2 = N * N,
                                   eg_loop(1, Kmax, N2, -LogN, One, -LogN, One, W, P)
                           end).

%% a = (a n^2 / k + b) / k with b = b n^2 / k^2 first: the B1 recurrences in fixed point
eg_loop(K, Kmax, _, _, _, U, V, W, P) when K > Kmax -> divide(norm(U, -W, P + 32), norm(V, -W, P + 32), P);
eg_loop(K, Kmax, N2, A, B, U, V, W, P) ->
    B2 = B * N2 div (K * K),
    A2 = (A * N2 div K + B2) div K,
    eg_loop(K + 1, Kmax, N2, A2, B2, U + A2, V + B2, W, P).

shift(M, S) when S >= 0 -> M bsl S;
shift(M, S) -> M bsr (-S).

%% zeta(k), k >= 2 an integer, by Borwein's algorithm 2 in exact rationals:
%%   zeta(s) = -1 / (d_n (1 - 2^(1-s))) sum_(j<n) (-1)^j (d_j - d_n) / (j + 1)^s,
%%   d_j = n sum_(i=0..j) (n + i - 1)! 4^i / ((n - i)! (2i)!), error about 3 (3 + sqrt 8)^-n
zeta(K, P) when is_integer(K), K >= 2 ->
    cached({zeta, K}, P, fun() ->
                                 N = ceil((P + 16) * 0.6931471805599453 / 1.7627471740390860) + 2,
                                 Ds = dseq(N),
                                 Dn = lists:last(Ds),
                                 Sum = lists:foldl(
                                         fun({J, Dj}, Acc) ->
                                                 S = case J rem 2 of 0 -> 1; 1 -> -1 end,
                                                 {A, B} = qadd(Dj, neg_q(Dn)),
                                                 qadd(Acc, {S * A, B * pow(J + 1, K)})
                                         end, {0, 1}, lists:zip(lists:seq(0, N - 1), lists:droplast(Ds))),
                                 % times -1 / (d_n (1 - 2^(1-k))) = -2^(k-1) / (d_n (2^(k-1) - 1))
                                 T = pow(2, K - 1),
                                 {Sn, Sd} = Sum, {Dnn, Dnd} = Dn,
                                 from(qnorm({-Sn * T * Dnd, Sd * Dnn * (T - 1)}), P)
                         end).

%% d_0 .. d_n as fractions
dseq(N) ->
    {_, Ds} = lists:foldl(fun(I, {Acc, L}) ->
                                  Acc2 = qadd(Acc, {fact(N + I - 1) * pow(4, I), fact(N - I) * fact(2 * I)}),
                                  {An, Ad} = Acc2,
                                  {Acc2, [qnorm({N * An, Ad}) | L]}
                          end, {{0, 1}, []}, lists:seq(0, N)),
    lists:reverse(Ds).

qadd({A, B}, {C, D}) ->
    N = A * D + C * B, M = B * D,
    G = gcd(abs(N), M),
    {N div G, M div G}.

gcd(A, 0) -> A;
gcd(A, B) -> gcd(B, A rem B).

fact(N) -> fact(N, 1).
fact(0, Acc) -> Acc;
fact(N, Acc) -> fact(N - 1, Acc * N).

%%% Elementary functions

%% exp(x): x = k log 2 + r, r halved s times, Taylor, squared s times, times 2^k
exp(X0, P) ->
    X = from(X0, P + 32),
    case is_zero(X) of
        true -> {mp, 1, 0};
        false ->
            W = P + 32 + max(0, mag(X)),
            L2 = ln2(W),
            K = round(to_float(divide(X, L2, 64))),
            R = sub(X, mul(from(K, W), L2, W), W),
            S = isqrt(W) + 1,
            Rs = ldexp(R, -S),
            T = taylor_exp(Rs, W),
            Sq = lists:foldl(fun(_, A) -> mul(A, A, W) end, T, lists:seq(1, S)),
            from(ldexp(Sq, K), P)
    end.

taylor_exp(R, W) -> taylor_exp(R, W, 1, {mp, 1, 0}, {mp, 1, 0}).
taylor_exp(R, W, K, Term, Acc) ->
    T2 = divide(mul(Term, R, W), from(K, W), W),
    case is_zero(T2) orelse mag(T2) < mag(Acc) - W - 2 of
        true -> Acc;
        false -> taylor_exp(R, W, K + 1, T2, add(Acc, T2, W))
    end.

%% log(x), x > 0: x = m 2^e, log m by Newton on exp (y += 2 (m - e^y) / (m + e^y)), plus e log 2
log({mp, M, _} = X0, P) when M > 0 -> log_(X0, P);
log(X0, P) when not is_tuple(X0) -> log_(from(X0, P + 32), P).

log_(X0, P) ->
    X = from(X0, P + 32),
    E = mag(X),
    Mx = ldexp(X, -E),               % in [1/2, 1)
    W = P + 32,
    Y0 = from(math:log(to_float(Mx)), W),
    Y = log_newton(Mx, Y0, W, 0),
    from(add(Y, mul(from(E, W), ln2(W), W), W), P).

log_newton(M, Y, W, I) ->
    Ey = exp(Y, W),
    D = divide(mul({mp, 2, 0}, sub(M, Ey, W), W), add(M, Ey, W), W),
    Y2 = add(Y, D, W),
    case is_zero(D) orelse mag(D) < mag(Y2) - W + 4 orelse I > 64 of
        true -> Y2;
        false -> log_newton(M, Y2, W, I + 1)
    end.

%% sin and cos: x reduced modulo pi/2 to r in [-pi/4, pi/4], quadrant q
sin(X, P) -> {S, C, Q} = sincos_r(X, P), quad(Q, S, C).
cos(X, P) -> {S, C, Q} = sincos_r(X, P), quad(Q + 1, S, C).

quad(Q, S, C) ->
    case Q band 3 of 0 -> S; 1 -> C; 2 -> neg(S); 3 -> neg(C) end.

sincos_r(X0, P) ->
    X = from(X0, P + 32),
    W = P + 32 + max(0, mag(X)),
    Pi2 = ldexp(pi(W), -1),
    Q = round(to_float(divide(X, Pi2, 64))),
    R = sub(X, mul(from(Q, W), Pi2, W), W),
    {from(taylor_sin(R, W), P), from(taylor_cos(R, W), P), Q}.

taylor_sin(R, W) -> series_alt(R, W, R, 1, R).
taylor_cos(R, W) -> series_alt(R, W, {mp, 1, 0}, 0, {mp, 1, 0}).

%% sum of (-1)^k r^(2k + s) / (2k + s)!, from the first term
series_alt(R, W, Term, S, Acc) -> series_alt(R, W, mul(R, R, W), Term, S, Acc, 1).
series_alt(R, W, R2, Term, S, Acc, K) ->
    T2 = neg(divide(mul(Term, R2, W), from((2 * K + S - 1) * (2 * K + S), W), W)),
    case is_zero(T2) orelse mag(T2) < mag(Acc) - W - 2 of
        true -> Acc;
        false -> series_alt(R, W, R2, T2, S, add(Acc, T2, W), K + 1)
    end.

%%% Gamma and polygamma

%% Bernoulli numbers B_0 .. B_N as fractions (cached), by the Akiyama-Tanigawa algorithm: for m = 0, 1, ..:
%% A[m] = 1/(m+1), then A[j-1] = j (A[j-1] - A[j]) for j = m .. 1; B_m = A[0] (B_1 = +1/2: only the even ones are used)
bernoulli(N) ->
    K = {?MODULE, bernoulli},
    case get(K) of
        L when is_list(L), length(L) > N -> lists:sublist(L, N + 1);
        _ ->
            {_, Bs} = lists:foldl(fun(M, {A, Acc}) ->
                                          A1 = A#{M => {1, M + 1}},
                                          A2 = lists:foldl(fun(J, Ax) ->
                                                                   {N1, D1} = qadd(maps:get(J - 1, Ax), neg_q(maps:get(J, Ax))),
                                                                   Ax#{J - 1 => qnorm({J * N1, D1})}
                                                           end, A1, lists:seq(M, 1, -1)),
                                          {A2, [maps:get(0, A2) | Acc]}
                                  end, {#{}, []}, lists:seq(0, N)),
            L = lists:reverse(Bs),
            put(K, L),
            L
    end.

neg_q({A, B}) -> {-A, B}.
qnorm({A, B}) -> G = gcd(abs(A), B), {A div G, B div G}.

%% Gamma(x), x real (not a pole)
gamma(X0, P) ->
    X = from(X0, P + 32),
    W = P + 32,
    case to_float(X) < 0.5 of
        true ->
            % reflection: Gamma(x) = pi / (sin(pi x) Gamma(1 - x))
            Pi = pi(W),
            from(divide(Pi, mul(sin(mul(Pi, X, W), W), gamma(sub({mp, 1, 0}, X, W), W), W), W), P);
        false ->
            {Z, Prod} = shift_up(X, W),
            from(divide(exp(lngamma_stirling(Z, W), W), Prod, W), P)
    end.

%% x shifted to z = x + n with n large enough for Stirling's series at W bits, and x (x+1) .. (x+n-1)
shift_up(X, W) ->
    Zmin = W * 0.35 + 10,
    N = max(0, ceil(Zmin - to_float(X))),
    lists:foldl(fun(I, {_, Pr}) ->
                        Xi = add(X, from(I, W), W),
                        {add(Xi, {mp, 1, 0}, W), mul(Pr, Xi, W)}
                end, {X, {mp, 1, 0}}, lists:seq(0, N - 1)).

%% log Gamma(z) = (z - 1/2) log z - z + log(2 pi)/2 + sum_k B_2k / (2k (2k - 1) z^(2k - 1))
lngamma_stirling(Z, W) ->
    L = log(Z, W),
    Base = add(sub(mul(sub(Z, from({1, 2}, W), W), L, W), Z, W), ldexp(log(ldexp(pi(W), 1), W), -1), W),
    Bs = bernoulli(2 * stirling_terms(Z, W)),
    Iz = divide({mp, 1, 0}, Z, W),
    Iz2 = mul(Iz, Iz, W),
    {S, _} = lists:foldl(fun(K, {Acc, Pw}) ->
                                 {Bn, Bd} = lists:nth(2 * K + 1, Bs),
                                 T = mul(from({Bn, Bd * 2 * K * (2 * K - 1)}, W), Pw, W),
                                 {add(Acc, T, W), mul(Pw, Iz2, W)}
                         end, {{mp, 0, 0}, Iz}, lists:seq(1, stirling_terms(Z, W))),
    add(Base, S, W).

stirling_terms(Z, W) -> max(4, min(200, ceil(W * 0.6931471805599453 / (2 * math:log(to_float(Z)))) + 2)).

%% psi^(m)(x): shifted up by the recurrence, then the asymptotic series
%%   psi(z) ~ log z - 1/(2z) - sum_k B_2k / (2k z^2k)
%%   psi^(m)(z) ~ (-1)^(m+1) [ (m-1)!/z^m + m!/(2 z^(m+1)) + sum_k B_2k (2k+m-1)! / ((2k)! z^(2k+m)) ]
polygamma(M, X0, P) when is_integer(M), M >= 0 ->
    W = P + 32,
    X = from(X0, W),
    N = max(0, ceil(W * 0.35 + 10 - to_float(X))),
    % psi^(m)(x) = psi^(m)(x + n) - (-1)^m m! sum_(i<n) 1/(x + i)^(m+1)
    Corr = lists:foldl(fun(I, Acc) ->
                               Xi = add(X, from(I, W), W),
                               add(Acc, divide({mp, 1, 0}, ipow(Xi, M + 1, W), W), W)
                       end, {mp, 0, 0}, lists:seq(0, N - 1)),
    Z = add(X, from(N, W), W),
    Asym = psi_asym(M, Z, W),
    Sgn = case M rem 2 of 0 -> 1; 1 -> -1 end,
    from(sub(Asym, mul(from(Sgn * fact(M), W), Corr, W), W), P).

psi_asym(0, Z, W) ->
    K = stirling_terms(Z, W),
    Bs = bernoulli(2 * K),
    Iz2 = divide({mp, 1, 0}, mul(Z, Z, W), W),
    {S, _} = lists:foldl(fun(J, {Acc, Pw}) ->
                                 {Bn, Bd} = lists:nth(2 * J + 1, Bs),
                                 {add(Acc, mul(from({Bn, Bd * 2 * J}, W), Pw, W), W), mul(Pw, Iz2, W)}
                         end, {{mp, 0, 0}, Iz2}, lists:seq(1, K)),
    sub(sub(log(Z, W), divide({mp, 1, 0}, ldexp(Z, 1), W), W), S, W);
psi_asym(M, Z, W) ->
    K = stirling_terms(Z, W),
    Bs = bernoulli(2 * K),
    Iz = divide({mp, 1, 0}, Z, W),
    T1 = mul(from(fact(M - 1), W), ipow(Iz, M, W), W),
    T2 = mul(from({fact(M), 2}, W), ipow(Iz, M + 1, W), W),
    Iz2 = mul(Iz, Iz, W),
    {S, _} = lists:foldl(fun(J, {Acc, Pw}) ->
                                 {Bn, Bd} = lists:nth(2 * J + 1, Bs),
                                 C = from({Bn * fact(2 * J + M - 1), Bd * fact(2 * J)}, W),
                                 {add(Acc, mul(C, Pw, W), W), mul(Pw, Iz2, W)}
                         end, {{mp, 0, 0}, mul(ipow(Iz, M, W), Iz2, W)}, lists:seq(1, K)),
    Sgn = case (M + 1) rem 2 of 0 -> 1; 1 -> -1 end,
    mul(from(Sgn, W), add(add(T1, T2, W), S, W), W).

ipow(_, 0, _) -> {mp, 1, 0};
ipow(X, K, W) when K rem 2 =:= 0 -> Y = ipow(X, K div 2, W), mul(Y, Y, W);
ipow(X, K, W) -> mul(X, ipow(X, K - 1, W), W).

%%% Bessel functions

%% J_nu(x), nu rational: sum_k (-1)^k (x/2)^(2k + nu) / (k! Gamma(k + nu + 1)); J_(-n) = (-1)^n J_n
besselj(Nu0, X0, P) ->
    {Nn, Nd} = rat(Nu0),
    case Nd =:= 1 andalso Nn < 0 of
        true ->
            V = besselj(-Nn, X0, P),
            case Nn rem 2 of 0 -> V; _ -> neg(V) end;
        false ->
            X = from(X0, P + 32),
            % guard bits for the cancellation: the largest term is about e^|x|
            W = P + 40 + ceil(abs(to_float(X)) * 1.4426950408889634),
            H = ldexp(from(X, W), -1),
            H2 = neg(mul(H, H, W)),
            Nu = from({Nn, Nd}, W),
            G0 = gamma(add(Nu, {mp, 1, 0}, W), W),
            T0 = divide(power(H, Nu, W), G0, W),
            S = series_sum(fun(K, T) -> divide(mul(T, H2, W), mul(from(K, W), add(Nu, from(K, W), W), W), W) end, T0, W),
            from(S, P)
    end.

%% x^a for x > 0 and a real, or a an integer
power(X, {mp, M, E} = A, W) ->
    case E >= 0 of
        true -> ipow_signed(X, M bsl E, W);
        false ->
            case is_zero(X) of
                true -> {mp, 0, 0};
                false -> exp(mul(A, log(X, W), W), W)
            end
    end.

ipow_signed(X, K, W) when K >= 0 -> ipow(X, K, W);
ipow_signed(X, K, W) -> divide({mp, 1, 0}, ipow(X, -K, W), W).

%% T_0 + T_1 + ..., T_k = Next(k, T_(k-1)), until the terms are negligible (past their peak)
series_sum(Next, T0, W) -> series_sum(Next, T0, W, 1, T0, mag(T0)).
series_sum(Next, T, W, K, Acc, Peak) ->
    T2 = Next(K, T),
    Peak2 = max(Peak, mag(T2)),
    case is_zero(T2) orelse (mag(T2) < Peak2 - W - 4 andalso mag(T2) < mag(Acc) - W - 4) of
        true -> Acc;
        false -> series_sum(Next, T2, W, K + 1, add(Acc, T2, W), Peak2)
    end.

%% Y_n(x), n an integer (the series of pq_series, evaluated):
%%   (2/pi) J_n(x) log(x/2) - (1/pi) sum_(k<n) (n-k-1)!/k! (x/2)^(2k-n)
%%   - (1/pi) sum_k (psi(k+1) + psi(n+k+1)) (-1)^k (x/2)^(2k+n) / (k! (n+k)!)
bessely(N, X0, P) when is_integer(N), N < 0 ->
    V = bessely(-N, X0, P),
    case N rem 2 of 0 -> V; _ -> neg(V) end;
bessely(N, X0, P) when is_integer(N) ->
    X = from(X0, P + 32),
    W = P + 48 + ceil(abs(to_float(X)) * 1.4426950408889634),
    H = ldexp(from(X, W), -1),
    Pi = pi(W),
    G = euler_gamma(W),
    J = besselj(N, X, W),
    T1 = mul(divide({mp, 2, 0}, Pi, W), mul(J, log(H, W), W), W),
    Hinv = divide({mp, 1, 0}, H, W),
    T2 = lists:foldl(fun(K, Acc) ->
                             Pw = ipow_signed(H, 2 * K - N, W),
                             add(Acc, mul(from({fact(N - K - 1), fact(K)}, W), Pw, W), W)
                     end, {mp, 0, 0}, lists:seq(0, N - 1)),
    _ = Hinv,
    Psi = fun(M) -> sub(from(harmonic(M), W), G, W) end,
    H2 = neg(mul(H, H, W)),
    % the k-th term: (psi(k+1) + psi(n+k+1)) t_k, t_k = (-1)^k (x/2)^(2k+n) / (k! (n+k)!)
    T0 = divide(ipow(H, N, W), from(fact(N), W), W),
    S = sum_y(Psi, H2, T0, N, W),
    from(sub(T1, divide(add(T2, S, W), Pi, W), W), P).

sum_y(Psi, H2, T0, N, W) -> sum_y(Psi, H2, N, W, 0, T0, mul(add(Psi(0), Psi(N), W), T0, W), mag(T0)).
sum_y(Psi, H2, N, W, K, T, Acc, Peak) ->
    K1 = K + 1,
    T2 = divide(mul(T, H2, W), from(K1 * (N + K1), W), W),
    Term = mul(add(Psi(K1), Psi(N + K1), W), T2, W),
    Peak2 = max(Peak, mag(T2)),
    case is_zero(T2) orelse (mag(T2) < Peak2 - W - 4 andalso mag(Term) < mag(Acc) - W - 4) of
        true -> Acc;
        false -> sum_y(Psi, H2, N, W, K1, T2, add(Acc, Term, W), Peak2)
    end.

harmonic(M) -> lists:foldl(fun(I, Acc) -> qadd(Acc, {1, I}) end, {0, 1}, lists:seq(1, M)).

%% the K-th positive zero of J_nu: McMahon's expansion for the start, then Newton with
%% J_nu'(x) = J_(nu-1)(x) - (nu/x) J_nu(x)
besselj_zero(Nu0, K, P) when is_integer(K), K >= 1 ->
    {Nn, Nd} = rat(Nu0),
    Nuf = Nn / Nd,
    Beta = (K + Nuf / 2 - 0.25) * math:pi(),
    Mu = 4 * Nuf * Nuf,
    X0 = Beta - (Mu - 1) / (8 * Beta) - 4 * (Mu - 1) * (7 * Mu - 31) / (3 * math:pow(8 * Beta, 3)),
    W = P + 16,
    zero_newton({Nn, Nd}, from(X0, W), W, 0).

zero_newton(Nu, X, W, I) ->
    J = besselj(Nu, X, W),
    Jm = besselj(qsub(Nu, {1, 1}), X, W),
    Dj = sub(Jm, mul(from(Nu, W), divide(J, X, W), W), W),
    D = divide(J, Dj, W),
    X2 = sub(X, D, W),
    case is_zero(D) orelse mag(D) < mag(X2) - W + 6 orelse I > 60 of
        true -> X2;
        false -> zero_newton(Nu, X2, W, I + 1)
    end.

qsub({A, B}, {C, D}) -> qadd({A, B}, {-C, D}).

rat(N) when is_integer(N) -> {N, 1};
rat({rat, N, D}) -> {N, D};
rat({N, D}) when is_integer(N), is_integer(D) -> {N, D}.

%%% Laguerre and 2F1

%% L_n^alpha(x) = sum_(k=0..n) (-1)^k binom(n + alpha, n - k) x^k / k!
laguerre(N, Alpha0, X0, P) ->
    W = P + 32 + 4 * N,
    A = from(Alpha0, W), X = from(X0, W),
    S = lists:foldl(fun(K, Acc) ->
                            B = lists:foldl(fun(I, Pr) -> mul(Pr, add(A, from(I, W), W), W) end, {mp, 1, 0}, lists:seq(K + 1, N)),
                            C = divide(B, from((case K rem 2 of 0 -> 1; 1 -> -1 end) * fact(N - K) * fact(K), W), W),
                            add(Acc, mul(C, ipow(X, K, W), W), W)
                    end, {mp, 0, 0}, lists:seq(0, N)),
    from(S, P).

%% 2F1(a, b; c; z) = sum_j (a)_j (b)_j / ((c)_j j!) z^j, |z| < 1
hyp2f1(A0, B0, C0, Z0, P) ->
    W = P + 32,
    A = from(A0, W), B = from(B0, W), C = from(C0, W), Z = from(Z0, W),
    S = series_sum(fun(J, T) ->
                           Jm = from(J - 1, W),
                           Num = mul(mul(add(A, Jm, W), add(B, Jm, W), W), Z, W),
                           divide(mul(T, Num, W), mul(add(C, Jm, W), from(J, W), W), W)
                   end, {mp, 1, 0}, W),
    from(S, P).

%%% Expressions

%% the value of an expression to D significant decimal digits: {ok, Digits} (a string), the symbols bound by
%% Bindings (numbers: integers, fractions, floats); the atoms of pq_series included
-spec evaluate(term(), map(), pos_integer()) -> {ok, string()} | {error, term()}.
evaluate(E, Bindings, D) ->
    P = ceil(D * 3.3219280948873622) + 16,
    try
        {ok, to_string(ev(E, Bindings, P), D)}
    catch
        throw:{mp, Why} -> {error, {symbolic_unsupported, #{evaluate => Why}}};
        error:badarith -> {error, {non_finite, arith}}
    end.

ev(N, _, _) when is_integer(N) -> {mp, N, 0};
ev({rat, N, D}, _, P) -> from({N, D}, P);
ev(F, _, P) when is_float(F) -> from(F, P);
ev({sym, X}, B, P) ->
    case maps:find(X, B) of
        {ok, V} -> from(V, P);
        error -> throw({mp, {unbound, X}})
    end;
ev({add, L}, B, P) -> lists:foldl(fun(T, Acc) -> add(Acc, ev(T, B, P), P) end, {mp, 0, 0}, L);
ev({mul, L}, B, P) -> lists:foldl(fun(T, Acc) -> mul(Acc, ev(T, B, P), P) end, {mp, 1, 0}, L);
ev({pow, X, K}, B, P) when is_integer(K) -> ipow_signed(ev(X, B, P), K, P);
ev({pow, X, Y}, B, P) -> power(ev(X, B, P), ev(Y, B, P), P);
ev({sqrt, X}, B, P) -> sqrt(ev(X, B, P), P);
ev({exp, X}, B, P) -> exp(ev(X, B, P), P);
ev({log, X}, B, P) -> log(ev(X, B, P), P);
ev({sin, X}, B, P) -> sin(ev(X, B, P), P);
ev({cos, X}, B, P) -> cos(ev(X, B, P), P);
ev({const, pi}, _, P) -> pi(P);
ev({const, euler_gamma}, _, P) -> euler_gamma(P);
ev({zeta, K}, _, P) -> zeta(K, P);
ev({gamma, X}, B, P) -> gamma(ev(X, B, P), P);
ev({polygamma, M, X}, B, P) -> polygamma(M, ev(X, B, P), P);
ev({besselj, Nu, X}, B, P) -> besselj(rational_order(Nu), ev(X, B, P), P);
ev({bessely, N, X}, B, P) when is_integer(N) -> bessely(N, ev(X, B, P), P);
ev({laguerre, N, A, X}, B, P) when is_integer(N) -> laguerre(N, ev(A, B, P), ev(X, B, P), P);
ev({hyp2f1, A, Bb, C, Z}, B, P) -> hyp2f1(ev(A, B, P), ev(Bb, B, P), ev(C, B, P), ev(Z, B, P), P);
ev(E, _, _) -> throw({mp, {unsupported, case E of T when is_tuple(T) -> element(1, T); _ -> E end}}).

rational_order(N) when is_integer(N) -> N;
rational_order({rat, _, _} = R) -> R;
rational_order(_) -> throw({mp, {order_must_be_rational, besselj}}).

%%% Integer relations: PSLQ

%% integers m (not all zero) with sum m_i x_i = 0 to the precision P, |m_i| =< MaxCoeff: {ok, Ms} or none. Ferguson
%% and Bailey's PSLQ (H. R. P. Ferguson, D. H. Bailey, S. Arno, Math. Comp. 68 (1999) 351), real version, gamma =
%% 2/sqrt(3); the relation is accepted only when it holds to P - 16 bits.
-spec pslq([term()], pos_integer(), pos_integer()) -> {ok, [integer()]} | none.
pslq(Xs0, P, MaxCoeff) ->
    W = P + 32,
    Xs = [from(X, W) || X <- Xs0],
    N = length(Xs),
    case [I || {I, X} <- lists:enumerate(Xs), is_zero(X)] of
        [I | _] -> {ok, [case J of I -> 1; _ -> 0 end || J <- lists:seq(1, N)]};
        [] ->
            % s_k = sqrt(sum_(j >= k) x_j^2), y = x / s_1, s = s / s_1
            Sq = [sqrt(lists:foldl(fun(X, A) -> add(A, mul(X, X, W), W) end, {mp, 0, 0}, lists:nthtail(K - 1, Xs)), W)
                  || K <- lists:seq(1, N)],
            S1 = hd(Sq),
            Y = [divide(X, S1, W) || X <- Xs],
            S = [divide(X, S1, W) || X <- Sq],
            St = list_to_tuple(S), Yt = list_to_tuple(Y),
            H = maps:from_list([{{I, J}, if
                                             J > I -> {mp, 0, 0};
                                             J =:= I -> divide(element(I + 1, St), element(I, St), W);
                                             true -> neg(divide(mul(element(I, Yt), element(J, Yt), W),
                                                                mul(element(J, St), element(J + 1, St), W), W))
                                         end} || I <- lists:seq(1, N), J <- lists:seq(1, N - 1)]),
            Id = maps:from_list([{{I, J}, case I of J -> 1; _ -> 0 end} || I <- lists:seq(1, N), J <- lists:seq(1, N)]),
            St0 = #{n => N, w => W, y => maps:from_list(lists:enumerate(Y)), h => H, a => Id, b => Id},
            pslq_loop(reduce(St0), Xs, P, MaxCoeff, 0)
    end.

pslq_loop(_, _, _, _, 5000) -> none;
pslq_loop(#{n := N, w := W, y := Y, h := H, b := B} = St, Xs, P, MaxCoeff, It) ->
    % a small y_j: column j of B is a relation
    Small = [J || J <- lists:seq(1, N), is_zero(maps:get(J, Y)) orelse mag(maps:get(J, Y)) < -(P - 16)],
    Col = fun(J) -> [maps:get({I, J}, B) || I <- lists:seq(1, N)] end,
    case [C || J <- Small, C <- [Col(J)], lists:max([abs(V) || V <- C]) =< MaxCoeff, holds(C, Xs, P, W)] of
        [Rel | _] -> {ok, normalise(Rel)};
        [] ->
            % a bound on any relation's norm: 1 / max |H_jj|; past MaxCoeff, give up
            Hmax = lists:max([abs(to_float(maps:get({J, J}, H))) || J <- lists:seq(1, N - 1)]),
            case Hmax > 0 andalso 1 / Hmax > MaxCoeff * math:sqrt(N) of
                true -> none;
                false -> pslq_loop(reduce(step(St)), Xs, P, MaxCoeff, It + 1)
            end
    end.

holds(C, Xs, P, W) ->
    R = lists:foldl(fun({M, X}, A) -> add(A, mul(from(M, W), X, W), W) end, {mp, 0, 0}, lists:zip(C, Xs)),
    Scale = lists:max([mag(X) || X <- Xs]),
    is_zero(R) orelse mag(R) < Scale - (P - 16).

normalise(C) ->
    case [V || V <- C, V =/= 0] of
        [F | _] when F < 0 -> [-V || V <- C];
        _ -> C
    end.

%% one PSLQ step: the m maximising gamma^m |H_mm|, rows m and m+1 exchanged, the corner restored
step(#{n := N, w := W, y := Y, h := H, a := A, b := B} = St) ->
    G = 1.1547005383792515,   % 2 / sqrt 3
    {_, M} = lists:max([{math:pow(G, I) * abs(to_float(maps:get({I, I}, H))), I} || I <- lists:seq(1, N - 1)]),
    SwapRow = fun(Mp, Cols) -> maps:merge(Mp, maps:from_list(
                                               lists:append([[{{M, K}, maps:get({M + 1, K}, Mp)}, {{M + 1, K}, maps:get({M, K}, Mp)}] || K <- Cols]))) end,
    Y2 = Y#{M => maps:get(M + 1, Y), M + 1 => maps:get(M, Y)},
    A2 = SwapRow(A, lists:seq(1, N)),
    B2 = maps:merge(B, maps:from_list(lists:append([[{{K, M}, maps:get({K, M + 1}, B)}, {{K, M + 1}, maps:get({K, M}, B)}] || K <- lists:seq(1, N)]))),
    H2 = SwapRow(H, lists:seq(1, N - 1)),
    H3 = case M =< N - 2 of
             false -> H2;
             true ->
                 Hmm = maps:get({M, M}, H2), Hm1 = maps:get({M, M + 1}, H2),
                 T0 = sqrt(add(mul(Hmm, Hmm, W), mul(Hm1, Hm1, W), W), W),
                 T1 = divide(Hmm, T0, W), T2 = divide(Hm1, T0, W),
                 lists:foldl(fun(I, Hx) ->
                                     T3 = maps:get({I, M}, Hx), T4 = maps:get({I, M + 1}, Hx),
                                     Hx#{{I, M} => add(mul(T1, T3, W), mul(T2, T4, W), W),
                                         {I, M + 1} => add(neg(mul(T2, T3, W)), mul(T1, T4, W), W)}
                             end, H2, lists:seq(M, N))
         end,
    St#{y => Y2, h => H3, a => A2, b => B2}.

%% Hermite reduction of H, applied to y, A and B
reduce(#{n := N, w := W} = St) ->
    lists:foldl(fun({I, J}, #{y := Y, h := H, a := A, b := B} = S) ->
                        Hjj = maps:get({J, J}, H),
                        case is_zero(Hjj) of
                            true -> S;
                            false ->
                                T = round_mp(divide(maps:get({I, J}, H), Hjj, W)),
                                case T of
                                    0 -> S;
                                    _ ->
                                        Tm = from(T, W),
                                        Y2 = Y#{J => add(maps:get(J, Y), mul(Tm, maps:get(I, Y), W), W)},
                                        H2 = lists:foldl(fun(K, Hx) -> Hx#{{I, K} => sub(maps:get({I, K}, Hx), mul(Tm, maps:get({J, K}, Hx), W), W)} end,
                                                         H, lists:seq(1, J)),
                                        A2 = lists:foldl(fun(K, Ax) -> Ax#{{I, K} => maps:get({I, K}, Ax) - T * maps:get({J, K}, Ax)} end,
                                                         A, lists:seq(1, N)),
                                        B2 = lists:foldl(fun(K, Bx) -> Bx#{{K, J} => maps:get({K, J}, Bx) + T * maps:get({K, I}, Bx)} end,
                                                         B, lists:seq(1, N)),
                                        S#{y => Y2, h => H2, a => A2, b => B2}
                                end
                        end
                end, St, [{I, J} || I <- lists:seq(2, N), J <- lists:seq(I - 1, 1, -1)]).

%% the nearest integer
round_mp({mp, M, E}) when E >= 0 -> M bsl E;
round_mp({mp, M, E}) ->
    S = -E,
    case M >= 0 of
        true -> (M + (1 bsl (S - 1))) bsr S;
        false -> -((-M + (1 bsl (S - 1))) bsr S)
    end.
