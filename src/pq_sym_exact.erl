%% Exact numbers of the symbolic engine (SY 3.1): rationals {N, D} in lowest terms (D > 0), coefficients that are
%% Gaussian rationals {q, Re, Im} or float pairs {f, Re, Im} (floats are contagious), squarefree parts for the
%% radicals, and the exact values of exp(i k pi / 12) (SY 4.4).
-module(pq_sym_exact).
-export([rat/1, rat/2, radd/2, rsub/2, rmul/2, rneg/1, rinv/1, rcmp/2, rfloor/1, rceil/1, is_int/1, to_float/1,
         c/1, cf/2, cq/2, cadd/2, cmul/2, cneg/1, cinv/1, cconj/1, cscale/2, is_zero/1, is_one/1, to_complex/1,
         is_exact/1, is_real/1, real_rational/1, csqrt_exact/1, cpow_int/2,
         squarefree/1, gcd/2, root_exact/2, pi_angle/1]).

-export_type([rational/0, coeff/0]).

-type rational() :: {integer(), pos_integer()}.
-type coeff() :: {q, rational(), rational()} | {f, float(), float()}.

-define(TRIAL_LIMIT, 1000000).

%%% Rationals

-spec rat(integer()) -> rational().
rat(N) when is_integer(N) -> {N, 1}.

-spec rat(integer(), integer()) -> rational().
rat(_, 0) -> error(badarith);
rat(N, D) when D < 0 -> rat(-N, -D);
rat(N, D) ->
    G = gcd(abs(N), D),
    {N div G, D div G}.

gcd(A, 0) -> abs(A);
gcd(A, B) -> gcd(B, A rem B).

radd({A, B}, {C, D}) -> rat(A * D + C * B, B * D).
rsub(X, Y) -> radd(X, rneg(Y)).
rmul({A, B}, {C, D}) -> rat(A * C, B * D).
rneg({A, B}) -> {-A, B}.
rinv({0, _}) -> error(badarith);
rinv({A, B}) -> rat(B, A).

%% -1, 0 or 1
rcmp({A, B}, {C, D}) ->
    L = A * D, R = C * B,
    if L < R -> -1; L > R -> 1; true -> 0 end.

rfloor({A, B}) when A >= 0 -> {A div B, 1};
rfloor({A, B}) -> {-((-A + B - 1) div B), 1}.

rceil(X) -> rneg(rfloor(rneg(X))).

is_int({_, 1}) -> true;
is_int(_) -> false.

to_float({A, B}) -> A / B.

%%% Coefficients

-spec c(integer() | rational() | float()) -> coeff().
c(N) when is_integer(N) -> {q, {N, 1}, {0, 1}};
c(X) when is_float(X) -> {f, X + 0.0, 0.0};
c({_, _} = R) -> {q, R, {0, 1}}.

%% A float pair, -0.0 turned into 0.0 so that equal values are equal terms.
cf(Re, Im) -> {f, float(Re) + 0.0, float(Im) + 0.0}.
cq(Re, Im) -> {q, Re, Im}.

to_complex({q, R, I}) -> {to_float(R), to_float(I)};
to_complex({f, R, I}) -> {R, I}.

cadd({q, A, B}, {q, C, D}) -> {q, radd(A, C), radd(B, D)};
cadd(X, Y) -> {A, B} = to_complex(X), {C, D} = to_complex(Y), cf(A + C, B + D).

cmul({q, A, B}, {q, C, D}) -> {q, rsub(rmul(A, C), rmul(B, D)), radd(rmul(A, D), rmul(B, C))};
cmul(X, Y) -> {A, B} = to_complex(X), {C, D} = to_complex(Y), cf(A * C - B * D, A * D + B * C).

cneg({q, A, B}) -> {q, rneg(A), rneg(B)};
cneg({f, A, B}) -> cf(-A, -B).

cinv({q, A, B}) ->
    N = radd(rmul(A, A), rmul(B, B)),
    {q, rmul(A, rinv(N)), rneg(rmul(B, rinv(N)))};
cinv({f, A, B}) ->
    N = A * A + B * B,
    N == 0 andalso error(badarith),
    cf(A / N, -B / N).

cconj({q, A, B}) -> {q, A, rneg(B)};
cconj({f, A, B}) -> cf(A, -B).

cscale(C, R) -> cmul(C, c(R)).

is_zero({q, {0, _}, {0, _}}) -> true;
is_zero({q, _, _}) -> false;
is_zero({f, A, B}) -> A == 0 andalso B == 0.

is_one({q, {1, 1}, {0, 1}}) -> true;
is_one(_) -> false.

is_exact({q, _, _}) -> true;
is_exact(_) -> false.

is_real({q, _, {0, _}}) -> true;
is_real({f, _, I}) -> I == 0;
is_real(_) -> false.

%% The rational value of a real exact coefficient, else false.
real_rational({q, R, {0, _}}) -> R;
real_rational(_) -> false.

%% C^N for an integer N, by repeated squaring.
cpow_int(_C, 0) -> c(1);
cpow_int(C, N) when N < 0 -> cinv(cpow_int(C, -N));
cpow_int(C, 1) -> C;
cpow_int(C, N) ->
    H = cpow_int(C, N div 2),
    H2 = cmul(H, H),
    case N rem 2 of 0 -> H2; 1 -> cmul(H2, C) end.

%% The principal square root of a real rational coefficient as {Coefficient, Radicand}: sqrt(C) = Coefficient *
%% sqrt(Radicand) with Radicand squarefree (1 when the root is rational); sqrt(-n) = i sqrt(n). false otherwise.
csqrt_exact({q, {0, _}, {0, _}}) -> {c(0), 1};
csqrt_exact({q, {N, D}, {0, _}}) ->
    {S, P} = squarefree(abs(N) * D),
    R = rat(S, D),
    case N > 0 of
        true -> {{q, R, {0, 1}}, P};
        false -> {{q, {0, 1}, R}, P}
    end;
csqrt_exact(_) -> false.

%% N = S^2 * P with P squarefree, by trial division up to 10^6; a larger factor stays in P (SY 3.1).
-spec squarefree(pos_integer()) -> {pos_integer(), pos_integer()}.
squarefree(N) when N >= 1 -> squarefree(N, 2, 1, 1).

squarefree(1, _F, S, P) -> {S, P};
squarefree(N, F, S, P) when F > ?TRIAL_LIMIT; F * F > N -> {S, P * N};
squarefree(N, F, S, P) ->
    case N rem F of
        0 ->
            {K, Rest} = strip(N, F, 0),
            squarefree(Rest, F + 1, S * pow(F, K div 2), P * pow(F, K rem 2));
        _ -> squarefree(N, F + 1, S, P)
    end.

strip(N, F, K) when N rem F =:= 0 -> strip(N div F, F, K + 1);
strip(N, _F, K) -> {K, N}.

pow(_, 0) -> 1;
pow(B, E) -> B * pow(B, E - 1).

%% The exact Q-th root of a positive rational, or false.
root_exact({N, D}, Q) when N > 0 ->
    case {iroot(N, Q), iroot(D, Q)} of
        {{ok, A}, {ok, B}} -> {ok, {A, B}};
        _ -> false
    end;
root_exact(_, _) -> false.

iroot(N, Q) ->
    A = round(math:pow(N, 1 / Q)),
    case [X || X <- [A - 1, A, A + 1], X > 0, pow(X, Q) =:= N] of
        [X | _] -> {ok, X};
        [] -> false
    end.

%% exp(i r pi) for r = k/12 with k in (-12, 12]: {Cos, Sin} as lists of {Coefficient, Radicand} pairs (radicand 1
%% means the number itself); false if r is not a multiple of 1/12.
pi_angle(R) ->
    case rmul(R, {12, 1}) of
        {K, 1} -> {trig(cos, K), trig(sin, K)};
        _ -> false
    end.

%% cos and sin of k pi/12 for k = 0..6 (first quadrant, SY 4.4)
quadrant(0) -> {[{{1, 1}, 1}], []};
quadrant(1) -> {[{{1, 4}, 6}, {{1, 4}, 2}], [{{1, 4}, 6}, {{-1, 4}, 2}]};
quadrant(2) -> {[{{1, 2}, 3}], [{{1, 2}, 1}]};
quadrant(3) -> {[{{1, 2}, 2}], [{{1, 2}, 2}]};
quadrant(4) -> {[{{1, 2}, 1}], [{{1, 2}, 3}]};
quadrant(5) -> {[{{1, 4}, 6}, {{-1, 4}, 2}], [{{1, 4}, 6}, {{1, 4}, 2}]};
quadrant(6) -> {[], [{{1, 1}, 1}]}.

trig(F, K) when K < 0 ->
    case F of
        cos -> trig(cos, -K);
        sin -> neg(trig(sin, -K))
    end;
trig(F, K) when K > 6 ->            % pi - x
    case F of
        cos -> neg(trig(cos, 12 - K));
        sin -> trig(sin, 12 - K)
    end;
trig(cos, K) -> element(1, quadrant(K));
trig(sin, K) -> element(2, quadrant(K)).

neg(L) -> [{rneg(C), P} || {C, P} <- L].
