%% Truncated series (decision record SY2, phase 2): an expression expanded in a variable x at 0, at a point or at
%% infinity, as a sum of terms c x^e log(x)^m, the exponents e rational, every coefficient c exact.
%%
%% - Coefficients are rational functions of the other symbols (pq_sym_ratfun). A part free of x that is not one
%%   (log 2, pi, Euler's gamma, Gamma(1/3), sin(a), ...) is an opaque atom: a placeholder symbol inside the
%%   coefficient, written back in the result.
%% - The functions: sums, products, powers (integer, rational, or free of x when the base starts at a constant),
%%   exp, log, sin, cos, sqrt; Bessel J_nu (nu rational) and Y_n (n an integer, with its logarithms), the Laguerre
%%   polynomials L_n^alpha, the hypergeometric function 2F1(a, b; c; z) and Gamma, the last at any rational point.
%% - A series is {Q, #{{K, M} => C}, Ord}: the terms C x^(K/Q) log(x)^M, known modulo x^(Ord/Q) (Ord an integer
%%   or inf). Each operation keeps the honest order; series/4 raises its working order until the result reaches the
%%   one asked for (cancellations lose order).
-module(pq_series).
-export([series/4, terms/4]).

-define(VARS, {?MODULE, vars}).
-define(ATOMS, {?MODULE, atoms}).

%% the series of E in X up to (not including) x^N: {ok, Sum} with Sum an expression, the terms sorted by exponent
%% then power of log(x). Opts: at => 0 (default) | A (an expression free of X) | infinity (then the result is in
%% powers of 1/x, written as negative powers of x).
-spec series(term(), atom(), integer(), map()) -> {ok, term()} | {error, term()}.
series(E, X, N, Opts) ->
    case terms(E, X, N, Opts) of
        {ok, {Ts, _}} -> {ok, to_sum(Ts, X)};
        Err -> Err
    end.

%% the terms [{Exponent, LogPower, Coefficient}] (Exponent an integer or {rat, P, Q}) and the order reached
-spec terms(term(), atom(), integer(), map()) -> {ok, {[{term(), non_neg_integer(), term()}], term()}} | {error, term()}.
terms(E0, X, N, Opts) ->
    At = maps:get(at, Opts, 0),
    {E, Inverted} = case At of
                        0 -> {E0, false};
                        infinity -> {subst(E0, X, {pow, {sym, X}, -1}), true};
                        A -> {subst(E0, X, {add, [{sym, X}, A]}), false}
                    end,
    put(?VARS, []),
    put(?ATOMS, #{}),
    try
        S = converge(E, X, N, N, 0),
        {Q, Map, Ord} = S,
        % at infinity, log(1/x) = -log(x): the sign (-1)^m
        Sign = fun(M, C) when Inverted, M rem 2 =:= 1 -> c_neg(C); (_, C) -> C end,
        Ts2 = [{exponent(K, Q, Inverted), M, c_to_ast(Sign(M, C))}
               || {{K, M}, C} <- lists:sort(fun({{K1, M1}, _}, {{K2, M2}, _}) -> {K1, M1} =< {K2, M2} end, maps:to_list(Map)),
                  K < N * Q, not c_zero(C)],
        {ok, {Ts2, order_ast(Ord, Q, Inverted)}}
    catch
        throw:{series, Why} -> {error, {symbolic_unsupported, #{series => Why}}}
    after
        erase(?VARS),
        erase(?ATOMS)
    end.

%% raise the working order W until the result is known to x^N
converge(_, _, N, W, 6) -> throw({series, {order_not_reached, N, W}});
converge(E, X, N, W, Tries) ->
    S = ser(E, X, W),
    {Q, _, Ord} = S,
    case Ord =:= inf orelse Ord >= N * Q of
        true -> S;
        false -> converge(E, X, N, W + max(1, N - Ord div Q) + 1, Tries + 1)
    end.

exponent(K, Q, Inv) ->
    {A, B} = qn({case Inv of true -> -K; false -> K end, Q}),
    case B of 1 -> A; _ -> {rat, A, B} end.

order_ast(inf, _, _) -> inf;
order_ast(Ord, Q, Inv) -> exponent(Ord, Q, Inv).

to_sum([], _) -> 0;
to_sum(Ts, X) ->
    {add, [{mul, [C] ++ [{pow, {sym, X}, Ex} || Ex =/= 0] ++ [{pow, {log, {sym, X}}, M} || M > 0]} || {Ex, M, C} <- Ts]}.

subst({sym, X}, X, R) -> R;
subst(T, X, R) when is_tuple(T) -> list_to_tuple([subst(E, X, R) || E <- tuple_to_list(T)]);
subst(L, X, R) when is_list(L) -> [subst(E, X, R) || E <- L];
subst(T, _, _) -> T.

%%% Series of the expression nodes, at working order W (x^W)

ser(E, X, W) ->
    case free_of(E, X) of
        true -> const(c_of(E));
        false -> ser_(E, X, W)
    end.

ser_({sym, X}, X, _) -> {1, #{{1, 0} => c_num(1)}, inf};
ser_({add, Ts}, X, W) -> lists:foldl(fun(T, Acc) -> s_add(Acc, ser(T, X, W)) end, zero(), Ts);
ser_({mul, Fs}, X, W) -> lists:foldl(fun(F, Acc) -> s_mul(Acc, ser(F, X, W), W) end, const(c_num(1)), Fs);
ser_({pow, B, K}, X, W) when is_integer(K), K >= 0 -> s_ipow(base(B, X, W, fun(V) -> W + max(0, 1 - K) * V end), K, W);
ser_({pow, B, K}, X, W) when is_integer(K) ->
    % B^-K: 1/B^|K| to x^W needs B^|K| = c x^(|K| v) (1 + U) to x^(W + 2 |K| v), v B's valuation
    {Q, _, _} = Bs = base(B, X, W, fun(V) -> W + (1 - K) * V end),
    s_inv(s_ipow(Bs, -K, W + 2 * (-K) * max(0, ceil_div(val(Bs), Q)) + 1), W);
ser_({pow, B, {rat, P, Q}}, X, W) -> s_pow(base(B, X, W, fun(V) -> W + max(0, 1 - P / Q) * V end), c_num({P, Q}), {P, Q}, W);
ser_({pow, B, Ex}, X, W) ->
    case free_of(Ex, X) of
        true -> s_pow(base(B, X, W, fun(V) -> W + V end), c_of(Ex), none, W);
        false -> s_exp(s_mul(ser(Ex, X, W), s_log(base(B, X, W, fun(V) -> W + V end), W), W), W)   % B^E = exp(E log B)
    end;
ser_({sqrt, B}, X, W) -> s_pow(base(B, X, W, fun(V) -> W + V / 2 end), c_num({1, 2}), {1, 2}, W);
ser_({exp, A}, X, W) -> s_exp(ser(A, X, W), W);
ser_({log, A}, X, W) -> s_log(base(A, X, W, fun(V) -> W + V end), W);
ser_({sin, A}, X, W) -> s_sincos(sin, ser(A, X, W), W);
ser_({cos, A}, X, W) -> s_sincos(cos, ser(A, X, W), W);
ser_({besselj, Nu, A}, X, W) -> s_besselj(rational(Nu), base(A, X, W, fun(V) -> W + V end), W);
ser_({bessely, N, A}, X, W) when is_integer(N) -> s_bessely(N, base(A, X, W, fun(V) -> W + (abs(N) + 1) * V end), W);
ser_({laguerre, N, Alpha, A}, X, W) when is_integer(N), N >= 0 -> s_laguerre(N, c_of(Alpha), ser(A, X, W), W);
ser_({hyp2f1, A, B, C, Z}, X, W) -> s_hyp2f1(c_of(A), c_of(B), c_of(C), ser(Z, X, W), W);
ser_({gamma, A}, X, W) -> s_gamma(base(A, X, W, fun(V) -> W + V end), W);
ser_(E, _, _) -> throw({series, {unsupported, kind(E)}}).

%% the series of an operand to the order its use needs: Need(V) given its valuation V (in powers of x), found from a
%% first expansion at W (raised while that has no terms: a small operand truncated away)
base(B, X, W, Need) -> base(B, X, W, Need, W, 0).

base(_, _, _, _, _, 8) -> throw({series, zero_operand});
base(B, X, W, Need, Wt, Tries) ->
    S = ser(B, X, Wt),
    {Q, M, _} = S,
    case map_size(M) of
        0 -> base(B, X, W, Need, 2 * Wt + 2, Tries + 1);
        _ ->
            W2 = max(W, ceil(Need(val(S) / Q))),
            case W2 > Wt of
                true -> ser(B, X, W2);
                false -> S
            end
    end.

kind(T) when is_tuple(T) -> element(1, T);
kind(T) -> T.

rational(K) when is_integer(K) -> {K, 1};
rational({rat, P, Q}) -> qn({P, Q});
rational(_) -> throw({series, {order_must_be_rational, besselj}}).

free_of({sym, X}, X) -> false;
free_of(T, X) when is_tuple(T) -> lists:all(fun(E) -> free_of(E, X) end, tuple_to_list(T));
free_of(L, X) when is_list(L) -> lists:all(fun(E) -> free_of(E, X) end, L);
free_of(_, _) -> true.

%%% Series arithmetic

zero() -> {1, #{}, inf}.
const(C) -> case c_zero(C) of true -> zero(); false -> {1, #{{0, 0} => C}, inf} end.

rescale({Q, M, O}, Q2) when Q =:= Q2 -> {Q, M, O};
rescale({Q, M, O}, Q2) ->
    F = Q2 div Q,
    {Q2, maps:from_list([{{K * F, L}, C} || {K, L} := C <- M]), case O of inf -> inf; _ -> O * F end}.

common(A, B) ->
    {Qa, _, _} = A, {Qb, _, _} = B,
    Q = Qa * Qb div pq_sym_exact:gcd(Qa, Qb),
    {rescale(A, Q), rescale(B, Q)}.

omin(inf, B) -> B;
omin(A, inf) -> A;
omin(A, B) -> min(A, B).

oadd(inf, _) -> inf;
oadd(_, inf) -> inf;
oadd(A, B) -> A + B.

%% drop the terms at or beyond the order, and those beyond the working order W
trunc({Q, M, O}, W) ->
    O2 = omin(O, W * Q),
    {Q, maps:filter(fun({K, _}, C) -> K < O2 andalso not c_zero(C) end, M), O2}.

s_add(A0, B0) ->
    {{Q, A, Oa}, {Q, B, Ob}} = common(A0, B0),
    M = maps:fold(fun(K, V, Acc) ->
                          S = c_add(maps:get(K, Acc, c_num(0)), V),
                          case c_zero(S) of true -> maps:remove(K, Acc); false -> Acc#{K => S} end
                  end, A, B),
    O = omin(Oa, Ob),
    {Q, maps:filter(fun({K, _}, _) -> O =:= inf orelse K < O end, M), O}.

s_neg({Q, M, O}) -> {Q, maps:map(fun(_, C) -> c_neg(C) end, M), O}.
s_scale(C, {Q, M, O}) -> {Q, maps:filter(fun(_, V) -> not c_zero(V) end, maps:map(fun(_, V) -> c_mul(C, V) end, M)), O}.

%% the lowest exponent with a nonzero term (Ord when there is none)
val({_, M, O}) ->
    case maps:keys(M) of
        [] -> O;
        Ks -> lists:min([K || {K, _} <- Ks])
    end.

s_mul(A0, B0, W) ->
    {{Q, A, Oa} = SA, {Q, B, Ob} = SB} = common(A0, B0),
    O = omin(oadd(Oa, val(SB)), oadd(Ob, val(SA))),
    Lim = omin(O, W * Q),
    M = lists:foldl(fun({K, V}, Acc) ->
                            S = c_add(maps:get(K, Acc, c_num(0)), V),
                            case c_zero(S) of true -> maps:remove(K, Acc); false -> Acc#{K => S} end
                    end, #{},
                    [{{Ka + Kb, Ma + Mb}, c_mul(Ca, Cb)} || {Ka, Ma} := Ca <- A, {Kb, Mb} := Cb <- B,
                                                            Lim =:= inf orelse Ka + Kb < Lim]),
    trunc({Q, M, O}, W).

s_ipow(_, 0, _) -> const(c_num(1));
s_ipow(A, 1, _) -> A;
s_ipow(A, K, W) ->
    H = s_ipow(A, K div 2, W),
    H2 = s_mul(H, H, W),
    case K rem 2 of 0 -> H2; 1 -> s_mul(H2, A, W) end.

%% A = C0 x^(K0/Q) (1 + U): {C0, K0, U}; the lowest exponent must carry no power of log(x)
lead(A, What) ->
    {Q, M, O} = A,
    case maps:keys(M) of
        [] -> throw({series, {zero_or_unknown_leading_term, What}});
        _ ->
            K0 = val(A),
            case [L || {K, L} <- maps:keys(M), K =:= K0] of
                [0] ->
                    C0 = maps:get({K0, 0}, M),
                    Ci = c_inv(C0),
                    U = {Q, maps:from_list([{{K - K0, L}, c_mul(Ci, C)} || {K, L} := C <- M, {K, L} =/= {K0, 0}]),
                         case O of inf -> inf; _ -> O - K0 end},
                    {C0, K0, U};
                _ -> throw({series, {logarithmic_leading_term, What}})
            end
    end.

%% sum_j coef(j) U^j for U of positive valuation, to the working order (U's exponents shifted by Shift in the
%% caller: the sum is needed to x^(W - Shift/Q))
compose(Coef, U, W, Shift) ->
    {Q, _, _} = U,
    Need = W * Q - Shift,   % in units of 1/Q
    V = val(U),
    V > 0 orelse V =:= inf orelse throw({series, {not_small, V}}),
    Wl = ceil_div(Need, Q) + 1,
    % in U's units throughout, so the caller's shift (in them) applies even when U is zero
    compose(Coef, U, Wl, Need, V, 0, rescale(const(Coef(0)), Q), rescale(const(c_num(1)), Q)).

compose(Coef, U, Wl, Need, V, J, Acc, Pow) ->
    Pow2 = s_mul(Pow, U, Wl),
    J2 = J + 1,
    case V =:= inf orelse J2 * V >= Need of
        true ->
            % the sum is complete below x^(Need/Q) only: the terms there and beyond are dropped
            {Q, M, O} = Acc,
            O2 = case V of inf -> O; _ -> omin(O, Need) end,
            {Q, maps:filter(fun({K, _}, _) -> O2 =:= inf orelse K < O2 end, M), O2};
        false -> compose(Coef, U, Wl, Need, V, J2, s_add(Acc, s_scale(Coef(J2), Pow2)), Pow2)
    end.

ceil_div(A, B) when A >= 0 -> (A + B - 1) div B;
ceil_div(A, B) -> -((-A) div B).

xpow({Q, M, O}, K0) -> {Q, maps:from_list([{{K + K0, L}, C} || {K, L} := C <- M]), oadd(O, K0)}.

s_inv(A, W) ->
    {C0, K0, U} = lead(A, inverse),
    S = compose(fun(J) -> c_num(case J rem 2 of 0 -> 1; 1 -> -1 end) end, U, W, -K0),
    trunc(xpow(s_scale(c_inv(C0), S), -K0), W).

%% A^R, R a coefficient (Rq its rational value when it has one): C0^R x^(K0 R) (1 + U)^R, (1 + U)^R by the binomial
%% series
s_pow(A, R, Rq, W) ->
    {C0, K0, U} = lead(A, power),
    {Q, _, _} = A,
    % the exponent K0/Q R as a multiple of 1/(Q D)
    {Shift, Ua} = case {K0, Rq} of
                      {0, _} -> {0, U};
                      {_, none} -> throw({series, {symbolic_power_of_a_vanishing_base, K0}});
                      {_, {P, D}} -> {K0 * P, rescale(U, Q * D)}
                  end,
    S = compose(fun(J) -> binom_c(R, J) end, Ua, W, Shift),
    trunc(xpow(s_scale(c_pow_atom(C0, R, Rq), S), Shift), W).

%% binom(R, J) = R (R - 1) .. (R - J + 1) / J!
binom_c(R, J) ->
    lists:foldl(fun(I, Acc) -> c_mul(Acc, c_mul(c_add(R, c_num(-(I - 1))), c_num({1, I}))) end, c_num(1), lists:seq(1, J)).

s_exp(A, W) ->
    {A0, U} = split_const(A, exp),
    S = compose(fun(J) -> c_num({1, fact(J)}) end, U, W, 0),
    trunc(s_scale(c_exp_atom(A0), S), W).

%% A = A0 + U, A0 the constant term, U of positive valuation (no negative powers, no log(x) at x^0)
split_const({Q, M, O}, What) ->
    lists:all(fun({K, L}) -> K > 0 orelse {K, L} =:= {0, 0} end, maps:keys(M))
        orelse throw({series, {singular_argument, What}}),
    A0 = maps:get({0, 0}, M, c_num(0)),
    {A0, {Q, maps:remove({0, 0}, M), O}}.

s_log(A, W) ->
    {C0, K0, U} = lead(A, log),
    {Q, _, _} = A,
    S = compose(fun(0) -> c_num(0); (J) -> c_num({case J rem 2 of 1 -> 1; 0 -> -1 end, J}) end, U, W, 0),
    LogC = case c_rational(C0) of
               {ok, {1, 1}} -> zero();
               _ -> const(c_atom({log, c_to_ast(C0)}))
           end,
    LogX = case K0 of 0 -> zero(); _ -> {Q, #{{0, 1} => c_num({K0, Q})}, inf} end,
    trunc(s_add(s_add(LogC, LogX), S), W).

s_sincos(F, A, W) ->
    {A0, U} = split_const(A, F),
    CosU = compose(fun(J) when J rem 2 =:= 0 -> c_num({case (J div 2) rem 2 of 0 -> 1; 1 -> -1 end, fact(J)});
                      (_) -> c_num(0) end, U, W, 0),
    SinU = compose(fun(J) when J rem 2 =:= 1 -> c_num({case (J div 2) rem 2 of 0 -> 1; 1 -> -1 end, fact(J)});
                      (_) -> c_num(0) end, U, W, 0),
    {SinA, CosA} = case c_zero(A0) of
                       true -> {c_num(0), c_num(1)};
                       false -> {c_atom({sin, c_to_ast(A0)}), c_atom({cos, c_to_ast(A0)})}
                   end,
    case F of
        sin -> s_add(s_scale(SinA, CosU), s_scale(CosA, SinU));
        cos -> s_add(s_scale(CosA, CosU), s_neg(s_scale(SinA, SinU)))
    end.

%%% Special functions

%% J_nu(A) = (A/2)^nu sum_k (-1)^k (A/2)^(2k) / (k! Gamma(k + nu + 1)), nu rational; J_(-n) = (-1)^n J_n
s_besselj({P, 1}, A, W) when P < 0 ->
    S = s_besselj({-P, 1}, A, W),
    case P rem 2 of 0 -> S; _ -> s_neg(S) end;
s_besselj({P, D} = Nu, A, W) ->
    val(A) > 0 orelse throw({series, {argument_not_small, besselj}}),
    H = s_scale(c_num({1, 2}), A),
    H2 = s_mul(H, H, W),
    {Gamma1, Rising} = gamma_shift(Nu),   % Gamma(k + nu + 1) = Gamma1 * Rising(k)
    Sum = compose(fun(K) -> c_mul(c_num({case K rem 2 of 0 -> 1; 1 -> -1 end, fact(K)}), c_inv(Rising(K))) end,
                  H2, W, nu_shift(A, Nu)),
    Pre = case Nu of
              {0, 1} -> const(c_num(1));
              _ -> s_pow(H, c_num({P, D}), {P, D}, W)
          end,
    trunc(s_scale(c_inv(Gamma1), s_mul(Pre, Sum, W)), W).

%% the shift for compose: the prefactor (A/2)^nu raises every exponent by nu times A's valuation
nu_shift(A, {P, D}) -> floor_div(val(A) * P, D).

floor_div(A, B) when A >= 0 -> A div B;
floor_div(A, B) -> -((-A + B - 1) div B).

%% Gamma(k + nu + 1) = G1 * R(k): G1 = Gamma(nu + 1) (exact for an integer nu >= 0, else an atom), R(k) the rising
%% factorial (nu + 1)_k
gamma_shift({P, D} = Nu) ->
    G1 = case D of
             1 when P >= 0 -> c_num(fact(P));
             _ -> c_atom({gamma, q_ast(qadd(Nu, {1, 1}))})
         end,
    {G1, fun(K) -> lists:foldl(fun(I, Acc) -> c_mul(Acc, c_num(qadd(Nu, {I, 1}))) end, c_num(1), lists:seq(1, K)) end}.

%% Y_n(A), n an integer:
%%   Y_n = (2/pi) J_n log(A/2) - (1/pi) sum_(k<n) (n-k-1)!/k! (A/2)^(2k-n)
%%         - (1/pi) sum_k (psi(k+1) + psi(n+k+1)) (-1)^k (A/2)^(2k+n) / (k! (n+k)!),  psi(m+1) = H_m - gamma
s_bessely(N, A, W) when N < 0 ->
    S = s_bessely(-N, A, W),
    case N rem 2 of 0 -> S; _ -> s_neg(S) end;
s_bessely(N, A, W) ->
    val(A) > 0 orelse throw({series, {argument_not_small, bessely}}),
    InvPi = c_inv(c_atom({const, pi})),
    EG = c_atom({const, euler_gamma}),
    H = s_scale(c_num({1, 2}), A),
    W2 = W + N + 2,   % the negative powers of the finite sum need the others to this much more
    J = s_besselj({N, 1}, A, W2),
    T1 = s_scale(c_mul(c_num(2), InvPi), s_mul(J, s_log(H, W2), W2)),
    Hinv = s_inv(H, W2),
    T2 = lists:foldl(fun(K, Acc) ->
                             s_add(Acc, s_scale(c_num({fact(N - K - 1), fact(K)}), s_ipow_signed(H, Hinv, 2 * K - N, W2)))
                     end, zero(), lists:seq(0, N - 1)),
    H2 = s_mul(H, H, W2),
    Psi = fun(M) -> c_add(c_num(harmonic(M)), c_neg(EG)) end,
    Sum = compose(fun(K) -> c_mul(c_add(Psi(K), Psi(N + K)), c_num({case K rem 2 of 0 -> 1; 1 -> -1 end, fact(K) * fact(N + K)})) end,
                  H2, W2, val(A) * N),
    T3 = s_mul(s_ipow(H, N, W2), Sum, W2),
    trunc(s_add(T1, s_neg(s_scale(InvPi, s_add(T2, T3)))), W).

s_ipow_signed(H, _, K, W) when K >= 0 -> s_ipow(H, K, W);
s_ipow_signed(_, Hinv, K, W) -> s_ipow(Hinv, -K, W).

%% L_n^alpha(A) = sum_(k=0..n) (-1)^k binom(n + alpha, n - k) A^k / k!, binom(n + alpha, n - k) =
%% (alpha + k + 1) .. (alpha + n) / (n - k)!
s_laguerre(N, Alpha, A, W) ->
    lists:foldl(fun(K, Acc) ->
                        B = lists:foldl(fun(I, P) -> c_mul(P, c_add(Alpha, c_num(I))) end, c_num(1), lists:seq(K + 1, N)),
                        C = c_mul(B, c_num({case K rem 2 of 0 -> 1; 1 -> -1 end, fact(N - K) * fact(K)})),
                        s_add(Acc, s_scale(C, s_ipow(A, K, W)))
                end, zero(), lists:seq(0, N)).

%% 2F1(a, b; c; Z) = sum_j (a)_j (b)_j / ((c)_j j!) Z^j, Z of positive valuation
s_hyp2f1(A, B, C, Z, W) ->
    Coef = fun(J) ->
                   lists:foldl(fun(I, Acc) ->
                                       Num = c_mul(c_add(A, c_num(I)), c_add(B, c_num(I))),
                                       c_mul(Acc, c_mul(Num, c_inv(c_mul(c_add(C, c_num(I)), c_num(I + 1)))))
                               end, c_num(1), lists:seq(0, J - 1))
           end,
    trunc(compose(Coef, Z, W, 0), W).

%% Gamma(a0 + U), a0 rational:
%%   a0 = n >= 1: Gamma(n) exp(sum_k psi^(k-1)(n) U^k / k!), psi(n) = H_(n-1) - gamma,
%%                psi^(k)(n) = (-1)^(k+1) k! (zeta(k+1) - H_(n-1)^(k+1));
%%   a0 = -m <= 0: Gamma(1 + (U - m - 1)) ... = Gamma(U + 1) / (U (U - 1) .. (U - m));
%%   otherwise: Gamma(a0) exp(sum_k psi^(k-1)(a0) U^k / k!) with Gamma(a0), psi^(k)(a0) atoms
s_gamma(A, W) ->
    {A0, U} = split_const(A, gamma),
    case c_rational(A0) of
        {ok, {Nn, 1}} when Nn >= 1 ->
            EG = c_atom({const, euler_gamma}),
            Coef = fun(1) -> c_add(c_num(harmonic(Nn - 1)), c_neg(EG));
                      (K) ->
                           Z = c_atom({zeta, K}),
                           Hk = c_num(harmonic(Nn - 1, K)),
                           c_mul(c_num({case K rem 2 of 0 -> 1; 1 -> -1 end * fact(K - 1), fact(K)}), c_add(Z, c_neg(Hk)))
                   end,
            Exp = compose(fun(0) -> c_num(0); (K) -> Coef(K) end, U, W, 0),
            trunc(s_scale(c_num(fact(Nn - 1)), s_exp(Exp, W)), W);
        {ok, {Mn, 1}} when Mn =< 0 ->
            M = -Mn,
            % Gamma(U - m) = Gamma(U + 1) / prod_(i=0..m) (U - i): U has positive valuation, so the product starts
            % at U and the inverse has negative powers
            G1 = s_gamma(s_add(U, const(c_num(1))), W + M + 2),
            Den = lists:foldl(fun(I, Acc) -> s_mul(Acc, s_add(U, const(c_num(-I))), W + M + 2) end, const(c_num(1)), lists:seq(0, M)),
            trunc(s_mul(G1, s_inv(Den, W + M + 2), W + M + 2), W);
        _ ->
            A0ast = c_to_ast(A0),
            Coef = fun(K) -> c_mul(c_atom({polygamma, K - 1, A0ast}), c_num({1, fact(K)})) end,
            Exp = compose(fun(0) -> c_num(0); (K) -> Coef(K) end, U, W, 0),
            trunc(s_scale(c_atom({gamma, A0ast}), s_exp(Exp, W)), W)
    end.

harmonic(M) -> harmonic(M, 1).
harmonic(M, P) -> lists:foldl(fun(I, Acc) -> qadd(Acc, {1, pow(I, P)}) end, {0, 1}, lists:seq(1, M)).

%%% Coefficients: rational functions over a growing list of symbols (the parameters, then atoms' placeholders)

c_num(K) when is_integer(K) -> c_num({K, 1});
c_num({P, Q}) ->
    {A, B} = qn({P, Q}),
    N = length(get(?VARS)),
    Z = lists:duplicate(N, 0),
    {case A of 0 -> #{}; _ -> #{Z => case B of 1 -> A; _ -> {A, B} end} end, #{Z => 1}}.

pad({Nu, De}) ->
    L = length(get(?VARS)),
    P = fun(M) -> maps:from_list([{K ++ lists:duplicate(L - length(K), 0), V} || K := V <- M]) end,
    {P(Nu), P(De)}.

c_add(A, B) -> pq_sym_ratfun:add(pad(A), pad(B)).
c_mul(A, B) -> pq_sym_ratfun:mul(pad(A), pad(B)).
c_neg(A) -> c_mul(c_num(-1), A).
c_inv(A) ->
    c_zero(A) andalso throw({series, division_by_zero}),
    pq_sym_ratfun:inv(pad(A)).
c_zero({Nu, _}) -> map_size(Nu) =:= 0.

c_rational(C) ->
    {Nu, De} = pad(C),
    L = length(get(?VARS)),
    Z = lists:duplicate(L, 0),
    case {maps:keys(Nu), De} of
        {[], _} -> {ok, {0, 1}};
        {[Z], #{Z := 1}} -> {ok, case maps:get(Z, Nu) of {A, B} -> {A, B}; A -> {A, 1} end};
        _ -> false
    end.

%% a symbol of the coefficients (a parameter, or an atom's placeholder)
c_var(Name) ->
    Vars = get(?VARS),
    Vars2 = case lists:member(Name, Vars) of true -> Vars; false -> Vars ++ [Name] end,
    put(?VARS, Vars2),
    I = index(Name, Vars2),
    L = length(Vars2),
    {#{[case J of I -> 1; _ -> 0 end || J <- lists:seq(1, L)] => 1}, #{lists:duplicate(L, 0) => 1}}.

index(X, L) -> length(lists:takewhile(fun(Y) -> Y =/= X end, L)) + 1.

%% an opaque atom: its placeholder symbol
c_atom(Ast) ->
    Atoms = get(?ATOMS),
    Name = case maps:find(Ast, Atoms) of
               {ok, Nm} -> Nm;
               error ->
                   Nm = list_to_atom("$series_atom_" ++ integer_to_list(map_size(Atoms) + 1)),
                   put(?ATOMS, Atoms#{Ast => Nm}),
                   Nm
           end,
    c_var(Name).

%% an expression free of x as a coefficient: rational operations exactly, anything else an atom
c_of(K) when is_integer(K) -> c_num(K);
c_of({rat, P, Q}) -> c_num({P, Q});
c_of({sym, Y}) -> c_var(Y);
c_of({add, Ts}) -> lists:foldl(fun(T, Acc) -> c_add(Acc, c_of(T)) end, c_num(0), Ts);
c_of({mul, Fs}) -> lists:foldl(fun(F, Acc) -> c_mul(Acc, c_of(F)) end, c_num(1), Fs);
c_of({pow, B, K}) when is_integer(K), K >= 0 -> lists:foldl(fun(_, Acc) -> c_mul(Acc, c_of(B)) end, c_num(1), lists:seq(1, K));
c_of({pow, B, K}) when is_integer(K) -> c_inv(c_of({pow, B, -K}));
c_of(E) -> c_atom(E).

c_exp_atom(A0) ->
    case c_zero(A0) of
        true -> c_num(1);
        false -> c_atom({exp, c_to_ast(A0)})
    end.

%% C0^R: exact for an integer R, or a rational C0 with an exact root; otherwise an atom
c_pow_atom(C0, R, Rq) ->
    case {Rq, c_rational(C0)} of
        {{P, 1}, _} when P >= 0 -> lists:foldl(fun(_, Acc) -> c_mul(Acc, C0) end, c_num(1), lists:seq(1, P));
        {{P, 1}, _} -> c_inv(c_pow_atom(C0, R, {-P, 1}));
        {_, {ok, {1, 1}}} -> c_num(1);
        {{P, D}, {ok, {A, B}}} when A > 0 ->
            case {root(A, D), root(B, D)} of
                {{ok, Ra}, {ok, Rb}} -> c_pow_atom(c_num({Ra, Rb}), R, {P, 1});
                _ -> c_atom({pow, q_ast({A, B}), q_ast({P, D})})
            end;
        _ -> c_atom({pow, c_to_ast(C0), c_to_ast(R)})
    end.

root(A, D) ->
    R = round(math:pow(A, 1 / D)),
    case [X || X <- [R - 1, R, R + 1], X >= 0, pow(X, D) =:= A] of
        [X | _] -> {ok, X};
        [] -> false
    end.

%% a coefficient as an expression: canonical, the atoms written back
c_to_ast(C) ->
    Vars = get(?VARS),
    {ok, Can} = pq_symbolic:canonicalise(pq_sym_ratfun:to_ast(pad(C), Vars)),
    Back = maps:from_list([{Nm, Ast} || Ast := Nm <- get(?ATOMS)]),
    unatom(Can, Back).

unatom({sym, Nm} = S, Back) ->
    case maps:find(Nm, Back) of
        {ok, Ast} -> Ast;
        error -> S
    end;
unatom(T, Back) when is_tuple(T) -> list_to_tuple([unatom(E, Back) || E <- tuple_to_list(T)]);
unatom(L, Back) when is_list(L) -> [unatom(E, Back) || E <- L];
unatom(T, _) -> T.

%%% Rationals

qn({P, Q}) when Q < 0 -> qn({-P, -Q});
qn({0, _}) -> {0, 1};
qn({P, Q}) -> G = pq_sym_exact:gcd(abs(P), Q), {P div G, Q div G}.
qadd({A, B}, {C, D}) -> qn({A * D + C * B, B * D}).
q_ast({P, 1}) -> P;
q_ast({P, Q}) -> {rat, P, Q}.

fact(0) -> 1;
fact(N) -> N * fact(N - 1).

pow(_, 0) -> 1;
pow(B, E) -> B * pow(B, E - 1).
