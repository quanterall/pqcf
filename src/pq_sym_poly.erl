%% Canonical sums of monomials (SY 3.3, 4.1, 4.4, 4.5).
%%
%%   poly()  :: #{mono() => coeff()}             coefficients never zero
%%   mono()  :: {Alg, Exp}
%%   Alg     :: [{atom_key(), rational()}]       sorted by key, exponents never zero
%%   Exp     :: [{mono(), coeff()}]              exp(Coeff * Base) factors, sorted by base, one per base
%%   atom_key() :: {sym, atom()} | pi | {rad, P} | {conj, atom()} | {base, poly()} | {cbase, poly()}
%%
%% {cbase, B}^e is conj(B^e) for a fractional e. The rule conj(B^e) = conj(B)^e of SY 4.2 fails on the branch
%% cut: for B = -8, conj((-8)^(1/3)) = 1 - sqrt(3) i while conj(-8)^(1/3) = 1 + sqrt(3) i. An integer exponent
%% turns it back into conj(B)^e.
%%
%% Products keep at most one radical atom (exponent 1), merge equal exponential bases, reduce exp(c pi) into
%% (-1, 1] times i pi and evaluate it exactly when 12 r is an integer, and expand Base(B)^k for small positive k.
%% Every limit is checked before the work: a product needs terms(A) * terms(B) =< max_poly_terms.
-module(pq_sym_poly).
-export([zero/0, one/0, const/1, symbol/2, atom/2, pi/0, is_const/1, const_value/1, single/1, m/1,
         add/2, sum/1, neg/1, scale/2, mul/3, pow_int/3, pow_rat/3, exp/2, conj/2,
         terms/1, pi_mono/0, one_mono/0, sorted/1, content/2, base_power/3, plain_terms/1]).

-import(pq_sym_exact, [c/1, cadd/2, cmul/2, cneg/1, cinv/1, cconj/1, is_zero/1, rat/2, radd/2, rmul/2, rneg/1]).

-define(ONE, {[], []}).
-define(PI, {[{pi, {1, 1}}], []}).
-define(RP, {rpoly, _, _}).
-define(IS_RP(X), (is_tuple(X) andalso element(1, X) =:= rpoly)).

one_mono() -> ?ONE.
pi_mono() -> ?PI.

zero() -> #{}.
one() -> #{?ONE => c(1)}.

const(C) ->
    case is_zero(C) of
        true -> #{};
        false -> #{?ONE => C}
    end.

pi() -> #{?PI => c(1)}.

%% x^1, or conj(z) as an atom
symbol(X, _Ctx) -> #{{[{{sym, X}, {1, 1}}], []} => c(1)}.
atom(Key, E) -> #{{[{Key, E}], []} => c(1)}.

is_const(?RP = H) -> case rp_terms(H) of 0 -> true; 1 -> is_const(m(H)); _ -> false end;
is_const(P) -> maps:size(P) =:= 0 orelse (maps:size(P) =:= 1 andalso maps:is_key(?ONE, P)).

const_value(?RP = H) -> const_value(m(H));
const_value(P) -> maps:get(?ONE, P, c(0)).

single(?RP = H) -> case rp_terms(H) of 1 -> single(m(H)); _ -> false end;
single(P) when map_size(P) =:= 1 -> [{M, C}] = maps:to_list(P), {ok, M, C};
single(_) -> false.

terms(?RP = H) -> rp_terms(H);
terms(P) -> maps:size(P).

%% the terms in canonical key order
sorted(P) -> lists:sort(maps:to_list(m(P))).

add(P, Q) when ?IS_RP(P); ?IS_RP(Q) ->
    case rp_add(P, Q) of
        {ok, R} -> R;
        false -> add(m(P), m(Q))
    end;
add(P, Q) when map_size(P) < map_size(Q) -> add(Q, P);
add(P, Q) -> maps:fold(fun(M, C, Acc) -> add_term(Acc, M, C) end, P, Q).

add_term(Acc, M, C) ->
    case Acc of
        #{M := C0} ->
            S = cadd(C0, C),
            case is_zero(S) of
                true -> maps:remove(M, Acc);
                false -> Acc#{M := S}
            end;
        _ ->
            case is_zero(C) of
                true -> Acc;
                false -> Acc#{M => C}
            end
    end.

sum(Ps) -> lists:foldl(fun add/2, #{}, Ps).

neg(?RP = H) -> scale(H, c(-1));
neg(P) -> maps:map(fun(_, C) -> cneg(C) end, P).

scale(?RP = H, C) ->
    case rp_scale(H, C) of
        {ok, R} -> R;
        false -> scale(m(H), C)
    end;
scale(P, C) ->
    case is_zero(C) of
        true -> #{};
        false -> maps:filter(fun(_, V) -> not is_zero(V) end, maps:map(fun(_, V) -> cmul(V, C) end, P))
    end.

%%% Products

mul(P, Q, Ctx) ->
    case rp_mul(P, Q, Ctx) of
        {ok, R} -> R;
        false -> mul_maps(m(P), m(Q), Ctx)
    end.

mul_maps(P, Q, Ctx) ->
    check_terms(maps:size(P) * maps:size(Q), Ctx),
    case kernel_mul(P, Q, Ctx) of
        {ok, R} -> R;
        false -> mul_terms(P, Q, Ctx)
    end.

mul_terms(P, Q, Ctx) ->
    maps:fold(
      fun(M1, C1, Acc1) ->
              maps:fold(
                fun(M2, C2, Acc2) ->
                        {CF, M, Extras} = mono_mul(M1, M2, Ctx),
                        C = cmul(cmul(C1, C2), CF),
                        case Extras of
                            [] -> add_term(Acc2, M, C);
                            _ ->
                                Term = lists:foldl(fun(X, T) -> mul(T, X, Ctx) end, #{M => C}, Extras),
                                add(Acc2, Term)
                        end
                end, Acc1, Q)
      end, #{}, P).

%% Plain polynomials (symbols, conj and pi with integer exponents, exact coefficients) multiply in Rust when both
%% have at least two terms and the product at least ?KERNEL_MIN term products (pq_nif:sym_poly_mul,
%% native/pq_nif/src/polymul.rs; decision record K1): the denominators scaled away, the exponents per variable, the
%% coefficients exact in i128. A monomial times a polynomial stays here: the conversion costs more than the product
%% (1 x 4096 terms 0.8x, 2 x 4096 1.7x). An overflow, exponents that don't pack into 64 bits, or no NIF loaded:
%% false, and the product is computed term by term here -- the same polynomial either way (D1). poly_kernel => false
%% in the context turns the kernel off.
-define(KERNEL_MIN, 64).

kernel_mul(P, Q, Ctx) when map_size(P) >= 2, map_size(Q) >= 2, map_size(P) * map_size(Q) >= ?KERNEL_MIN ->
    case maps:get(poly_kernel, Ctx, true) andalso kernel_vars(P, kernel_vars(Q, #{})) of
        false -> false;
        VarSet ->
            Vars = lists:sort(maps:keys(VarSet)),
            case {kernel_in(P, Vars), kernel_in(Q, Vars)} of
                {{EP, CP, LP}, {EQ, CQ, LQ}} when LP * LQ < 1 bsl 127 ->
                    try pq_nif:sym_poly_mul(Vars, EP, CP, EQ, CQ, LP * LQ) of
                        {ok, R} -> {ok, R};
                        _ -> false
                    catch
                        error:nif_not_loaded -> false;
                        error:undef -> false
                    end;
                _ -> false
            end
    end;
kernel_mul(_, _, _) -> false.

%% P^N for a plain polynomial of at least two terms in Rust (pq_nif:sym_poly_pow, polymul.rs; K1), when the result
%% may have at least ?KERNEL_MIN terms: repeated products in the packed layout of P^N's exponent ranges, exact in
%% i128; the coefficients over L^N. Refused (overflow, packing, L^N beyond i128, no NIF): false, and the
%% multinomial enumeration below -- the same polynomial.
kernel_pow(P, N, Ctx) ->
    Bound = binom(N + map_size(P) - 1, map_size(P) - 1),
    check_terms(Bound, Ctx),
    Flag = maps:get(rpoly, Ctx, false),
    case Bound >= ?KERNEL_MIN andalso maps:get(poly_kernel, Ctx, true) andalso kernel_vars(P, #{}) of
        false -> false;
        VarSet ->
            Vars = lists:sort(maps:keys(VarSet)),
            case kernel_in(P, Vars) of
                {E, C, L} ->
                    D = ipow(L, N),
                    case D < 1 bsl 127 of
                        true when Flag -> rp_pow(rp_from(P), N, Ctx);
                        true -> nif_pow(Vars, E, C, N, D);
                        false -> false
                    end;
                false -> false
            end
    end.

nif_pow(Vars, E, C, N, D) ->
    try pq_nif:sym_poly_pow(Vars, E, C, N, D) of
        {ok, R} -> {ok, R};
        _ -> false
    catch
        error:nif_not_loaded -> false;
        error:undef -> false
    end.

ipow(_, 0) -> 1;
ipow(B, K) -> B * ipow(B, K - 1).

%% the atoms of a plain polynomial with integer exponents and exact coefficients, or false. An exponential factor
%% exp(k i x) of a symbol x, k an integer, counts as the atom {expi, X} to the k: such factors only ever merge by
%% adding k (normalize/3 drops a zero one), exactly as an atom's exponents do.
kernel_vars(_P, false) -> false;
kernel_vars(P, Acc0) ->
    try maps:fold(fun({Alg, Exp}, {q, _, _}, Acc) ->
                          lists:foldl(fun({K, {_, 1}}, A) when K =:= pi; element(1, K) =:= sym; element(1, K) =:= conj;
                                                                 element(1, K) =:= expi ->
                                              A#{K => []};
                                         (_, _) -> throw(not_plain)
                                      end, Acc, Alg ++ expi(Exp));
                     (_, _, _) -> throw(not_plain)
                  end, Acc0, P)
    catch throw:not_plain -> false
    end.

expi(Exp) ->
    [case F of
         {{[{{sym, X}, {1, 1}}], []}, {q, {0, 1}, {K, 1}}} -> {{expi, X}, {K, 1}};
         _ -> throw(not_plain)
     end || F <- Exp].

%% {Exponents, Coefficients, L}: the coefficients times L, their common denominator, as integers; false if one
%% of them is outside i128
kernel_in(P, Vars) ->
    Terms = maps:to_list(P),
    L = lists:foldl(fun({_, {q, {_, D1}, {_, D2}}}, A) -> lcm(lcm(A, D1), D2) end, 1, Terms),
    Coef = lists:append([[N1 * (L div D1), N2 * (L div D2)] || {_, {q, {N1, D1}, {N2, D2}}} <- Terms]),
    case lists:all(fun(X) -> X >= -(1 bsl 127) andalso X < 1 bsl 127 end, Coef) of
        true ->
            Keys = fun(Alg, []) -> Alg; (Alg, Exp) -> lists:sort(Alg ++ expi(Exp)) end,
            {lists:append([exps(Vars, Keys(Alg, Exp)) || {{Alg, Exp}, _} <- Terms]), Coef, L};
        false -> false
    end.

exps([], _) -> [];
exps([V | Vs], [{V, {E, 1}} | Alg]) -> [E | exps(Vs, Alg)];
exps([_ | Vs], Alg) -> [0 | exps(Vs, Alg)].

lcm(A, B) -> A div pq_sym_exact:gcd(A, B) * B.

check_terms(N, Ctx) ->
    Max = maps:get(max_poly_terms, Ctx, 100000),
    N > Max andalso throw({symbolic_limit, #{reason => poly_terms, terms => N, max => Max}}),
    ok.

%% plain monomials (symbols, conj, pi; no radical, opaque base or exponential: 31x of the work of a product of
%% such polynomials went into normalize, comparison with FORM, 2 Nov 2026) only merge their exponents
mono_mul({A1, []}, {A2, []}, Ctx) ->
    case plain(A1) andalso plain(A2) of
        true -> {c(1), {[{K, E} || {K, E} <- merge(A1, A2, fun pq_sym_exact:radd/2), E =/= {0, 1}], []}, []};
        false -> normalize(merge(A1, A2, fun pq_sym_exact:radd/2), [], Ctx)
    end;
mono_mul({A1, E1}, {A2, E2}, Ctx) ->
    normalize(merge(A1, A2, fun pq_sym_exact:radd/2), merge(E1, E2, fun pq_sym_exact:cadd/2), Ctx).

plain(Alg) -> lists:all(fun({{K, _}, _}) -> K =:= sym orelse K =:= conj; ({pi, _}) -> true; (_) -> false end, Alg).

%% merge two key-sorted lists, combining equal keys
merge([], L, _) -> L;
merge(L, [], _) -> L;
merge([{K, V1} | T1], [{K, V2} | T2], F) -> [{K, F(V1, V2)} | merge(T1, T2, F)];
merge([{K1, _} = H1 | T1], [{K2, _} | _] = L2, F) when K1 < K2 -> [H1 | merge(T1, L2, F)];
merge(L1, [H2 | T2], F) -> [H2 | merge(L1, T2, F)].

%% {CoefficientFactor, Mono, ExtraPolys}: the monomial with its atoms and exponentials in canonical form; the
%% coefficient factor and the extra polynomials multiply it.
normalize(Alg0, Exp0, Ctx) ->
    Alg1 = [{K, E} || {K, E} <- Alg0, E =/= {0, 1}],
    {Rads, Alg2} = lists:partition(fun({{rad, _}, _}) -> true; (_) -> false end, Alg1),
    {CR, RadAtoms} = radicals(Rads),
    MaxExp = maps:get(max_expand_power, Ctx, 8),
    {Alg3, Extras1} = lists:foldr(
                        fun({{base, B}, {K, 1}} = A, {As, Xs}) ->
                                %% an integer power of an opaque base goes back through the rules: a monomial's
                                %% power, a small power of a sum expanded, or the sum made primitive (sign included)
                                case map_size(B) =:= 1 orelse (K > 0 andalso K =< MaxExp) of
                                    true -> {As, [pow_int(B, K, Ctx) | Xs]};
                                    false ->
                                        case content(B, true) of
                                            {{1, 1}, _} -> {[A | As], Xs};
                                            _ -> {As, [base_power(B, {K, 1}, Ctx) | Xs]}
                                        end
                                end;
                           ({{cbase, B}, {K, 1}}, {As, Xs}) ->
                                {As, [pow_int(conj(B, Ctx), K, Ctx) | Xs]};
                           (A, {As, Xs}) -> {[A | As], Xs}
                        end, {[], []}, Alg2),
    Alg = lists:sort(Alg3 ++ RadAtoms),
    {CE, Exp, Extras2} = exponentials(Exp0),
    {cmul(CR, CE), {Alg, Exp}, Extras1 ++ Extras2}.

%% R(p)^e = p^(e div 2) R(p)^(e rem 2), then the product of the remaining radicals: S R(squarefree part)
radicals([]) -> {c(1), []};
radicals(Rads) ->
    {C, Prod} = lists:foldl(
                  fun({{rad, P}, {E, 1}}, {C0, Pr}) ->
                          Half = floor_div(E, 2),
                          C1 = cmul(C0, pq_sym_exact:cpow_int(c(P), Half)),
                          case E - 2 * Half of
                              1 -> {C1, Pr * P};
                              0 -> {C1, Pr}
                          end
                  end, {c(1), 1}, Rads),
    {S, Q} = pq_sym_exact:squarefree(Prod),
    {cmul(C, c(S)), [{{rad, Q}, {1, 1}} || Q > 1]}.

floor_div(A, B) when A >= 0 -> A div B;
floor_div(A, B) -> -((-A + B - 1) div B).

%% exp(c pi): Im c reduced into (-1, 1] and evaluated exactly when 12 Im c is an integer; exp(float) folded into the
%% coefficient; zero exponents dropped.
exponentials(Exp0) ->
    lists:foldr(
      fun({_B, C}, Acc) when element(1, C) =:= q, element(2, C) =:= {0, 1}, element(3, C) =:= {0, 1} -> Acc;
         ({?ONE, {f, Re, Im}}, {CF, Es, Xs}) ->
              M = math:exp(Re),
              {cmul(CF, pq_sym_exact:cf(M * math:cos(Im), M * math:sin(Im))), Es, Xs};
         ({?PI, {q, _A, {0, 1}}} = F, {CF, Es, Xs}) -> {CF, [F | Es], Xs};    % exp(a pi), a real: kept
         ({?PI, {q, A, R}}, {CF, Es, Xs}) ->
              R1 = reduce_turn(R),
              case pq_sym_exact:pi_angle(R1) of
                  {Cos, Sin} ->
                      Val = add(trig_poly(Cos, c(1)), trig_poly(Sin, {q, {0, 1}, {1, 1}})),
                      Es1 = [{?PI, {q, A, {0, 1}}} || A =/= {0, 1}] ++ Es,
                      {CF, Es1, [Val | Xs]};
                  false -> {CF, [{?PI, {q, A, R1}} | Es], Xs}
              end;
         ({_B, C} = F, {CF, Es, Xs}) ->
              case is_zero(C) of
                  true -> {CF, Es, Xs};
                  false -> {CF, [F | Es], Xs}
              end
      end, {c(1), [], []}, Exp0).

%% r into (-1, 1]: r - 2 ceil((r - 1) / 2)
reduce_turn(R) ->
    K = pq_sym_exact:rceil(rmul(radd(R, {-1, 1}), {1, 2})),
    radd(R, rmul({-2, 1}, K)).

%% sum of Coefficient * sqrt(Radicand), times Unit
trig_poly(Parts, Unit) ->
    sum([scale(case P of 1 -> one(); _ -> atom({rad, P}, {1, 1}) end, cmul(Unit, c(Q))) || {Q, P} <- Parts]).

%%% Powers

pow_int(_P, 0, _Ctx) -> one();
pow_int(?RP = H, N, Ctx) ->
    case N > 0 andalso rp_pow(H, N, Ctx) of
        {ok, R} -> R;
        _ -> pow_int(m(H), N, Ctx)
    end;
pow_int(P, N, Ctx) when N < 0 ->
    case single(P) of
        {ok, M, C} -> pow_int(inverse(M, C, Ctx), -N, Ctx);
        false when map_size(P) =:= 0 -> throw({non_finite, division_by_zero});
        false -> base_power(P, {N, 1}, Ctx)
    end;
pow_int(P, N, Ctx) ->
    case single(P) of
        %% a plain monomial (symbols, conj, pi): its exponents times N and its coefficient to the N, directly --
        %% square-and-multiply did log2 N products of it (found by the comparison with FORM's published
        %% benchmarks, 5 Nov 2026: every term of FORM's output is a product of such powers)
        {ok, {A, []}, C} when A =/= [] ->
            case plain(A) of
                true -> #{{[{X, rmul(E, {N, 1})} || {X, E} <- A], []} => cpow(C, N)};
                false -> pow_int_sum(P, N, Ctx)
            end;
        _ -> pow_int_sum(P, N, Ctx)
    end.

pow_int_sum(P, N, Ctx) ->
    Max = maps:get(max_expand_power, Ctx, 8),
    case maps:size(P) > 1 andalso N > Max of
        true -> base_power(P, {N, 1}, Ctx);
        false ->
            case maps:size(P) > 1 andalso lists:all(fun({A, E}) -> E =:= [] andalso plain(A) end, maps:keys(P)) of
                true ->
                    case kernel_pow(P, N, Ctx) of
                        {ok, R} -> R;
                        false -> multinomial(P, N, Ctx)
                    end;
                false -> square_and_multiply(P, N, Ctx)
            end
    end.

%% (sum c_i M_i)^N = sum over k_1 + ... + k_m = N of N! / prod k_i! prod (c_i M_i)^k_i, for plain monomials: one
%% step per term of the result instead of squaring the half power (the size check is the result's bound,
%% C(N + m - 1, m - 1), not the square of the half power's terms)
multinomial(P, N, Ctx) ->
    Terms = maps:to_list(P),
    check_terms(binom(N + length(Terms) - 1, length(Terms) - 1), Ctx),
    Fact = fun F(0) -> 1; F(K) -> K * F(K - 1) end,
    FN = Fact(N),
    Pow = [fun(K) -> {cpow(C, K), [{X, rmul(E, {K, 1})} || {X, E} <- A]} end || {{A, []}, C} <- Terms],
    distribute(Pow, N, {FN, c(1), []}, Fact, #{}).

%% every split of the remaining exponent over the remaining terms; Acc = {multinomial numerator over the chosen
%% factorials so far, coefficient, exponents}
distribute([Last], K, {Num, C, A}, Fact, Map) ->
    {Ck, Ak} = Last(K),
    Coef = cmul(cmul(C, Ck), c(Num div Fact(K))),
    add_term(Map, {[{X, E} || {X, E} <- merge(A, Ak, fun pq_sym_exact:radd/2), E =/= {0, 1}], []}, Coef);
distribute([F | Rest], K, {Num, C, A}, Fact, Map) ->
    lists:foldl(fun(J, M) ->
                        {Cj, Aj} = F(J),
                        distribute(Rest, K - J, {Num div Fact(J), cmul(C, Cj), merge(A, Aj, fun pq_sym_exact:radd/2)},
                                   Fact, M)
                end, Map, lists:seq(0, K)).

cpow(_C, 0) -> c(1);
cpow(C, K) -> pq_sym_exact:cpow_int(C, K).

binom(N, K) -> binom(N, K, 1, 1).
binom(_N, 0, A, B) -> A div B;
binom(N, K, A, B) -> binom(N - 1, K - 1, A * N, B * K).

square_and_multiply(P, 1, _Ctx) -> P;
square_and_multiply(P, N, Ctx) ->
    H = square_and_multiply(P, N div 2, Ctx),
    H2 = mul(H, H, Ctx),
    case N rem 2 of 0 -> H2; 1 -> mul(H2, P, Ctx) end.

%% 1 / (C M): the coefficient inverted, every exponent negated; R(p)^-1 becomes R(p)/p in normalize.
inverse({Alg, Exp}, C, Ctx) ->
    {CF, M, Extras} = normalize([{K, rneg(E)} || {K, E} <- Alg], [{B, cneg(X)} || {B, X} <- Exp], Ctx),
    lists:foldl(fun(X, T) -> mul(T, X, Ctx) end, #{M => cmul(cinv(C), CF)}, Extras).

%% P^(N/D), D > 1 (SY 4.5). A monomial splits when every factor allows it (a positive or nonnegative symbol, pi,
%% an exponential of a real multiple of such a base, and a real rational coefficient with D = 2 or an exact
%% D-th root); otherwise the whole of it becomes Base(B)^(N/D).
pow_rat(P, {N, 1}, Ctx) -> pow_int(P, N, Ctx);
pow_rat(?RP = H, E, Ctx) -> pow_rat(m(H), E, Ctx);
pow_rat(P, {_, _} = E, _Ctx) when map_size(P) =:= 0 ->
    case pq_sym_exact:rcmp(E, {0, 1}) of
        1 -> #{};
        _ -> throw({non_finite, division_by_zero})
    end;
pow_rat(P, {N, D} = E, Ctx) ->
    case single(P) of
        {ok, {Alg, Exp} = M, C} ->
            case splittable(M, Ctx) andalso coeff_root(C, D) of
                {ok, Root} when M =:= ?ONE -> pow_int(Root, N, Ctx);
                {ok, Root} ->
                    Alg1 = [{K, rmul(X, E)} || {K, X} <- Alg],
                    Exp1 = [{B, cmul(X, c(E))} || {B, X} <- Exp],
                    {CF, M1, Extras} = normalize(Alg1, Exp1, Ctx),
                    Mono = lists:foldl(fun(X, T) -> mul(T, X, Ctx) end, #{M1 => CF}, Extras),
                    mul(pow_int(Root, N, Ctx), Mono, Ctx);
                false -> base_power(P, E, Ctx)
            end;
        false -> base_power(P, E, Ctx)
    end.

splittable({Alg, Exp}, Ctx) ->
    lists:all(fun({K, _}) -> positive_atom(K, Ctx) end, Alg) andalso
        lists:all(fun({{BA, BE}, X}) ->
                          pq_sym_exact:is_real(X) andalso BE =:= [] andalso
                              lists:all(fun({K, _}) -> positive_atom(K, Ctx) end, BA)
                  end, Exp).

positive_atom(pi, _) -> true;
positive_atom({sym, X}, Ctx) -> lists:member(domain(X, Ctx), [positive, nonnegative]);
positive_atom(_, _) -> false.

domain(X, Ctx) -> maps:get(X, maps:get(domains, Ctx, #{}), real).

%% C^(1/D) as a polynomial, for a real rational C: D = 2 through radicals (sqrt(-n) = i sqrt(n)); D > 2 only for
%% a positive exact root.
coeff_root(C, 1) -> {ok, const(C)};
coeff_root(C, 2) ->
    case pq_sym_exact:csqrt_exact(C) of
        {K, 1} -> {ok, const(K)};
        {K, R} -> {ok, scale(atom({rad, R}, {1, 1}), K)};
        false -> false
    end;
coeff_root(C, D) ->
    case pq_sym_exact:real_rational(C) of
        false -> false;
        R ->
            case pq_sym_exact:root_exact(R, D) of
                {ok, Root} -> {ok, const(c(Root))};
                false -> false
            end
    end.

%% G^E * Base(B')^E with P = G B' (SY 3.4): G is the positive numeric content, and also the sign for an integer E
%% (moving a sign out of a fractional power would change the branch). G is taken out only when G^E is exact.
base_power(?RP = H, E, Ctx) -> base_power(m(H), E, Ctx);
base_power(P, E, Ctx) ->
    {G, Prim} = content(P, pq_sym_exact:is_int(E)),
    {N, D} = E,
    case coeff_root(c(G), D) of
        {ok, Root} when Prim =:= #{?ONE => {q, {1, 1}, {0, 1}}} -> pow_int(Root, N, Ctx);
        {ok, Root} -> mul(pow_int(Root, N, Ctx), atom({base, Prim}, E), Ctx);
        false -> atom({base, P}, E)
    end.

%% {G, P / G}: G the gcd of the numerators over the lcm of the denominators of every rational part; for Signed,
%% also the sign that makes the first term's first non-zero part positive. Float coefficients: G = 1.
content(?RP = H, Signed) -> content(m(H), Signed);
content(P, Signed) ->
    Cs = [C || {_, C} <- sorted(P)],
    case lists:all(fun pq_sym_exact:is_exact/1, Cs) of
        false -> {{1, 1}, P};
        true ->
            Parts = lists:append([[R, I] || {q, R, I} <- Cs]),
            Num = lists:foldl(fun({A, _}, G) -> pq_sym_exact:gcd(A, G) end, 0, Parts),
            Den = lists:foldl(fun({_, B}, L) -> L * B div pq_sym_exact:gcd(L, B) end, 1, Parts),
            Sign = case Signed andalso first_negative(hd(Cs)) of true -> -1; false -> 1 end,
            G = rat(Sign * Num, Den),
            {G, scale(P, c(pq_sym_exact:rinv(G)))}
    end.

first_negative({q, {A, _}, {B, _}}) -> A < 0 orelse (A =:= 0 andalso B < 0).

%%% exp and conj

%% exp of a canonical polynomial: one factor exp(c b) per term
exp(?RP = H, Ctx) -> exp(m(H), Ctx);
exp(P, Ctx) ->
    lists:foldl(fun({M, C}, Acc) ->
                        {CF, M1, Extras} = normalize([], [{M, C}], Ctx),
                        Factor = lists:foldl(fun(X, T) -> mul(T, X, Ctx) end, #{M1 => CF}, Extras),
                        mul(Acc, Factor, Ctx)
                end, one(), sorted(P)).

%% termwise conjugate: coefficients; z <-> z* for complex symbols; exponents and bases conjugated
conj(?RP = H, Ctx) -> conj(m(H), Ctx);
conj(P, Ctx) ->
    maps:fold(fun(M, C, Acc) -> add(Acc, conj_term(M, cconj(C), Ctx)) end, #{}, P).

conj_term(M, C, Ctx) ->
    {Alg, Exp} = conj_mono(M, Ctx),
    {CF, M1, Extras} = normalize(lists:sort(Alg), lists:sort(Exp), Ctx),
    lists:foldl(fun(X, T) -> mul(T, X, Ctx) end, #{M1 => cmul(C, CF)}, Extras).

conj_mono({Alg, Exp}, Ctx) ->
    {[conj_atom(A, Ctx) || A <- Alg], [{conj_base(B, Ctx), cconj(X)} || {B, X} <- Exp]}.

conj_base(B, Ctx) ->
    {BA, BE} = conj_mono(B, Ctx),
    {lists:sort(BA), lists:sort(BE)}.

conj_atom({{sym, X}, E} = A, Ctx) ->
    case domain(X, Ctx) of
        complex -> {{conj, X}, E};
        _ -> A
    end;
conj_atom({{conj, X}, E}, _Ctx) -> {{sym, X}, E};
conj_atom({{base, B}, E}, _Ctx) -> {{cbase, B}, E};
conj_atom({{cbase, B}, E}, _Ctx) -> {{base, B}, E};
conj_atom(A, _Ctx) -> A.

%%% Polynomials kept in Rust (decision record K2)
%%
%% Under the context flag rpoly (set by pq_symbolic for the calls whose result goes straight to pq_sym_norm:to_ast:
%% canonicalise, permanent), a product or power that the kernel computes stays in Rust as {rpoly, Vars, Resource}:
%% Vars the sorted atom keys, the resource the terms (native/pq_nif/src/rpoly.rs). add, sum, neg, scale, mul,
%% pow_int and terms work on it in Rust, converting a plain map operand; every other function of this module turns
%% it back into the map first (m/1), so nothing outside sees one. A refusal (overflow, packing) also falls back to
%% the maps: the same polynomial either way (D1). Its memory is outside the process heap, so each new one is
%% checked against sym_max_heap_words (terms x (variables + 4) words).

m({rpoly, Vars, R}) -> pq_nif:rpoly_to_map(R, Vars);
m(P) -> P.

rp_terms({rpoly, _, R}) -> element(1, pq_nif:rpoly_info(R)).

%% a plain map as a handle (a handle as itself), or false
rp_from(?RP = H) -> H;
rp_from(P) ->
    case kernel_vars(P, #{}) of
        false -> false;
        VarSet ->
            Vars = lists:sort(maps:keys(VarSet)),
            case kernel_in(P, Vars) of
                {E, C, L} when L < 1 bsl 127 -> {rpoly, Vars, pq_nif:rpoly_new(length(Vars), E, C, L)};
                _ -> false
            end
    end.

size_of(?RP = H) -> rp_terms(H);
size_of(P) -> map_size(P).

rp_mul(P, Q, Ctx) ->
    Use = ?IS_RP(P) orelse ?IS_RP(Q) orelse
          (maps:get(rpoly, Ctx, false) andalso map_size(P) >= 2 andalso map_size(Q) >= 2
           andalso map_size(P) * map_size(Q) >= ?KERNEL_MIN),
    case Use andalso {rp_from(P), rp_from(Q)} of
        {{rpoly, VA, A}, {rpoly, VB, B}} ->
            check_terms(size_of(P) * size_of(Q), Ctx),
            {Vars, MA, MB} = union(VA, VB),
            rp_new(pq_nif:rpoly_mul(A, MA, B, MB, length(Vars)), Vars, Ctx);
        _ -> false
    end.

rp_add(P, Q) ->
    case {rp_from(P), rp_from(Q)} of
        {{rpoly, VA, A}, {rpoly, VB, B}} ->
            {Vars, MA, MB} = union(VA, VB),
            rp_new(pq_nif:rpoly_add(A, MA, B, MB, length(Vars)), Vars, none);
        _ -> false
    end.

%% P^N as pow_int would expand it: more than one term and N beyond max_expand_power stays a power (false here)
rp_pow(false, _N, _Ctx) -> false;
rp_pow({rpoly, Vars, R} = H, N, Ctx) ->
    T = rp_terms(H),
    case T > 1 andalso N =< maps:get(max_expand_power, Ctx, 8) of
        true ->
            check_terms(binom(N + T - 1, T - 1), Ctx),
            rp_new(pq_nif:rpoly_pow(R, N), Vars, Ctx);
        false -> false
    end.

%% times an exact coefficient (Re + i Im) / Den
rp_scale({rpoly, Vars, R}, {q, {Rn, Rd}, {In, Id}}) ->
    Den = lcm(Rd, Id),
    Re = Rn * (Den div Rd), Im = In * (Den div Id),
    case lists:all(fun(X) -> abs(X) < 1 bsl 126 end, [Re, Im, Den]) of
        true -> rp_new(pq_nif:rpoly_scale(R, Re, Im, Den), Vars, none);
        false -> false
    end;
rp_scale(_, _) -> false.

rp_new({ok, R}, Vars, none) -> {ok, {rpoly, Vars, R}};
rp_new({ok, R}, Vars, Ctx) ->
    {T, NV} = pq_nif:rpoly_info(R),
    check_terms(T, Ctx),
    Words = T * (NV + 4),
    Max = maps:get(sym_max_heap_words, Ctx, infinity),
    Words > Max andalso throw({symbolic_limit, #{reason => heap}}),
    {ok, {rpoly, Vars, R}};
rp_new(_, _, _) -> false.

%% the sorted union of two variable lists and where each list's variables go in it
union(VA, VB) ->
    Vars = lists:umerge(VA, VB),
    Pos = maps:from_list([{V, I} || {I, V} <- lists:enumerate(0, Vars)]),
    {Vars, [maps:get(V, Pos) || V <- VA], [maps:get(V, Pos) || V <- VB]}.

%% Polynomials of symbols only (integer exponents, exact coefficients) as kernel data over their common variables:
%% {Vars, [{Exps, Coef, L}]} (the coefficients times L), or false (pq_sym_slp's order search scores with it, H4)
plain_terms(Ps0) ->
    Ps = [m(P) || P <- Ps0],
    case lists:foldl(fun(P, Acc) -> kernel_vars(P, Acc) end, #{}, Ps) of
        false -> false;
        VarSet ->
            Vars = lists:sort(maps:keys(VarSet)),
            Data = [kernel_in(P, Vars) || P <- Ps],
            case lists:all(fun({sym, _}) -> true; (_) -> false end, Vars) andalso not lists:member(false, Data)
                 andalso lists:all(fun({_, _, L}) -> L < 1 bsl 127 end, Data) of
                true -> {Vars, Data};
                false -> false
            end
    end.
