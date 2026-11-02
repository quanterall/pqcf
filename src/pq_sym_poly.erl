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
-export([zero/0, one/0, const/1, symbol/2, atom/2, pi/0, is_const/1, const_value/1, single/1,
         add/2, sum/1, neg/1, scale/2, mul/3, pow_int/3, pow_rat/3, exp/2, conj/2,
         terms/1, pi_mono/0, one_mono/0, sorted/1, content/2, base_power/3]).

-import(pq_sym_exact, [c/1, cadd/2, cmul/2, cneg/1, cinv/1, cconj/1, is_zero/1, rat/2, radd/2, rmul/2, rneg/1]).

-define(ONE, {[], []}).
-define(PI, {[{pi, {1, 1}}], []}).

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

is_const(P) -> maps:size(P) =:= 0 orelse (maps:size(P) =:= 1 andalso maps:is_key(?ONE, P)).

const_value(P) -> maps:get(?ONE, P, c(0)).

single(P) when map_size(P) =:= 1 -> [{M, C}] = maps:to_list(P), {ok, M, C};
single(_) -> false.

terms(P) -> maps:size(P).

%% the terms in canonical key order
sorted(P) -> lists:sort(maps:to_list(P)).

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

neg(P) -> maps:map(fun(_, C) -> cneg(C) end, P).

scale(P, C) ->
    case is_zero(C) of
        true -> #{};
        false -> maps:filter(fun(_, V) -> not is_zero(V) end, maps:map(fun(_, V) -> cmul(V, C) end, P))
    end.

%%% Products

mul(P, Q, Ctx) ->
    check_terms(maps:size(P) * maps:size(Q), Ctx),
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
pow_int(P, N, Ctx) when N < 0 ->
    case single(P) of
        {ok, M, C} -> pow_int(inverse(M, C, Ctx), -N, Ctx);
        false when map_size(P) =:= 0 -> throw({non_finite, division_by_zero});
        false -> base_power(P, {N, 1}, Ctx)
    end;
pow_int(P, N, Ctx) ->
    Max = maps:get(max_expand_power, Ctx, 8),
    case maps:size(P) > 1 andalso N > Max of
        true -> base_power(P, {N, 1}, Ctx);
        false ->
            case maps:size(P) > 1 andalso lists:all(fun({A, E}) -> E =:= [] andalso plain(A) end, maps:keys(P)) of
                true -> multinomial(P, N, Ctx);
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
exp(P, Ctx) ->
    lists:foldl(fun({M, C}, Acc) ->
                        {CF, M1, Extras} = normalize([], [{M, C}], Ctx),
                        Factor = lists:foldl(fun(X, T) -> mul(T, X, Ctx) end, #{M1 => CF}, Extras),
                        mul(Acc, Factor, Ctx)
                end, one(), sorted(P)).

%% termwise conjugate: coefficients; z <-> z* for complex symbols; exponents and bases conjugated
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
