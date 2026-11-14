%% Factorisation of multivariate polynomials over the rationals into irreducibles (decision record SY2).
%%
%% - Square-free decomposition first (pq_sym_gcd:squarefree/1); then each square-free part.
%% - One variable: Zassenhaus's method. Factors modulo a small prime (pq_nif:factor_mod/3: distinct-degree, then
%%   Cantor-Zassenhaus equal-degree splitting), lifted by quadratic Hensel lifting on a factor tree to beyond Mignotte's bound, then
%%   recombined by trial division (von zur Gathen and Gerhard, "Modern Computer Algebra", 14-15).
%% - Several variables: Wang's method. The other variables at an integer point where the degree and the
%%   square-freeness survive, the univariate image factored, the leading coefficients fixed by giving each factor
%%   L, the leading coefficient of the whole (and the polynomial L^(r-1) f), then lifted one variable at a time
%%   modulo a prime beyond the coefficients' bound (P. S. Wang, Math. Comp. 32 (1978) 1215; Geddes, Czapor, Labahn
%%   6.4-6.8). If that fails at every point tried, a two-factor split is searched over subsets of the image's
%%   factors.
%%
%% Every factor found is checked: the product of the factors is the input. A factor is irreducible because its image
%% at the point, of the same degree in the main variable, is.
%%
%% Polynomials are pq_sym_gcd's: #{Exps => Coef}. Univariate ones inside are dense lists, the constant first.
-module(pq_sym_factor).
-export([factor/1, factor_univariate/1, result/3]).

-type poly() :: pq_sym_gcd:poly().

%% P = C * prod F_i^m_i, the F_i irreducible over Q, primitive with integer coefficients and a positive leading
%% coefficient: {C, [{F_i, m_i}]}, sorted by multiplicity, then by the factors
-spec factor(poly()) -> {{integer(), pos_integer()}, [{poly(), pos_integer()}]}.
factor(P) ->
    {C, Sqf} = pq_sym_gcd:squarefree(P),
    try
        Fs = [{F, M} || {S, M} <- Sqf, F <- factor_sqf(S)],
        {C, lists:sort(fun({F1, M1}, {F2, M2}) -> {M1, F1} =< {M2, F2} end, Fs)}
    after
        [erase(K) || {?MODULE, udio, _, _} = K <- get_keys()]
    end.

%% a factorisation {C, [{F, M}]} over the symbols Vars as pq_symbolic returns it: {Constant, [{Factor, M}]}, the
%% factors canonical, equal ones merged, sorted by multiplicity, then by the factor
-spec result({integer(), pos_integer()}, [{poly(), pos_integer()}], [atom()]) -> {term(), [{term(), pos_integer()}]}.
result({N, D}, Fs, Vars) ->
    Merged = maps:to_list(lists:foldl(fun({F, M}, Acc) -> Acc#{F => maps:get(F, Acc, 0) + M} end, #{}, Fs)),
    Asts = [{M, element(2, pq_symbolic:canonicalise(pq_sym_gcd:to_ast(F, Vars)))} || {F, M} <- Merged],
    {case D of 1 -> N; _ -> {rat, N, D} end, [{F, M} || {M, F} <- lists:sort(Asts)]}.

%% the irreducible factors of a primitive, square-free, non-constant integer polynomial
factor_sqf(S) ->
    N = length(hd(maps:keys(S))),
    Present = [V || V <- lists:seq(1, N), degree(S, V) > 0],
    case Present of
        [V] -> [from_dense(D, V, N) || D <- factor_univariate(to_dense(S, V))];
        _ -> factor_multi(S, N, Present)
    end.

%% Every variable is tried as the main one: a content in one splits S at once; an irreducible image proves S
%% irreducible; otherwise the main variable is the one whose images have the fewest factors (then the fewest terms in
%% its leading coefficient, as the lifted factors carry the whole of it)
factor_multi(S, N, Present) ->
    Try = fun(V) ->
                  Cont = pq_sym_gcd:content(S, V),
                  case is_const(Cont) of
                      false -> {split, Cont};
                      true ->
                          Order = [V | Present -- [V]],
                          F = project(S, Order),
                          Rng = rand:seed_s(exsss, {erlang:phash2(F), 17, 4}),
                          case points(F, lc1(F), length(Order), Rng, 10, 30, []) of
                              [] -> none;
                              Cands -> {ok, lists:min([R || {R, _, _} <- Cands]), map_size(lc_in(S, V)), V, Order, F, Cands}
                          end
                  end
          end,
    Results = try_vars(Try, Present, []),
    case Results of
        {split, Cont} ->
            {ok, Q} = pq_sym_gcd:divide(S, Cont),
            factor_sqf(prim(Cont)) ++ factor_sqf(prim(Q));
        irreducible -> [S];
        [] -> error({pq_sym_factor, no_good_point});
        Oks ->
            {ok, _, _, _, Order, F, Cands} = lists:min(Oks),
            {_R, A, Us} = hd(lists:keysort(1, Cands)),
            % normalised in the original variables: the leading term of the reordered polynomial is another one
            [prim(unproject(X, Order, N)) || X <- lift_all(F, lc1(F), A, Us, Cands)]
    end.

try_vars(_, [], Acc) -> Acc;
try_vars(Try, [V | Vs], Acc) ->
    case Try(V) of
        {split, _} = Sp -> Sp;
        {ok, 1, _, _, _, _, _} -> irreducible;
        none -> try_vars(Try, Vs, Acc);
        Ok -> try_vars(Try, Vs, [Ok | Acc])
    end.

%%% Several variables: Wang's method

%% a factor found by the two-factor search, factored further with every variable tried again
factor_mv(F) -> factor_sqf(prim(F)).

%% up to 3 good points [{R, A, UnivariateFactors}], drawn from [-Range, Range]
points(_, _, _, _, _, 0, Acc) -> Acc;
points(_, _, _, _, _, _, Acc) when length(Acc) >= 3 -> Acc;
points(F, L, N, Rng, Range, Tries, Acc) ->
    {A, Rng2} = lists:foldl(fun(_, {Xs, R0}) -> {X, R1} = rand:uniform_s(2 * Range + 1, R0), {Xs ++ [X - Range - 1], R1} end,
                            {[], Rng}, lists:seq(2, N)),
    Pts = lists:zip(lists:seq(2, N), A),
    U = eval_all(F, Pts),
    Ok = eval_all(L, Pts) =/= #{} andalso degree(U, 1) =:= degree(F, 1) andalso squarefree_u(U),
    case Ok of
        false -> points(F, L, N, Rng2, Range + 2, Tries - 1, Acc);
        true ->
            Us = factor_univariate(element(2, zprim(to_dense(U, 1)))),
            points(F, L, N, Rng2, Range, Tries - 1, [{length(Us), Pts, Us} | Acc])
    end.

%% gcd(U, U') = 1 modulo a large prime not dividing the leading coefficient proves it over Q; otherwise the point is
%% treated as bad (rarely wrongly: an unlucky prime only costs another point)
squarefree_u(U) ->
    P = pq_ff:prime(1),
    D = to_dense(U, 1),
    case lists:last(D) rem P of
        0 -> false;
        _ ->
            T = fun(L) -> [{[K - 1], md(C, P)} || {K, C} <- lists:enumerate(L), md(C, P) =/= 0] end,
            Du = tl([0 | [K * C || {K, C} <- lists:enumerate(tl(D))]]),
            case pq_nif:mgcd_mod(P, 1, T(D), T(Du)) of
                {ok, [{[0], _}]} -> true;
                _ -> false
            end
    end.

%% the lifting at the chosen point; then, if it fails, at the others; then two-factor splits
lift_all(F, L, A, Us, Cands) ->
    Try = [{A, Us} | [{A2, U2} || {_, A2, U2} <- Cands, A2 =/= A, length(U2) =:= length(Us)]],
    case first_ok([fun() -> safe_wang(F, L, A2, U2) end || {A2, U2} <- Try]) of
        {ok, Fs} -> Fs;
        error -> split2(F, L, A, Us)
    end.

first_ok([]) -> error;
first_ok([T | Ts]) ->
    case T() of
        {ok, _} = R -> R;
        error -> first_ok(Ts)
    end.

%% a failed lifting (a composite modulus, a factor that does not divide) is a failed point
safe_wang(F, L, Pts, Us) ->
    try wang(F, L, Pts, Us)
    catch error:_ -> error
    end.

%% Wang's lifting of all r factors: {ok, Factors} or error. The factors given L each (lifting L^(r-1) F) are unique
%% integer polynomials: lifted modulo primes near 2^63 in Rust (pq_nif:mfactor_lift/6) and brought back by Chinese
%% remaindering until they multiply back to F; modulo one prime beyond their bound in Erlang when the exponents do
%% not pack into 64 bits.
wang(F, L, Pts, Us) ->
    R = length(Us),
    N = length(hd(maps:keys(F))),
    FF = pq_sym_gcd:mul(mpow(L, R - 1, none), F),
    La = maps:get(lists:duplicate(N, 0), eval_all(L, Pts)),
    % each univariate factor scaled to the leading coefficient L(a)
    case lists:all(fun(U) -> La rem lists:last(U) =:= 0 end, Us) of
        false -> error;
        true ->
            U0 = [[C * (La div lists:last(U)) || C <- U] || U <- Us],
            MaxE = lists:max([X || E <- maps:keys(FF), X <- E]),
            case N =< 16 andalso (MaxE * R + 1) < (1 bsl (64 div N)) of
                true -> lift_primes(F, FF, L, Pts, U0, La, N, 1, none, 0);
                false -> wang_big(F, FF, L, Pts, U0, N)
            end
    end.

lift_primes(_, _, _, _, _, _, _, K, _, _) when K > 60 -> error;
lift_primes(_, _, _, _, _, _, _, _, none, Fails) when Fails >= 3 -> error;   % a bad point: no prime lifts
lift_primes(F, FF, L, Pts, U0, La, N, K, State, Fails) ->
    P = pq_ff:prime(K),
    case La rem P of
        0 -> lift_primes(F, FF, L, Pts, U0, La, N, K + 1, State, Fails);
        _ ->
            Arg = fun(M) -> [{E, md(C, P)} || E := C <- M] end,
            case pq_nif:mfactor_lift(P, N, Arg(FF), [md(A, P) || {_, A} <- Pts], [[md(C, P) || C <- U] || U <- U0], Arg(L)) of
                failed -> lift_primes(F, FF, L, Pts, U0, La, N, K + 1, State, Fails + 1);
                {ok, Ts} ->
                    Imgs = [maps:from_list(T) || T <- Ts],
                    {Hs, M} = case State of
                                  none -> {Imgs, P};
                                  {H0, M0} -> {[crt_map(H, I, M0, P) || {H, I} <- lists:zip(H0, Imgs)], M0 * P}
                              end,
                    Syms = [msym(H, M) || H <- Hs],
                    Small = lists:all(fun(X) -> maps:fold(fun(_, V, Mx) -> max(abs(V), Mx) end, 0, X) bsl 20 < M end, Syms),
                    Fs = [prim(pp1(X)) || X <- Syms],
                    case Small andalso check(F, Fs) of
                        true -> {ok, Fs};
                        false -> lift_primes(F, FF, L, Pts, U0, La, N, K + 1, {Hs, M}, Fails)
                    end
            end
    end.

crt_map(H, I, M, P) ->
    maps:filter(fun(_, V) -> V =/= 0 end,
                maps:from_list([{E, element(1, pq_ff:crt({maps:get(E, H, 0), M}, {maps:get(E, I, 0), P}))}
                                || E <- lists:usort(maps:keys(H) ++ maps:keys(I))])).

%% the lifting in Erlang, modulo one prime beyond the factors' bound
wang_big(F, FF, L, Pts, U0, N) ->
    P = big_prime(FF),
    case hensel(mmod(FF, P), Pts, [mmod(from_dense(U, 1, N), P) || U <- U0], mmod(L, P), P) of
        error -> error;
        {ok, Lifted} ->
            Fs = [prim(pp1(msym(X, P))) || X <- Lifted],
            case check(F, Fs) of
                true -> {ok, Fs};
                false -> error
            end
    end.

%% a factorisation into two at the point, over subsets of the image's factors (the smaller side first): used when
%% the image has more factors than the polynomial
split2(F, L, Pts, Us) ->
    R = length(Us),
    Subsets = [S || K <- lists:seq(1, R div 2), S <- combinations(K, lists:seq(1, R))],
    case search2(F, L, Pts, Us, Subsets) of
        none -> [prim(F)];
        {G, H} -> factor_mv(G) ++ factor_mv(H)
    end.

search2(_, _, _, _, []) -> none;
search2(F, L, Pts, Us, [S | Rest]) ->
    In = [lists:nth(I, Us) || I <- S],
    Out = [lists:nth(I, Us) || I <- lists:seq(1, length(Us)), not lists:member(I, S)],
    Two = [zmul_all(In), zmul_all(Out)],
    case safe_wang(F, L, Pts, Two) of
        {ok, [G, H]} -> {G, H};
        _ -> search2(F, L, Pts, Us, Rest)
    end.

check(F, Fs) ->
    Prod = lists:foldl(fun pq_sym_gcd:mul/2, hd(Fs), tl(Fs)),
    prim(Prod) =:= prim(F) andalso lists:all(fun(X) -> degree(X, 1) > 0 end, Fs).

%% Multivariate Hensel lifting (Geddes, Czapor, Labahn, algorithm 6.4) of F = prod U_i at x_j = a_j (j = 2..n) to
%% F, modulo P, every factor's leading coefficient in x1 forced to Lc evaluated at the variables not yet lifted
hensel(F, Pts, Us, Lc, P) ->
    N = length(hd(maps:keys(F))),
    % A_j: F with x_(j+1) .. x_n at the point, for j = 1 .. n
    As = lists:foldl(fun({J, Aj}, [Top | _] = Acc) -> [meval(Top, J, Aj, P) | Acc] end, [F], lists:reverse(Pts)),
    D = lists:max([degree(F, J) || J <- lists:seq(2, N)]),
    Lifted = lists:foldl(
               fun(_, error) -> error;
                  ({J, Aj}, U) ->
                       Aj_poly = lists:nth(J, As),
                       Later = [{K, X} || {K, X} <- Pts, K > J],
                       LcJ = lists:foldl(fun({K, X}, Acc) -> meval(Acc, K, X, P) end, Lc, Later),
                       Unew = [replace_lc(Ui, LcJ, P) || Ui <- U],
                       Earlier = [{K, X} || {K, X} <- Pts, K < J],
                       lift_var(Aj_poly, J, Aj, U, Unew, Earlier, D, P)
               end, Us, Pts),
    case Lifted of
        error -> error;
        _ ->
            case msub(F, mprod(Lifted, P), P) of
                E when map_size(E) =:= 0 -> {ok, Lifted};
                _ -> error
            end
    end.

lift_var(A, J, Aj, Uold, U, Earlier, D, P) ->
    lift_var(A, J, Aj, Uold, U, Earlier, D, P, 1, degree_or0(A, J), msub(A, mprod(U, P), P)).

lift_var(_, _, _, _, U, _, _, _, _, _, E) when map_size(E) =:= 0 -> U;
lift_var(_, _, _, _, U, _, _, _, K, Dj, _) when K > Dj -> U;
lift_var(A, J, Aj, Uold, U, Earlier, D, P, K, Dj, E) ->
    C = taylor(E, J, Aj, K, P),
    case map_size(C) of
        0 -> lift_var(A, J, Aj, Uold, U, Earlier, D, P, K + 1, Dj, E);
        _ ->
            case diophant(Uold, C, Earlier, D, P) of
                error -> error;
                Ds ->
                    Mono = shiftpow(J, Aj, K, length(hd(maps:keys(A))), P),
                    U2 = [madd(Ui, mmul(Di, Mono, P), P) || {Ui, Di} <- lists:zip(U, Ds)],
                    lift_var(A, J, Aj, Uold, U2, Earlier, D, P, K + 1, Dj, msub(A, mprod(U2, P), P))
            end
    end.

%% sigma_i with sum sigma_i prod_(k /= i) A_k = C, deg_x1 sigma_i < deg_x1 A_i, in x1 and the variables of Pts
%% (Geddes, Czapor, Labahn, algorithm 6.2), modulo P
diophant(As, C, [], _, P) -> udiophant(As, C, P);
diophant(As, C, Pts, D, P) ->
    {J, Aj} = lists:last(Pts),
    Rest = lists:droplast(Pts),
    Anew = [meval(A, J, Aj, P) || A <- As],
    Cnew = meval(C, J, Aj, P),
    case diophant(Anew, Cnew, Rest, D, P) of
        error -> error;
        Sigma ->
            Bs = cofactors(As, P),
            E = msub(C, msum([mmul(S, B, P) || {S, B} <- lists:zip(Sigma, Bs)], P), P),
            dio_loop(Anew, Rest, D, P, J, Aj, Bs, Sigma, E, 1)
    end.

dio_loop(_, _, _, _, _, _, _, Sigma, E, _) when map_size(E) =:= 0 -> Sigma;
dio_loop(_, _, D, _, _, _, _, Sigma, _, M) when M > D -> Sigma;
dio_loop(Anew, Rest, D, P, J, Aj, Bs, Sigma, E, M) ->
    Cm = taylor(E, J, Aj, M, P),
    case map_size(Cm) of
        0 -> dio_loop(Anew, Rest, D, P, J, Aj, Bs, Sigma, E, M + 1);
        _ ->
            case diophant(Anew, Cm, Rest, D, P) of
                error -> error;
                Ds0 ->
                    Mono = shiftpow(J, Aj, M, length(hd(maps:keys(E))), P),
                    Ds = [mmul(X, Mono, P) || X <- Ds0],
                    Sigma2 = [madd(S, X, P) || {S, X} <- lists:zip(Sigma, Ds)],
                    E2 = msub(E, msum([mmul(X, B, P) || {X, B} <- lists:zip(Ds, Bs)], P), P),
                    dio_loop(Anew, Rest, D, P, J, Aj, Bs, Sigma2, E2, M + 1)
            end
    end.

%% the univariate case: sigma_i = C s_i rem A_i, where s_i is the inverse of prod_(k /= i) A_k modulo A_i
udiophant(As, C, P) ->
    N = length(hd(maps:keys(C) ++ maps:keys(hd(As)))),
    Ad = [to_dense(A, 1) || A <- As],
    Cd = to_dense(C, 1),
    Key = {?MODULE, udio, Ad, P},
    Ss = case get(Key) of
             undefined ->
                 S = [inverse_mod(upmod(uprod([B || {K, B} <- lists:enumerate(Ad), K =/= I], P), P), lists:nth(I, Ad), P)
                      || I <- lists:seq(1, length(Ad))],
                 put(Key, S),
                 S;
             S -> S
         end,
    case lists:member(error, Ss) of
        true -> error;
        false -> [from_dense(element(2, udivrem(umul(Cd, Si, P), Ai, P)), 1, N) || {Si, Ai} <- lists:zip(Ss, Ad)]
    end.

%% B^-1 mod A (both mod P), or error when they are not coprime
inverse_mod(B, A, P) ->
    {G, S, _} = uegcd(B, A, P),
    case G of
        [_] -> element(2, udivrem(uscale(S, inv(hd(G), P), P), A, P));
        _ -> error
    end.

cofactors(As, P) -> [mprod([B || {K, B} <- lists:enumerate(As), K =/= I], P) || I <- lists:seq(1, length(As))].

replace_lc(U, Lc, P) ->
    D = degree(U, 1),
    Rest = maps:filter(fun([E | _], _) -> E =/= D end, U),
    madd(Rest, maps:from_list([{[D | T], C} || [_ | T] := C <- Lc]), P).

%% the K-th coefficient of E in powers of (x_J - A): sum over terms of C binom(e, K) A^(e - K)
taylor(E, J, A, K, P) ->
    maps:fold(fun(Ex, C, Acc) ->
                      case lists:nth(J, Ex) of
                          X when X < K -> Acc;
                          X ->
                              case C * binom(X, K) rem P * powm(A, X - K, P) rem P of
                                  0 -> Acc;
                                  V -> acc(setnth(J, Ex, 0), V, Acc, P)
                              end
                      end
              end, #{}, E).

%% (x_J - A)^K
shiftpow(J, A, K, N, P) ->
    maps:filter(fun(_, C) -> C =/= 0 end,
                maps:from_list([{setnth(J, lists:duplicate(N, 0), I), md(binom(K, I) * powm(-A, K - I, P), P)} || I <- lists:seq(0, K)])).

%%% One variable: Zassenhaus's method

%% the irreducible factors over Z of a primitive, square-free integer polynomial (dense, constant first) of degree
%% at least 1, each primitive with a positive leading coefficient
-spec factor_univariate([integer()]) -> [[integer()]].
factor_univariate(F0) ->
    F = positive(F0),
    case length(F) - 1 of
        1 -> [F];
        _ ->
            case F of
                [0 | T] -> [[0, 1] | factor_univariate(T)];   % a factor x
                _ -> zassenhaus(F)
            end
    end.

zassenhaus(F) ->
    Lc = lists:last(F),
    P = choose_prime(F),
    {_, Gs} = pq_nif:factor_mod(P, umonic(upmod(F, P), P), false),
    case Gs of
        [_] -> [F];
        _ ->
            B = mignotte(F),
            Target = 2 * abs(Lc) * B + 1,
            M = modulus(P, Target),
            Lifted = tree_lift(F, Gs, P, M),
            recombine(F, Lifted, M, 1)
    end.

%% P^(2^k) >= Target
modulus(P, Target) when P >= Target -> P;
modulus(P, Target) -> modulus(P * P, Target).

%% the prime (of the first ones that keep the degree and square-freeness) whose factorisation has the fewest factors
choose_prime(F) -> choose_prime(F, small_primes(), []).

choose_prime(_, [], Acc) -> element(2, lists:min(Acc));
choose_prime(_, _, Acc) when length(Acc) >= 3 -> element(2, lists:min(Acc));
choose_prime(F, [P | Ps], Acc) ->
    case lists:last(F) rem P of
        0 -> choose_prime(F, Ps, Acc);
        _ ->
            case pq_nif:factor_mod(P, umonic(upmod(F, P), P), true) of
                {0, _} -> choose_prime(F, Ps, Acc);   % not square-free mod P
                {N, _} -> choose_prime(F, Ps, [{N, P} | Acc])
            end
    end.

small_primes() ->
    Key = {?MODULE, small_primes},
    try persistent_term:get(Key)
    catch error:badarg ->
            L = [X || X <- lists:seq(3, 3000, 2), lists:all(fun(D) -> X rem D =/= 0 end, lists:seq(3, isqrt(X), 2))],
            persistent_term:put(Key, L),
            L
    end.

%% quadratic Hensel lifting on a factor tree: F = lc(F) prod Gs mod P, the Gs monic; the factors monic mod M
tree_lift(F, [_], _, M) -> [umonic(upmod(F, M), M)];
tree_lift(F, Gs, P, M) ->
    {G1, G2} = lists:split(length(Gs) div 2, Gs),
    G0 = uscale(uprod(G1, P), md(lists:last(F), P), P),
    H0 = uprod(G2, P),
    {One, S0, T0} = uegcd(G0, H0, P),
    I = inv(hd(One), P),
    {G, H} = hensel_steps(F, G0, H0, uscale(S0, I, P), uscale(T0, I, P), P, M),
    tree_lift(G, G1, P, M) ++ tree_lift(H, G2, P, M).

%% von zur Gathen and Gerhard, algorithm 15.10, repeated from modulus Mod until M
hensel_steps(_, G, H, _, _, Mod, M) when Mod >= M -> {G, H};
hensel_steps(F, G, H, S, T, Mod, M) ->
    M2 = Mod * Mod,
    E = usub(upmod(F, M2), umul(G, H, M2), M2),
    {Q, R} = udivrem(umul(S, E, M2), H, M2),
    G2 = uadd(uadd(G, umul(T, E, M2), M2), umul(Q, G, M2), M2),
    H2 = uadd(H, R, M2),
    B = usub(uadd(umul(S, G2, M2), umul(T, H2, M2), M2), [1], M2),
    {C, Dd} = udivrem(umul(S, B, M2), H2, M2),
    S2 = usub(S, Dd, M2),
    T2 = usub(usub(T, umul(T, B, M2), M2), umul(C, G2, M2), M2),
    hensel_steps(F, G2, H2, S2, T2, M2, M).

%% Zassenhaus recombination: subsets of K lifted factors, times lc(F), symmetric mod M, tried by exact division
recombine(F, Gs, _, K) when 2 * K > length(Gs) -> [positive(F)];
recombine(F, Gs, M, K) ->
    Lc = lists:last(F),
    case find_factor(F, Lc, Gs, M, combinations(K, lists:seq(1, length(Gs)))) of
        none -> recombine(F, Gs, M, K + 1);
        {Fac, Q, S} ->
            Rest = [G || {I, G} <- lists:enumerate(Gs), not lists:member(I, S)],
            [Fac | case Rest of [] -> []; _ -> recombine(Q, Rest, M, K) end]
    end.

find_factor(_, _, _, _, []) -> none;
find_factor(F, Lc, Gs, M, [S | Rest]) ->
    G = usym(uscale(uprod([lists:nth(I, Gs) || I <- S], M), md(Lc, M), M), M),
    {_, Pp} = zprim(G),
    case zdiv(F, Pp) of
        {ok, Q} -> {positive(Pp), Q, S};
        false -> find_factor(F, Lc, Gs, M, Rest)
    end.

mignotte(F) ->
    N = length(F) - 1,
    (isqrt(N + 1) + 1) * (1 bsl N) * lists:max([abs(C) || C <- F]).

%%% Dense univariate polynomials mod P (constant first, no trailing zeros)

md(X, P) -> ((X rem P) + P) rem P.
utrim(L) -> lists:reverse(lists:dropwhile(fun(X) -> X =:= 0 end, lists:reverse(L))).
upmod(A, P) -> utrim([md(X, P) || X <- A]).
usym(A, P) -> [case X > P div 2 of true -> X - P; false -> X end || X <- A].
uadd(A, B, P) -> utrim(zip_with(fun(X, Y) -> md(X + Y, P) end, A, B)).
usub(A, B, P) -> utrim(zip_with(fun(X, Y) -> md(X - Y, P) end, A, B)).
uscale(A, K, P) -> utrim([X * K rem P || X <- A]).

zip_with(F, A, B) ->
    N = max(length(A), length(B)),
    [F(X, Y) || {X, Y} <- lists:zip(A ++ lists:duplicate(N - length(A), 0), B ++ lists:duplicate(N - length(B), 0))].

umul([], _, _) -> [];
umul(_, [], _) -> [];
umul(A, B, P) ->
    Bt = list_to_tuple(B),
    Nb = tuple_size(Bt),
    N = length(A) + Nb - 1,
    At = list_to_tuple(A),
    utrim([lists:foldl(fun(I, S) -> S + element(I + 1, At) * element(K - I + 1, Bt) end, 0,
                       lists:seq(max(0, K - Nb + 1), min(K, tuple_size(At) - 1))) rem P || K <- lists:seq(0, N - 1)]).

uprod(L, P) -> lists:foldl(fun(X, Acc) -> umul(Acc, X, P) end, [1], L).

%% (Q, R) with A = Q B + R, deg R < deg B; B's leading coefficient invertible mod P
udivrem(A, B, P) ->
    Db = length(B) - 1,
    Li = inv(lists:last(B), P),
    udivrem_(lists:reverse(A), lists:reverse(B), Db, Li, P, []).

udivrem_(Ra, _, Db, _, _, Q) when length(Ra) - 1 < Db -> {utrim(Q), utrim(lists:reverse(Ra))};  % Q: the constant first
udivrem_([Lead | _] = Ra, Rb, Db, Li, P, Q) ->
    C = Lead * Li rem P,
    Sub = [md(X - C * Y, P) || {X, Y} <- lists:zip(lists:sublist(Ra, length(Rb)), Rb)],
    Ra2 = tl(Sub ++ lists:nthtail(length(Rb), Ra)),
    udivrem_(Ra2, Rb, Db, Li, P, [C | Q]).

umonic([], _) -> [];
umonic(A, P) -> uscale(A, inv(lists:last(A), P), P).

%% {G, S, T}: S A + T B = G, G the gcd (not normalised)
uegcd(A, B, P) -> uegcd(A, B, [1], [], [], [1], P).
uegcd(R0, [], S0, _, T0, _, _) -> {R0, S0, T0};
uegcd(R0, R1, S0, S1, T0, T1, P) ->
    {Q, R2} = udivrem(R0, R1, P),
    uegcd(R1, R2, S1, usub(S0, umul(Q, S1, P), P), T1, usub(T0, umul(Q, T1, P), P), P).

%%% Dense univariate polynomials over Z

positive(F) -> case lists:last(F) < 0 of true -> [-X || X <- F]; false -> F end.

zprim(F) ->
    G = lists:foldl(fun(X, A) -> pq_sym_exact:gcd(abs(X), A) end, 0, F),
    {G, positive([X div G || X <- F])}.

zmul_all(L) -> lists:foldl(fun(X, Acc) -> zmul(Acc, X) end, [1], L).
zmul(A, B) ->
    N = length(A) + length(B) - 1,
    At = list_to_tuple(A), Bt = list_to_tuple(B),
    [lists:foldl(fun(I, S) -> S + element(I + 1, At) * element(K - I + 1, Bt) end, 0,
                 lists:seq(max(0, K - tuple_size(Bt) + 1), min(K, tuple_size(At) - 1))) || K <- lists:seq(0, N - 1)].

%% {ok, A / B} when B divides A exactly over Z, else false
zdiv(A, B) when length(A) < length(B) -> false;
zdiv(A, B) -> zdiv_(lists:reverse(A), lists:reverse(B), []).
zdiv_(Ra, Rb, Q) when length(Ra) < length(Rb) ->
    case lists:all(fun(X) -> X =:= 0 end, Ra) of
        true -> {ok, Q};
        false -> false
    end;
zdiv_([Lead | _] = Ra, [Lb | _] = Rb, Q) ->
    case Lead rem Lb of
        0 ->
            C = Lead div Lb,
            Sub = [X - C * Y || {X, Y} <- lists:zip(lists:sublist(Ra, length(Rb)), Rb)],
            zdiv_(tl(Sub ++ lists:nthtail(length(Rb), Ra)), Rb, [C | Q]);
        _ -> false
    end.

%%% Sparse polynomials mod P (pq_sym_gcd's maps, coefficients in [0, P))

mnorm(A, P) -> maps:filter(fun(_, C) -> C =/= 0 end, maps:map(fun(_, C) -> md(C, P) end, A)).
mmod(A, P) -> mnorm(A, P).
msym(A, P) -> maps:map(fun(_, C) when C > P div 2 -> C - P; (_, C) -> C end, A).
madd(A, B, P) when map_size(A) < map_size(B) -> madd(B, A, P);
madd(A, B, P) -> maps:fold(fun(K, V, Acc) -> acc(K, V, Acc, P) end, A, B).

%% Acc + V x^K mod P, a zero sum removed
acc(K, V, Acc, P) ->
    case Acc of
        #{K := W} ->
            case (W + V) rem P of
                0 -> maps:remove(K, Acc);
                S -> Acc#{K := S}
            end;
        _ -> Acc#{K => V}
    end.
msub(A, B, P) -> madd(A, maps:map(fun(_, C) -> P - C end, B), P).
msum(L, P) -> lists:foldl(fun(X, Acc) -> madd(Acc, X, P) end, #{}, L).
mmul(A, B, P) ->
    Bl = maps:to_list(B),
    maps:fold(fun(Ka, Va, Acc0) ->
                      lists:foldl(fun({Kb, Vb}, Acc) -> acc(addv(Ka, Kb), Va * Vb rem P, Acc, P) end, Acc0, Bl)
              end, #{}, A).

addv([X | Xs], [Y | Ys]) -> [X + Y | addv(Xs, Ys)];
addv([], []) -> [].
mprod([H | T], P) -> lists:foldl(fun(X, Acc) -> mmul(Acc, X, P) end, H, T).
mpow(A, K, none) -> lists:foldl(fun(_, Acc) -> pq_sym_gcd:mul(Acc, A) end, #{lists:duplicate(length(hd(maps:keys(A))), 0) => 1}, lists:seq(1, K)).

%% x_J = A, the exponent of x_J then 0 (P = none: over Z)
meval(F, J, A, P) ->
    maps:fold(fun(E, C, Acc) ->
                      K = setnth(J, E, 0),
                      V = case P of none -> C * pow(A, lists:nth(J, E)); _ -> C * powm(A, lists:nth(J, E), P) rem P end,
                      S = maps:get(K, Acc, 0) + V,
                      S2 = case P of none -> S; _ -> md(S, P) end,
                      case S2 of 0 -> maps:remove(K, Acc); _ -> Acc#{K => S2} end
              end, #{}, F).

eval_all(F, Pts) -> lists:foldl(fun({J, A}, Acc) -> meval(Acc, J, A, none) end, F, Pts).

%% the leading coefficient in variable V (its terms, V's exponent kept)
lc_in(F, V) ->
    D = degree(F, V),
    maps:filter(fun(E, _) -> lists:nth(V, E) =:= D end, F).

%% the leading coefficient in x1, as a polynomial in the others
lc1(F) ->
    D = degree(F, 1),
    maps:from_list([{[0 | T], C} || [E | T] := C <- F, E =:= D]).

%% the primitive part in x1 (the content in x1 divided out), over Z
pp1(F) ->
    C = pq_sym_gcd:content(F, 1),
    case is_const(C) of
        true -> F;
        false -> {ok, Q} = pq_sym_gcd:divide(F, C), Q
    end.

%% the smallest prime above 2 B, B a bound for the coefficients of F's factors (Mignotte's, in all variables)
big_prime(F) ->
    N = lists:sum([degree(F, J) || J <- lists:seq(1, length(hd(maps:keys(F))))]),
    Lc = abs(element(2, lists:max([{E, C} || E := C <- F]))),
    Max = lists:max([abs(C) || _ := C <- F]),
    B = (isqrt(N + 1) + 1) * (1 bsl N) * Max * Lc,
    next_prime(2 * B + 1).

next_prime(N) when N rem 2 =:= 0 -> next_prime(N + 1);
next_prime(N) ->
    case probable_prime(N) of
        true -> N;
        false -> next_prime(N + 2)
    end.

%% Miller-Rabin with the first 20 primes as bases: deterministic far beyond 2^64 in practice; a composite would
%% only make an inverse fail (an error, not a wrong factor, since every factorisation is checked)
probable_prime(N) when N < 3 -> N =:= 2;
probable_prime(N) ->
    Bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71],
    case [B || B <- Bases, N rem B =:= 0] of
        [] ->
            {S, D} = split(N - 1, 0),
            lists:all(fun(A) -> witness(A, S, D, N) end, Bases);
        _ -> lists:member(N, Bases)
    end.

split(D, S) when D band 1 =:= 0 -> split(D bsr 1, S + 1);
split(D, S) -> {S, D}.

witness(A, S, D, N) ->
    X = powm(A, D, N),
    X =:= 1 orelse X =:= N - 1 orelse squares(X, S - 1, N).

squares(_, 0, _) -> false;
squares(X, K, N) ->
    Y = X * X rem N,
    Y =:= N - 1 orelse (Y =/= 1 andalso squares(Y, K - 1, N)).

%%% Helpers

inv(A, P) -> pq_ff:inv(A, P).

powm(_, 0, _) -> 1;
powm(B, E, M) ->
    B2 = md(B, M),
    case E band 1 of
        1 -> B2 * powm(B2 * B2 rem M, E bsr 1, M) rem M;
        0 -> powm(B2 * B2 rem M, E bsr 1, M)
    end.

pow(_, 0) -> 1;
pow(B, E) -> B * pow(B, E - 1).

binom(N, K) -> binom(N, K, 1, 1).
binom(_, 0, Num, Den) -> Num div Den;
binom(N, K, Num, Den) -> binom(N - 1, K - 1, Num * N, Den * K).

isqrt(0) -> 0;
isqrt(N) -> isqrt(N, N).
isqrt(N, X) ->
    Y = (X + N div X) div 2,
    case Y >= X of true -> X; false -> isqrt(N, Y) end.

combinations(0, _) -> [[]];
combinations(_, []) -> [];
combinations(K, [H | T]) -> [[H | C] || C <- combinations(K - 1, T)] ++ combinations(K, T).

degree(P, V) -> lists:max([lists:nth(V, E) || E <- maps:keys(P)]).
degree_or0(P, _) when map_size(P) =:= 0 -> 0;
degree_or0(P, V) -> degree(P, V).

is_const(P) -> map_size(P) =:= 1 andalso lists:all(fun(X) -> X =:= 0 end, hd(maps:keys(P))).

prim(P) -> element(2, pq_sym_gcd:primitive(P)).

setnth(1, [_ | T], X) -> [X | T];
setnth(K, [H | T], X) -> [H | setnth(K - 1, T, X)].

to_dense(P, V) ->
    D = degree(P, V),
    [maps:get(K, maps:from_list([{lists:nth(V, E), C} || E := C <- P]), 0) || K <- lists:seq(0, D)].

from_dense(L, V, N) ->
    maps:from_list([{setnth(V, lists:duplicate(N, 0), K - 1), C} || {K, C} <- lists:enumerate(L), C =/= 0]).

%% the variables reordered: position I of the result is variable lists:nth(I, Order); variables not in Order dropped
project(P, Order) -> maps:from_list([{[lists:nth(V, E) || V <- Order], C} || E := C <- P]).

unproject(P, Order, N) ->
    maps:from_list([{lists:foldl(fun({I, V}, Acc) -> setnth(V, Acc, lists:nth(I, E)) end, lists:duplicate(N, 0),
                                 lists:enumerate(Order)), C} || E := C <- P]).
