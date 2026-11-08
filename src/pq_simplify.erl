%% Simplifying exact results by verified rewriting (decision record R1, option 3): a rewrite is proposed -- here,
%% factors from an alphabet of candidates (the variables, small linear forms, and any the caller proposes: a person,
%% a heuristic, or a trained model) -- and kept only when it is exact. Exactness is by construction (exact division
%% over the rationals), and the result is checked again against the input at random points (verify/4): a proposal
%% can make the result smaller, never wrong.
%%
%% Polynomials are pq_recon's maps (rational/3, linear/4); the result is an AST: content * prod f_i^e_i * rest, the
%% rest the cofactor no candidate divides, kept expanded.
-module(pq_simplify).
-export([rational/4, polynomial/3, alphabet/1, adds/1, verify/4, format/1, divide/2, discover/2, gcd_poly/2]).

-type poly() :: #{[integer()] => {integer(), pos_integer()}}.
-export_type([poly/0]).

%% N / D (pq_recon maps over Vars) as {Ast, #{size_in := Before, size_out := After, form := Form}}: the smallest of
%% three forms by additions (adds/1), the more structured among equals:
%%   factored: content * prod f^e * rest over prod g^e * rest, factors from the alphabet;
%%   collected: the factors of D in the collect variables taken out, then a sum over the monomials of N in them,
%%     each with its coefficient (a rational function of the other variables) factored and cancelled against the
%%     rest of D -- for a reduction, a sum over monomials in s, t with coefficients factored in d.
%%   partial: partial fractions (apart/5) in some of the variables, each term's numerator factored.
%% Opts: alphabet => [Ast] (the candidates; default alphabet(Vars)), extra => [Ast] (proposals added to it),
%% collect => [Var] (default: every variable but d, when d is one), apart => [[Var]] (the variable sets to try partial
%% fractions in; default: all of them, d alone, all but d), form => expanded | factored | collected | partial (that
%% form, when it exists, instead of the smallest).
-spec rational(map(), map(), [atom()], map()) -> {term(), map()}.
rational(N, _D, _Vars, _Opts) when map_size(N) =:= 0 -> {0, #{size_in => 0, size_out => 0, form => expanded}};
rational(N, D, Vars, Opts) ->
    Alpha = candidates(Vars, Opts),
    One = const(1, Vars),
    {Cn, Fn, Rn} = factor(from_recon(N, Vars), Alpha),
    {Cd, Fd, Rd} = factor(from_recon(D, Vars), Alpha),
    %% factors on both sides cancel (a reconstruction's N / D is already reduced; another input may not be)
    Fn2 = [{F, E - min(E, mult(F, Fd))} || {F, E} <- Fn],
    Fd2 = [{F, E - min(E, mult(F, Fn))} || {F, E} <- Fd],
    C = qdiv(Cn, Cd),
    Ast = {mul, [rat_ast(C)]
                ++ [{pow, poly_ast(F, Vars), E} || {F, E} <- Fn2, E > 0]
                ++ [poly_ast(Rn, Vars) || Rn =/= One]
                ++ [{pow, poly_ast(F, Vars), -E} || {F, E} <- Fd2, E > 0]
                ++ [{pow, poly_ast(Rd, Vars), -1} || Rd =/= One]},
    Expanded = {mul, [poly_ast(from_recon(N, Vars), Vars), {pow, poly_ast(from_recon(D, Vars), Vars), -1}]},
    Collect = maps:get(collect, Opts, case lists:member(d, Vars) of true -> Vars -- [d]; false -> [] end),
    Forms = [{expanded, Expanded}, {factored, Ast}]
            ++ [{collected, collected(from_recon(N, Vars), {Cd, Fd, Rd}, Vars, Collect, Alpha)} || Collect =/= [], Collect =/= Vars]
            ++ [{partial, P} || X <- maps:get(apart, Opts, apart_sets(Vars)),
                                P <- [apart(from_recon(N, Vars), {Cd, Fd, Rd}, Vars, X, Alpha)], P =/= none],
    %% the fewest additions; among equals the more structured (factored, collected, partial fractions, expanded)
    Rank = #{factored => 0, collected => 1, partial => 2, expanded => 3},
    %% form => F: that form when there is one (partial fractions for an expansion in d, say), whatever its size
    Wanted = [X || X = {F, _} <- Forms, F =:= maps:get(form, Opts, none)],
    [{_, _, Form, Best} | _] = case Wanted of
                                   [{F, A} | _] -> [{adds(A), 0, F, A}];
                                   [] -> lists:sort([{adds(A), maps:get(F, Rank), F, A} || {F, A} <- Forms])
                               end,
    {Best, #{size_in => adds(Expanded), size_out => adds(Best), form => Form}}.

%% the collected form
collected(N, {Cd, Fd, Rd}, Vars, Collect, Alpha) ->
    Mask = [lists:member(V, Collect) || V <- Vars],
    Kin = fun(E) -> [case M of true -> X; false -> 0 end || {M, X} <- lists:zip(Mask, E)] end,
    Rest = fun(E) -> [case M of true -> 0; false -> X end || {M, X} <- lists:zip(Mask, E)] end,
    Free = fun(P) -> lists:all(fun(E) -> Kin(E) =:= [0 || _ <- E] end, maps:keys(P)) end,
    One = const(1, Vars),
    %% D's factors free of the collect variables stay with each term; the others are taken out
    {Fin, Fout} = lists:partition(fun({F, _}) -> Free(F) end, Fd),
    {RdIn, RdOut} = case Free(Rd) of true -> {Rd, One}; false -> {One, Rd} end,
    Groups = maps:groups_from_list(fun({E, _}) -> Kin(E) end, fun({E, C}) -> {Rest(E), C} end, maps:to_list(N)),
    Terms = [begin
                 {Cc, Fc, Rc} = factor(maps:from_list(Cs), Alpha),
                 %% cancel against D's free factors and its free rest
                 FcN = [{F, E - min(E, mult(F, Fin))} || {F, E} <- Fc],
                 FcD = [{F, E - min(E, mult(F, Fc))} || {F, E} <- Fin],
                 {RcN, RdN} = case RdIn =/= One andalso divide(Rc, RdIn) of
                                  {ok, Q} -> {Q, One};
                                  _ -> {Rc, RdIn}
                              end,
                 Mono = [{pow, {sym, V}, X} || {V, X} <- lists:zip(Vars, K), X > 0],
                 {Mono, {mul, [rat_ast(qdiv(Cc, Cd))]
                              ++ [{pow, poly_ast(F, Vars), E} || {F, E} <- FcN, E > 0]
                              ++ [poly_ast(RcN, Vars) || RcN =/= One]
                              ++ Mono
                              ++ [{pow, poly_ast(F, Vars), -E} || {F, E} <- FcD, E > 0]
                              ++ [{pow, poly_ast(RdN, Vars), -1} || RdN =/= One]}}
             end || K := Cs <- Groups],
    Outside = [{pow, poly_ast(F, Vars), -E} || {F, E} <- Fout] ++ [{pow, poly_ast(RdOut, Vars), -1} || RdOut =/= One],
    {mul, Outside ++ [{add, [T || {_, T} <- lists:sort(Terms)]}]}.

%%% Factors outside the alphabet: found by exact greatest common divisors, never guessed. A factor that is
%%% repeated in a polynomial divides its derivatives too (gcd(P, dP/dx)); a factor shared by several of a result's
%%% numerators and denominators divides their gcd. The polynomials' rests (what the alphabet left) are refined into
%%% a coprime set this way; the non-trivial pieces become candidates, and factor/2 still keeps one only where it
%%% divides exactly.

%% candidates (ASTs, for the extra option) from the coefficients [{N, D}] (pq_recon maps over Vars) of one result
-spec discover([{map(), map()}], [atom()]) -> [term()].
discover(Coeffs, Vars) ->
    Alpha = candidates(Vars, #{}),
    One = const(1, Vars),
    Rests = lists:usort([R || {N, D} <- Coeffs, P <- [N, D], map_size(P) > 0,
                              {_, _, R} <- [factor(from_recon(P, Vars), Alpha)], R =/= One, total_degree(R) > 1]),
    Pieces = refine(Rests ++ [G || R <- Rests, G <- repeated(R, Vars)], Vars),
    Found = [P || P <- Pieces, P =/= One, total_degree(P) > 0, not lists:member(P, Rests) orelse shared(P, Rests)],
    [poly_ast(P, Vars) || P <- lists:usort(Found)].

shared(P, Rests) -> length([R || R <- Rests, divide(R, P) =/= false]) > 1.

total_degree(P) -> lists:max([lists:sum(E) || E <- maps:keys(P)]).

%% the repeated part of P: gcd(P, dP/dx) for each variable x it depends on
repeated(P, Vars) ->
    lists:usort([G || I <- lists:seq(1, length(Vars)), Dp <- [deriv(P, I)], map_size(Dp) > 0,
                      G <- [gcd_poly(P, Dp)], total_degree(G) > 0]).

deriv(P, I) ->
    maps:fold(fun(E, C, Acc) ->
                      case lists:nth(I, E) of
                          0 -> Acc;
                          K -> padd(Acc, #{setnth(I, E, K - 1) => qmul({K, 1}, C)})
                      end
              end, #{}, P).

setnth(1, [_ | T], V) -> [V | T];
setnth(N, [H | T], V) -> [H | setnth(N - 1, T, V)].

%% a coprime refinement: while two pieces share a non-trivial gcd G, replace them by G and the cofactors
refine(Ps, Vars) -> refine(lists:usort([element(2, primitive(P)) || P <- Ps]), Vars, 50).
refine(Ps, _, 0) -> Ps;
refine(Ps, Vars, Budget) ->
    One = const(1, Vars),
    Pairs = [{A, B, G} || A <- Ps, B <- Ps, A < B, G <- [gcd_poly(A, B)], G =/= One],
    case Pairs of
        [] -> Ps;
        [{A, B, G} | _] ->
            {ok, Qa} = divide(A, G), {ok, Qb} = divide(B, G),
            New = [element(2, primitive(Q)) || Q <- [G, Qa, Qb], total_degree(Q) > 0],
            refine(lists:usort((Ps -- [A, B]) ++ New), Vars, Budget - 1)
    end.

%% the greatest common divisor of two polynomials over Q (primitive, positive leading coefficient): recursive in the
%% variables, by primitive pseudo-remainder sequences in the first variable either depends on
-spec gcd_poly(poly(), poly()) -> poly().
gcd_poly(A, B) when map_size(A) =:= 0 -> element(2, primitive(B));
gcd_poly(A, B) when map_size(B) =:= 0 -> element(2, primitive(A));
gcd_poly(A, B) ->
    N = length(hd(maps:keys(A))),
    case [I || I <- lists:seq(1, N), degree(A, I) > 0 orelse degree(B, I) > 0] of
        [] -> const(1, nvars(N));
        [V | _] ->
            Ca = content(A, V), Cb = content(B, V),
            {ok, Pa} = divide(A, Ca), {ok, Pb} = divide(B, Cb),
            C = gcd_poly(Ca, Cb),
            G = prs(Pa, Pb, V),
            element(2, primitive(pmul(C, G)))
    end.

degree(P, I) -> lists:max([lists:nth(I, E) || E <- maps:keys(P)]).

%% the coefficients of P in variable V (polynomials free of it)
coeffs(P, V) ->
    maps:values(maps:groups_from_list(fun({E, _}) -> lists:nth(V, E) end,
                                      fun({E, C}) -> {setnth(V, E, 0), C} end, maps:to_list(P))).

content(P, V) ->
    Cs = [maps:from_list(L) || L <- coeffs(P, V)],
    lists:foldl(fun gcd_poly/2, hd(Cs), tl(Cs)).

%% the gcd of two polynomials primitive in V
prs(A, _, V) when map_size(A) =:= 0 -> error({prs, zero, V});
prs(A, B, V) when map_size(B) =:= 0 -> pp(A, V);
prs(A, B, V) ->
    case degree(A, V) < degree(B, V) of
        true -> prs(B, A, V);
        false ->
            case map_size(B) =:= 0 of
                true -> pp(A, V);
                false ->
                    case degree(B, V) of
                        0 -> const(1, nvars(length(hd(maps:keys(A)))));
                        _ -> prs(B, pp(prem(A, B, V), V), V)
                    end
            end
    end.

pp(P, _) when map_size(P) =:= 0 -> P;
pp(P, V) -> {ok, Q} = divide(P, content(P, V)), element(2, primitive(Q)).

%% the pseudo-remainder of A by B in V
prem(A, B, V) ->
    Db = degree(B, V),
    Lb = lc(B, V),
    prem_loop(A, B, V, Db, Lb).

prem_loop(R, _, _, _, _) when map_size(R) =:= 0 -> R;
prem_loop(R, B, V, Db, Lb) ->
    case degree(R, V) of
        Dr when Dr < Db -> R;
        Dr ->
            Lr = lc(R, V),
            Shift = maps:fold(fun(E, C, Acc) -> Acc#{setnth(V, E, lists:nth(V, E) + Dr - Db) => C} end, #{}, B),
            prem_loop(psub(pmul(Lb, R), pmul(Lr, Shift)), B, V, Db, Lb)
    end.

lc(P, V) ->
    D = degree(P, V),
    maps:from_list([{setnth(V, E, 0), C} || E := C <- P, lists:nth(V, E) =:= D]).

%%% Multivariate partial fractions (Leinartas's decomposition, for denominators whose factors are affine in the
%%% variables X: as in MultivariateApart). Two moves, each an exact identity:
%%%   - the factors' linear parts dependent, sum l_i L_i = c with c free of X and not 0: 1 = (sum l_i L_i / c)^M;
%%%   - with c = 0: L_k = sum m_i L_i, so 1 = (sum m_i L_i / L_k)^M.
%%% With M = sum (e_i - 1) + 1 over the factors involved, every term of the expansion has one of them cancelled
%%% entirely (pigeonhole), so the number of distinct factors falls until the linear parts are independent. Then each
%%% numerator is reduced modulo the denominator: in coordinates where the factors are variables, the powers of a
%%% factor in the numerator cancel against the denominator's. Factors free of X (and c's) are constants here, kept
%%% outside each term.

apart_sets(Vars) ->
    lists:usort([Vars] ++ [[d] || lists:member(d, Vars), Vars =/= [d]] ++ [Vars -- [d] || lists:member(d, Vars), Vars =/= [d]]).

%% N / D, D = Cd * prod Fd * Rd, in partial fractions in the variables X: an AST, or none (a factor of D that is not
%% affine in X, or none that depends on X)
apart(N, {Cd, Fd, Rd}, Vars, X, Alpha) ->
    In = [V || V <- X, lists:member(V, Vars)] =/= [],
    Mask = [lists:member(V, X) || V <- Vars],
    One = const(1, Vars),
    Den = Fd ++ [{Rd, 1} || Rd =/= One],
    case In andalso classify(Den, Mask) of
        {ok, Lin, Out} when Lin =/= [] ->
            Terms = reduce_terms([{pscale(qdiv({1, 1}, Cd), N), maps:from_list(Lin), maps:from_list(Out)}], Vars, Mask, Alpha, []),
            Merged = merge_terms(Terms, Vars),
            case Merged of
                [] -> 0;
                _ -> apart_ast(Merged, Vars, X, Alpha)
            end;
        _ -> none
    end.

%% {ok, Affine, Free} or error: the factors affine in X (each with its multiplicity), the ones free of it
classify(Den, Mask) ->
    lists:foldl(fun(_, error) -> error;
                   ({F, E}, {ok, Lin, Out}) ->
                        case xdegree(F, Mask) of
                            0 -> {ok, Lin, [{F, E} | Out]};
                            1 -> case affine(F, Mask) of
                                     true -> {ok, [{F, E} | Lin], Out};
                                     false -> error
                                 end;
                            _ -> error
                        end
                end, {ok, [], []}, Den).

xdegree(P, Mask) -> lists:max([0 | [lists:sum([X || {true, X} <- lists:zip(Mask, E)]) || E <- maps:keys(P)]]).

%% degree one in X, its linear part with constant coefficients
affine(P, Mask) ->
    lists:all(fun(E) -> lists:sum([X || {true, X} <- lists:zip(Mask, E)]) =:= 0
                            orelse lists:sum([X || {false, X} <- lists:zip(Mask, E)]) =:= 0 end, maps:keys(P)).

%% an affine factor's linear part (a coefficient per variable of X) and its constant (a polynomial in the others)
linpart(F, Mask) ->
    Xs = [I || {I, true} <- lists:enumerate(Mask)],
    Vec = [maps:get([case J of I -> 1; _ -> 0 end || J <- lists:seq(1, length(Mask))], F, {0, 1}) || I <- Xs],
    C = maps:filter(fun(E, _) -> lists:sum([X || {true, X} <- lists:zip(Mask, E)]) =:= 0 end, F),
    {Vec, C}.

%% the terms {Num, Lin, Out} until every term's affine factors have independent linear parts, then their
%% numerators reduced
reduce_terms([], _, _, _, Done) -> Done;
reduce_terms([T = {Num, Lin, Out} | More], Vars, Mask, Alpha, Done) ->
    Fs = [F || F := E <- Lin, E > 0],
    case map_size(Num) =:= 0 of
        true -> reduce_terms(More, Vars, Mask, Alpha, Done);
        false ->
            case nullvec([element(1, linpart(F, Mask)) || F <- Fs]) of
                none -> reduce_terms(More, Vars, Mask, Alpha, numerator_reduce(T, Vars, Mask) ++ Done);
                Lambda ->
                    I = [{F, L} || {F, L} <- lists:zip(Fs, Lambda), L =/= {0, 1}],
                    C = lists:foldl(fun({F, L}, Acc) -> padd(Acc, pscale(L, element(2, linpart(F, Mask)))) end, #{}, I),
                    New = case map_size(C) of
                              0 ->
                                  %% L_k = sum m_i L_i: 1 = (sum m_i L_i / L_k)^M
                                  [{Fk, Lk} | Rest] = lists:reverse(I),
                                  Mu = [{F, qdiv(qsub({0, 1}, L), Lk)} || {F, L} <- Rest],
                                  M = lists:sum([maps:get(F, Lin) - 1 || {F, _} <- Mu]) + 1,
                                  [{pscale(Co, Num), Lin#{Fk => maps:get(Fk, Lin) + M}, Out, B} || {Co, B} <- multinomial(Mu, M)];
                              _ ->
                                  %% 1 = (sum l_i L_i / c)^M, c free of X: a constant outside the term
                                  M = lists:sum([maps:get(F, Lin) - 1 || {F, _} <- I]) + 1,
                                  {Cc, Cf, Cr} = factor(C, Alpha),
                                  Out2 = lists:foldl(fun({F, E}, O) -> O#{F => maps:get(F, O, 0) + E * M} end, Out,
                                                     Cf ++ [{Cr, 1} || Cr =/= const(1, Vars)]),
                                  Inv = qpow(qdiv({1, 1}, Cc), M),
                                  [{pscale(qmul(Inv, Co), Num), Lin, Out2, B} || {Co, B} <- multinomial(I, M)]
                          end,
                    %% each term: the factors' powers in B cancel against the denominator's, the excess to the numerator
                    Next = [begin
                                {Lin3, Num3} = maps:fold(fun(F, Bf, {L, Nm}) ->
                                                                 E = maps:get(F, L),
                                                                 case Bf >= E of
                                                                     true -> {maps:remove(F, L), pmul(Nm, ppow(F, Bf - E, Vars))};
                                                                     false -> {L#{F => E - Bf}, Nm}
                                                                 end
                                                         end, {Lin2, Num2}, B),
                                {Num3, Lin3, Out3}
                            end || {Num2, Lin2, Out3, B} <- New],
                    reduce_terms(Next ++ More, Vars, Mask, Alpha, Done)
            end
    end.

%% (sum c_i x_i)^M as [{Coefficient, #{x_i => k_i}}]
multinomial(Cs, M) ->
    [{lists:foldl(fun({{_, C}, K}, Acc) -> qmul(Acc, qpow(C, K)) end, {multicoef(M, Ks), 1}, lists:zip(Cs, Ks)),
      maps:from_list([{F, K} || {{F, _}, K} <- lists:zip(Cs, Ks), K > 0])}
     || Ks <- compositions(M, length(Cs))].

compositions(M, 1) -> [[M]];
compositions(M, N) -> [[K | R] || K <- lists:seq(0, M), R <- compositions(M - K, N - 1)].

multicoef(M, Ks) -> lists:foldl(fun(K, Acc) -> Acc div fact(K) end, fact(M), Ks).
fact(0) -> 1;
fact(N) -> N * fact(N - 1).

%% a nonzero l with sum l_i v_i = 0, or none (Gaussian elimination on [v_i | e_i])
nullvec([]) -> none;
nullvec(Vs) ->
    N = length(Vs),
    Rows = [V ++ [case J of I -> {1, 1}; _ -> {0, 1} end || J <- lists:seq(1, N)] || {I, V} <- lists:enumerate(Vs)],
    nullvec(Rows, length(hd(Vs)), 1).

nullvec(Rows, W, Col) when Col > W ->
    case [R || R <- Rows, lists:all(fun(X) -> X =:= {0, 1} end, lists:sublist(R, W))] of
        [R | _] -> lists:nthtail(W, R);
        [] -> none
    end;
nullvec(Rows, W, Col) ->
    case lists:partition(fun(R) -> lists:nth(Col, R) =/= {0, 1} end, Rows) of
        {[], _} -> nullvec(Rows, W, Col + 1);
        {[P | Ps], Zs} ->
            Pc = lists:nth(Col, P),
            Elim = [[qsub(X, qmul(qdiv(lists:nth(Col, R), Pc), Y)) || {X, Y} <- lists:zip(R, P)] || R <- Ps],
            %% the pivot row leaves: the rest no longer involves this column
            nullvec(Elim ++ Zs, W, Col + 1)
    end.

%% a term whose affine factors are independent: its numerator rewritten in coordinates y_i = L_i (the pivot
%% variables solved for), the powers of y_i cancelled against the denominator; the terms by their new denominators
numerator_reduce({Num, Lin, Out}, Vars, Mask) ->
    Fs = [F || F := E <- Lin, E > 0],
    case Fs of
        [] -> [{Num, Lin, Out}];
        _ ->
            N = length(Vars), R = length(Fs),
            Xs = [I || {I, true} <- lists:enumerate(Mask)],
            Lp = [linpart(F, Mask) || F <- Fs],
            %% the linear parts over all of Vars (0 outside X), and the pivot columns
            Full = [[case lists:member(J, Xs) of
                         true -> lists:nth(index(J, Xs), Vec);
                         false -> {0, 1}
                     end || J <- lists:seq(1, N)] || {Vec, _} <- Lp],
            Piv = pivots(Full, N),
            Inv = inverse([[lists:nth(P, Row) || P <- Piv] || Row <- Full]),
            Ext = fun(P) -> maps:fold(fun(E, C, Acc) -> Acc#{E ++ [0 || _ <- Fs] => C} end, #{}, P) end,
            Y = fun(I) -> #{[0 || _ <- Vars] ++ [case J of I -> 1; _ -> 0 end || J <- lists:seq(1, R)] => {1, 1}} end,
            XV = fun(J) -> #{[case K of J -> 1; _ -> 0 end || K <- lists:seq(1, N)] ++ [0 || _ <- Fs] => {1, 1}} end,
            %% L_i - c_i - (non-pivot part) = sum over pivots: x_P = Inv (y - c - M_rest x_rest)
            Rhs = [lists:foldl(fun(J, Acc) ->
                                       case lists:member(J, Piv) of
                                           true -> Acc;
                                           false -> psub(Acc, pscale(lists:nth(J, Row), XV(J)))
                                       end
                               end, psub(Y(I), Ext(C)), lists:seq(1, N))
                   || {I, Row, {_, C}} <- lists:zip3(lists:seq(1, R), Full, Lp)],
            Sub = maps:from_list([{P, lists:foldl(fun({Q, Rh}, Acc) -> padd(Acc, pscale(Q, Rh)) end, #{}, lists:zip(InvRow, Rhs))}
                                  || {P, InvRow} <- lists:zip(Piv, Inv)]),
            Num2 = substitute(Num, Sub, N, R),
            %% each monomial: y^b over prod y_i^e_i
            Es = [maps:get(F, Lin) || F <- Fs],
            Groups = maps:groups_from_list(
                       fun({Ex, _}) -> [max(0, E - B) || {E, B} <- lists:zip(Es, lists:nthtail(N, Ex))] end,
                       fun({Ex, C}) ->
                               {Xe, Ye} = lists:split(N, Ex),
                               lists:foldl(fun({F, E, B}, Acc) when B > E -> pmul(Acc, ppow(F, B - E, Vars));
                                              (_, Acc) -> Acc end,
                                           #{Xe => C}, lists:zip3(Fs, Es, Ye))
                       end, maps:to_list(Num2)),
            [{lists:foldl(fun padd/2, #{}, Ps),
              maps:from_list([{F, K} || {F, K} <- lists:zip(Fs, Key), K > 0]), Out} || Key := Ps <- Groups]
    end.

index(X, L) -> length(lists:takewhile(fun(Y) -> Y =/= X end, L)) + 1.

%% the pivot columns of independent rows
pivots(Rows, W) -> pivots(Rows, W, 1, []).
pivots([], _, _, Acc) -> lists:reverse(Acc);
pivots(Rows, W, Col, Acc) ->
    case lists:partition(fun(R) -> lists:nth(Col, R) =/= {0, 1} end, Rows) of
        {[], _} -> pivots(Rows, W, Col + 1, Acc);
        {[P | Ps], Zs} ->
            Pc = lists:nth(Col, P),
            Elim = [[qsub(X, qmul(qdiv(lists:nth(Col, R), Pc), Y)) || {X, Y} <- lists:zip(R, P)] || R <- Ps],
            pivots(Elim ++ Zs, W, Col + 1, [Col | Acc])
    end.

%% the inverse of a square matrix over Q (Gauss-Jordan; the rows are independent)
inverse(M) ->
    N = length(M),
    Aug = [Row ++ [case J of I -> {1, 1}; _ -> {0, 1} end || J <- lists:seq(1, N)] || {I, Row} <- lists:enumerate(M)],
    Done = lists:foldl(fun(Col, Rows) ->
                               {Before, [P0 | After]} = lists:split(Col - 1, Rows),
                               {Swap, Rest} = case lists:nth(Col, P0) of
                                                  {0, 1} ->
                                                      Q = hd([R || R <- After, lists:nth(Col, R) =/= {0, 1}]),
                                                      {Q, [case R of Q -> P0; _ -> R end || R <- After]};
                                                  _ -> {P0, After}
                                              end,
                               Pc = lists:nth(Col, Swap),
                               P = [qdiv(X, Pc) || X <- Swap],
                               Elim = fun(R) -> [qsub(X, qmul(lists:nth(Col, R), Y)) || {X, Y} <- lists:zip(R, P)] end,
                               [Elim(R) || R <- Before] ++ [P] ++ [Elim(R) || R <- Rest]
                       end, Aug, lists:seq(1, N)),
    [lists:nthtail(N, R) || R <- Done].

%% P (over Vars) with the pivot variables replaced (Sub: index => polynomial over Vars ++ R coordinates)
substitute(P, Sub, N, R) ->
    maps:fold(fun(E, C, Acc) ->
                      Mono = lists:foldl(fun({J, K}, M) ->
                                                 case Sub of
                                                     #{J := S} when K > 0 -> pmul(M, ppow(S, K, nvars(N + R)));
                                                     _ -> M
                                                 end
                                         end, #{[case maps:is_key(J, Sub) of true -> 0; false -> K end
                                                 || {J, K} <- lists:enumerate(E)] ++ [0 || _ <- lists:seq(1, R)] => C},
                                         lists:enumerate(E)),
                      padd(Acc, Mono)
              end, #{}, P).

nvars(N) -> [x || _ <- lists:seq(1, N)].

%% the terms with the same affine denominator summed (over the common multiple of their constants), the
%% constants cancelled where they divide
merge_terms(Terms, Vars) ->
    Groups = maps:groups_from_list(fun({_, Lin, _}) -> Lin end, fun({Num, _, Out}) -> {Num, Out} end, Terms),
    lists:sort([{Num, Lin, Out} || Lin := Ts <- Groups, {Num, Out} <- [sum_over(Ts, Vars)], map_size(Num) > 0]).

sum_over(Ts, _Vars) ->
    Lcm = lists:foldl(fun({_, O}, Acc) -> maps:fold(fun(F, E, A) -> A#{F => max(E, maps:get(F, A, 0))} end, Acc, O) end, #{}, Ts),
    Num = lists:foldl(fun({Nm, O}, Acc) ->
                              padd(Acc, maps:fold(fun(F, E, P) -> mpow(P, F, E - maps:get(F, O, 0)) end, Nm, Lcm))
                      end, #{}, Ts),
    %% cancel
    maps:fold(fun(F, E, {Nm, O}) ->
                      {K, Nm2} = divide_upto(Nm, F, E, 0),
                      {Nm2, case E - K of 0 -> maps:remove(F, O); R -> O#{F => R} end}
              end, {Num, Lcm}, Lcm).

mpow(P, _, 0) -> P;
mpow(P, F, K) -> mpow(pmul(P, F), F, K - 1).

divide_upto(P, _, 0, K) -> {K, P};
divide_upto(P, F, E, K) ->
    case map_size(P) > 0 andalso divide(P, F) of
        {ok, Q} -> divide_upto(Q, F, E - 1, K + 1);
        _ -> {K, P}
    end.

%% the sum as an AST: the constants every term shares taken out; each term's numerator over its constants factored,
%% or (the smaller) in partial fractions again in the other variables
apart_ast(Terms, Vars, X, Alpha) ->
    Common = lists:foldl(fun({_, _, O}, Acc) -> maps:filtermap(fun(F, E) -> case maps:get(F, O, 0) of 0 -> false; E2 -> {true, min(E, E2)} end end, Acc) end,
                         element(3, hd(Terms)), tl(Terms)),
    One = const(1, Vars),
    Others = Vars -- X,
    Asts = [begin
                {C, Fs, R} = factor(Num, Alpha),
                Own = [{F, E - maps:get(F, Common, 0)} || F := E <- Out, E > maps:get(F, Common, 0)],
                Factored = {mul, [rat_ast(C)] ++ [{pow, poly_ast(F, Vars), E} || {F, E} <- Fs]
                                 ++ [poly_ast(R, Vars) || R =/= One]
                                 ++ [{pow, poly_ast(F, Vars), -E} || {F, E} <- lists:sort(maps:to_list(Lin) ++ Own)]},
                Nested = case Own =/= [] andalso Others =/= [] andalso apart(Num, {{1, 1}, Own, One}, Vars, Others, Alpha) of
                             false -> none;
                             none -> none;
                             P -> {mul, [{pow, poly_ast(F, Vars), -E} || {F, E} <- lists:sort(maps:to_list(Lin))] ++ [P]}
                         end,
                case Nested =/= none andalso adds(Nested) < adds(Factored) of
                    true -> Nested;
                    false -> Factored
                end
            end || {Num, Lin, Out} <- Terms],
    Sum = case Asts of [A] -> A; _ -> {add, Asts} end,
    {mul, [{pow, poly_ast(F, Vars), -E} || {F, E} <- lists:sort(maps:to_list(Common))] ++ [Sum]}.

%% an AST as text, its structure kept (pq_symbolic:pretty/1 canonicalises, which expands products of sums)
-spec format(term()) -> iolist().
format(N) when is_integer(N) -> integer_to_list(N);
format({rat, N, D}) -> [integer_to_list(N), "/", integer_to_list(D)];
format({sym, X}) -> atom_to_list(X);
format({add, [E0 | Es]}) ->
    [format(E0) | [case negated(E) of
                       {true, P} -> [" - ", format(P)];
                       false -> [" + ", format(E)]
                   end || E <- Es]];
format({mul, Es}) ->
    case [E || E <- Es, E =/= 1] of
        [] -> "1";
        Fs -> lists:join("*", [paren(E) || E <- Fs])
    end;
format({pow, B, 1}) -> paren(B);
format({pow, B, E}) -> [paren(B), "^", case E < 0 of true -> ["(", integer_to_list(E), ")"]; false -> integer_to_list(E) end].

%% a term written with a leading minus: {true, its negation}
negated(N) when is_integer(N), N < 0 -> {true, -N};
negated({rat, N, D}) when N < 0 -> {true, {rat, -N, D}};
negated({mul, [C | Rest]}) ->
    case negated(C) of
        {true, 1} -> {true, {mul, Rest}};
        {true, P} -> {true, {mul, [P | Rest]}};
        false -> false
    end;
negated(_) -> false.

paren({add, [_, _ | _]} = A) -> ["(", format(A), ")"];
paren({rat, _, _} = R) -> ["(", format(R), ")"];
paren(N) when is_integer(N), N < 0 -> ["(", integer_to_list(N), ")"];
paren({mul, _} = A) ->
    case negated(A) of {true, _} -> ["(", format(A), ")"]; false -> format(A) end;
paren(A) -> format(A).

%% a polynomial (pq_recon map) as {Ast, Stats}
-spec polynomial(map(), [atom()], map()) -> {term(), map()}.
polynomial(P, Vars, Opts) ->
    rational(P, #{{[], []} => {q, {1, 1}, {0, 1}}}, Vars, Opts).

%% The default candidates: every variable x; a x + b for a = 1..3, |b| =< 12, coprime; x + y, x - y, x + 2y and
%% 2x + y for every pair (as in the denominators of reductions: d - 4, 2d - 7, s + t, ...)
-spec alphabet([atom()]) -> [term()].
alphabet(Vars) ->
    X = fun(V) -> {sym, V} end,
    [X(V) || V <- Vars]
    ++ [{add, [{mul, [A, X(V)]}, B]} || V <- Vars, A <- [1, 2, 3], B <- lists:seq(-12, 12), B =/= 0, gcd(A, abs(B)) =:= 1]
    ++ lists:append([[{add, [X(V), X(W)]}, {add, [X(V), {mul, [-1, X(W)]}]}, {add, [X(V), {mul, [2, X(W)]}]},
                      {add, [{mul, [2, X(V)]}, X(W)]}]
                     || {I, V} <- lists:enumerate(Vars), {J, W} <- lists:enumerate(Vars), I < J]).

%% the size the rewriting reduces: the additions written in an AST (a sum of k terms has k - 1), so that a product
%% of linear factors counts one per factor, and its expansion one per term
-spec adds(term()) -> non_neg_integer().
adds({add, Es}) -> length(Es) - 1 + lists:sum([adds(E) || E <- Es]);
adds({mul, Es}) -> lists:sum([adds(E) || E <- Es]);
adds({pow, B, _}) -> adds(B);
adds(_) -> 0.

%% Ast = N / D at Points random integer points (exact rational arithmetic; points where D vanishes are skipped):
%% a rewrite that changed the function is caught with probability 1 - (degree / range) per point
-spec verify(term(), map(), map(), [atom()]) -> boolean().
verify(Ast, N, D, Vars) ->
    S0 = rand:seed_s(exsss, {7, 11, 13}),
    {Points, _} = lists:mapfoldl(fun(_, S) ->
                                         lists:mapfoldl(fun(V, Sx) -> {R, Sy} = rand:uniform_s(1 bsl 40, Sx), {{V, R}, Sy} end, S, Vars)
                                 end, S0, [1, 2, 3, 4]),
    Orig = {mul, [pq_sym_norm:to_ast(N), {pow, pq_sym_norm:to_ast(D), -1}]},
    lists:all(fun(Pt) ->
                      Bind = maps:from_list(Pt),
                      case pq_symbolic:evaluate(pq_sym_norm:to_ast(D), Bind, #{exact => true}) of
                          {ok, 0} -> true;
                          _ -> pq_symbolic:evaluate(Ast, Bind, #{exact => true}) =:= pq_symbolic:evaluate(Orig, Bind, #{exact => true})
                      end
              end, Points).

%%% Factoring by trial division

%% {Content, [{Factor, Multiplicity}], Rest}: P (not 0) = Content * prod Factor^Multiplicity * Rest, the factors
%% primitive with a positive leading coefficient, Rest primitive
factor(P, Alpha) ->
    {C, P1} = primitive(P),
    {Fs, Rest} = lists:foldl(fun(F, {Acc, R}) ->
                                     case divide_all(R, F, 0) of
                                         {0, _} -> {Acc, R};
                                         {E, R2} -> {[{F, E} | Acc], R2}
                                     end
                             end, {[], P1}, Alpha),
    {C2, Rest2} = primitive(Rest),
    {qmul(C, C2), lists:reverse(Fs), Rest2}.

divide_all(P, F, E) ->
    case map_size(F) > 0 andalso divide(P, F) of
        {ok, Q} when map_size(Q) > 0 -> divide_all(Q, F, E + 1);
        _ -> {E, P}
    end.

%% exact division over the rationals: {ok, P / F} or false (polynomials as maps Exponents => {N, D}, the form
%% pq_recon works in). With a lexicographic order, F divides P only if F's
%% leading term divides P's; subtracting keeps the remainder divisible, so the division is exact iff it never stops.
-spec divide(poly(), poly()) -> {ok, poly()} | false.
divide(P, F) ->
    {Ef, Cf} = lead(F),
    divide(P, F, Ef, Cf, #{}).

divide(R, _, _, _, Q) when map_size(R) =:= 0 -> {ok, Q};
divide(R, F, Ef, Cf, Q) ->
    {Er, Cr} = lead(R),
    case lists:all(fun({A, B}) -> A >= B end, lists:zip(Er, Ef)) of
        false -> false;
        true ->
            Et = [A - B || {A, B} <- lists:zip(Er, Ef)],
            Ct = qdiv(Cr, Cf),
            R2 = maps:fold(fun(E, C, Acc) ->
                                   K = [A + B || {A, B} <- lists:zip(E, Et)],
                                   V = qsub(maps:get(K, Acc, {0, 1}), qmul(Ct, C)),
                                   case V of {0, _} -> maps:remove(K, Acc); _ -> Acc#{K => V} end
                           end, R, F),
            divide(R2, F, Ef, Cf, Q#{Et => Ct})
    end.

lead(P) -> E = lists:max(maps:keys(P)), {E, maps:get(E, P)}.

%% {Content, P / Content}: integer coefficients without a common factor, the leading one positive
primitive(P) when map_size(P) =:= 0 -> {{0, 1}, P};
primitive(P) ->
    Den = lists:foldl(fun({_, D}, L) -> L * D div gcd(L, D) end, 1, maps:values(P)),
    Ints = maps:map(fun(_, {N, D}) -> N * (Den div D) end, P),
    G0 = lists:foldl(fun(N, G) -> gcd(abs(N), G) end, 0, maps:values(Ints)),
    {_, Lc} = lead(Ints),
    G = case Lc < 0 of true -> -G0; false -> G0 end,
    {qnorm({G, Den}), maps:map(fun(_, N) -> {N div G, 1} end, Ints)}.

%%% Conversions

candidates(Vars, Opts) ->
    Asts = maps:get(alphabet, Opts, alphabet(Vars)) ++ maps:get(extra, Opts, []),
    Polys = [element(2, primitive(from_ast(A, Vars))) || A <- Asts],
    lists:usort([F || F <- Polys, map_size(F) > 0, F =/= const(1, Vars)]).

from_recon(M, Vars) ->
    Index = maps:from_list(lists:zip(Vars, lists:seq(1, length(Vars)))),
    maps:from_list([{exps(Mono, Index, length(Vars)), qnorm(R)} || {Mono, []} := {q, R, {0, 1}} <- M]).

exps(Mono, Index, N) ->
    E = maps:from_list([{maps:get(X, Index), K} || {{sym, X}, {K, 1}} <- Mono]),
    [maps:get(I, E, 0) || I <- lists:seq(1, N)].

from_ast(Ast, Vars) ->
    {ok, C} = pq_symbolic:canonicalise(Ast),
    poly_of(C, Vars).

%% a canonical AST that is a polynomial in Vars, as a poly()
poly_of(N, Vars) when is_integer(N) -> const(N, Vars);
poly_of({rat, N, D}, Vars) -> maps:map(fun(_, _) -> qnorm({N, D}) end, const(1, Vars));
poly_of({sym, X}, Vars) -> #{[case V of X -> 1; _ -> 0 end || V <- Vars] => {1, 1}};
poly_of({add, Es}, Vars) -> lists:foldl(fun(E, Acc) -> padd(Acc, poly_of(E, Vars)) end, #{}, Es);
poly_of({mul, Es}, Vars) -> lists:foldl(fun(E, Acc) -> pmul(Acc, poly_of(E, Vars)) end, const(1, Vars), Es);
poly_of({pow, B, K}, Vars) when is_integer(K), K >= 0 ->
    lists:foldl(fun(_, Acc) -> pmul(Acc, poly_of(B, Vars)) end, const(1, Vars), lists:seq(1, K)).

const(0, _) -> #{};
const(N, Vars) -> #{[0 || _ <- Vars] => {N, 1}}.

psub(A, B) -> padd(A, pscale({-1, 1}, B)).
pscale({0, _}, _) -> #{};
pscale(C, P) -> maps:map(fun(_, V) -> qmul(C, V) end, P).
ppow(P, K, Vars) -> lists:foldl(fun(_, Acc) -> pmul(Acc, P) end, const(1, Vars), lists:seq(1, K)).
qpow(_, 0) -> {1, 1};
qpow(Q, K) -> qmul(Q, qpow(Q, K - 1)).

padd(A, B) ->
    maps:filter(fun(_, {N, _}) -> N =/= 0 end,
                maps:fold(fun(K, V, Acc) -> Acc#{K => qadd(maps:get(K, Acc, {0, 1}), V)} end, A, B)).
pmul(A, B) ->
    lists:foldl(fun(T, Acc) -> padd(Acc, T) end, #{},
                [#{[X + Y || {X, Y} <- lists:zip(Ka, Kb)] => qmul(Va, Vb)} || Ka := Va <- A, Kb := Vb <- B]).

poly_ast(P, Vars) ->
    Terms = [case [{pow, {sym, V}, E} || {V, E} <- lists:zip(Vars, K), E > 0] of
                 [] -> rat_ast(C);
                 Mono when C =:= {1, 1} -> {mul, Mono};
                 Mono -> {mul, [rat_ast(C) | Mono]}
             end || {K, C} <- lists:reverse(lists:sort(maps:to_list(P)))],
    case Terms of [T] -> T; _ -> {add, Terms} end.

rat_ast({N, 1}) -> N;
rat_ast({N, D}) -> {rat, N, D}.

mult(F, Fs) -> case lists:keyfind(F, 1, Fs) of {_, E} -> E; false -> 0 end.

-spec qnorm({integer(), integer()}) -> {integer(), pos_integer()}.
qnorm({N, D}) when D < 0 -> qnorm({-N, -D});
qnorm({0, _}) -> {0, 1};
qnorm({N, D}) -> G = gcd(abs(N), D), {N div G, D div G}.
qadd({A, B}, {C, D}) -> qnorm({A * D + C * B, B * D}).
qsub({A, B}, {C, D}) -> qnorm({A * D - C * B, B * D}).
qmul({A, B}, {C, D}) -> qnorm({A * C, B * D}).
qdiv({A, B}, {C, D}) -> qnorm({A * D, B * C}).

gcd(A, 0) -> A;
gcd(A, B) -> gcd(B, A rem B).
