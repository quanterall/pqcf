%% Simplifying exact results by verified rewriting (decision record R1, option 3): a rewrite is proposed -- here,
%% factors from an alphabet of candidates (the variables, small linear forms, and any the caller proposes: a person,
%% a heuristic, or a trained model) -- and kept only when it is exact. Exactness is by construction (exact division
%% over the rationals), and the result is checked again against the input at random points (verify/4): a proposal
%% can make the result smaller, never wrong.
%%
%% Polynomials are pq_recon's maps (rational/3, linear/4); the result is an AST: content * prod f_i^e_i * rest, the
%% rest the cofactor no candidate divides, kept expanded.
-module(pq_simplify).
-export([rational/4, polynomial/3, alphabet/1, adds/1, verify/4, format/1, divide/2]).

-type poly() :: #{[integer()] => {integer(), pos_integer()}}.
-export_type([poly/0]).

%% N / D (pq_recon maps over Vars) as {Ast, #{size_in := Before, size_out := After, form := Form}}: the smallest of
%% three forms by additions (adds/1), the more structured among equals:
%%   factored: content * prod f^e * rest over prod g^e * rest, factors from the alphabet;
%%   collected: the factors of D in the collect variables taken out, then a sum over the monomials of N in them,
%%     each with its coefficient (a rational function of the other variables) factored and cancelled against the
%%     rest of D -- for a reduction, a sum over monomials in s, t with coefficients factored in d.
%% Opts: alphabet, extra (as for the candidates) and collect => [Var] (default: every variable but d, when d is one). Opts: alphabet => [Ast] (the
%% candidates; default alphabet(Vars)), extra => [Ast] (proposals added to it).
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
            ++ [{collected, collected(from_recon(N, Vars), {Cd, Fd, Rd}, Vars, Collect, Alpha)} || Collect =/= [], Collect =/= Vars],
    %% the fewest additions; among equals the more structured (factored, then collected, then expanded)
    Rank = #{factored => 0, collected => 1, expanded => 2},
    [{_, _, Form, Best} | _] = lists:sort([{adds(A), maps:get(F, Rank), F, A} || {F, A} <- Forms]),
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
