%% The greatest common divisor of multivariate polynomials over the rationals (decision record SY2), by the modular
%% method: the gcd modulo primes near 2^63 (pq_nif:mgcd_mod/4, Brown's algorithm in native/pq_nif/src/mgcd.rs),
%% scaled to the gcd of the leading coefficients, brought back to the integers by Chinese remaindering, and accepted
%% only when it divides both inputs exactly.
%%
%% A polynomial is a map #{Exps => Coef} over a fixed list of variables: Exps a list of non-negative exponents, one
%% per variable; Coef an integer or a fraction {N, D} (D > 0); no zero coefficients. Terms are ordered by their
%% exponents, lexicographically, the first variable most significant. A gcd is returned primitive, with integer
%% coefficients and a positive leading coefficient: unique, since over Q a gcd is defined up to a constant factor.
-module(pq_sym_gcd).
-export([gcd/2, divide/2, primitive/1, from_ast/2, to_ast/2, squarefree/1, mul/2, derivative/2, content/2]).

-export_type([poly/0]).
-type poly() :: #{[non_neg_integer()] => integer() | {integer(), pos_integer()}}.

-define(MAX_PRIMES, 2000).

%% the gcd of A and B: primitive, integer coefficients, positive leading coefficient (the gcd of 0 and 0 is 0)
-spec gcd(poly(), poly()) -> poly().
gcd(A, B) when map_size(A) =:= 0 -> element(2, primitive(B));
gcd(A, B) when map_size(B) =:= 0 -> element(2, primitive(A));
gcd(A0, B0) ->
    {_, A} = primitive(A0),
    {_, B} = primitive(B0),
    N = length(hd(maps:keys(A))),
    {_, La} = lead(A),
    {_, Lb} = lead(B),
    case is_const(A) orelse is_const(B) of
        true -> one(N);
        false ->
            Gamma = pq_sym_exact:gcd(abs(La), abs(Lb)),
            St = #{a => A, b => B, n => N, la => La, lb => Lb, gamma => Gamma,
                   ta => maps:to_list(A), tb => maps:to_list(B)},
            primes(1, St, none)
    end.

%% State: none, or {LeadingExps, #{Exps => Residue}, Modulus, Symmetric} for the primes combined so far
primes(K, _, _) when K > ?MAX_PRIMES -> error({pq_sym_gcd, too_many_primes});
primes(K, #{n := N, la := La, lb := Lb, gamma := Gamma, ta := Ta, tb := Tb} = St, State) ->
    P = pq_ff:prime(K),
    case La rem P =:= 0 orelse Lb rem P =:= 0 of
        true -> primes(K + 1, St, State);
        false ->
            case pq_nif:mgcd_mod(P, N, [{E, mod(C, P)} || {E, C} <- Ta], [{E, mod(C, P)} || {E, C} <- Tb]) of
                failed -> primes(K + 1, St, State);
                {ok, [{E0, _} | _] = G} ->
                    case lists:all(fun(X) -> X =:= 0 end, E0) of
                        true -> one(N);
                        false ->
                            Img = maps:from_list([{E, C * Gamma rem P} || {E, C} <- G]),
                            combine(K, St, State, E0, Img, P)
                    end
            end
    end.

combine(K, St, State, Lm, Img, P) ->
    case State of
        {Cur, _, _, _} when Lm > Cur -> primes(K + 1, St, State);  % an unlucky prime: skip it
        {Cur, H, M, Sym} when Lm =:= Cur ->
            Keys = lists:usort(maps:keys(H) ++ maps:keys(Img)),
            H2 = maps:from_list([{E, element(1, pq_ff:crt({maps:get(E, H, 0), M}, {maps:get(E, Img, 0), P}))}
                                 || E <- Keys]),
            M2 = M * P,
            Sym2 = symmetric(H2, M2),
            try_accept(K, St, {Lm, H2, M2, Sym2}, Sym2 =:= Sym);
        _ -> try_accept(K, St, {Lm, Img, P, symmetric(Img, P)}, false)  % the first prime, or a better one
    end.

%% A candidate that divides both inputs is the gcd: it divides the gcd, and its leading monomial is that of the gcd
%% modulo a prime, which is at least the gcd's. So it is tried as soon as the residues look settled (unchanged by
%% the last prime, or far below the modulus), without waiting for a confirming prime.
try_accept(K, St, {_, _, M, Sym} = State, Same) ->
    Small = maps:fold(fun(_, V, Mx) -> max(abs(V), Mx) end, 0, Sym) bsl 20 < M,
    case (Same orelse Small) andalso accept(Sym, St) of
        {ok, G} -> G;
        _ -> primes(K + 1, St, State)
    end.

symmetric(H, M) ->
    Half = M div 2,
    maps:filter(fun(_, V) -> V =/= 0 end, maps:map(fun(_, V) when V > Half -> V - M; (_, V) -> V end, H)).

accept(Sym, #{a := A, b := B}) ->
    {_, G} = primitive(Sym),
    case divide(A, G) =/= false andalso divide(B, G) =/= false of
        true -> {ok, G};
        false -> false
    end.

%% {ok, A / G} when G divides A exactly (over Q), else false
-spec divide(poly(), poly()) -> {ok, poly()} | false.
divide(_, G) when map_size(G) =:= 0 -> error(badarith);
divide(A, _) when map_size(A) =:= 0 -> {ok, #{}};
divide(A0, G0) ->
    {Ca, A} = primitive(A0),
    {Cg, G} = primitive(G0),
    {Lm, Lc} = lead(G),
    Gt = maps:to_list(G),
    case idiv(gb_trees:from_orddict(lists:sort(maps:to_list(A))), Lm, Lc, Gt, []) of
        false -> false;
        Q -> {ok, scale(qdiv(Ca, Cg), maps:from_list(Q))}
    end.

%% division of the integer polynomial R by the primitive G (leading term Lc x^Lm): the quotient terms, or false
%% when a leading term of the rest is not divisible (G does not divide over Q: Gauss's lemma, G is primitive)
idiv(R, Lm, Lc, Gt, Q) ->
    case gb_trees:is_empty(R) of
        true -> Q;
        false ->
            {E, C} = gb_trees:largest(R),
            case C rem Lc =:= 0 andalso lists:all(fun({X, Y}) -> X >= Y end, lists:zip(E, Lm)) of
                false -> false;
                true ->
                    S = [X - Y || {X, Y} <- lists:zip(E, Lm)],
                    T = C div Lc,
                    R2 = lists:foldl(fun({Eg, Cg}, Acc) ->
                                             Key = [X + Y || {X, Y} <- lists:zip(Eg, S)],
                                             case gb_trees:lookup(Key, Acc) of
                                                 none -> gb_trees:insert(Key, -T * Cg, Acc);
                                                 {value, V} when V =:= T * Cg -> gb_trees:delete(Key, Acc);
                                                 {value, V} -> gb_trees:update(Key, V - T * Cg, Acc)
                                             end
                                     end, R, Gt),
                    idiv(R2, Lm, Lc, Gt, [{S, T} | Q])
            end
    end.

%% {Content, P / Content}: the content a fraction (positive when P's leading coefficient is), the rest with integer
%% coefficients, coprime, the leading one positive
-spec primitive(poly()) -> {{integer(), pos_integer()}, poly()}.
primitive(P) when map_size(P) =:= 0 -> {{0, 1}, P};
primitive(P) ->
    case maps:fold(fun(_, C, Ok) -> Ok andalso is_integer(C) end, true, P) of
        true ->
            G0 = maps:fold(fun(_, N, G) -> pq_sym_exact:gcd(abs(N), G) end, 0, P),
            {_, L} = lead(P),
            G = case L < 0 of true -> -G0; false -> G0 end,
            {{G, 1}, case G of 1 -> P; _ -> maps:map(fun(_, N) -> N div G end, P) end};
        false -> primitive_q(P)
    end.

primitive_q(P0) ->
    P = maps:map(fun(_, C) when is_integer(C) -> {C, 1}; (_, {_, _} = C) -> C end, P0),
    Den = maps:fold(fun(_, {_, D}, L) -> L * D div pq_sym_exact:gcd(L, D) end, 1, P),
    Ints = maps:map(fun(_, {N, D}) -> N * (Den div D) end, P),
    G0 = maps:fold(fun(_, N, G) -> pq_sym_exact:gcd(abs(N), G) end, 0, Ints),
    {_, L} = lead(Ints),
    G = case L < 0 of true -> -G0; false -> G0 end,
    {qnorm({G, Den}), maps:map(fun(_, N) -> N div G end, Ints)}.

scale({N, 1}, P) -> maps:map(fun(_, C) -> C * N end, P);
scale(Q, P) -> maps:map(fun(_, C) -> case qnorm(qmul(Q, {C, 1})) of {X, 1} -> X; R -> R end end, P).

lead(P) ->
    E = lists:max(maps:keys(P)),
    {E, maps:get(E, P)}.

is_const(P) -> map_size(P) =:= 1 andalso lists:all(fun(X) -> X =:= 0 end, hd(maps:keys(P))).

one(N) -> #{lists:duplicate(N, 0) => 1}.

mod(A, M) -> ((A rem M) + M) rem M.

qmul({A, B}, {C, D}) -> {A * C, B * D}.
qdiv({A, B}, {C, D}) -> qnorm({A * D, B * C}).
qnorm({N, D}) when D < 0 -> qnorm({-N, -D});
qnorm({0, _}) -> {0, 1};
qnorm({N, D}) ->
    case pq_sym_exact:gcd(abs(N), D) of
        1 -> {N, D};
        G -> {N div G, D div G}
    end.

%%% Square-free decomposition

%% P = C * prod F_i^i, the F_i square-free, pairwise coprime, primitive with a positive leading coefficient, each
%% non-constant: {C, [{F_i, i}]} (i ascending). Yun's algorithm in the first variable P depends on, its content in that
%% variable decomposed recursively (D. Y. Y. Yun, SYMSAC 1976; Geddes, Czapor, Labahn 8.2).
-spec squarefree(poly()) -> {{integer(), pos_integer()}, [{poly(), pos_integer()}]}.
squarefree(P0) when map_size(P0) =:= 0 -> {{0, 1}, []};
squarefree(P0) ->
    {C, P} = primitive(P0),
    % the monomial content first: x_i^e contributes x_i with multiplicity e
    N = length(hd(maps:keys(P))),
    Low = [lists:min([lists:nth(I, E) || E <- maps:keys(P)]) || I <- lists:seq(1, N)],
    Rest = case lists:all(fun(X) -> X =:= 0 end, Low) of
               true -> P;
               false -> maps:from_list([{[X - Y || {X, Y} <- lists:zip(E, Low)], V} || E := V <- P])
           end,
    Mono = [{#{[case J of I -> 1; _ -> 0 end || J <- lists:seq(1, N)] => 1}, E} || {I, E} <- lists:enumerate(Low), E > 0],
    {C, lists:keysort(2, merge(Mono, sqf(Rest)))}.

sqf(P) ->
    N = length(hd(maps:keys(P))),
    case [V || V <- lists:seq(1, N), degree(P, V) > 0] of
        [] -> [];
        [V | _] ->
            Cont = content(P, V),
            {ok, Pp} = divide(P, Cont),
            merge(sqf(Cont), yun(Pp, V))
    end.

%% Yun's algorithm for P primitive in V
yun(P, V) ->
    D = derivative(P, V),
    G = gcd(P, D),
    {ok, C} = divide(P, G),
    {ok, Dg} = divide(D, G),
    yun(C, sub(Dg, derivative(C, V)), V, 1, []).

yun(C, Dd, V, I, Acc) ->
    case is_const(C) of
        true -> lists:reverse(Acc);
        false ->
            A = gcd(C, Dd),
            {ok, C2} = divide(C, A),
            {ok, Da} = divide(Dd, A),
            Acc2 = case is_const(A) of true -> Acc; false -> [{A, I} | Acc] end,
            yun(C2, sub(Da, derivative(C2, V)), V, I + 1, Acc2)
    end.

%% the factors of two decompositions of coprime polynomials, multiplied where the multiplicities agree
merge(As, Bs) ->
    Ms = lists:usort([M || {_, M} <- As ++ Bs]),
    [{element(2, primitive(lists:foldl(fun mul/2, hd(Fs), tl(Fs)))), M}
     || M <- Ms, Fs <- [[F || {F, M2} <- As ++ Bs, M2 =:= M]]].

%% the gcd of P's coefficients in variable V
-spec content(poly(), pos_integer()) -> poly().
content(P, V) ->
    Groups = maps:groups_from_list(fun({E, _}) -> lists:nth(V, E) end,
                                   fun({E, C}) -> {setnth(V, E, 0), C} end, maps:to_list(P)),
    Cs = lists:sort(fun(X, Y) -> length(X) =< length(Y) end, maps:values(Groups)),
    lists:foldl(fun(L, G) ->
                        case is_const(G) of
                            true -> G;   % 1: done
                            false -> gcd(G, maps:from_list(L))
                        end
                end, element(2, primitive(maps:from_list(hd(Cs)))), tl(Cs)).

degree(P, V) -> lists:max([lists:nth(V, E) || E <- maps:keys(P)]).

setnth(1, [_ | T], X) -> [X | T];
setnth(K, [H | T], X) -> [H | setnth(K - 1, T, X)].

%% dP/dV
-spec derivative(poly(), pos_integer()) -> poly().
derivative(P, V) ->
    maps:from_list([{setnth(V, E, K - 1), cmul(K, C)} || E := C <- P, K <- [lists:nth(V, E)], K > 0]).

%% A * B
-spec mul(poly(), poly()) -> poly().
mul(A, B) ->
    maps:filter(fun(_, C) -> C =/= 0 end,
                lists:foldl(fun({K, V}, Acc) -> Acc#{K => cadd(maps:get(K, Acc, 0), V)} end, #{},
                            [{[X + Y || {X, Y} <- lists:zip(Ka, Kb)], cmul(Va, Vb)} || Ka := Va <- A, Kb := Vb <- B])).

sub(A, B) ->
    maps:filter(fun(_, C) -> C =/= 0 end,
                maps:fold(fun(K, V, Acc) -> Acc#{K => cadd(maps:get(K, Acc, 0), cmul(-1, V))} end, A, B)).

%% coefficients: integers, or fractions {N, D} (D > 1) in lowest terms
cmul(A, B) when is_integer(A), is_integer(B) -> A * B;
cmul(A, B) -> cnorm(qmul(q(A), q(B))).
cadd(A, B) when is_integer(A), is_integer(B) -> A + B;
cadd(A, B) -> cnorm(qadd(q(A), q(B))).
q(A) when is_integer(A) -> {A, 1};
q(Q) -> Q.
cnorm(Q) -> case qnorm(Q) of {N, 1} -> N; R -> R end.

%%% Conversions from and to the engine's canonical forms

%% a canonical polynomial (pq_symbolic) over the symbols Vars as a poly(); anything else is refused
-spec from_ast(term(), [atom()]) -> poly().
from_ast(Ast, Vars) ->
    Index = maps:from_list(lists:zip(Vars, lists:seq(1, length(Vars)))),
    Terms = case Ast of
                {add, Ts} -> Ts;
                _ -> [Ast]
            end,
    lists:foldl(fun(T, Acc) ->
                        {E, C} = term(T, Index, length(Vars)),
                        case qnorm(qadd(maps:get(E, Acc, {0, 1}), C)) of
                            {0, _} -> maps:remove(E, Acc);
                            {N, 1} -> Acc#{E => N};
                            Q -> Acc#{E => Q}
                        end
                end, #{}, [T || T <- Terms, T =/= 0]).

term(T, Index, N) ->
    Fs = case T of {mul, L} -> L; _ -> [T] end,
    {Exps, C} = lists:foldl(fun(F, {Ex, Q}) -> factor(F, Index, Ex, Q) end, {#{}, {1, 1}}, Fs),
    {[maps:get(I, Exps, 0) || I <- lists:seq(1, N)], C}.

factor(K, _, Ex, Q) when is_integer(K) -> {Ex, qmul(Q, {K, 1})};
factor({rat, A, B}, _, Ex, Q) -> {Ex, qmul(Q, {A, B})};
factor({sym, X}, Index, Ex, Q) -> factor({pow, {sym, X}, 1}, Index, Ex, Q);
factor({pow, {sym, X}, K}, Index, Ex, Q) when is_integer(K), K > 0 ->
    I = maps:get(X, Index),
    {Ex#{I => maps:get(I, Ex, 0) + K}, Q};
factor(_, _, _, _) -> throw({symbolic_unsupported, gcd}).

%% a poly() over the symbols Vars as an expression (not yet canonical)
-spec to_ast(poly(), [atom()]) -> term().
to_ast(P, _) when map_size(P) =:= 0 -> 0;
to_ast(P, Vars) ->
    {add, [{mul, [case C of {A, B} -> {rat, A, B}; _ -> C end
                  | [{pow, {sym, V}, K} || {V, K} <- lists:zip(Vars, E), K > 0]]} || E := C <- P]}.

qadd({A, B}, {C, D}) -> {A * D + C * B, B * D}.
