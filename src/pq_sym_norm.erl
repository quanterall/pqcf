%% The normaliser (SY 4.2): one bottom-up pass from sym_expr() to the canonical form, and back to the canonical
%% AST shape of SY 4.1. Matrices are expanded entrywise (SY 4.6, mode `expand`); {identity, N} stays structural
%% until it is combined. Throws {Category, Detail}; pq_symbolic catches.
-module(pq_sym_norm).
-export([norm/2, to_ast/1, value_to_ast/1, coeff_to_ast/1]).
-export_type([value/0]).

-import(pq_sym_exact, [c/1]).

-type value() :: {scalar, map()} | {mat, pos_integer(), pos_integer(), [[map()]]} | {identity, pos_integer()}.

-spec norm(term(), map()) -> value().
norm(N, _Ctx) when is_integer(N) -> {scalar, pq_sym_poly:const(c(N))};
norm(X, _Ctx) when is_float(X) -> {scalar, pq_sym_poly:const(c(X))};
norm({rat, N, D}, _Ctx) when is_integer(N), is_integer(D), D > 0 ->
    {scalar, pq_sym_poly:const(c(pq_sym_exact:rat(N, D)))};
norm({complex, Re, Im}, _Ctx) ->
    {scalar, pq_sym_poly:const(complex_coeff(Re, Im))};
norm({const, pi}, _Ctx) -> {scalar, pq_sym_poly:pi()};
norm({const, e}, Ctx) -> {scalar, pq_sym_poly:exp(pq_sym_poly:one(), Ctx)};
norm({const, i}, _Ctx) -> {scalar, pq_sym_poly:const(pq_sym_exact:cq({0, 1}, {1, 1}))};
norm({sym, X}, _Ctx) when is_atom(X) -> {scalar, pq_sym_poly:atom({sym, X}, {1, 1})};
norm({add, Args}, Ctx) when is_list(Args) ->
    {scalar, pq_sym_poly:sum([scalar(A, Ctx) || A <- Args])};
norm({mul, Args}, Ctx) when is_list(Args) ->
    {scalar, lists:foldl(fun(A, Acc) -> pq_sym_poly:mul(Acc, scalar(A, Ctx), Ctx) end, pq_sym_poly:one(), Args)};
norm({pow, B, E}, Ctx) -> {scalar, power(scalar(B, Ctx), E, Ctx)};
norm({exp, A}, Ctx) -> {scalar, pq_sym_poly:exp(scalar(A, Ctx), Ctx)};
norm({sin, A}, Ctx) -> {scalar, trig(sin, scalar(A, Ctx), Ctx)};
norm({cos, A}, Ctx) -> {scalar, trig(cos, scalar(A, Ctx), Ctx)};
norm({conj, A}, Ctx) ->
    case norm(A, Ctx) of
        {scalar, P} -> {scalar, pq_sym_poly:conj(P, Ctx)};
        _ -> throw({bad_shape, {conj, matrix_operand}})
    end;
norm({matrix, [First | _] = Rows}, Ctx) when is_list(First), First =/= [] ->
    C = length(First),
    lists:all(fun(R) -> is_list(R) andalso length(R) =:= C end, Rows)
        orelse throw({bad_shape, {matrix, ragged}}),
    {mat, length(Rows), C, [[scalar(E, Ctx) || E <- R] || R <- Rows]};
norm({identity, N}, _Ctx) when is_integer(N), N > 0 -> {identity, N};
norm({madd, [_ | _] = Ms}, Ctx) ->
    [First | Rest] = [matrix(M, Ctx) || M <- Ms],
    lists:foldl(fun(M, Acc) -> entrywise2(fun pq_sym_poly:add/2, Acc, M, madd) end, First, Rest);
norm({smul, S, M}, Ctx) ->
    P = scalar(S, Ctx),
    entrywise(fun(E) -> pq_sym_poly:mul(P, E, Ctx) end, matrix(M, Ctx));
norm({matmul, A, B}, Ctx) -> matmul(matrix(A, Ctx), matrix(B, Ctx), Ctx);
norm({kron, A, B}, Ctx) -> kron(matrix(A, Ctx), matrix(B, Ctx), Ctx);
norm({transpose, M}, Ctx) ->
    case norm(M, Ctx) of
        {identity, _} = I -> I;
        {mat, _, _, _} = X -> transpose(X);
        _ -> throw({bad_shape, {transpose, scalar_operand}})
    end;
norm({dagger, M}, Ctx) ->
    case norm(M, Ctx) of
        {identity, _} = I -> I;
        {mat, _, _, _} = X -> entrywise(fun(E) -> pq_sym_poly:conj(E, Ctx) end, transpose(X));
        _ -> throw({bad_shape, {dagger, scalar_operand}})
    end;
norm({subs, E, Map}, Ctx) when is_map(Map) -> norm(pq_symbolic:substitute(E, Map), Ctx);
norm(Other, _Ctx) -> throw({bad_expr, Other}).

complex_coeff(Re, Im) ->
    case {exact(Re), exact(Im)} of
        {{ok, R}, {ok, I}} -> pq_sym_exact:cq(R, I);
        _ -> pq_sym_exact:cf(num(Re), num(Im))
    end.

exact(N) when is_integer(N) -> {ok, {N, 1}};
exact({rat, N, D}) when is_integer(N), is_integer(D), D > 0 -> {ok, pq_sym_exact:rat(N, D)};
exact(X) when is_float(X) -> float;
exact(X) -> throw({bad_expr, X}).

num(N) when is_integer(N) -> float(N);
num(X) when is_float(X) -> X;
num({rat, N, D}) -> N / D.

scalar(E, Ctx) ->
    case norm(E, Ctx) of
        {scalar, P} -> P;
        _ -> throw({bad_shape, {scalar_expected, E}})
    end.

matrix(E, Ctx) ->
    case norm(E, Ctx) of
        {scalar, _} -> throw({bad_shape, {matrix_expected, E}});
        M -> M
    end.

%% exponents: integer, {rat, N, D} or a float (a float exponent only of a constant)
power(P, N, Ctx) when is_integer(N) -> pq_sym_poly:pow_int(P, N, Ctx);
power(P, {rat, N, D}, Ctx) when is_integer(N), is_integer(D), D > 0 ->
    pq_sym_poly:pow_rat(P, pq_sym_exact:rat(N, D), Ctx);
power(P, X, _Ctx) when is_float(X) ->
    case pq_sym_poly:is_const(P) of
        true ->
            {Re, Im} = pq_sym_exact:to_complex(pq_sym_poly:const_value(P)),
            {R, I} = pq_sym_eval:cpow({Re, Im}, {X, 0.0}),
            pq_sym_poly:const(pq_sym_exact:cf(R, I));
        false -> throw({symbolic_unsupported, float_exponent})
    end;
power(_P, E, _Ctx) -> throw({bad_expr, {pow, E}}).

%% sin x = (exp(i x) - exp(-i x)) / (2 i), cos x = (exp(i x) + exp(-i x)) / 2
trig(F, P, Ctx) ->
    I = pq_sym_exact:cq({0, 1}, {1, 1}),
    Plus = pq_sym_poly:exp(pq_sym_poly:scale(P, I), Ctx),
    Minus = pq_sym_poly:exp(pq_sym_poly:scale(P, pq_sym_exact:cneg(I)), Ctx),
    case F of
        cos -> pq_sym_poly:scale(pq_sym_poly:add(Plus, Minus), c({1, 2}));
        sin -> pq_sym_poly:scale(pq_sym_poly:add(Plus, pq_sym_poly:neg(Minus)), pq_sym_exact:cq({0, 1}, {-1, 2}))
    end.

%%% Matrices

expand({identity, N}) ->
    {mat, N, N, [[case I =:= J of true -> pq_sym_poly:one(); false -> pq_sym_poly:zero() end
                  || J <- lists:seq(1, N)] || I <- lists:seq(1, N)]};
expand(M) -> M.

entrywise(F, M) ->
    {mat, R, C, Rows} = expand(M),
    {mat, R, C, [[F(E) || E <- Row] || Row <- Rows]}.

entrywise2(F, A, B, Op) ->
    {mat, R, C, RA} = expand(A),
    case expand(B) of
        {mat, R, C, RB} -> {mat, R, C, lists:zipwith(fun(X, Y) -> lists:zipwith(F, X, Y) end, RA, RB)};
        {mat, R2, C2, _} -> throw({bad_shape, {Op, {R, C}, {R2, C2}}})
    end.

transpose(M) ->
    {mat, R, C, Rows} = expand(M),
    {mat, C, R, transpose_rows(Rows)}.

transpose_rows([[] | _]) -> [];
transpose_rows(Rows) -> [[hd(R) || R <- Rows] | transpose_rows([tl(R) || R <- Rows])].

matmul({identity, N}, B, _Ctx) -> check_inner(N, rows(B)), B;
matmul(A, {identity, N}, _Ctx) -> check_inner(cols(A), N), A;
matmul(A, B, Ctx) ->
    {mat, R, K, RA} = expand(A),
    {mat, K2, C, RB} = expand(B),
    check_inner(K, K2),
    Cols = transpose_rows(RB),
    {mat, R, C, [[pq_sym_poly:sum(lists:zipwith(fun(X, Y) -> pq_sym_poly:mul(X, Y, Ctx) end, Row, Col))
                  || Col <- Cols] || Row <- RA]}.

check_inner(K, K) -> ok;
check_inner(K1, K2) -> throw({bad_shape, {matmul, K1, K2}}).

rows({identity, N}) -> N;
rows({mat, R, _, _}) -> R.
cols({mat, _, C, _}) -> C.

kron({identity, N}, {identity, M}, _Ctx) -> {identity, N * M};
kron(A, B, Ctx) ->
    {mat, RA, CA, XA} = expand(A),
    {mat, RB, CB, XB} = expand(B),
    {mat, RA * RB, CA * CB,
     [lists:append([[pq_sym_poly:mul(X, Y, Ctx) || Y <- RowB] || X <- RowA]) || RowA <- XA, RowB <- XB]}.

%%% Back to the canonical AST (SY 4.1)

value_to_ast({scalar, P}) -> to_ast(P);
value_to_ast({identity, N}) -> {identity, N};
value_to_ast({mat, _, _, Rows}) -> {matrix, [[to_ast(E) || E <- Row] || Row <- Rows]}.

%% {add, Terms} in print order (constants first, then by total degree, then by symbol names), or one term; a term
%% is a number or {mul, [Coeff | Factors]}.
to_ast({rpoly, Vars, R} = H) ->
    %% a polynomial kept in Rust (K2): its terms in print order straight from there, when every atom is a symbol
    case lists:all(fun({sym, _}) -> true; (_) -> false end, Vars) of
        true -> pq_nif:rpoly_to_ast(R, Vars);
        false -> to_ast(pq_sym_poly:m(H))
    end;
to_ast(P) ->
    L = maps:to_list(P),
    Terms = case byte_keys(L) of
                %% print_key ends with the monomial, so the order needs no presorting
                false -> [{M, C} || {_, M, C} <- lists:sort([{print_key(M), M, C} || {M, C} <- L])];
                Keyed -> [MC || {_, MC} <- lists:sort(Keyed)]
            end,
    case [term_ast(M, C) || {M, C} <- Terms] of
        [] -> 0;
        [T] -> T;
        Ts -> {add, Ts}
    end.

print_key({Alg, Exp} = M) ->
    Deg = lists:sum([pq_sym_exact:to_float(E) || {{sym, _}, E} <- Alg]),
    Names = [X || {{sym, X}, _} <- Alg],
    {M =/= pq_sym_poly:one_mono(), Exp =/= [], Deg, Names, M}.

%% The same order as print_key's, as binaries compared byte by byte, when every monomial is a product of symbols
%% with integer exponents (no exponential, pi, conj or radical): 591k terms sorted in 0.5 s instead of 2.4 s
%% (decision record K1). [{Key, {M, C}}], or false for the general path.
%%
%%   <<NotOne:8, Degree + 2^62:64, Name_1, 1, ..., Name_k, 1, 0, E_1 + 2^62:64, ..., E_k + 2^62:64>>
%%
%% NotOne orders the constant first (false < true); no exponential, so print_key's second element is always false;
%% the degree is an integer, which orders as print_key's float sum does; the names are atoms' UTF-8 text, which
%% orders as the atoms do, each closed by byte 1 and the list by byte 0 (a shorter list of the same prefix is
%% smaller), so a name may not contain bytes 0 or 1; with equal names the monomials differ only in exponents, in
%% the same positions.
-define(OFF, (1 bsl 62)).

byte_keys(L) ->
    try
        {Keyed, _} = lists:mapfoldl(fun byte_key/2, #{}, L),
        Keyed
    catch throw:general -> false
    end.

byte_key({{Alg, []}, _} = MC, Names0) ->
    {Ns, Es, Deg, Names} =
        lists:foldr(fun({{sym, X}, {E, 1}}, {N, Ex, D, Cache}) when is_integer(E), abs(E) < ?OFF ->
                            {B, Cache2} = name(X, Cache),
                            {[B, 1 | N], [<<(E + ?OFF):64>> | Ex], D + E, Cache2};
                       (_, _) -> throw(general)
                    end, {[], [], 0, Names0}, Alg),
    abs(Deg) < ?OFF orelse throw(general),
    NotOne = case Alg of [] -> 0; _ -> 1 end,
    {{iolist_to_binary([NotOne, <<(Deg + ?OFF):64>>, Ns, 0 | Es]), MC}, Names};
byte_key(_, _) -> throw(general).

name(X, Cache) ->
    case Cache of
        #{X := B} -> {B, Cache};
        _ ->
            B = atom_to_binary(X, utf8),
            binary:match(B, [<<0>>, <<1>>]) =:= nomatch orelse throw(general),
            {B, Cache#{X => B}}
    end.

term_ast(M, C) ->
    case M =:= pq_sym_poly:one_mono() of
        true -> coeff_to_ast(C);
        false -> {mul, [coeff_to_ast(C) | factors(M)]}
    end.

factors({Alg, Exp}) ->
    [atom_ast(K, E) || {K, E} <- Alg] ++ [{exp, exp_arg(B, X)} || {B, X} <- Exp].

exp_arg(B, X) ->
    case B =:= pq_sym_poly:one_mono() of
        true -> coeff_to_ast(X);
        false -> {mul, [coeff_to_ast(X) | factors(B)]}
    end.

atom_ast({sym, X}, {1, 1}) -> {sym, X};
atom_ast({sym, X}, E) -> {pow, {sym, X}, rat_ast(E)};
atom_ast(pi, {1, 1}) -> {const, pi};
atom_ast(pi, E) -> {pow, {const, pi}, rat_ast(E)};
atom_ast({rad, P}, {1, 1}) -> {pow, P, {rat, 1, 2}};
atom_ast({conj, X}, {1, 1}) -> {conj, {sym, X}};
atom_ast({conj, X}, E) -> {pow, {conj, {sym, X}}, rat_ast(E)};
atom_ast({base, B}, E) -> {pow, to_ast(B), rat_ast(E)};
atom_ast({cbase, B}, E) -> {conj, {pow, to_ast(B), rat_ast(E)}}.

rat_ast({N, 1}) -> N;
rat_ast({N, D}) -> {rat, N, D}.

coeff_to_ast({q, R, {0, 1}}) -> rat_ast(R);
coeff_to_ast({q, R, I}) -> {complex, rat_ast(R), rat_ast(I)};
coeff_to_ast({f, R, I}) when I == 0 -> R;
coeff_to_ast({f, R, I}) -> {complex, R, I}.
