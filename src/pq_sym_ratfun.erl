%% Rational functions in normal form (decision record SY2): N / D with N and D coprime polynomials over the
%% rationals, D primitive with integer coefficients and a positive leading coefficient. Equal rational functions have
%% equal normal forms, which the engine's canonical form alone does not give (a sum in a denominator is an opaque atom,
%% SY 3.4).
%%
%% Read from a canonical expression: sums of products of numbers, symbols to integer powers and powers of sums (the
%% opaque atoms, recursively), combined over the lcm of the denominators and reduced by gcds (pq_sym_gcd).
-module(pq_sym_ratfun).
-export([from_ast/2, to_ast/2, equal/2]).

-type poly() :: pq_sym_gcd:poly().
-type ratfun() :: {poly(), poly()}.

%% the normal form of a canonical expression over the symbols Vars; throws {symbolic_unsupported, rational} for
%% anything that is not a rational function of them with rational coefficients
-spec from_ast(term(), [atom()]) -> ratfun().
from_ast(Ast, Vars) ->
    Index = maps:from_list(lists:zip(Vars, lists:seq(1, length(Vars)))),
    rf(Ast, Index, length(Vars)).

rf(K, _, N) when is_integer(K) -> norm({const(K, N), one(N)});
rf({rat, A, B}, _, N) -> norm({const({A, B}, N), one(N)});
rf({sym, X}, Index, N) -> {mono(Index, X, 1, N), one(N)};
rf({pow, B, K}, Index, N) when is_integer(K) ->
    {Nb, Db} = rf(B, Index, N),
    case K >= 0 of
        true -> {ppow(Nb, K, N), ppow(Db, K, N)};
        false -> norm({ppow(Db, -K, N), ppow(Nb, -K, N)})
    end;
rf({add, Ts}, Index, N) -> lists:foldl(fun(T, Acc) -> add(Acc, rf(T, Index, N)) end, {#{}, one(N)}, Ts);
rf({mul, Fs}, Index, N) -> lists:foldl(fun(F, Acc) -> mul(Acc, rf(F, Index, N)) end, {one(N), one(N)}, Fs);
rf(_, _, _) -> throw({symbolic_unsupported, rational}).

%% N1/D1 + N2/D2 over lcm(D1, D2)
add({N1, _}, B) when map_size(N1) =:= 0 -> B;
add(A, {N2, _}) when map_size(N2) =:= 0 -> A;
add({N1, D1}, {N2, D2}) ->
    G = pq_sym_gcd:gcd(D1, D2),
    {ok, A1} = pq_sym_gcd:divide(D1, G),
    {ok, A2} = pq_sym_gcd:divide(D2, G),
    norm({padd(pq_sym_gcd:mul(N1, A2), pq_sym_gcd:mul(N2, A1)), pq_sym_gcd:mul(A1, D2)}).

%% (N1 N2) / (D1 D2), the cross gcds taken out first
mul({N1, _} = A, _) when map_size(N1) =:= 0 -> A;
mul(_, {N2, _} = B) when map_size(N2) =:= 0 -> B;
mul({N1, D1}, {N2, D2}) ->
    {A1, B2} = cancel(N1, D2),
    {A2, B1} = cancel(N2, D1),
    norm({pq_sym_gcd:mul(A1, A2), pq_sym_gcd:mul(B1, B2)}).

cancel(N, D) ->
    G = pq_sym_gcd:gcd(N, D),
    {ok, N2} = pq_sym_gcd:divide(N, G),
    {ok, D2} = pq_sym_gcd:divide(D, G),
    {N2, D2}.

%% gcd(N, D) divided out; D primitive with a positive leading coefficient, its content moved into N
norm({N, D}) when map_size(N) =:= 0 -> {#{}, one(length(hd(maps:keys(D))))};
norm({N, D}) ->
    {N2, D2} = cancel(N, D),
    {C, Dp} = pq_sym_gcd:primitive(D2),
    {scale(qinv(C), N2), Dp}.

%% equal rational functions? (N1 D2 = N2 D1)
-spec equal(ratfun(), ratfun()) -> boolean().
equal({N1, D1}, {N2, D2}) -> {N1, D1} =:= {N2, D2} orelse padd(pq_sym_gcd:mul(N1, D2), scale({-1, 1}, pq_sym_gcd:mul(N2, D1))) =:= #{}.

%% the normal form as an expression (not yet canonical): N, or N * D^-1
-spec to_ast(ratfun(), [atom()]) -> term().
to_ast({N, D}, Vars) ->
    case D =:= one(length(Vars)) of
        true -> pq_sym_gcd:to_ast(N, Vars);
        false -> {mul, [pq_sym_gcd:to_ast(N, Vars), {pow, pq_sym_gcd:to_ast(D, Vars), -1}]}
    end.

%%% Polynomials with rational coefficients (integers, or {N, D} with D > 1)

one(N) -> #{lists:duplicate(N, 0) => 1}.
const(0, _) -> #{};
const({A, B}, N) -> #{lists:duplicate(N, 0) => cnorm({A, B})};
const(K, N) -> #{lists:duplicate(N, 0) => K}.

mono(Index, X, K, N) ->
    I = maps:get(X, Index),
    #{[case J of I -> K; _ -> 0 end || J <- lists:seq(1, N)] => 1}.

ppow(P, K, N) -> lists:foldl(fun(_, Acc) -> pq_sym_gcd:mul(Acc, P) end, one(N), lists:seq(1, K)).

padd(A, B) ->
    maps:fold(fun(K, V, Acc) ->
                      case cadd(maps:get(K, Acc, 0), V) of
                          0 -> maps:remove(K, Acc);
                          S -> Acc#{K => S}
                      end
              end, A, B).

scale({1, 1}, P) -> P;
scale(Q, P) -> maps:map(fun(_, C) -> cmul(Q, C) end, P).

qinv({A, B}) when A < 0 -> {-B, -A};
qinv({A, B}) -> {B, A}.

q(A) when is_integer(A) -> {A, 1};
q(Q) -> Q.
cmul(A, B) -> cnorm(qm(q(A), q(B))).
cadd(A, B) when is_integer(A), is_integer(B) -> A + B;
cadd(A, B) -> {X, Y} = q(A), {Z, W} = q(B), cnorm({X * W + Z * Y, Y * W}).
qm({A, B}, {C, D}) -> {A * C, B * D}.
cnorm({N, D}) when D < 0 -> cnorm({-N, -D});
cnorm({0, _}) -> 0;
cnorm({N, D}) ->
    G = pq_sym_exact:gcd(abs(N), D),
    case D div G of 1 -> N div G; D2 -> {N div G, D2} end.
