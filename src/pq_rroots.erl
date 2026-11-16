%% Real roots of a polynomial with rational coefficients in an interval, to any precision: the squarefree part, Sturm's
%% sequence in exact rationals to isolate the roots, bisection in exact rationals to refine them. Meant for small
%% degrees (the singular points of one-loop integrands are roots of polynomials of degree up to 4).
%%
%% A polynomial is a list of coefficients [c0, c1, .., cn] (integers or fractions {N, D}); roots/4 returns the distinct
%% real roots strictly inside (Lo, Hi), sorted, as pq_mp numbers within 2^-P of the root.
-module(pq_rroots).
-export([roots/4, mul/2, add/2, sub/2, scale/2, eval/2]).

roots(Poly0, Lo0, Hi0, P) ->
    Poly = trim([q(C) || C <- Poly0]),
    {Lo, Hi} = {q(Lo0), q(Hi0)},
    case degree(Poly) of
        D when D =< 0 -> [];
        _ ->
            Sq = squarefree(Poly),
            Chain = chain(Sq),
            Tol = {1, 1 bsl (P + 8)},
            lists:sort(fun(A, B) -> pq_mp:cmp(A, B) =< 0 end,
                       [refine(Sq, Chain, Iv, Tol, P) || Iv <- isolate(Sq, Chain, Lo, Hi)])
    end.

%% the number of distinct roots in the open interval (A, B): Sturm counts (A, B]
count(Poly, Chain, A, B) ->
    variations(Chain, A) - variations(Chain, B) - case sgn(eval(Poly, B)) of 0 -> 1; _ -> 0 end.

%% isolating intervals: {exact, X}, or {A, B} with exactly one root in the open interval
isolate(Poly, Chain, Lo, Hi) ->
    case count(Poly, Chain, Lo, Hi) of
        0 -> [];
        1 -> [{Lo, Hi}];
        _ ->
            Mid = qmul(qadd(Lo, Hi), {1, 2}),
            Here = case sgn(eval(Poly, Mid)) of 0 -> [{exact, Mid}]; _ -> [] end,
            Here ++ isolate(Poly, Chain, Lo, Mid) ++ isolate(Poly, Chain, Mid, Hi)
    end.

refine(_, _, {exact, X}, _, P) -> pq_mp:from(X, P);
refine(Poly, Chain, {A, B}, Tol, P) ->
    M = qmul(qadd(A, B), {1, 2}),
    case qcmp(qsub(B, A), Tol) =< 0 of
        true -> pq_mp:from(M, P);
        false ->
            case sgn(eval(Poly, M)) of
                0 -> pq_mp:from(M, P);
                _ ->
                    case count(Poly, Chain, A, M) of
                        1 -> refine(Poly, Chain, {A, M}, Tol, P);
                        _ -> refine(Poly, Chain, {M, B}, Tol, P)
                    end
            end
    end.

%% Sturm's sequence: p, p', then the negated remainders
chain(P) -> sturm([deriv(P), P]).
sturm([A, B | _] = Acc) ->
    case degree(A) of
        D when D =< 0 -> lists:reverse(Acc);
        _ ->
            R = trim(neg(rem_(B, A))),
            case R of
                [] -> lists:reverse(Acc);
                _ -> sturm([R | Acc])
            end
    end.

variations(Chain, X) ->
    Signs = [S || S <- [sgn(eval(C, X)) || C <- Chain], S =/= 0],
    case Signs of [] -> 0; _ -> length([x || {S1, S2} <- lists:zip(lists:droplast(Signs), tl(Signs)), S1 =/= S2]) end.

squarefree(P) ->
    G = gcd(P, deriv(P)),
    case degree(G) of
        0 -> P;
        _ -> {Q, _} = divmod(P, G), Q
    end.

gcd(A, []) -> monic(A);
gcd(A, B) -> gcd(B, trim(rem_(A, B))).

monic(P) -> L = lists:last(P), [qdiv(C, L) || C <- P].

%%% Polynomials: [c0, .., cn] of exact rationals {N, D}

degree(P) -> length(trim(P)) - 1.

trim(P) -> lists:reverse(lists:dropwhile(fun(C) -> sgn(C) =:= 0 end, lists:reverse(P))).

deriv([]) -> [];
deriv([_ | T]) -> trim([qmul(C, {K, 1}) || {K, C} <- lists:enumerate(T)]).

eval(P, X) -> lists:foldr(fun(C, Acc) -> qadd(qmul(Acc, X), C) end, {0, 1}, P).

add(A, B) -> trim(zipl(fun qadd/2, [q(X) || X <- A], [q(X) || X <- B])).
sub(A, B) -> add(A, neg([q(X) || X <- B])).
neg(P) -> [qmul(C, {-1, 1}) || C <- P].
scale(P, K) -> trim([qmul(q(C), q(K)) || C <- P]).
mul([], _) -> [];
mul(_, []) -> [];
mul(A0, B0) ->
    A = [q(X) || X <- A0], B = [q(X) || X <- B0],
    trim(lists:foldl(fun({I, Ai}, Acc) -> add(Acc, lists:duplicate(I - 1, {0, 1}) ++ [qmul(Ai, Bj) || Bj <- B]) end,
                     [], lists:enumerate(A))).

zipl(F, [A | As], [B | Bs]) -> [F(A, B) | zipl(F, As, Bs)];
zipl(_, [], Bs) -> Bs;
zipl(_, As, []) -> As.

%% the remainder and quotient of A by B (B non-zero)
rem_(A, B) -> {_, R} = divmod(A, B), R.
divmod(A, B) ->
    Db = degree(B), Lb = lists:last(trim(B)),
    divmod(trim(A), trim(B), Db, Lb, []).
divmod(A, B, Db, Lb, Q) ->
    case A =:= [] orelse degree(A) < Db of
        true -> {qpoly(Q), A};
        false ->
            C = qdiv(lists:last(A), Lb),
            Shift = degree(A) - Db,
            Sub = lists:duplicate(Shift, {0, 1}) ++ [qmul(C, X) || X <- B],
            divmod(trim(add(A, neg(Sub))), B, Db, Lb, [{Shift, C} | Q])
    end.

qpoly([]) -> [];
qpoly(Terms) ->
    N = lists:max([K || {K, _} <- Terms]),
    [lists:foldl(fun({K2, C}, Acc) when K2 =:= K -> qadd(Acc, C); (_, Acc) -> Acc end, {0, 1}, Terms) || K <- lists:seq(0, N)].

%%% Exact rationals

q(N) when is_integer(N) -> {N, 1};
q({N, D}) when is_integer(N), is_integer(D), D =/= 0 -> norm(N, D).
norm(N, D) when D < 0 -> norm(-N, -D);
norm(N, D) -> G = max(1, gcd_i(abs(N), D)), {N div G, D div G}.
gcd_i(A, 0) -> A;
gcd_i(A, B) -> gcd_i(B, A rem B).
qadd({A, B}, {C, D}) -> norm(A * D + C * B, B * D).
qsub(X, {C, D}) -> qadd(X, {-C, D}).
qmul({A, B}, {C, D}) -> norm(A * C, B * D).
qdiv({A, B}, {C, D}) when C =/= 0 -> norm(A * D, B * C).
qcmp(X, Y) -> sgn(qsub(X, Y)).
sgn({0, _}) -> 0;
sgn({N, _}) when N > 0 -> 1;
sgn(_) -> -1.
