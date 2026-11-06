%% Polynomials kept in Rust (pq_sym_poly's {rpoly, Vars, Resource}; decision record K2): canonicalise with them
%% gives exactly the canonical form of the maps-only path (poly_kernel => false), on random expressions that mix
%% them with everything they don't cover; every handle operation equals its map operation; their memory counts
%% against sym_max_heap_words.
-module(pq_rpoly_tests).
-include_lib("eunit/include/eunit.hrl").

-define(LIMITS, #{max_poly_terms => 5000000, max_expand_power => 16, max_expr_nodes => 50000000,
                  sym_timeout_ms => 120000, sym_max_heap_words => 2000000000}).

syms() -> [{sym, x}, {sym, y}, {sym, z}, {sym, w}].

leaf() ->
    case rand:uniform(10) of
        1 -> rand:uniform(7) - 4;
        2 -> {rat, rand:uniform(9) - 5, 1 + rand:uniform(4)};
        3 -> {complex, rand:uniform(5) - 3, {rat, rand:uniform(3), 2}};
        4 -> {const, pi};
        5 -> {const, i};
        _ -> lists:nth(rand:uniform(4), syms())
    end.

%% a random expression; sums of several terms raised to powers make the large products the handles are for
expr(0) -> leaf();
expr(D) ->
    case rand:uniform(14) of
        N when N =< 3 -> {add, [expr(D - 1) || _ <- lists:seq(1, 2 + rand:uniform(3))]};
        N when N =< 6 -> {mul, [expr(D - 1) || _ <- lists:seq(1, 2 + rand:uniform(2))]};
        N when N =< 9 -> {pow, {add, [expr(D - 1) || _ <- lists:seq(1, 3 + rand:uniform(3))]}, 2 + rand:uniform(4)};
        10 -> {conj, expr(D - 1)};
        11 -> {exp, {mul, [{const, i}, lists:nth(rand:uniform(4), syms())]}};
        12 -> {pow, lists:nth(rand:uniform(4), syms()), {rat, 1, 2}};
        13 -> {pow, {add, [expr(D - 1), 1]}, -1};
        _ -> leaf()
    end.

same(E) ->
    A = pq_symbolic:canonicalise(E, ?LIMITS),
    B = pq_symbolic:canonicalise(E, maps:put(poly_kernel, false, ?LIMITS)),
    ?assertEqual(B, A).

random_expressions_test_() ->
    {timeout, 600, fun() ->
                           rand:seed(exsss, {11, 22, 33}),
                           [same(expr(1 + K rem 2)) || K <- lists:seq(1, 200)]
                   end}.

%% the benchmarks' shapes: Fateman's f (f + 1), a sparse product, a power of a sum, with rationals and i
benchmarks_test_() ->
    {timeout, 300, fun() ->
                           F = {pow, {add, [1 | syms()]}, 6},
                           G = {pow, {add, [{rat, 1, 3}, {mul, [{const, i}, {sym, x}]}, {pow, {sym, y}, 3}, {pow, {sym, z}, -2}]}, 5},
                           [same(E) || E <- [{mul, [F, {add, [F, 1]}]}, {mul, [F, G]}, {pow, {add, [{sym, x}, {mul, [2, {sym, y}]}, {rat, 1, 2}]}, 12},
                                             {add, [{mul, [F, G]}, {mul, [-1, F, G]}]}, {mul, [F, {pow, F, -1}]}]]
                   end}.

%% a matrix product of polynomial entries, and the result through matrices
matrix_test() ->
    [X, Y, Z, W] = syms(),
    Row = fun(A, B) -> [{pow, {add, [A, B, 1]}, 3}, {pow, {add, [A, {mul, [-1, B]}, {const, i}]}, 3}] end,
    M = {matrix, [Row(X, Y), Row(Z, W)]},
    same({matmul, M, M}),
    same({kron, M, {matrix, [[X, Y], [Z, W]]}}).

%% the handle operations against their map operations
rand_poly(N) ->
    Ctx = #{max_poly_terms => 1000000},
    {scalar, P} = pq_sym_norm:norm({add, [{mul, [{rat, rand:uniform(9) - 5, rand:uniform(3)} | [{pow, lists:nth(rand:uniform(4), syms()), rand:uniform(4) - 1} || _ <- [1, 2]]]}
                                          || _ <- lists:seq(1, N)]}, Ctx),
    P.

operations_test() ->
    rand:seed(exsss, {1, 1, 2}),
    On = #{max_poly_terms => 1000000, max_expand_power => 16, rpoly => true},
    Off = #{max_poly_terms => 1000000, max_expand_power => 16, poly_kernel => false},
    M = fun pq_sym_poly:m/1,
    [begin
         P = rand_poly(9 + K), Q = rand_poly(12),
         H = pq_sym_poly:mul(P, Q, On),
         ?assertMatch({rpoly, _, _}, H),
         Ref = pq_sym_poly:mul(P, Q, Off),
         ?assertEqual(Ref, M(H)),
         ?assertEqual(pq_sym_poly:terms(Ref), pq_sym_poly:terms(H)),
         ?assertEqual(pq_sym_poly:add(Ref, P), M(pq_sym_poly:add(H, P))),
         ?assertEqual(pq_sym_poly:add(P, Ref), M(pq_sym_poly:add(P, H))),
         ?assertEqual(#{}, M(pq_sym_poly:add(H, pq_sym_poly:neg(H)))),
         ?assertEqual(pq_sym_poly:scale(Ref, {q, {-3, 7}, {2, 5}}), M(pq_sym_poly:scale(H, {q, {-3, 7}, {2, 5}}))),
         ?assertEqual(pq_sym_poly:scale(Ref, {f, 0.5, 0.0}), M(pq_sym_poly:scale(H, {f, 0.5, 0.0}))),
         ?assertEqual(pq_sym_poly:mul(Ref, Ref, Off), M(pq_sym_poly:mul(H, H, On))),
         ?assertEqual(pq_sym_poly:mul(Ref, P, Off), M(pq_sym_poly:mul(H, P, On))),
         ?assertEqual(pq_sym_poly:pow_int(P, 3, Off), M(pq_sym_poly:pow_int(P, 3, On))),
         ?assertEqual(pq_sym_poly:pow_int(Ref, 2, Off), M(pq_sym_poly:pow_int(H, 2, On))),
         ?assertEqual(pq_sym_poly:conj(Ref, Off), pq_sym_poly:conj(H, On)),
         ?assertEqual(pq_sym_norm:to_ast(Ref), pq_sym_norm:to_ast(H)),
         ?assertEqual(pq_sym_poly:is_const(Ref), pq_sym_poly:is_const(H))
     end || K <- lists:seq(1, 20)].

%% a handle's memory counts against sym_max_heap_words
memory_limit_test() ->
    F = {pow, {add, [1 | syms()]}, 10},
    ?assertMatch({error, {symbolic_limit, #{reason := heap}}},
                 pq_symbolic:canonicalise({mul, [F, {add, [F, 1]}]}, maps:put(sym_max_heap_words, 50000, ?LIMITS))).
