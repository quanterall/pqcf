%% The Rust kernel for products of plain polynomials (pq_sym_poly:mul/3, native/pq_nif/src/polymul.rs; decision
%% record K1): the same polynomial as the product term by term in Erlang, term for term (D1), and the Erlang path
%% whenever the kernel refuses.
-module(pq_poly_kernel_tests).
-include_lib("eunit/include/eunit.hrl").

-define(ON, #{max_poly_terms => 10000000, poly_kernel => true}).
-define(OFF, #{max_poly_terms => 10000000, poly_kernel => false}).

atoms() -> [pi, {conj, z}, {sym, x}, {sym, y}, {sym, w}].

%% a random plain polynomial: N terms over the atoms, exponents in Lo..Hi, coefficients integer, rational or Gaussian
random_poly(N, Lo, Hi) ->
    Terms = [{lists:sort([{K, {E, 1}} || K <- atoms(), E <- [Lo + rand:uniform(Hi - Lo + 1) - 1], E =/= 0]), coeff()}
             || _ <- lists:seq(1, N)],
    pq_sym_poly:sum([#{{Alg, []} => C} || {Alg, C} <- Terms]).

coeff() ->
    R = pq_sym_exact:rat(rand:uniform(41) - 21, rand:uniform(3)),
    I = case rand:uniform(3) of 1 -> pq_sym_exact:rat(rand:uniform(9) - 5, rand:uniform(4)); _ -> {0, 1} end,
    case {R, I} of {{0, 1}, {0, 1}} -> {q, {1, 1}, {0, 1}}; _ -> {q, R, I} end.

both(P, Q) -> {pq_sym_poly:mul(P, Q, ?ON), pq_sym_poly:mul(P, Q, ?OFF)}.

random_products_test() ->
    rand:seed(exsss, {7, 11, 13}),
    lists:foreach(fun(K) ->
                          {Lo, Hi} = case K rem 3 of 0 -> {-3, 3}; _ -> {0, 6} end,
                          P = random_poly(5 + K rem 23, Lo, Hi),
                          Q = random_poly(9 + K rem 17, Lo, Hi),
                          {On, Off} = both(P, Q),
                          ?assertEqual(Off, On)
                  end, lists:seq(1, 120)).

%% the kernel is really used: the NIF multiplies (x + y / 2)^2 directly, coefficients over L in lowest terms
nif_used_test() ->
    X = {sym, x}, Y = {sym, y},
    ?assertEqual({ok, #{{[{X, {2, 1}}], []} => {q, {1, 1}, {0, 1}}, {[{X, {1, 1}}, {Y, {1, 1}}], []} => {q, {1, 1}, {0, 1}},
                        {[{Y, {2, 1}}], []} => {q, {1, 4}, {0, 1}}}},
                 pq_nif:sym_poly_mul([X, Y], [1, 0, 0, 1], [2, 0, 1, 0], [1, 0, 0, 1], [2, 0, 1, 0], 4)).

%% coefficients beyond i128 (scaled or in the product), exponents that don't pack, a radical, a float, an
%% exponential: the Erlang path, the same result
refusals_test() ->
    rand:seed(exsss, {1, 2, 3}),
    Base = random_poly(12, 0, 4),
    X = pq_sym_poly:symbol(x, #{}),
    Big = pq_sym_poly:scale(Base, pq_sym_exact:c(1 bsl 100)),
    Huge = pq_sym_poly:scale(Base, pq_sym_exact:c(1 bsl 130)),
    Wide = pq_sym_poly:add(Base, pq_sym_poly:pow_int(X, 1 bsl 40, #{})),
    Rad = pq_sym_poly:add(Base, pq_sym_poly:atom({rad, 2}, {1, 1})),
    Flt = pq_sym_poly:add(Base, pq_sym_poly:const(pq_sym_exact:c(0.5))),
    Frac = pq_sym_poly:add(Base, pq_sym_poly:atom({sym, x}, {1, 2})),
    [begin {On, Off} = both(P, Q), ?assertEqual(Off, On) end
     || {P, Q} <- [{Big, Big}, {Huge, Base}, {Wide, Wide}, {Rad, Base}, {Flt, Base}, {Frac, Base}]].

%% small products (below the threshold) and the empty polynomial
edges_test() ->
    rand:seed(exsss, {4, 5, 6}),
    P = random_poly(3, 0, 3),
    ?assertEqual(#{}, pq_sym_poly:mul(P, #{}, ?ON)),
    {On, Off} = both(P, P),
    ?assertEqual(Off, On),
    Big = random_poly(40, 0, 5),
    Neg = pq_sym_poly:neg(Big),
    {On2, Off2} = both(pq_sym_poly:add(Big, pq_sym_poly:one()), pq_sym_poly:add(Neg, pq_sym_poly:one())),
    ?assertEqual(Off2, On2).

%% end to end: Fateman's f (f + 1) and a sparse product, with and without the kernel, the same canonical form
end_to_end_test() ->
    F = {pow, {add, [1, {sym, x}, {sym, y}, {sym, z}, {sym, t}]}, 6},
    G = {pow, {add, [{rat, 1, 3}, {mul, [{complex, 0, 1}, {sym, x}]}, {pow, {sym, y}, 3}, {pow, {sym, t}, -2}]}, 5},
    [begin
         {ok, A} = pq_symbolic:canonicalise(E, #{max_poly_terms => 1000000}),
         {ok, B} = pq_symbolic:canonicalise(E, #{max_poly_terms => 1000000, poly_kernel => false}),
         ?assertEqual(B, A)
     end || E <- [{mul, [F, {add, [F, 1]}]}, {mul, [G, F]}]].

%% powers: the kernel (pq_nif:sym_poly_pow) against the multinomial enumeration, term for term
powers_test() ->
    rand:seed(exsss, {9, 8, 7}),
    Ctx = fun(K) -> #{max_poly_terms => 10000000, max_expand_power => 64, poly_kernel => K} end,
    lists:foreach(fun(K) ->
                          {Lo, Hi} = case K rem 3 of 0 -> {-2, 2}; _ -> {0, 3} end,
                          P = random_poly(2 + K rem 7, Lo, Hi),
                          N = 2 + K rem 9,
                          ?assertEqual(pq_sym_poly:pow_int(P, N, Ctx(false)), pq_sym_poly:pow_int(P, N, Ctx(true)))
                  end, lists:seq(1, 60)).

%% the power kernel is really used, and refuses what it can't hold: a coefficient or L^N beyond i128
power_refusals_test() ->
    X = {sym, x}, Y = {sym, y},
    ?assertEqual({ok, #{{[{X, {2, 1}}], []} => {q, {1, 1}, {0, 1}}, {[{X, {1, 1}}, {Y, {1, 1}}], []} => {q, {1, 1}, {0, 1}},
                        {[{Y, {2, 1}}], []} => {q, {1, 4}, {0, 1}}}},
                 pq_nif:sym_poly_pow([X, Y], [1, 0, 0, 1], [2, 0, 1, 0], 2, 4)),
    Ctx = fun(K) -> #{max_poly_terms => 10000000, max_expand_power => 64, poly_kernel => K} end,
    Base = pq_sym_poly:add(pq_sym_poly:symbol(x, #{}), pq_sym_poly:add(pq_sym_poly:symbol(y, #{}), pq_sym_poly:one())),
    Big = pq_sym_poly:scale(Base, pq_sym_exact:c(1 bsl 40)),
    Thirds = pq_sym_poly:scale(Base, pq_sym_exact:c({1, 3})),
    [?assertEqual(pq_sym_poly:pow_int(P, N, Ctx(false)), pq_sym_poly:pow_int(P, N, Ctx(true)))
     || {P, N} <- [{Big, 5}, {Thirds, 60}, {Thirds, 12}]].
