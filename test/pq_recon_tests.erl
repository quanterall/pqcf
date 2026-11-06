%% Polynomials reconstructed from their values over prime fields (pq_recon; decision record R1, phase 2): exactly the
%% canonical form, from programs that compute them however they like.
-module(pq_recon_tests).
-include_lib("eunit/include/eunit.hrl").

%% reconstruct E from its compiled program; the canonical form of the result, and E's own canonical form
both(E, Vars) -> both(E, Vars, #{}).
both(E, Vars, Opts) ->
    {ok, Prog} = pq_symbolic:compile([E], Vars, #{exact_constants => true, optimize => maps:get(optimize, Opts, none)}),
    {ok, Poly} = pq_recon:polynomial(Prog, Vars, Opts),
    {ok, Want} = pq_symbolic:canonicalise(E),
    {pq_sym_norm:to_ast(Poly), Want}.

vars() -> [w, x, y, z].

random_poly() ->
    {add, [{mul, [{rat, rand:uniform(41) - 21, rand:uniform(9)} |
                  [{pow, {sym, lists:nth(rand:uniform(4), vars())}, rand:uniform(5)} || _ <- lists:seq(1, rand:uniform(3))]]}
           || _ <- lists:seq(1, 1 + rand:uniform(25))]}.

random_polynomials_test_() ->
    {timeout, 300, fun() ->
                           rand:seed(exsss, {3, 3, 3}),
                           [begin
                                {Got, Want} = both(random_poly(), vars(), #{optimize => case K rem 2 of 0 -> none; 1 -> horner_cse end}),
                                ?assertEqual(Want, Got)
                            end || K <- lists:seq(1, 25)]
                   end}.

%% a program that divides: (x^2 - y^2) / (x - y) is x + y; a product of sums expanded only by the reconstruction
cancellation_test() ->
    X = {sym, x}, Y = {sym, y},
    E = {mul, [{add, [{pow, X, 2}, {mul, [-1, {pow, Y, 2}]}]}, {pow, {add, [X, {mul, [-1, Y]}]}, -1}]},
    {ok, Prog} = pq_symbolic:compile([E], [x, y], #{exact_constants => true}),
    {ok, Poly} = pq_recon:polynomial(Prog, [x, y], #{}),
    ?assertEqual({add, [{mul, [1, X]}, {mul, [1, Y]}]}, pq_sym_norm:to_ast(Poly)).

%% coefficients far beyond one prime: several primes, Chinese remaindering
big_coefficients_test() ->
    E = {add, [{mul, [{rat, (1 bsl 200) + 7, 3}, {sym, x}, {pow, {sym, y}, 3}]},
               {mul, [{rat, -(1 bsl 150), (1 bsl 61) - 1}, {pow, {sym, z}, 2}]}, {rat, 5, 7}]},
    {Got, Want} = both(E, [x, y, z]),
    ?assertEqual(Want, Got).

%% the zero polynomial, a constant, one variable
edges_test() ->
    [begin {Got, Want} = both(E, [x, y]), ?assertEqual(Want, Got) end
     || E <- [{add, [{sym, x}, {mul, [-1, {sym, x}]}]}, {rat, -7, 3}, {pow, {add, [{sym, x}, 1]}, 6}]].

%% a resultant, res(3, 2): the Sylvester determinant written out, reconstructed from its expanded form's program
resultant_test_() ->
    {timeout, 120, fun() ->
                           {ok, E} = pq_symbolic:parse(<<"(a0*b2^2 - a1*b1*b2 + a2*b1^2)*a0 + (a3*b1 - a2*b2)*a1*b0 "
                                                         "+ a3^2*b0^3 - 2*a1*a3*b0*b2 + a2^2*b0^2*b2 - a3*a0*b1^3">>),
                           Vars = pq_symbolic:free_symbols(E),
                           [begin {Got, Want} = both(E, Vars, #{optimize => O}), ?assertEqual(Want, Got) end
                            || O <- [none, horner_cse]]
                   end}.

%% the same polynomial for another seed, and the same result twice
determinism_test() ->
    rand:seed(exsss, {8, 8, 8}),
    E = random_poly(),
    {ok, Prog} = pq_symbolic:compile([E], vars(), #{exact_constants => true}),
    A = pq_recon:polynomial(Prog, vars(), #{seed => 1}),
    ?assertEqual(A, pq_recon:polynomial(Prog, vars(), #{seed => 1})),
    ?assertEqual(A, pq_recon:polynomial(Prog, vars(), #{seed => 42})).
