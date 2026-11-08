%% Reducing Feynman integrals (pq_ibp, pq_ibp_family; decision record R1): families, zero sectors, symmetries,
%% reductions checked against textbook results, several targets at once, differential equations.
-module(pq_ibp_tests).
-include_lib("eunit/include/eunit.hrl").

d() -> {sym, d}.
s() -> {sym, s}.
t() -> {sym, t}.

%% N / D = the AST as functions: equal exactly at random integer points (the canonical form doesn't cancel
%% fractions; two different rational functions agree at a random point with probability about degree / range)
same({N, D}, Ast) ->
    Lhs = {mul, [pq_sym_norm:to_ast(N), {pow, pq_sym_norm:to_ast(D), -1}]},
    [?assertEqual(pq_symbolic:evaluate(Ast, Pt, #{exact => true}), pq_symbolic:evaluate(Lhs, Pt, #{exact => true}))
     || Pt <- [#{d => 7919, s => 104729, t => -1299709}, #{d => -611953, s => 15485863, t => 32452843}]].

%% a family is data, validated: its propagators must span the scalar products
family_test() ->
    Box = pq_ibp_family:example(box),
    ?assertEqual(4, pq_ibp_family:size(Box)),
    ?assertEqual([d, s, t], pq_ibp_family:vars(Box)),
    ?assertEqual({mul, [{rat, 1, 2}, s()]}, pq_ibp_family:kin(Box, p2, p1)),
    ?assertEqual(0, pq_ibp_family:kin(Box, p1, p1)),
    ?assertError({bad_family, {not_a_basis, 2, 1}},
                 pq_ibp_family:new(#{loops => [k], ext => [p], vars => [d, s], props => [{1, [{1, k}]}]})),
    ?assertError({bad_family, props},
                 pq_ibp_family:new(#{loops => [k], ext => [p], vars => [d, s], props => [{2, [{1, k}]}, {1, [{1, k}, {1, p}]}]})).

%% scaleless sectors are zero (the massless tadpole; the box needs an s- or a t-cut); symmetric integrals are named
%% alike (the kite's two sunsets)
zero_and_symmetry_test() ->
    Box = pq_ibp_family:example(box),
    ?assert(pq_ibp_family:zero(Box, {1, 0, 0, 0})),
    ?assert(pq_ibp_family:zero(Box, {1, 1, 0, 0})),
    ?assertNot(pq_ibp_family:zero(Box, {1, 0, 1, 0})),
    Kite = pq_ibp_family:example(kite),
    ?assertEqual(pq_ibp_family:canon(Kite, {1, 0, 0, 1, 1}), pq_ibp_family:canon(Kite, {0, 1, 1, 0, 1})),
    ?assert(pq_ibp_family:zero(Kite, {0, 1, 0, 1, 0})).

%% the one-loop box: I(2,1,1,1) = -(d - 5)/s Box + 4 (d - 5)(d - 3) / (s t^2 (d - 6)) Bubble_t; a zero target maps to [],
%% a master to itself
box_test_() ->
    {timeout, 120, fun() ->
        Box = pq_ibp_family:example(box),
        {ok, R, #{seeds := {dots, _, _}, masters := M}} = pq_ibp:reduce(Box, [{2, 1, 1, 1}, {1, 1, 1, 1}, {1, 0, 0, 0}], #{}),
        ?assert(M >= 2),
        ?assertEqual([], maps:get({1, 0, 0, 0}, R)),
        [{{1, 1, 1, 1}, One}] = maps:get({1, 1, 1, 1}, R),
        same(One, 1),
        Coeffs = maps:from_list(maps:get({2, 1, 1, 1}, R)),
        ?assertEqual([{0, 1, 0, 1}, {1, 1, 1, 1}], lists:sort(maps:keys(Coeffs))),
        same(maps:get({1, 1, 1, 1}, Coeffs), {mul, [-1, {add, [d(), -5]}, {pow, s(), -1}]}),
        same(maps:get({0, 1, 0, 1}, Coeffs),
             {mul, [4, {add, [d(), -5]}, {add, [d(), -3]}, {pow, s(), -1}, {pow, t(), -2}, {pow, {add, [d(), -6]}, -1}]})
    end}.

%% several targets at once give what each gives alone; explicit seeds as well as automatic ones
targets_test_() ->
    {timeout, 120, fun() ->
        Kite = pq_ibp_family:example(kite),
        {ok, Both, _} = pq_ibp:reduce(Kite, [{1, 1, 1, 1, 1}, {2, 1, 1, 1, 1}], #{}),
        {ok, One, _} = pq_ibp:reduce(Kite, [{2, 1, 1, 1, 1}], #{seeds => {dots, 1, 0}}),
        ?assertEqual(maps:get({2, 1, 1, 1, 1}, Both), maps:get({2, 1, 1, 1, 1}, One)),
        ?assertEqual(2, length(maps:get({1, 1, 1, 1, 1}, Both))),
        ?assertMatch({error, {bad_targets, [{1, 1}]}}, pq_ibp:reduce(Kite, [{1, 1}], #{}))
    end}.

%% the box's masters' differential equation in t: the t-bubble's row is (d - 4)/(2t), the s-bubble's 0
differential_equation_test_() ->
    {timeout, 120, fun() ->
        Box = pq_ibp_family:example(box),
        Masters = [pq_ibp_family:canon(Box, X) || X <- [{1, 1, 1, 1}, {1, 0, 1, 0}, {0, 1, 0, 1}]],
        {ok, [RowBox, RowS, RowT]} = pq_ibp:differential_equation(Box, Masters, t, #{}),
        [same(E, 0) || E <- RowS],
        same(lists:nth(1, RowT), 0),
        same(lists:nth(3, RowT), {mul, [{rat, 1, 2}, {add, [d(), -4]}, {pow, t(), -1}]}),
        same(lists:nth(1, RowBox), {mul, [{rat, 1, 2}, {add, [{mul, [d(), s()]}, {mul, [-6, s()]}, {mul, [-2, t()]}]},
                                          {pow, t(), -1}, {pow, {add, [s(), t()]}, -1}]}),
        ?assertMatch({error, {masters_incomplete, _}}, pq_ibp:differential_equation(Box, [hd(Masters)], t, #{}))
    end}.

%% massive propagators: the tadpole's I(2) = (d - 2) / (2 m^2) I(1), and its differential equation in m^2
massive_tadpole_test_() ->
    {timeout, 120, fun() ->
        M2 = {sym, msq},
        Tad = pq_ibp_family:new(#{loops => [k], ext => [], vars => [d, msq], props => [{1, [{1, k}], M2}]}),
        ?assertNot(pq_ibp_family:zero(Tad, {1})),
        {ok, #{{2} := [{{1}, C}]}, _} = pq_ibp:reduce(Tad, [{2}], #{}),
        Expected = {mul, [{rat, 1, 2}, {add, [d(), -2]}, {pow, M2, -1}]},
        same(C, Expected),
        {ok, [[E]]} = pq_ibp:differential_equation(Tad, [{1}], msq, #{}),
        same(E, Expected)
    end}.

%% a bubble with one massive line: the massless tadpole is scaleless, the massive one isn't
massive_bubble_test() ->
    Bub = pq_ibp_family:new(#{loops => [k], ext => [p], vars => [d, s, msq], kin => #{{p, p} => s()},
                              props => [{1, [{1, k}], {sym, msq}}, {1, [{1, k}, {1, p}]}]}),
    ?assert(pq_ibp_family:zero(Bub, {0, 1})),
    ?assertNot(pq_ibp_family:zero(Bub, {1, 0})),
    ?assertNot(pq_ibp_family:zero(Bub, {1, 1})).

%% three loops: the massless banana (four lines, five numerators); by dimensional analysis its master obeys
%% d I / ds = (3d/2 - 4) / s I
three_loop_banana_test_() ->
    {timeout, 300, fun() ->
        K = fun(L) -> L end,
        Fam = pq_ibp_family:new(#{loops => [k1, k2, k3], ext => [p], vars => [d, s], kin => #{{p, p} => s()},
                                  top => [1, 2, 3, 4],
                                  props => [{1, K([{1, k1}])}, {1, [{1, k2}]}, {1, [{1, k3}]},
                                            {1, [{1, p}, {-1, k1}, {-1, k2}, {-1, k3}]},
                                            {1, [{1, k1}, {-1, k2}]}, {1, [{1, k1}, {-1, k3}]}, {1, [{1, k2}, {-1, k3}]},
                                            {1, [{1, k1}, {-1, p}]}, {1, [{1, k2}, {-1, p}]}]}),
        Banana = {1, 1, 1, 1, 0, 0, 0, 0, 0},
        ?assertNot(pq_ibp_family:zero(Fam, Banana)),
        ?assert(pq_ibp_family:zero(Fam, {1, 1, 1, 0, 0, 0, 0, 0, 0})),
        {ok, [[E]]} = pq_ibp:differential_equation(Fam, [pq_ibp_family:canon(Fam, Banana)], s, #{}),
        same(E, {mul, [{add, [{mul, [{rat, 3, 2}, d()]}, -4]}, {pow, s(), -1}]})
    end}.

