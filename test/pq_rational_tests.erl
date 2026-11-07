%% Rational functions reconstructed from their values over prime fields (pq_recon:rational/3; decision record R1,
%% phase 4): N / D equal to the input as functions (cross-multiplied exactly), reduced, D normalised.
-module(pq_rational_tests).
-include_lib("eunit/include/eunit.hrl").

rat(E, Vars) -> rat(E, Vars, #{}).
rat(E, Vars, Opts) ->
    {ok, Prog} = pq_symbolic:compile([E], Vars, #{exact_constants => true}),
    pq_recon:rational(Prog, Vars, Opts).

%% N / D = Nin / Din exactly: N Din - Nin D is 0
same(Num, Den, Nin, Din) ->
    N = pq_sym_norm:to_ast(Num), D = pq_sym_norm:to_ast(Den),
    ?assertEqual({ok, 0}, pq_symbolic:canonicalise({add, [{mul, [N, Din]}, {mul, [-1, Nin, D]}]})).

x() -> {sym, x}.
y() -> {sym, y}.
z() -> {sym, z}.

simple_test() ->
    Nin = {add, [x(), {mul, [2, y()]}]}, Din = {add, [x(), {mul, [-1, y()]}, 3]},
    {ok, {N, D}} = rat({mul, [Nin, {pow, Din, -1}]}, [x, y]),
    same(N, D, Nin, Din),
    %% D normalised, degrees minimal
    ?assertEqual(3, maps:size(D)),
    ?assertEqual(2, maps:size(N)).

%% the denominator vanishes at the origin: the shift is what makes it work
pole_at_origin_test() ->
    Din = {mul, [x(), y()]},
    {ok, {N, D}} = rat({pow, Din, -1}, [x, y]),
    same(N, D, 1, Din),
    ?assertEqual(1, maps:size(D)).

%% a sum of fractions combined into one: 1/(x + y) + 1/(x - y) = 2x / (x^2 - y^2)
sum_of_fractions_test() ->
    E = {add, [{pow, {add, [x(), y()]}, -1}, {pow, {add, [x(), {mul, [-1, y()]}]}, -1}]},
    {ok, {N, D}} = rat(E, [x, y]),
    same(N, D, {mul, [2, x()]}, {add, [{pow, x(), 2}, {mul, [-1, {pow, y(), 2}]}]}),
    ?assertEqual(1, maps:size(N)),
    ?assertEqual(2, maps:size(D)).

%% a fraction that is a polynomial: D is 1
cancellation_test() ->
    E = {mul, [{add, [{pow, x(), 2}, {mul, [-1, {pow, y(), 2}]}]}, {pow, {add, [x(), {mul, [-1, y()]}]}, -1}]},
    {ok, {N, D}} = rat(E, [x, y]),
    ?assertEqual(#{{[], []} => {q, {1, 1}, {0, 1}}}, D),
    same(N, D, {add, [x(), y()]}, 1).

one_variable_test() ->
    Nin = {add, [{pow, x(), 3}, {rat, -1, 2}]}, Din = {add, [{pow, x(), 2}, {mul, [3, x()]}, {rat, 5, 7}]},
    {ok, {N, D}} = rat({mul, [Nin, {pow, Din, -1}]}, [x]),
    same(N, D, Nin, Din).

%% random rational functions in three variables, rational coefficients, possibly sharing no factor
random_test_() ->
    {timeout, 300, fun() ->
                           rand:seed(exsss, {6, 6, 6}),
                           Vars = [x, y, z],
                           Poly = fun(Terms) ->
                                          {add, [{mul, [{rat, rand:uniform(19) - 10, rand:uniform(5)} |
                                                        [{pow, lists:nth(rand:uniform(3), [x(), y(), z()]), rand:uniform(3) - 1}
                                                         || _ <- [1, 2]]]} || _ <- lists:seq(1, Terms)]}
                                  end,
                           [begin
                                Nin = Poly(1 + rand:uniform(6)), Din = Poly(1 + rand:uniform(5)),
                                case pq_symbolic:canonicalise(Din) of
                                    {ok, 0} -> ok;
                                    _ ->
                                        {ok, {N, D}} = rat({mul, [Nin, {pow, Din, -1}]}, Vars),
                                        same(N, D, Nin, Din)
                                end
                            end || _ <- lists:seq(1, 12)]
                   end}.

%% deterministic, and the same for another seed
determinism_test() ->
    E = {add, [{pow, {add, [x(), z()]}, -1}, {mul, [y(), {pow, {add, [x(), y(), 1]}, -2}]}]},
    A = rat(E, [x, y, z], #{seed => 1}),
    ?assertMatch({ok, _}, A),
    ?assertEqual(A, rat(E, [x, y, z], #{seed => 1})),
    ?assertEqual(A, rat(E, [x, y, z], #{seed => 9})).

%% the zero function: 0 / 1
zero_test() ->
    {ok, {N, D}} = rat({add, [{pow, {add, [x(), y()]}, -1}, {mul, [-1, {pow, {add, [y(), x()]}, -1}]}, {mul, [0, z()]}]}, [x, y, z]),
    ?assertEqual(#{}, N),
    ?assertEqual(#{{[], []} => {q, {1, 1}, {0, 1}}}, D).
