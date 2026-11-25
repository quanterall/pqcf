%% pq_lean: the Lean statements PQCF writes for its results (no Lean needed; lean-check's own tests and the
%% integration run inside the lean-check image check them with Lean).
-module(pq_lean_tests).
-include_lib("eunit/include/eunit.hrl").

src(R) -> {ok, #{source := S}} = R, S.

to_lean_test() ->
    ?assertEqual(<<"((x ^ 2) + (-1))">>, unicode:characters_to_binary(pq_lean:to_lean({add, [{pow, {sym, x}, 2}, -1]}))),
    ?assertEqual(<<"(((3 : ℚ) / 4) * «end»)"/utf8>>, unicode:characters_to_binary(pq_lean:to_lean({mul, [{rat, 3, 4}, {sym, 'end'}]}))),
    ?assertError({pq_lean, {not_a_polynomial, _}}, pq_lean:to_lean({pow, {sym, x}, -1})).

factor_test() ->
    {ok, E} = pq_symbolic:parse(<<"x^2 - y^2">>),
    S = src(pq_lean:factor(E, #{checker => none})),
    ?assertMatch({_, _}, binary:match(S, <<"theorem pqcf_factor (x y : ℚ) :"/utf8>>)),
    ?assertMatch({_, _}, binary:match(S, <<":= by\n  ring">>)).

gcd_test() ->
    {ok, A} = pq_symbolic:parse(<<"x^2 - 1">>),
    {ok, B} = pq_symbolic:parse(<<"x^2 + 2*x + 1">>),
    S = src(pq_lean:gcd(A, B, #{checker => none})),
    ?assertMatch({_, _}, binary:match(S, <<"theorem pqcf_gcd (x : ℚ) :"/utf8>>)),
    ?assertMatch({_, _}, binary:match(S, <<" ∧ "/utf8>>)).

linear_test() ->
    % 2 x + y = 3, x - y = 0: x = y = 1
    S = src(pq_lean:linear([{[2, 1], 3}, {[1, -1], 0}], [1, 1], #{checker => none})),
    ?assertEqual(<<"theorem pqcf_linear :\n    (((2 * 1) + (1 * 1)) : ℚ) = 3 ∧ (((1 * 1) + ((-1) * 1)) : ℚ) = 0 := by\n  norm_num\n"/utf8>>, S).
