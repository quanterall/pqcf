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
    % integer coefficients: a reflective certificate (PqcfCert, decided by the kernel)
    S = src(pq_lean:factor(E, #{checker => none})),
    ?assertMatch({_, _}, binary:match(S, <<"-- variables: v 0 = x, v 1 = y">>)),
    ?assertMatch({_, _}, binary:match(S, <<"theorem pqcf_factor : ∀ ρ : ℕ → ℚ,"/utf8>>)),
    ?assertMatch({_, _}, binary:match(S, <<"PqcfCert.cert _ _ (by decide +kernel)">>)),
    % on request, the statement over ℚ with ring
    T = src(pq_lean:factor(E, #{checker => none, method => tactic})),
    ?assertMatch({_, _}, binary:match(T, <<"theorem pqcf_factor (x y : ℚ) :"/utf8>>)),
    ?assertMatch({_, _}, binary:match(T, <<":= by\n  ring">>)).

%% the E terms: balanced trees of sums and products
cert_term_test() ->
    I = #{x => 0, y => 1},
    ?assertEqual(<<".add (.add (.v 0) (.c 2)) (.add (.mul (.c (-3)) (.v 1)) (.pow (.v 0) 4))">>,
                 iolist_to_binary(pq_lean:cert_term({add, [{sym, x}, 2, {mul, [-3, {sym, y}]}, {pow, {sym, x}, 4}]}, I))).

gcd_test() ->
    {ok, A} = pq_symbolic:parse(<<"x^2 - 1">>),
    {ok, B} = pq_symbolic:parse(<<"x^2 + 2*x + 1">>),
    S = src(pq_lean:gcd(A, B, #{checker => none})),
    ?assertMatch({_, _}, binary:match(S, <<"theorem pqcf_gcd : ∀ ρ : ℕ → ℚ,"/utf8>>)),
    % two identities, two certificates
    ?assertMatch({_, _}, binary:match(S, <<"⟨PqcfCert.cert _ _ (by decide +kernel) ρ, PqcfCert.cert _ _ (by decide +kernel) ρ⟩"/utf8>>)).

linear_test() ->
    % 2 x + y = 3, x - y = 0: x = y = 1
    S = src(pq_lean:linear([{[2, 1], 3}, {[1, -1], 0}], [1, 1], #{checker => none})),
    ?assertEqual(<<"set_option maxHeartbeats 0 in\ntheorem pqcf_linear :\n    (((2 * 1) + (1 * 1)) : ℚ) = 3 ∧ (((1 * 1) + ((-1) * 1)) : ℚ) = 0 := by\n  norm_num\n"/utf8>>, S).
