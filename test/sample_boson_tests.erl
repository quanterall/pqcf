%% FS 8.4 sample_boson through the NIF: determinism (D1 across splits), the HOM dip, errors.
-module(sample_boson_tests).
-include_lib("eunit/include/eunit.hrl").

hom() ->
    H = 1 / math:sqrt(2),
    {ok, U} = pq_core:pack({cmat, [[H, -H], [H, H]]}),
    U.

rows({imat, _, _, _} = M) -> {ok, {imat, R}} = pq_core:unpack(M), R.

hom_test_() ->
    [?_assertNot(lists:member([1, 1], rows(element(1, element(2,
        pq_nif:sample_boson(hom(), [1, 1], 5000, #{seed => S, algorithm => A})))))) ||
        S <- [1, 2], A <- [clifford_clifford, brute_force]].

split_test_() ->
    [fun() ->
         {ok, {All, #{shots := 30}}} = pq_nif:sample_boson(hom(), [2, 1], 30, #{seed => 9, stream_start => 5, algorithm => A}),
         {ok, {P1, _}} = pq_nif:sample_boson(hom(), [2, 1], 12, #{seed => 9, stream_start => 5, algorithm => A}),
         {ok, {P2, _}} = pq_nif:sample_boson(hom(), [2, 1], 18, #{seed => 9, stream_start => 17, algorithm => A}),
         ?assertEqual(rows(All), rows(P1) ++ rows(P2))
     end || A <- [clifford_clifford, brute_force]].

errors_test() ->
    ?assertMatch({error, {bad_shape, _}}, pq_nif:sample_boson(hom(), [1, 1, 0], 1, #{})),
    ?assertMatch({error, {dimension_limit, #{patterns := 3, max := 2}}},
                 pq_nif:sample_boson(hom(), [1, 1], 1, #{algorithm => brute_force, max_enumeration => 2})),
    ?assertError(badarg, pq_nif:sample_boson(hom(), [1, 1], 1, #{algorithm => metropolis})),
    {ok, {Empty, _}} = pq_nif:sample_boson(hom(), [0, 0], 3, #{}),
    ?assertEqual([[0, 0], [0, 0], [0, 0]], rows(Empty)),
    T = pq_nif:cancel_token_new(),
    ok = pq_nif:cancel_token_set(T),
    ?assertEqual({error, cancelled}, pq_nif:sample_boson(hom(), [1, 1], 100, #{cancel => T})),
    %% a non-unitary matrix makes the brute-force probabilities miss 1: an internal error
    {ok, Bad} = pq_core:pack({cmat, [[1, 1], [1, 1]]}),
    ?assertError({pq_bug, {probability_sum, _}}, pq_nif:sample_boson(Bad, [1, 1], 1, #{algorithm => brute_force})).
