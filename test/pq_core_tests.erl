-module(pq_core_tests).
-include_lib("proper/include/proper.hrl").
-include_lib("eunit/include/eunit.hrl").

%% FS 8.2: the exact bytes of a 1 x 2 cmat and an imat
layout_test() ->
    {ok, {cmat, 1, 2, B}} = pq_core:pack({cmat, [[{complex, 1, -2}, 3]]}),
    ?assertEqual(<<1.0:64/float-little, -2.0:64/float-little, 3.0:64/float-little, 0.0:64/float-little>>, B),
    ?assertEqual({ok, {imat, 1, 2, <<7:16/little, 65535:16/little>>}}, pq_core:pack({imat, [[7, 65535]]})),
    ?assertEqual({ok, {rvec, 1, <<0.5:64/float-little>>}}, pq_core:pack({rvec, [{rat, 1, 2}]})).

errors_test() ->
    ?assertEqual({error, {photon_number_overflow, 65536}}, pq_core:pack({imat, [[65536]]})),
    ?assertMatch({error, {bad_shape, _}}, pq_core:pack({cmat, [[1, 2], [3]]})),
    ?assertMatch({error, {bad_shape, _}}, pq_core:pack({cmat, []})),
    ?assertMatch({error, {bad_shape, _}}, pq_core:pack({rvec, [a]})),
    ?assertMatch({error, {bad_shape, #{bytes := 8, expected := 16}}}, pq_core:validate({cvec, 1, <<0:64>>})),
    ?assertEqual(ok, pq_core:validate({imat, 0, 3, <<>>})).

non_finite_test() ->
    NaN = <<0:48, 16#F8, 16#7F>>,                     % a quiet NaN, little-endian
    ?assertEqual({error, {non_finite, unpack}}, pq_core:unpack({rvec, 1, NaN})),
    ?assertEqual({error, {non_finite, unpack}}, pq_core:unpack({cvec, 1, <<0.0:64/float-little, NaN/binary>>})).

config_test() ->
    ?assertEqual(2, pq_core:get(hbar)),
    ?assertEqual(5, pq_core:get(max_photons, #{max_photons => 5})),
    ?assertEqual(16, pq_core:get(max_photons, #{})),
    ?assertEqual(maps:keys(pq_core:defaults()), maps:keys(pq_core:config())).

taxonomy_test() ->
    ?assert(pq_core:is_reason({bad_shape, x})),
    ?assert(pq_core:is_reason(cancelled)),
    ?assert(pq_core:is_reason({non_finite, input})),
    ?assertNot(pq_core:is_reason({bad_count, 1})),
    ?assertEqual(warning, maps:get(truncation, pq_core:categories())).

%% FS 12.3: pack and unpack round trip
finite() -> ?LET(X, oneof([integer(-1000000, 1000000), float()]), float(X)).

matrix(Gen) ->
    ?LET({R, C}, {integer(1, 6), integer(1, 6)}, vector(R, vector(C, Gen))).

prop_roundtrip() ->
    ?FORALL({Kind, Rows}, oneof([{cmat, matrix(?LET({A, B}, {finite(), finite()}, {complex, A, B}))},
                                 {rmat, matrix(finite())},
                                 {imat, matrix(integer(0, 65535))}]),
            begin
                {ok, P} = pq_core:pack({Kind, Rows}),
                ok =:= pq_core:validate(P) andalso {ok, {Kind, Rows}} =:= pq_core:unpack(P)
            end).

prop_roundtrip_vec() ->
    ?FORALL(V, non_empty(list(finite())),
            begin
                {ok, P} = pq_core:pack({rvec, V}),
                {ok, {rvec, V}} =:= pq_core:unpack(P)
            end).

proper_test_() ->
    {timeout, 60, [?_assert(proper:quickcheck(prop_roundtrip(), [{numtests, 500}, quiet])),
                   ?_assert(proper:quickcheck(prop_roundtrip_vec(), [{numtests, 500}, quiet]))]}.
