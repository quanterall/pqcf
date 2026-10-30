-module(pq_nif_tests).
-include_lib("eunit/include/eunit.hrl").

fixture(Name) ->
    Path = filename:join([filename:dirname(code:which(?MODULE)), "..", "..", "..", "..", "test", "fixtures", Name]),
    {ok, Bin} = file:read_file(case filelib:is_regular(Path) of
                                   true -> Path;
                                   false -> filename:join(["test", "fixtures", Name])
                               end),
    json:decode(Bin).

version_test() ->
    {0, 1, _, Triple} = pq_nif:nif_version(),
    ?assert(is_binary(Triple) orelse is_list(Triple)),
    ?assertEqual(ok, pq_nif:check_version()).

cancel_token_test() ->
    T = pq_nif:cancel_token_new(),
    ?assertEqual(ok, pq_nif:cancel_token_set(T)).

%% GT-R01: the raw stream equals the independent reference of test/oracle/chacha8_ref.py
gt_r01_test_() ->
    #{<<"cases">> := Cases} = fixture("gt_r01.json"),
    [{lists:flatten(io_lib:format("seed ~p stream ~p", [S, T])),
      ?_assertEqual({ok, V}, pq_nif:rng_u64(S, T, length(V)))}
     || #{<<"seed">> := S, <<"stream">> := T, <<"values">> := V} <- Cases].

%% a stream's prefix does not depend on how many values are asked for
prefix_test() ->
    {ok, Long} = pq_nif:rng_u64(9, 3, 100),
    {ok, Short} = pq_nif:rng_u64(9, 3, 17),
    ?assertEqual(Short, lists:sublist(Long, 17)).

count_limit_test() ->
    ?assertEqual({error, {dimension_limit, #{count => 65537, max => 65536}}}, pq_nif:rng_u64(1, 1, 65537)).

bad_arguments_test() ->
    ?assertError(badarg, pq_nif:rng_u64(-1, 0, 1)),
    ?assertError(badarg, pq_nif:rng_u64(1 bsl 64, 0, 1)).
