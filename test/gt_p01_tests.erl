%% GT-P01 (FS 12.2): the permanent family against thewalrus 0.22.0 (test/oracle/gt_p01.py), with tolerances
%% relative to the scale S of AN 2.4 (amendment A-1, NFR-06: 1e-10).
-module(gt_p01_tests).
-include_lib("eunit/include/eunit.hrl").

-define(TOL, 1.0e-10).

fixture() ->
    {ok, Bin} = file:read_file(filename:join([code:lib_dir(pq), "..", "..", "..", "..", "test", "fixtures",
                                              "gt_p01.json"])),
    json:decode(Bin).

cmat(Rows) ->
    {ok, M} = pq_core:pack({cmat, [[{complex, Re, Im} || [Re, Im] <- Row] || Row <- Rows]}),
    M.

close({Re, Im}, [WRe, WIm], S) ->
    abs(Re - WRe) =< ?TOL * max(S, 1.0) andalso abs(Im - WIm) =< ?TOL * max(S, 1.0).

permanents_test_() ->
    #{<<"permanents">> := Cases} = fixture(),
    [{binary_to_list(K) ++ " n=" ++ integer_to_list(length(M)),
      fun() ->
          {ok, {V, #{scale := S}}} = pq_nif:permanent(cmat(M), #{}),
          ?assert(close(V, W, S))
      end} || #{<<"kind">> := K, <<"matrix">> := M, <<"perm">> := W} <- Cases].

repeated_test_() ->
    #{<<"repeated">> := Cases} = fixture(),
    [fun() ->
         {ok, {V, #{scale := S}}} = pq_nif:permanent_repeated(cmat(M), R, C, #{}),
         ?assert(close(V, W, S))
     end || #{<<"matrix">> := M, <<"rows">> := R, <<"cols">> := C, <<"perm">> := W} <- Cases].

probabilities_test_() ->
    #{<<"probabilities">> := Cases} = fixture(),
    [fun() ->
         {ok, {P, #{permanent_scale := S}}} = pq_nif:bs_probability(cmat(U), In, Out, #{}),
         ?assert(abs(P - W) =< ?TOL * max(S * S, 1.0e-6))
     end || #{<<"unitary">> := U, <<"input">> := In, <<"output">> := Out, <<"probability">> := W} <- Cases].

%% the check bites: one entry moved by 1e-6 is caught
power_test() ->
    #{<<"permanents">> := Cases} = fixture(),
    [#{<<"matrix">> := [[[Re, Im] | Row] | Rows], <<"perm">> := W} | _] =
        [C || #{<<"kind">> := <<"random">>, <<"matrix">> := M} = C <- Cases, length(M) =:= 8],
    {ok, {V, #{scale := S}}} = pq_nif:permanent(cmat([[[Re + 1.0e-6, Im] | Row] | Rows]), #{}),
    ?assertNot(close(V, W, S)).

%% D1 from Erlang: threads change no bit (n = 19, four blocks)
threads_test() ->
    M = cmat([[[math:sin(I * 19 + J), math:cos(I * J + 0.5)] || J <- lists:seq(1, 19)] || I <- lists:seq(1, 19)]),
    {ok, One} = pq_nif:permanent(M, #{}),
    {ok, Four} = pq_nif:permanent(M, #{threads => 4}),
    ?assertMatch({_, #{blocks := 4}}, One),
    ?assert(One =:= Four).

cancel_test() ->
    T = pq_nif:cancel_token_new(),
    ok = pq_nif:cancel_token_set(T),
    M = cmat([[[1.0 / (I + J), 0.0] || J <- lists:seq(1, 20)] || I <- lists:seq(1, 20)]),
    ?assertEqual({error, cancelled}, pq_nif:permanent(M, #{cancel => T})).

errors_test() ->
    ?assertMatch({error, {bad_shape, _}}, pq_nif:permanent({cmat, 2, 3, <<0:(6 * 128)>>}, #{})),
    ?assertMatch({error, {bad_shape, _}}, pq_nif:permanent({cmat, 2, 2, <<0:64>>}, #{})),
    ?assertMatch({error, {bad_shape, _}}, pq_nif:permanent_repeated(cmat([[[1, 0]]]), [1], [2], #{})),
    NaN = <<0:48, 16#F8, 16#7F>>,
    ?assertEqual({error, {non_finite, input}}, pq_nif:permanent({cmat, 1, 1, <<NaN/binary, 0:64>>}, #{})),
    ?assertMatch({error, {dimension_limit, _}}, pq_nif:permanent(cmat(lists:duplicate(52, lists:duplicate(52, [1, 0]))), #{})).
