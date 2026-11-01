%% Phase 3 API review against FS 7 and 9.4-9.6: the three-argument facade, pq_sampling:probabilities/3 and
%% postselect/2, job handles (notify, status, cancel), partial results on timeout, auto_throttle, telemetry events.
-module(phase3_api_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

%% 8 photons in 12 modes: long enough to watch a job run
big() ->
    {ok, C} = pq_circuit:new(12, {fock, lists:duplicate(8, 1) ++ lists:duplicate(4, 0)}),
    C#circuit{instrs = lists:append([[{beamsplitter, J, J + 1, 0.3 + 0.05 * J, 0.1 * L} || J <- lists:seq(L rem 2, 10, 2)]
                                     || L <- lists:seq(0, 11)]) ++ [{measure_pnr, lists:seq(0, 11)}]}.

hom() ->
    {ok, C} = pq_circuit:new(2, {fock, [1, 1]}),
    C#circuit{instrs = [{beamsplitter, 0, 1, math:pi() / 4, 0}, {measure_pnr, [0, 1]}]}.

api_test_() ->
    {setup, fun() -> {ok, _} = application:ensure_all_started(telemetry), {ok, _} = application:ensure_all_started(pq) end,
     fun(_) -> ok end,
     [{"three-argument facade: the backend first", fun facade/0},
      {"pq_sampling:probabilities/3 and postselect/2", fun sampling/0},
      {"notify, status and cancel", {timeout, 60, fun cancel/0}},
      {"timeout with partial results", {timeout, 60, fun partial/0}},
      {"auto_throttle lowers max_inflight instead of refusing", fun throttle/0},
      {"telemetry events of a job", fun telemetry_events/0},
      {"FS 7 helpers: two_mode_squeezing/3, from_xpxp/1", fun helpers/0}]}.

facade() ->
    C = hom(),
    {ok, A} = pq_api:run(boson, C, #{shots => 200, seed => 3}),
    {ok, B} = pq_api:run(C, #{shots => 200, seed => 3, backend => boson}),
    ?assertEqual(A, B),
    {ok, A2} = pq_api:sample(fock, C, #{shots => 50, seed => 1}),
    ?assertEqual(fock, maps:get(backend, A2#result.metadata, fock)),
    {ok, #fock_dm{}, _} = pq_api:evolve(fock_dm, C#circuit{instrs = [{beamsplitter, 0, 1, 0.3, 0}]}, #{}),
    ?assertEqual({error, {unsupported_input, shots_required}}, pq_api:sample(fock, C, #{})).

sampling() ->
    C = hom(),
    {ok, P, _} = pq_sampling:probabilities(boson, C, #{}),
    {ok, P2, _} = pq_api:probabilities(C, #{backend => boson}),
    ?assertEqual(P, P2),
    {ok, R} = pq_api:run(boson, C, #{shots => 1000, seed => 9}),
    {ok, R1} = pq_sampling:postselect(R, {eq, 0, 2}),
    ?assert(lists:all(fun(S) -> S =:= [2, 0] end, R1#result.samples)),
    N = length(R1#result.samples),
    ?assert(N > 400 andalso N < 600),
    ?assertEqual(N, maps:get(shots_accepted, R1#result.metadata)),
    ?assert(abs(maps:get(acceptance_rate, R1#result.metadata) - N / 1000) < 1.0e-12),
    %% the same as post-selecting during the run in attempts mode: the streams do not move
    {ok, R2} = pq_api:run(boson, C, #{shots => 1000, seed => 9}),
    {ok, R3} = pq_api:run(fock, C, #{shots => 1000, seed => 9, postselect => {eq, 0, 2}}),
    {ok, R4} = pq_api:run(fock, C, #{shots => 1000, seed => 9}),
    {ok, R5} = pq_sampling:postselect(R4, {eq, 0, 2}),
    ?assertEqual(R3#result.samples, R5#result.samples),
    ?assertEqual(R#result.samples, R2#result.samples),
    ?assertMatch({error, {bad_term, _}}, pq_sampling:postselect(R, {eq, 7, 0})).

cancel() ->
    Self = self(),
    spawn_link(fun() -> Self ! {res, pq_api:run(boson, big(), #{shots => 200000, seed => 5, chunk => 200, notify => Self})} end),
    Job = receive {pq_job, J} -> J after 10000 -> error(no_job) end,
    timer:sleep(50),
    {ok, #{state := running, chunks := 1000, done := D, nodes := [_]}} = pq_parallel:status(Job),
    ?assert(D < 1000),
    T0 = erlang:monotonic_time(millisecond),
    ok = pq_parallel:cancel(Job),
    ?assertEqual({error, cancelled}, receive {res, R} -> R after 10000 -> error(no_result) end),
    %% FS 12.5: cancellation latency at most 100 ms
    ?assert(erlang:monotonic_time(millisecond) - T0 =< 100),
    timer:sleep(20),
    ?assertEqual({ok, #{state => done}}, pq_parallel:status(Job)).

partial() ->
    %% a timeout that cuts the job off with part of its chunks done: a fixed fraction of one measured run was
    %% flaky (a slow first run let the cut-off one finish), so the timeout grows from 2 ms until a run is cut off
    %% with some chunks done
    Opts = #{shots => 200000, seed => 5, chunk => 200},
    {ok, Full} = pq_api:run(boson, big(), Opts),
    {S, I} = cut_off(Opts, 2),
    %% the partial samples are the full run's rows at those shot indices (D1)
    Rows = list_to_tuple(Full#result.samples),
    ?assertEqual([element(K + 1, Rows) || K <- I], S),
    ?assertEqual(length(S), length(I)),
    ?assert(length(S) > 0 andalso length(S) < 200000),
    ?assertEqual(lists:usort(I), I),
    ?assert(lists:all(fun(Row) -> lists:sum(Row) =:= 8 end, S)),
    ?assertEqual({error, timeout}, pq_api:run(boson, big(), Opts#{timeout => 1})).

cut_off(_Opts, T) when T > 60000 -> error(never_cut_off_with_chunks_done);
cut_off(Opts, T) ->
    case pq_api:run(boson, big(), Opts#{timeout => T, partial_results => true}) of
        {error, {timeout, #{partial := #{samples := [_ | _] = S, shot_indices := I}}}} -> {S, I};
        _ -> cut_off(Opts, 2 * T)
    end.

throttle() ->
    Limit = pq_core:get(spec_cache_bytes) + 10 * 1024 * 1024,
    Opts = #{shots => 2000, seed => 4, chunk => 1, max_inflight => 100000, max_job_memory => Limit},
    ?assertMatch({error, {dimension_limit, _, Limit}}, pq_api:run(boson, big(), Opts)),
    {ok, R} = pq_api:run(boson, big(), Opts#{auto_throttle => true}),
    {ok, R0} = pq_api:run(boson, big(), #{shots => 2000, seed => 4}),
    ?assertEqual(R0#result.samples, R#result.samples).

telemetry_events() ->
    Self = self(),
    Events = [[pq, job, start], [pq, job, stop], [pq, chunk, start], [pq, chunk, stop]],
    ok = telemetry:attach_many(guide_handler, Events, fun(E, M, Meta, _) -> Self ! {ev, E, M, Meta} end, []),
    try
        {ok, _} = pq_api:run(boson, hom(), #{shots => 400, seed => 1, chunk => 100}),
        Got = collect([]),
        Kinds = [E || {E, _, _} <- Got],
        ?assertEqual(1, length([x || [pq, job, start] <- Kinds])),
        ?assertEqual(4, length([x || [pq, chunk, start] <- Kinds])),
        ?assertEqual(4, length([x || [pq, chunk, stop] <- Kinds])),
        [{_, #{duration := Dur}, #{result := ok, chunks := 4}}] = [X || {[pq, job, stop], _, _} = X <- Got],
        ?assert(Dur >= 0)
    after
        telemetry:detach(guide_handler)
    end.

collect(Acc) ->
    receive {ev, E, M, Meta} -> collect([{E, M, Meta} | Acc])
    after 200 -> lists:reverse(Acc)
    end.

helpers() ->
    ?assertEqual(pq_circuit:two_mode_squeezing(0, 1, 0.5, 0), pq_circuit:two_mode_squeezing(0, 1, 0.5)),
    S = pq_state:squeezed_vacuum(2, [0.3, 0.1], 2),
    ?assertEqual(pq_state:from_xpxp(pq_state:to_xpxp(S), 2), pq_state:from_xpxp(pq_state:to_xpxp(S))).
