%% Phase 1, boson sampling end to end through the job server (FS 9, FR-P01/P03/P04/P05, FS 12.5): D1 across
%% workers and chunk sizes, retry on a worker crash, failure after max_retries, cancellation and timeout.
-module(phase1_jobs_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

setup() -> {ok, _} = application:ensure_all_started(pq), ok.

c() ->
    {ok, C} = pq_circuit:new(4, {fock, [1, 1, 1, 0]}),
    C#circuit{instrs = [{beamsplitter, 0, 1, 0.7, 0.3}, {beamsplitter, 1, 2, 0.4, -1.0}, {phaseshift, 2, 0.5},
                        {beamsplitter, 2, 3, 1.1, 0.2}, {beamsplitter, 0, 1, 0.3, 0}, {measure_pnr, [0, 1, 2, 3]}]}.

samples(Opts) ->
    {ok, #result{samples = S}} = pq_sampling:sample(boson, c(), Opts#{seed => 77}),
    S.

jobs_test_() ->
    {setup, fun setup/0, fun(_) -> ok end,
     [{"1 worker = 8 workers = other chunk sizes = no job (D1)",
       ?_test(begin
                  One = samples(#{shots => 3000, max_inflight => 1}),
                  ?assertEqual(3000, length(One)),
                  ?assertEqual(One, samples(#{shots => 3000, max_inflight => 8})),
                  ?assertEqual(One, samples(#{shots => 3000, max_inflight => 3, chunk => 7})),
                  ?assertEqual(One, samples(#{shots => 3000, parallel => false})),
                  ?assertEqual(lists:sublist(One, 100), samples(#{shots => 100, max_inflight => 4}))
              end)},
      {"a crashed worker is retried with identical output",
       ?_test(begin
                  One = samples(#{shots => 1000, max_inflight => 4}),
                  ?assertEqual(One, samples(#{shots => 1000, max_inflight => 4, crash => [{2, 0}, {5, 0}, {5, 1}]}))
              end)},
      {"after max_retries the chunk fails",
       ?_assertMatch({error, {chunk_failed, 1, injected_crash}},
                     pq_sampling:sample(boson, c(), #{shots => 500, seed => 1, max_inflight => 2, chunk => 100,
                                                      max_retries => 1, crash => [{1, 0}, {1, 1}]}))},
      {"probabilities match the exact distribution (brute force, one chunk)",
       ?_test(begin
                  {ok, #result{samples = S}} = pq_sampling:sample(boson, c(), #{shots => 20000, seed => 5,
                                                                                algorithm => brute_force}),
                  U = pq_sampling:passive_unitary((c())#circuit{instrs = lists:droplast((c())#circuit.instrs)}),
                  {ok, {P, _}} = pq_nif:bs_probability(U, [1, 1, 1, 0], [1, 1, 1, 0], #{}),
                  F = length([x || [1, 1, 1, 0] <- S]) / 20000,
                  ?assert(abs(F - P) < 4 * math:sqrt(P * (1 - P) / 20000))
              end)},
      {"timeout",
       ?_test(begin
                  {ok, C0} = pq_circuit:new(24, {fock, lists:duplicate(14, 1) ++ lists:duplicate(10, 0)}),
                  C = C0#circuit{instrs = [{beamsplitter, J, J + 1, 0.7, 0.1} || J <- lists:seq(0, 22)] ++
                                          [{measure_pnr, lists:seq(0, 23)}]},
                  ?assertEqual({error, timeout},
                               pq_sampling:sample(boson, C, #{shots => 100000, seed => 1, max_inflight => 2,
                                                              timeout => 200}))
              end)},
      {"no stray workers or jobs",
       ?_test(begin
                  timer:sleep(100),
                  ?assertEqual([], supervisor:which_children(pq_worker_sup)),
                  ?assertEqual([], supervisor:which_children(pq_job_sup))
              end)}]}.
