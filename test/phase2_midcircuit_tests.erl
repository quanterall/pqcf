%% Phase 2: mid-circuit measurements (FR-G08, FR-M06, FR-M07) and post-selection (FR-M05) through pq_simulator.
-module(phase2_midcircuit_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

setup() -> {ok, _} = application:ensure_all_started(pq), ok.

c(N, In, Instrs) -> {ok, C} = pq_circuit:new(N, In), C#circuit{instrs = Instrs}.

cov(S) ->
    N = length(S),
    lists:sum([A * B || [A, B] <- S]) / N.

all_test_() ->
    {setup, fun setup/0, fun(_) -> ok end,
     [{"homodyne mid-circuit: the joint distribution of terminal homodyne, D1 across workers",
       ?_test(begin
                  R = 0.4,
                  Tms = pq_circuit:two_mode_squeezing(0, 1, R, 0),
                  Term = c(2, vacuum, Tms ++ [{measure_homodyne, 0, 0}, {measure_homodyne, 1, 0}]),
                  Mid = c(2, vacuum, Tms ++ [{measure_homodyne, 0, 0}, {phaseshift, 1, 0}, {measure_homodyne, 1, 0}]),
                  {ok, #result{samples = St}} = pq_simulator:run(pq_gaussian, Term, #{shots => 40000, seed => 3}),
                  {ok, #result{samples = Sm, measured = [0, 1], kind = homodyne}} =
                      pq_simulator:run(pq_gaussian, Mid, #{shots => 40000, seed => 3}),
                  Want = math:sinh(2 * R),
                  ?assert(abs(abs(cov(St)) - Want) < 0.03),
                  ?assert(abs(cov(Sm) - cov(St)) < 0.04),
                  {ok, #result{samples = S1}} = pq_simulator:run(pq_gaussian, Mid, #{shots => 300, seed => 9, max_inflight => 1}),
                  {ok, #result{samples = S4}} = pq_simulator:run(pq_gaussian, Mid, #{shots => 300, seed => 9, max_inflight => 4}),
                  {ok, #result{samples = S0}} = pq_simulator:run(pq_gaussian, Mid, #{shots => 300, seed => 9, parallel => false}),
                  ?assertEqual(S1, S4),
                  ?assertEqual(S1, S0)
              end)},
      {"Fock heralding: zero photons in mode 0 leave two in mode 1",
       ?_test(begin
                  H = math:pi() / 4,
                  C = c(3, {fock, [1, 1, 0]}, [{beamsplitter, 0, 1, H, 0}, {measure_pnr, [0]},
                                               {beamsplitter, 1, 2, H, 0}, {measure_pnr, [1, 2]}]),
                  {ok, #result{kind = pnr, measured = [0, 1, 2], samples = S}} =
                      pq_simulator:run(pq_fock, C, #{shots => 20000, seed => 5}),
                  ?assertEqual([], [X || [_, 1, _] = X <- S, hd(X) =:= 1]),       % HOM: never one photon in mode 0
                  Zero = [[A, B] || [0, A, B] <- S],
                  F = fun(P) -> length([x || X <- Zero, X =:= P]) / length(Zero) end,
                  ?assert(abs(F([2, 0]) - 0.25) < 0.02 andalso abs(F([1, 1]) - 0.5) < 0.02),
                  ?assertEqual([], [X || [2, A, B] = X <- S, A + B > 0])
              end)},
      {"post-selection: attempts and accepted",
       ?_test(begin
                  C = c(2, {fock, [1, 1]}, [{beamsplitter, 0, 1, 0.3, 0}, {measure_pnr, [0, 1]}]),
                  P = {eq, 0, 1},
                  {ok, #result{samples = S, metadata = M}} =
                      pq_simulator:run(pq_fock, C, #{shots => 2000, seed => 1, postselect => P}),
                  ?assert(lists:all(fun([A, _]) -> A =:= 1 end, S)),
                  #{attempts := 2000, acceptance_rate := Rate} = M,
                  ?assertEqual(length(S) / 2000, Rate),
                  {ok, #result{samples = Sa, metadata = Ma}} =
                      pq_simulator:run(pq_fock, C, #{shots => 500, seed => 1, postselect => P, shots_mode => accepted}),
                  ?assertEqual(500, length(Sa)),
                  ?assert(maps:get(attempts, Ma) > 500),
                  %% the accepted samples are the first accepted attempts: the same as filtering in attempts mode
                  ?assertEqual(lists:sublist(S, 500), Sa),
                  ?assertMatch({error, {postselection_exhausted, #{attempts := 100, accepted := 0}}},
                               pq_simulator:run(pq_fock, C, #{shots => 10, seed => 1, postselect => {eq, 0, 5},
                                                              shots_mode => accepted, max_attempts => 100})),
                  {ok, #result{samples = Sn}} =
                      pq_simulator:run(pq_fock, C, #{shots => 300, seed => 2,
                                                     postselect => {'and', [{total_eq, [0, 1], 2}, {'not', {eq, 1, 2}}]}}),
                  ?assert(lists:all(fun([_, B]) -> B =/= 2 end, Sn))
              end)},
      {"placement: a measured mode is not used again; homodyne needs a mode left",
       ?_test(begin
                  ?assertMatch({error, {measured_mode_reuse, 0, 2}},
                               pq_simulator:run(pq_gaussian, c(2, vacuum, [{measure_homodyne, 0, 0}, {phaseshift, 1, 0},
                                                                           {squeezing, 0, 1, 0}]), #{shots => 1})),
                  %% two modes: measuring mode 0 mid-circuit leaves mode 1 (and mode 0 reset to vacuum)
                  ?assertMatch({ok, #result{measured = [0, 1]}},
                               pq_simulator:run(pq_gaussian, c(2, vacuum, [{squeezing, 1, 0.1, 0}, {measure_homodyne, 0, 0},
                                                                           {phaseshift, 1, 0}, {measure_homodyne, 1, 0}]),
                                                #{shots => 2, parallel => false}))
              end)}]}.
