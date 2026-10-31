%% GT-H02 (FS 12.2): displaced-state GBS and the torontonian's sign. Coherent state P(0) = exp(-|alpha|^2),
%% P(1) = |alpha|^2 exp(-|alpha|^2) (our displacement r e^(i phi) has |alpha| = r); squeezed vacuum click probability
%% 1 - 1/cosh r; threshold probabilities with odd and even numbers of clicked modes equal to sums of photon-counting
%% probabilities (the torontonian against the hafnian).
-module(gt_h02_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

c(N, Is) -> {ok, C} = pq_circuit:new(N), C#circuit{instrs = Is}.

coherent_test() ->
    [begin
         {ok, P, _} = pq_api:probabilities(c(1, [{displacement, 0, R, Phi}, {measure_pnr, [0]}]), #{}),
         A2 = R * R,
         ?assert(abs(maps:get([0], P) - math:exp(-A2)) < 1.0e-14),
         ?assert(abs(maps:get([1], P) - A2 * math:exp(-A2)) < 1.0e-14),
         ?assert(abs(maps:get([3], P) - A2 * A2 * A2 * math:exp(-A2) / 6) < 1.0e-14)
     end || {R, Phi} <- [{0.3, 0.0}, {0.8, 1.2}, {1.5, -2.0}]].

squeezed_click_test() ->
    [begin
         {ok, P, _} = pq_api:probabilities(c(1, [{squeezing, 0, R, 0.4}, {measure_threshold, [0]}]), #{}),
         ?assert(abs(maps:get([1], P) - (1 - 1 / math:cosh(R))) < 1.0e-14),
         ?assert(abs(maps:get([0], P) - 1 / math:cosh(R)) < 1.0e-14)
     end || R <- [0.1, 0.5, 1.2]].

%% 1, 2 and 3 clicked modes of 3: the torontonian (with its (-1)^(N - |Z|) signs) against the sum of the hafnian
%% probabilities of every photon pattern with those modes occupied and the others empty, up to the truncation
odd_even_test() ->
    Is = [{squeezing, 0, 0.25, 0}, {squeezing, 1, 0.2, 0.7}, {squeezing, 2, 0.3, -0.4},
          {beamsplitter, 0, 1, 0.6, 0.3}, {beamsplitter, 1, 2, 0.9, -0.2}, {beamsplitter, 0, 1, 0.4, 0.1}],
    {ok, T, _} = pq_api:probabilities(c(3, Is ++ [{measure_threshold, [0, 1, 2]}]), #{}),
    {ok, Pn, #{truncated_mass := Tm}} = pq_api:probabilities(c(3, Is ++ [{measure_pnr, [0, 1, 2]}]), #{max_photons => 16}),
    io:format(user, "~ntruncated mass ~.3g~n", [Tm]),
    ?assert(Tm < 1.0e-6),
    [begin
         Sum = lists:sum([V || K := V <- Pn, [min(X, 1) || X <- K] =:= Clicks]),
         %% the missing patterns can only add to the photon-counting sum, by at most the truncated mass
         D = maps:get(Clicks, T) - Sum,
         ?assert(D > -1.0e-13 andalso D < Tm + 1.0e-13),
         ?assert(maps:get(Clicks, T) > 1.0e-4)
     end || Clicks <- [[1, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 1], [1, 1, 1]]],
    ?assert(abs(maps:get([0, 0, 0], T) - maps:get([0, 0, 0], Pn)) < 1.0e-14),
    ?assert(abs(lists:sum(maps:values(T)) - 1) < 1.0e-12).
