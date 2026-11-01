%% Phase 4, item 2: feed-forward (FS 15 item 1). {feed_forward, Gate, #{Symbol => Mode}} binds the symbols to the
%% shot's outcome on each mode; {conditional, Predicate, Gates} applies the gates when the predicate holds.
-module(phase4_feed_forward_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

c(N, In, Is) -> {ok, C} = pq_circuit:new(N, In), C#circuit{instrs = Is}.

-define(T0, 0.6).
-define(T1, 1.1).

%% a photon in mode 0 meets a beamsplitter; when mode 0 is found empty, the photon (now in mode 1) meets a second
%% beamsplitter towards mode 2
exact() ->
    C0 = math:cos(?T0), S0 = math:sin(?T0), C1 = math:cos(?T1), S1 = math:sin(?T1),
    #{[1, 0, 0] => C0 * C0, [0, 1, 0] => S0 * S0 * C1 * C1, [0, 0, 1] => S0 * S0 * S1 * S1}.

conditional() ->
    c(3, {fock, [1, 0, 0]}, [{beamsplitter, 0, 1, ?T0, 0}, {measure_pnr, [0]},
                             {conditional, {eq, 0, 0}, [{beamsplitter, 1, 2, ?T1, 0}]}, {measure_pnr, [1, 2]}]).

%% the same with the angle a function of the outcome k: T1 (1 - k)
fed() ->
    Angle = {mul, [?T1, {add, [1, {mul, [-1, {sym, k}]}]}]},
    c(3, {fock, [1, 0, 0]}, [{beamsplitter, 0, 1, ?T0, 0}, {measure_pnr, [0]},
                             {feed_forward, {beamsplitter, 1, 2, Angle, 0}, #{k => 0}}, {measure_pnr, [1, 2]}]).

ff_test_() ->
    {setup, fun() -> {ok, _} = application:ensure_all_started(pq) end, fun(_) -> ok end,
     [{timeout, 120, {"conditional gate: the exact branch distribution", fun() -> fits(conditional()) end}},
      {timeout, 120, {"feed-forward parameter: the same distribution", fun() -> fits(fed()) end}},
      {timeout, 120, {"Gaussian feed-forward: the corrected variance", fun gaussian/0}},
      {timeout, 120, {"the same samples for any chunking (D1)", fun determinism/0}},
      {"validation and refusals", fun refusals/0}]}.

fits(C) ->
    {ok, #result{measured = [0, 1, 2], samples = S}} = pq_api:sample(fock, C, #{shots => 200000, seed => 11}),
    {P, Info} = pq_chi2:p_value(exact(), pq_chi2:counts(S)),
    io:format(user, "~n  feed-forward fit: p = ~.3g ~p", [P, Info]),
    ?assert(P > 1.0e-4),
    %% without the conditional the photon would never reach mode 2
    ?assert(lists:any(fun(X) -> X =:= [0, 0, 1] end, S)).

%% two-mode squeezed vacuum; x measured on mode 0, mode 1 displaced by g x0 (displacement r = g x0 / sqrt(2 hbar)
%% moves x by g x0), then x measured on mode 1: Var = S11 + 2 g S01 + g^2 S00, the optimal g = -S01 / S00 gives the
%% conditional variance S11 - S01^2 / S00
gaussian() ->
    {ok, C0} = pq_circuit:new(2),
    Tms = pq_circuit:two_mode_squeezing(0, 1, 0.8),
    {ok, St, _} = pq_api:evolve(C0#circuit{instrs = Tms}, #{}),
    {_, Cov} = pq_state:mean_cov(St),
    [S00, S01, S11] = [lists:nth(1, lists:nth(1, Cov)), lists:nth(2, lists:nth(1, Cov)), lists:nth(2, lists:nth(2, Cov))],
    G = -S01 / S00,
    Disp = {mul, [G / 2.0, {sym, x0}]},                       % sqrt(2 hbar) = 2 at hbar = 2
    C = C0#circuit{instrs = Tms ++ [{measure_homodyne, 0, 0}, {feed_forward, {displacement, 1, Disp, 0}, #{x0 => 0}},
                                    {measure_homodyne, 1, 0}]},
    gaussian = pq_api:backend_for(C, #{}),
    N = 200000,
    {ok, #result{kind = homodyne, samples = S}} = pq_api:sample(C, #{shots => N, seed => 5}),
    X1 = [B || [_, B] <- S],
    Mean = lists:sum(X1) / N,
    Var = lists:sum([(X - Mean) * (X - Mean) || X <- X1]) / (N - 1),
    Want = S11 + 2 * G * S01 + G * G * S00,
    io:format(user, "~n  Gaussian feed-forward: variance ~.5f, closed form ~.5f (uncorrected ~.4f)", [Var, Want, S11]),
    ?assert(abs(Var - Want) < 6 * Want * math:sqrt(2 / N)),
    ?assert(abs(Mean) < 6 * math:sqrt(Want / N)),
    ?assert(Want < S11 / 2).

determinism() ->
    {ok, A} = pq_api:sample(fock, fed(), #{shots => 3000, seed => 9, max_inflight => 1}),
    {ok, B} = pq_api:sample(fock, fed(), #{shots => 3000, seed => 9, max_inflight => 8, chunk => 37}),
    ?assertEqual(A#result.samples, B#result.samples).

refusals() ->
    %% the outcome must come from an earlier measurement
    ?assertMatch({error, {unsupported_input, #{reason := not_measured_before}}},
                 pq_circuit:validate(c(2, {fock, [1, 0]}, [{feed_forward, {phaseshift, 1, {sym, k}}, #{k => 0}},
                                                          {measure_pnr, [0, 1]}]))),
    ?assertMatch({error, {unsupported_input, #{reason := not_measured_before}}},
                 pq_circuit:validate(c(2, {fock, [1, 0]}, [{measure_pnr, [0]}, {conditional, {eq, 1, 0}, [{phaseshift, 1, 0.1}]},
                                                          {measure_pnr, [1]}]))),
    %% only gates inside, and a predicate that is data
    ?assertMatch({error, {unsupported_instruction, measure_pnr, _}},
                 pq_circuit:validate(c(2, {fock, [1, 0]}, [{measure_pnr, [0]}, {conditional, {eq, 0, 0}, [{measure_pnr, [1]}]}]))),
    ?assertMatch({error, {bad_term, #{predicate := _}}},
                 pq_circuit:validate(c(2, {fock, [1, 0]}, [{measure_pnr, [0]},
                                                          {conditional, {local_fun, fun(_) -> true end}, [{phaseshift, 1, 0.1}]},
                                                          {measure_pnr, [1]}]))),
    %% the outcome symbols are not the caller's to bind, and do not count as unbound
    ?assertEqual([], pq_circuit:free_symbols(fed())),
    %% exact probabilities and the state are not defined without the outcomes
    ?assertMatch({error, {unsupported_input, #{feed_forward := true}}}, pq_api:probabilities(fed(), #{})),
    ?assertMatch({error, {unsupported_input, #{feed_forward := true}}}, pq_api:evolve(fed(), #{backend => fock})),
    %% a kind the backend lacks is refused inside feed-forward too
    ?assertMatch({error, {unsupported_instruction, kerr, gaussian}},
                 pq_api:sample(gaussian, c(2, vacuum, [{squeezing, 0, 0.3, 0}, {measure_homodyne, 0, 0},
                                                       {conditional, {ge, 0, 0.0}, [{kerr, 1, 0.1}]},
                                                       {measure_homodyne, 1, 0}]), #{shots => 10})).

%% JSON: both instructions round-trip (the outcome symbols made atoms by the same checked path as the others)
json_test() ->
    [begin
         {ok, J} = pq_circuit:to_json(C),
         {ok, Back} = pq_circuit:from_json(J),
         ?assertEqual(C#circuit.instrs, canonical(Back#circuit.instrs, C#circuit.instrs))
     end || C <- [conditional(), fed()]],
    Pred = {'and', [{'not', {in, 0, [1, 2]}}, {total_ge, [0], 0}, {le, 0, 3}]},
    C2 = c(2, {fock, [1, 0]}, [{measure_pnr, [0]}, {conditional, Pred, [{phaseshift, 1, 0.5}]}, {measure_pnr, [1]}]),
    {ok, J2} = pq_circuit:to_json(C2),
    {ok, B2} = pq_circuit:from_json(J2),
    ?assertEqual(C2#circuit.instrs, B2#circuit.instrs),
    %% an outcome symbol is checked like any symbol name before it becomes an atom
    Bad = binary:replace(element(2, pq_circuit:to_json(fed())), <<"\"k\"">>, <<"\"1k\"">>, [global]),
    ?assertMatch({error, {json_symbol, _}}, pq_circuit:from_json(Bad)).

%% parameters come back canonicalised: compare them by value
canonical(Back, Orig) ->
    [case {B, O} of
         {{feed_forward, {beamsplitter, A1, B1, T, P}, M}, {feed_forward, {beamsplitter, _, _, T0, _}, _}} ->
             {ok, Eq} = pq_symbolic:equiv(T, T0, #{}),
             true = Eq =:= true,
             {feed_forward, {beamsplitter, A1, B1, T0, P}, M};
         _ -> B
     end || {B, O} <- lists:zip(Back, Orig)].
