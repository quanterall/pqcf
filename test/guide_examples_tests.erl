%% Every example of docs/guide.md, run as written: if one of these fails, the guide is wrong.
-module(guide_examples_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

close(A, B) -> abs(A - B) < 1.0e-9.

guide_test_() ->
    {setup, fun() -> {ok, _} = application:ensure_all_started(pq) end, fun(_) -> ok end,
     [fun hom/0, fun determinism/0, fun squeezed/0, fun parameters/0, fun gradient/0, fun symbolic/0, fun json/0,
      fun loss/0, fun postselection/0, fun jobs/0, fun symbolic_gaussian/0, fun feed_forward/0, fun errors/0]}.

%% 1. Hong-Ou-Mandel
hom() ->
    {ok, C0} = pq_api:circuit(2, {fock, [1, 1]}),
    C1 = pq_circuit:add(C0, {beamsplitter, 0, 1, math:pi() / 4, 0}),
    C = pq_circuit:add(C1, {measure_pnr, [0, 1]}),
    fock = pq_api:backend_for(C, #{}),
    {ok, P, _Meta} = pq_api:probabilities(C, #{}),
    ?assert(close(maps:get([2, 0], P), 0.5)),
    ?assert(close(maps:get([0, 2], P), 0.5)),
    ?assert(maps:get([1, 1], P, 0.0) < 1.0e-12),
    {ok, #result{kind = pnr, measured = [0, 1], samples = S}} = pq_api:sample(C, #{shots => 1000, seed => 42}),
    ?assertEqual(1000, length(S)),
    ?assertEqual([], [X || [1, 1] = X <- S]).

%% 2. Reproducible sampling: the same seed gives the same samples for any number of workers and chunk size
determinism() ->
    {ok, C0} = pq_api:circuit(4, {fock, [1, 1, 1, 0]}),
    C = lists:foldl(fun(I, Acc) -> pq_circuit:add(Acc, I) end, C0,
                    [{beamsplitter, 0, 1, 0.7, 0.1}, {beamsplitter, 1, 2, 0.4, 0.0}, {beamsplitter, 2, 3, 1.1, 0.3},
                     {measure_pnr, [0, 1, 2, 3]}]),
    {ok, A} = pq_api:sample(C, #{backend => boson, shots => 2000, seed => 7, max_inflight => 1}),
    {ok, B} = pq_api:sample(C, #{backend => boson, shots => 2000, seed => 7, max_inflight => 8, chunk => 37}),
    ?assertEqual(A#result.samples, B#result.samples),
    ?assertEqual(7, maps:get(seed, A#result.metadata)).

%% 3. Squeezed vacuum on the Gaussian backend
squeezed() ->
    R = 0.5,
    {ok, C0} = pq_api:circuit(1),
    C = pq_circuit:add(pq_circuit:add(C0, {squeezing, 0, R, 0}), {measure_pnr, [0]}),
    gaussian = pq_api:backend_for(C, #{}),
    {ok, P, #{truncated_mass := T}} = pq_api:probabilities(C, #{}),
    ?assert(close(maps:get([0], P), 1 / math:cosh(R))),
    ?assertEqual(0.0, maps:get([1], P)),
    ?assert(T < 1.0e-6),
    {ok, S, _Warnings} = pq_api:evolve(C0#circuit{instrs = [{squeezing, 0, R, 0}]}, #{}),
    {ok, N} = pq_stats:mean_photon_number(S, [0]),
    ?assert(close(N, math:pow(math:sinh(R), 2))).

%% 4. Parameters: a symbol, bound per run
parameters() ->
    {ok, Th} = pq_symbolic:parse(<<"theta">>),
    {ok, C0} = pq_api:circuit(2, {fock, [1, 0]}),
    C = pq_circuit:add(pq_circuit:add(C0, {beamsplitter, 0, 1, Th, 0}), {measure_pnr, [0, 1]}),
    ?assertEqual([theta], pq_circuit:free_symbols(C)),
    ?assertMatch({error, {unbound_symbols, [theta]}}, pq_api:probabilities(C, #{})),
    {ok, P, _} = pq_api:probabilities(C, #{bindings => #{theta => math:pi() / 3}}),
    ?assert(close(maps:get([1, 0], P), 0.25)),
    {ok, [R0, R1]} = pq_api:run_batch(C, [#{theta => 0}, #{theta => math:pi() / 2}], #{shots => 10, seed => 1}),
    ?assertEqual(lists:duplicate(10, [1, 0]), R0#result.samples),
    ?assertEqual(lists:duplicate(10, [0, 1]), R1#result.samples).

%% 5. Gradients
gradient() ->
    {ok, Th} = pq_symbolic:parse(<<"theta">>),
    {ok, C0} = pq_api:circuit(2, {fock, [1, 0]}),
    C = pq_circuit:add(C0, {beamsplitter, 0, 1, Th, 0}),
    Obs = {probability, [0, 1], [1, 0]},
    T = 0.3,
    {ok, #{theta := Gn}} = pq_gradient:gradient(C, Obs, #{theta => T}, #{}),
    {ok, #{theta := Gs}} = pq_gradient:gradient(C, Obs, #{theta => T}, #{mode => symbolic}),
    ?assert(abs(Gn + math:sin(2 * T)) < 1.0e-6),
    ?assert(close(Gs, -math:sin(2 * T))),
    {ok, #{theta := E}} = pq_gradient:symbolic_gradient(C, Obs, [theta]),
    {ok, V} = pq_symbolic:evaluate(E, #{theta => T}),
    ?assert(close(re(V), -math:sin(2 * T))).

re({complex, Re, _Im}) -> Re;
re(X) when is_number(X) -> X.

%% 6. The symbolic engine
symbolic() ->
    {ok, E} = pq_symbolic:parse(<<"cos(x)^2 + sin(x)^2">>),
    {ok, One} = pq_symbolic:simplify(E),
    ?assertEqual(<<"1">>, pq_symbolic:pretty(One)),
    {ok, F} = pq_symbolic:parse(<<"x^3 * exp(i*y)">>),
    {ok, Dx} = pq_symbolic:differentiate(F, x),
    {ok, Prog} = pq_symbolic:compile([Dx], [x, y]),
    {ok, [V]} = pq_symbolic:eval_slp(Prog, #{x => 2.0, y => 0.0}),
    ?assert(close(re(V), 12.0)),
    {ok, A} = pq_symbolic:parse(<<"sqrt(x^2)">>),
    {ok, X} = pq_symbolic:parse(<<"x">>),
    ?assertMatch({ok, {false, _}}, pq_symbolic:equiv(A, X, #{numeric => 8})).

%% 7. JSON
json() ->
    {ok, Th} = pq_symbolic:parse(<<"theta">>),
    {ok, C0} = pq_api:circuit(2, {fock, [1, 1]}),
    C = pq_circuit:add(pq_circuit:add(C0, {beamsplitter, 0, 1, Th, 0}), {measure_pnr, [0, 1]}),
    {ok, Json} = pq_circuit:to_json(C),
    ?assert(is_binary(Json)),
    {ok, Back} = pq_circuit:from_json(Json),
    ?assertEqual([theta], pq_circuit:free_symbols(Back)),
    {ok, P1, _} = pq_api:probabilities(C, #{bindings => #{theta => 0.4}}),
    {ok, P2, _} = pq_api:probabilities(Back, #{bindings => #{theta => 0.4}}),
    ?assert(lists:all(fun(K) -> close(maps:get(K, P1), maps:get(K, P2)) end, maps:keys(P1))).

%% 8. Loss
loss() ->
    {ok, C0} = pq_api:circuit(1, {fock, [2]}),
    C = pq_circuit:add(pq_circuit:add(C0, {loss, 0, 0.5}), {measure_pnr, [0]}),
    fock_dm = pq_api:backend_for(C, #{}),
    {ok, P, _} = pq_api:probabilities(C, #{}),
    ?assert(close(maps:get([0], P), 0.25)),
    ?assert(close(maps:get([1], P), 0.5)),
    ?assert(close(maps:get([2], P), 0.25)),
    {ok, #result{samples = S}} = pq_api:sample(C, #{shots => 500, seed => 3, noise => trajectory}),
    ?assertEqual(500, length(S)).

%% 9. Errors are values
errors() ->
    ?assertMatch({error, {bad_shape, _}}, pq_api:circuit(0)),
    {ok, C0} = pq_api:circuit(1),
    ?assertEqual({error, {unsupported_input, shots_required}}, pq_api:sample(pq_circuit:add(C0, {measure_pnr, [0]}), #{})).

%% 10. Post-selection (heralding)
postselection() ->
    {ok, C0} = pq_api:circuit(2, {fock, [1, 0]}),
    C = pq_circuit:add(pq_circuit:add(C0, {beamsplitter, 0, 1, math:pi() / 4, 0}), {measure_pnr, [0, 1]}),
    {ok, #result{samples = S, metadata = M}} =
        pq_api:sample(C, #{shots => 100, seed => 5, postselect => {eq, 1, 0}, shots_mode => accepted}),
    ?assertEqual(lists:duplicate(100, [1, 0]), S),
    ?assertEqual(100, maps:get(shots_accepted, M)).

%% 13. Jobs
jobs() ->
    {ok, C0} = pq_api:circuit(12, {fock, lists:duplicate(8, 1) ++ lists:duplicate(4, 0)}),
    C = lists:foldl(fun(I, Acc) -> pq_circuit:add(Acc, I) end, C0,
                    [{beamsplitter, J, J + 1, 0.3 + 0.05 * J, 0.1 * L} || L <- lists:seq(0, 11), J <- lists:seq(L rem 2, 10, 2)]
                    ++ [{measure_pnr, lists:seq(0, 11)}]),
    Self = self(),
    spawn(fun() -> Self ! {done, pq_api:run(C, #{backend => boson, shots => 1000000, seed => 1, notify => Self})} end),
    Job = receive {pq_job, J} -> J end,
    {ok, #{state := running, done := Done, chunks := Total}} = pq_parallel:status(Job),
    ?assert(Done < Total),
    ok = pq_parallel:cancel(Job),
    ?assertEqual({error, cancelled}, receive {done, R} -> R end).

%% 11. Symbolic Gaussian states
symbolic_gaussian() ->
    {ok, R} = pq_symbolic:parse(<<"r">>),
    {ok, C0} = pq_api:circuit(1),
    C = pq_circuit:add(C0, {squeezing, 0, R, 0}),
    {ok, S} = pq_gaussian:symbolic_evolve(C, #{}),
    {matrix, [[Vx, _], [_, _]]} = S#gaussian_state.cov,
    <<"exp(-2*r)">> = pq_symbolic:pretty(Vx),
    {ok, L} = pq_gaussian:lower(S, #{r => 0.5}),
    {_, [[X, _], _]} = pq_state:mean_cov(L),
    ?assert(close(X, math:exp(-1.0))),
    {ok, #{r := G}} = pq_gradient:gradient(C, {mean_photon_number, [0]}, #{r => 0.5}, #{mode => symbolic}),
    ?assert(close(G, math:sinh(1.0))).

%% 10. Feed-forward
feed_forward() ->
    {ok, C0} = pq_api:circuit(3, {fock, [1, 0, 0]}),
    C = lists:foldl(fun(I, Acc) -> pq_circuit:add(Acc, I) end, C0,
                    [{beamsplitter, 0, 1, 0.6, 0}, {measure_pnr, [0]},
                     {conditional, {eq, 0, 0}, [{beamsplitter, 1, 2, 1.1, 0}]},
                     {feed_forward, {phaseshift, 2, {mul, [{const, pi}, {sym, k}]}}, #{k => 0}},
                     {measure_pnr, [1, 2]}]),
    {ok, #result{measured = [0, 1, 2], samples = S}} = pq_api:sample(C, #{shots => 1000, seed => 1}),
    ?assertEqual(1000, length(S)),
    ?assert(lists:member([0, 0, 1], S)).

