%% FR-G04, SY 6.2: symbolic Gaussian evolution. Lowered symbolic states equal the numeric backend's; closed forms;
%% symbolic gradients of the mean photon number through squeezing, displacement and loss (SY 6.4).
-module(phase3_symbolic_gaussian_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

s(X) -> {sym, X}.

circuit(N, Is, Domains) ->
    {ok, C} = pq_circuit:new(N),
    C#circuit{instrs = Is, symbols = Domains}.

max_diff(A, B) ->
    {Ma, Ca} = pq_state:mean_cov(A),
    {Mb, Cb} = pq_state:mean_cov(B),
    lists:max([abs(X - Y) || {X, Y} <- lists:zip(Ma ++ lists:append(Ca), Mb ++ lists:append(Cb))]).

%% every supported instruction, symbolic parameters
circuits() ->
    {ok, U} = pq_core:pack({cmat, [[{complex, 0.6, 0.0}, {complex, 0.0, 0.8}], [{complex, 0.0, 0.8}, {complex, 0.6, 0.0}]]}),
    [circuit(1, [{squeezing, 0, s(r), 0}], #{}),
     circuit(1, [{squeezing, 0, s(r), s(p)}, {displacement, 0, s(d), s(q)}, {phaseshift, 0, s(t)}], #{}),
     circuit(2, [{squeezing, 0, s(r), 0}, {squeezing, 1, {mul, [-1, s(r)]}, s(p)}, {beamsplitter, 0, 1, s(t), s(q)},
                 {loss, 0, s(eta)}, {displacement, 1, s(d), 0}], #{eta => nonnegative}),
     circuit(3, [{displacement, 0, s(d), s(p)}, {squeezing, 1, s(r), 0}, {beamsplitter, 0, 1, s(t), 0},
                 {unitary, [1, 2], U}, {beamsplitter, 1, 2, 0.3, s(q)}, {loss, 2, s(eta)}, {phaseshift, 2, s(t)}],
             #{eta => nonnegative})].

binding() -> #{r => 0.37, p => 0.9, d => 0.55, q => -0.4, t => 1.1, eta => 0.7}.

lowered_equals_numeric_test_() ->
    [?_test(begin
                {ok, S} = pq_gaussian:symbolic_evolve(C, #{}),
                B = maps:with(pq_circuit:free_symbols(C), binding()),
                {ok, L} = pq_gaussian:lower(S, B),
                {ok, Num, _} = pq_api:evolve(C, #{bindings => B, backend => gaussian}),
                ?assert(max_diff(L, Num) < 1.0e-12)
            end) || C <- circuits()].

squeezed_vacuum_closed_form_test() ->
    {ok, #gaussian_state{mean = {matrix, M}, cov = {matrix, [[Sxx, Sxp], [Spx, Spp]]}}} =
        pq_gaussian:symbolic_evolve(circuit(1, [{squeezing, 0, s(r), 0}], #{}), #{}),
    ?assertEqual([[0], [0]], M),
    {ok, E2} = pq_symbolic:canonicalise({exp, {mul, [-2, s(r)]}}),
    {ok, E2p} = pq_symbolic:canonicalise({exp, {mul, [2, s(r)]}}),
    ?assertEqual({E2, E2p, 0, 0}, {Sxx, Spp, Sxp, Spx}).

%% n = eta (sinh^2 s + d^2) on mode 0 (hbar 2: the displacement moves x by 2 d)
gradient_test() ->
    C = circuit(1, [{squeezing, 0, s(sq), 0}, {displacement, 0, s(d), s(p)}, {loss, 0, s(eta)}], #{eta => nonnegative}),
    Obs = {mean_photon_number, [0]},
    B = #{sq => 0.4, d => 0.3, p => 0.8, eta => 0.6},
    {ok, Gs} = pq_gradient:gradient(C, Obs, B, #{mode => symbolic}),
    {ok, Gn} = pq_gradient:gradient(C, Obs, B, #{}),
    Exact = #{sq => 0.6 * math:sinh(0.8), d => 0.6 * 2 * 0.3, p => 0.0,
              eta => math:pow(math:sinh(0.4), 2) + 0.09},
    [?assert(abs(maps:get(K, Gs) - V) < 1.0e-12) || K := V <- Exact],
    [?assert(abs(maps:get(K, Gn) - V) < 1.0e-6) || K := V <- Exact],
    %% the expression itself, and a weighted sum
    {ok, #{eta := E}} = pq_gradient:symbolic_gradient(C, Obs, [eta]),
    ?assertEqual([d, sq], pq_symbolic:free_symbols(E)),
    {ok, #{d := Wd}} = pq_gradient:gradient(C, {weighted_sum, [{3, Obs}]}, B, #{mode => symbolic}),
    ?assert(abs(Wd - 3 * 0.36) < 1.0e-12).

two_mode_gradient_test() ->
    %% two-mode squeezing then a beamsplitter: symbolic and numeric gradients agree on every symbol
    {ok, C0} = pq_circuit:new(2),
    C = C0#circuit{instrs = pq_circuit:two_mode_squeezing(0, 1, s(r)) ++ [{beamsplitter, 0, 1, s(t), s(q)}, {loss, 1, s(eta)}],
                   symbols = #{eta => nonnegative}},
    Obs = {mean_photon_number, [1]},
    B = #{r => 0.5, t => 0.3, q => 0.2, eta => 0.8},
    {ok, Gs} = pq_gradient:gradient(C, Obs, B, #{mode => symbolic}),
    {ok, Gn} = pq_gradient:gradient(C, Obs, B, #{}),
    [?assert(abs(maps:get(K, Gs) - maps:get(K, Gn)) < 1.0e-6) || K <- [r, t, q, eta]].

conversions_test() ->
    S = pq_state:squeezed_vacuum(2, [0.3, 0.1], 2),
    Sym = pq_gaussian:to_symbolic(S),
    {ok, Back} = pq_gaussian:lower(Sym, #{}),
    ?assert(max_diff(S, Back) < 1.0e-15),
    %% a symbolic Gaussian input evolves further
    {ok, Mid} = pq_gaussian:symbolic_evolve(circuit(2, [{squeezing, 0, s(r), 0}], #{}), #{}),
    {ok, C1} = pq_circuit:new(2, {gaussian, Mid}),
    {ok, End} = pq_gaussian:symbolic_evolve(C1#circuit{instrs = [{beamsplitter, 0, 1, s(t), 0}]}, #{}),
    {ok, L} = pq_gaussian:lower(End, #{r => 0.2, t => 0.5}),
    {ok, Num, _} = pq_api:evolve(circuit(2, [{squeezing, 0, 0.2, 0}, {beamsplitter, 0, 1, 0.5, 0}], #{}), #{}),
    ?assert(max_diff(L, Num) < 1.0e-12),
    %% sampling needs a lowered state; unbound symbols are reported
    ?assertMatch({error, {unbound_symbols, _}}, pq_gaussian:lower(End, #{r => 0.2})).

limits_test() ->
    ?assertMatch({error, {symbolic_limit, #{max_symbolic_modes := 8}}},
                 pq_gaussian:symbolic_evolve(circuit(9, [{squeezing, 0, s(r), 0}], #{}), #{})),
    ?assertMatch({error, {unsupported_instruction, {kerr, symbolic_gaussian}}},
                 pq_gaussian:symbolic_evolve(circuit(1, [{kerr, 0, s(k)}], #{}), #{})),
    ?assertMatch({error, {unsupported_instruction, _}},
                 pq_gaussian:symbolic_evolve(circuit(2, [{measure_homodyne, 0, 0}, {squeezing, 1, s(r), 0}], #{}), #{})),
    %% terminal measurements are ignored
    {ok, _} = pq_gaussian:symbolic_evolve(circuit(1, [{squeezing, 0, s(r), 0}, {measure_pnr, [0]}], #{}), #{}),
    ?assertMatch({error, {differentiate_unsupported, _}},
                 pq_gradient:symbolic_gradient(circuit(1, [{squeezing, 0, s(r), 0}], #{}), {probability, [0], [0]}, [r])).
