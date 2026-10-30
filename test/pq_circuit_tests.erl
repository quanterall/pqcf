-module(pq_circuit_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

-define(C, pq_circuit).

c(N, Instrs) -> {ok, C} = ?C:new(N), C#circuit{instrs = Instrs}.

new_and_edit_test() ->
    ?assertMatch({error, {bad_shape, _}}, ?C:new(0)),
    ?assertMatch({error, {unsupported_input, _}}, ?C:new(2, {fock, [1]})),
    {ok, C0} = ?C:new(2, {fock, [1, 1]}),
    C1 = ?C:add(C0, {phaseshift, 0, 1.0}),
    {ok, C2} = ?C:insert(C1, 0, {beamsplitter, 0, 1, {sym, t}, 0}),
    ?assertMatch([{beamsplitter, 0, 1, _, _}, {phaseshift, 0, 1.0}], C2#circuit.instrs),
    {ok, C3} = ?C:replace(C2, 1, {kerr, 1, 0.5}),
    ?assertMatch([_, {kerr, 1, 0.5}], C3#circuit.instrs),
    {ok, C4} = ?C:delete(C3, 0),
    ?assertEqual([{kerr, 1, 0.5}], C4#circuit.instrs),
    ?assertEqual({error, {index_out_of_range, 5}}, ?C:insert(C4, 5, {kerr, 0, 1})),
    ?assertEqual({error, {index_out_of_range, 1}}, ?C:delete(C4, 1)),
    ?assertEqual(C1, C1#circuit{}),                    % functional: the old circuit is unchanged
    ?assertEqual([{phaseshift, 0, 1.0}], C1#circuit.instrs).

validate_test() ->
    ?assertEqual(ok, ?C:validate(c(3, [{beamsplitter, 0, 1, {sym, t}, 0}, {squeezing, 2, 0.3, 0},
                                       {loss, 1, {rat, 9, 10}}, {measure_pnr, [0, 1]}, {measure_homodyne, 2, 0}]))),
    ?assertEqual({error, {mode_out_of_range, #{mode => 3, instr_index => 0}}}, ?C:validate(c(3, [{kerr, 3, 1}]))),
    ?assertEqual({error, {duplicate_modes, #{modes => [1, 1], instr_index => 0}}},
                 ?C:validate(c(3, [{beamsplitter, 1, 1, 0, 0}]))),
    ?assertEqual({error, {measured_mode_reuse, 0, 1}}, ?C:validate(c(2, [{measure_pnr, [0]}, {phaseshift, 0, 1}]))),
    ?assertMatch({error, {domain_violation, #{param := eta, value := 1.5}}}, ?C:validate(c(1, [{loss, 0, {rat, 3, 2}}]))),
    ?assertEqual(ok, ?C:validate(c(1, [{loss, 0, {sym, eta}}]))),
    ?assertMatch({error, {bad_shape, _}}, ?C:validate(c(1, [{phaseshift, 0, {matrix, [[1]]}}]))),
    ?assertMatch({error, {unsupported_instruction, {teleport, 0}, _}}, ?C:validate(c(1, [{teleport, 0}]))),
    ?assertEqual({error, {unsupported_instruction, nope, custom}}, ?C:validate(c(1, [{custom, nope, #{}}]))),
    %% a default makes eta numeric, so its range is checked
    C = c(1, [{loss, 0, {sym, eta}}]),
    ?assertMatch({error, {domain_violation, _}}, ?C:validate(C#circuit{params = #{eta => 2}})).

unitary_test() ->
    H = 1 / math:sqrt(2),
    {ok, U} = pq_core:pack({cmat, [[H, -H], [H, H]]}),
    {ok, Bad} = pq_core:pack({cmat, [[1, 1], [0, 1]]}),
    ?assertEqual(ok, ?C:validate(c(2, [{unitary, [0, 1], U}]))),
    ?assertMatch({error, {not_unitary, #{deviation := 1.0}}}, ?C:validate(c(2, [{unitary, [0, 1], Bad}]))),
    ?assertMatch({error, {bad_shape, _}}, ?C:validate(c(3, [{unitary, [0, 1, 2], U}]))),
    ?assertEqual(ok, ?C:validate(c(2, [{unitary, [1, 0], {matrix, [[{cos, {sym, a}}, 0], [0, 1]]}}]))).

input_test() ->
    G = #gaussian_state{n_modes = 1, hbar = 2, mean = element(2, pq_core:pack({rvec, [0, 0]})),
                        cov = element(2, pq_core:pack({rmat, [[1, 0], [0, 1]]}))},
    {ok, C} = ?C:new(1, {gaussian, G}),
    ?assertEqual(ok, ?C:validate(C)),
    ?assertEqual({error, {hbar_mismatch, 1, 2}}, ?C:validate(C#circuit{hbar = 1})),
    ?assertMatch({error, {bad_shape, _}}, ?C:new(2, {gaussian, G})),
    F = #fock_state{n_modes = 2, cutoff = 2, amps = element(2, pq_core:pack({cvec, lists:duplicate(6, 0)}))},
    ?assertMatch({ok, _}, ?C:new(2, {fock_state, F})),
    ?assertMatch({error, {bad_shape, _}}, ?C:new(2, {fock_state, F#fock_state{cutoff = 3}})).

binding_test() ->
    C0 = c(2, [{beamsplitter, 0, 1, {sym, theta}, {sym, phi}}]),
    C = C0#circuit{params = #{theta => {mul, [{rat, 1, 4}, {const, pi}]}}},
    ?assertEqual({error, {unbound_symbols, [phi]}}, ?C:effective_binding(C, #{})),
    {ok, B} = ?C:effective_binding(C, #{phi => 0}),
    ?assert(abs(maps:get(theta, B) - math:pi() / 4) < 1.0e-15),
    {ok, B2} = ?C:effective_binding(C, #{phi => 0, theta => 1}),
    ?assertEqual(1.0, maps:get(theta, B2)),
    ?assertEqual([phi], ?C:free_symbols(C)),
    ?assertEqual([phi, theta], ?C:all_symbols(C)).

compose_test() ->
    A = (c(2, [{phaseshift, 0, {sym, a}}]))#circuit{params = #{a => 1}},
    B = (c(2, [{phaseshift, 1, {sym, a}}]))#circuit{params = #{a => 2}},
    ?assertEqual({error, {param_conflict, a}}, ?C:compose(A, B)),
    {ok, AB} = ?C:compose(A, B, #{rename => #{a => b}}),
    ?assertEqual([{phaseshift, 0, {sym, a}}, {phaseshift, 1, {sym, b}}], AB#circuit.instrs),
    ?assertEqual(#{a => 1, b => 2}, AB#circuit.params),
    Small = c(1, [{squeezing, 0, 0.5, 0}]),
    ?assertMatch({error, {bad_shape, _}}, ?C:compose(A, Small)),
    {ok, Emb} = ?C:compose(A, Small, #{modes => #{0 => 1}}),
    ?assertEqual({squeezing, 1, 0.5, 0}, lists:last(Emb#circuit.instrs)),
    ?assertMatch({error, {mode_out_of_range, _}}, ?C:compose(A, Small, #{modes => #{0 => 5}})),
    ?assertEqual({error, {hbar_mismatch, 1, 2}}, ?C:compose(A#circuit{hbar = 1}, A#circuit{hbar = 2})),
    ?assertEqual({error, {param_conflict, a}},
                 ?C:compose(A#circuit{symbols = #{a => real}}, A#circuit{symbols = #{a => positive}})),
    %% equal defaults written differently are not a conflict
    ?assertMatch({ok, _}, ?C:compose(A#circuit{params = #{a => {add, [1, 1]}}}, A#circuit{params = #{a => 2}})).

binary_test() ->
    C = (c(2, [{beamsplitter, 0, 1, {sym, t}, 0}, {measure_pnr, [0, 1]}]))#circuit{params = #{t => 1}},
    ?assertEqual({ok, C}, ?C:from_binary(?C:to_binary(C))),
    ?assertEqual({error, {bad_term, {version, 99}}}, ?C:from_binary(term_to_binary({pq_circuit, 99, C}))),
    ?assertEqual({error, {bad_term, not_a_circuit}},
                 ?C:from_binary(term_to_binary({pq_circuit, 1, C#circuit{meta = #{f => fun() -> ok end}}}))),
    ?assertEqual({error, {bad_term, not_a_circuit}}, ?C:from_binary(term_to_binary({pq_circuit, 1, {circuit, 1}}))),
    ?assertEqual({error, {bad_term, decode}}, ?C:from_binary(<<131, 100, 0, 5, "zzzzz">>)),
    %% an atom the node has never seen is refused, not created
    Fresh = <<131, 119, 20, "pqcf_never_an_atom_x">>,
    ?assertEqual({error, {bad_term, decode}}, ?C:from_binary(Fresh)),
    ?assertMatch({error, {mode_out_of_range, _}}, ?C:from_binary(term_to_binary({pq_circuit, 1, c(1, [{kerr, 4, 1}])}))).

%% FR-C11 against SY 10.2: the MZI's single-photon matrix gives P(stay in mode 0) = (1 - cos Theta)/2
mzi_test() ->
    Th = {sym, theta},
    Instrs = ?C:mzi(0, 1, Th, {sym, phi}),
    ?assertEqual(ok, ?C:validate(c(2, Instrs))),
    U = lists:foldl(fun(I, Acc) -> {matmul, single_photon(I), Acc} end, {identity, 2}, Instrs),
    {ok, {matrix, [[U00, _], _]}} = pq_symbolic:canonicalise(U),
    {ok, P} = pq_symbolic:canonicalise({mul, [U00, {conj, U00}]}),
    ?assertEqual(pq_symbolic:canonicalise({add, [{rat, 1, 2}, {mul, [{rat, -1, 2}, {cos, Th}]}]}), {ok, P}),
    ?assertEqual(4, length(?C:two_mode_squeezing(0, 1, 0.5, 0))),
    ?assertEqual(ok, ?C:validate(c(2, ?C:two_mode_squeezing(0, 1, {sym, r}, 0)))).

%% the single-photon matrices of FS 3.3 on two modes
single_photon({beamsplitter, 0, 1, T, P}) ->
    {matrix, [[{cos, T}, {mul, [-1, {exp, {mul, [{complex, 0, -1}, P]}}, {sin, T}]}],
              [{mul, [{exp, {mul, [{complex, 0, 1}, P]}}, {sin, T}]}, {cos, T}]]};
single_photon({phaseshift, 0, P}) -> {matrix, [[{exp, {mul, [{complex, 0, 1}, P]}}, 0], [0, 1]]}.
