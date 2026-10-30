%% FR-S06, SY 6.1, P10: the symbolic unitary of passive circuits and is_unitary.
-module(phase2_symbolic_unitary_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

c(N, Instrs) -> {ok, C} = pq_circuit:new(N), C#circuit{instrs = Instrs}.

mzi_test() ->
    Th = {sym, theta},
    {ok, {matrix, [[U00, _], _]} = U} = pq_symbolic:circuit_unitary(c(2, pq_circuit:mzi(0, 1, Th, {sym, phi}))),
    ?assertEqual({ok, true}, pq_symbolic:is_unitary(U)),
    {ok, P} = pq_symbolic:canonicalise({mul, [U00, {conj, U00}]}),
    ?assertEqual(pq_symbolic:canonicalise({add, [{rat, 1, 2}, {mul, [{rat, -1, 2}, {cos, Th}]}]}), {ok, P}),
    %% evaluated, it equals the numeric unitary of the lowered circuit
    {ok, Num} = pq_symbolic:evaluate(U, #{theta => 0.7, phi => -0.3}),
    {ok, L} = pq_circuit:lower(c(2, pq_circuit:mzi(0, 1, 0.7, -0.3)), #{}),
    {ok, {cmat, A}} = pq_core:unpack(Num),
    {ok, {cmat, B}} = pq_core:unpack(pq_sampling:passive_unitary(L)),
    [?assert(abs(Ar - Br) + abs(Ai - Bi) < 1.0e-14) || {{complex, Ar, Ai}, {complex, Br, Bi}} <- lists:zip(lists:append(A), lists:append(B))].

%% P10: random circuits of beamsplitters and phase shifts with symbolic angles are unitary; one angle's sign changed
%% in one entry, or a gate dropped from one row, is not
p10_test_() ->
    {timeout, 120,
     [?_test(begin
                 rand:seed(exsss, {K, 5, 7}),
                 N = 3,
                 Instrs = [case rand:uniform(2) of
                               1 -> A = rand:uniform(N) - 1, {beamsplitter, A, (A + 1) rem N, {sym, list_to_atom("t" ++ integer_to_list(J))}, {sym, list_to_atom("p" ++ integer_to_list(J))}};
                               2 -> {phaseshift, rand:uniform(N) - 1, {sym, list_to_atom("f" ++ integer_to_list(J))}}
                           end || J <- lists:seq(1, 5)],
                 {ok, {matrix, Rows} = U} = pq_symbolic:circuit_unitary(c(N, Instrs)),
                 ?assertEqual({ok, true}, pq_symbolic:is_unitary(U)),
                 %% break it: double one non-zero entry (a negated entry alone in its row and column is still unitary)
                 [[E | Rest] | More] = Rows,
                 Broken = {matrix, [[{mul, [2, E]} | Rest] | More]},
                 case E of
                     0 -> ok;
                     _ -> ?assertEqual({ok, false}, pq_symbolic:is_unitary(Broken))
                 end
             end) || K <- lists:seq(1, 6)]}.

limits_test() ->
    ?assertMatch({error, {unsupported_instruction, {squeezing, symbolic_unitary}}},
                 pq_symbolic:circuit_unitary(c(1, [{squeezing, 0, 1, 0}]))),
    ?assertMatch({error, {symbolic_limit, _}},
                 pq_symbolic:circuit_unitary(c(9, [{phaseshift, 0, 1}]))),
    %% an opaque base: undecidable
    ?assertEqual({ok, unknown}, pq_symbolic:is_unitary({matrix, [[{pow, {add, [1, {sym, x}]}, {rat, 1, 2}}]]})).
