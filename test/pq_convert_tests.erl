%% The Perceval and Piquasso converters: our unitary equals the unitary of the converted components under the
%% conventions asserted in test/oracle/converter_conventions.py, and round trips preserve the circuit's action.
-module(pq_convert_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

close(A, B) ->
    {ok, {cmat, X}} = pq_core:unpack(A),
    {ok, {cmat, Y}} = pq_core:unpack(B),
    lists:all(fun({{complex, Ar, Ai}, {complex, Br, Bi}}) -> abs(Ar - Br) + abs(Ai - Bi) < 1.0e-13 end,
              lists:zip(lists:append(X), lists:append(Y))).

passive() ->
    {ok, U} = pq_core:pack({cmat, [[{complex, 0.6, 0.0}, {complex, 0.0, 0.8}], [{complex, 0.0, 0.8}, {complex, 0.6, 0.0}]]}),
    {ok, C} = pq_circuit:new(4, {fock, [1, 1, 0, 1]}),
    C#circuit{instrs = [{beamsplitter, 0, 1, 0.3, 1.1}, {beamsplitter, 2, 1, 0.7, -0.4}, {beamsplitter, 0, 3, 0.5, 0.2},
                        {phaseshift, 2, 0.9}, {unitary, [3, 1], U}, {beamsplitter, 2, 3, 1.2, 0}, {measure_pnr, [0, 1, 2, 3]}]}.

perceval_test() ->
    C = passive(),
    {ok, #{n_modes := 4, input := [1, 1, 0, 1], components := Cs}} = pq_convert:to_perceval(C, #{}),
    ?assert(lists:any(fun({K, _, _}) -> K =:= 'BS.Rx' end, Cs)),
    {ok, L} = pq_circuit:lower(C, #{}),
    Ours = pq_sampling:passive_unitary(L#circuit{instrs = lists:droplast(L#circuit.instrs)}),
    ?assert(close(Ours, pq_convert:perceval_unitary(4, Cs))),
    {ok, Back} = pq_convert:from_perceval(#{n_modes => 4, components => Cs, input => [1, 1, 0, 1]}),
    ?assert(close(Ours, pq_sampling:passive_unitary(Back))),
    %% a general BS.Rx with four phases comes back as phases around B(T/2, pi/2)
    G = [{'BS.Rx', [0, 1], #{theta => 0.7, phi_tl => 0.1, phi_bl => 0.2, phi_tr => 0.3, phi_br => 0.4}}],
    {ok, {cmat, [[A, B], [Cc, D]]}} = pq_core:unpack(pq_convert:perceval_unitary(2, G)),
    Ct = math:cos(0.35), St = math:sin(0.35),
    Ex = fun(Ph, M) -> {complex, M * math:cos(Ph), M * math:sin(Ph)} end,
    [?assert(abs(Xr - Yr) + abs(Xi - Yi) < 1.0e-14) || {{complex, Xr, Xi}, {complex, Yr, Yi}} <-
        [{A, Ex(0.4, Ct)}, {B, Ex(0.5 + math:pi() / 2, St)}, {Cc, Ex(0.5 + math:pi() / 2, St)}, {D, Ex(0.6, Ct)}]],
    {ok, Sq} = pq_circuit:new(1),
    ?assertEqual({error, {unsupported_instruction, squeezing, perceval}},
                 pq_convert:to_perceval(Sq#circuit{instrs = [{squeezing, 0, 1, 0}]}, #{})).

piquasso_test() ->
    G = pq_state:squeezed_vacuum(2, [0.3, 0], 2),
    {ok, C0} = pq_circuit:new(2, {gaussian, G}),
    C = C0#circuit{instrs = [{beamsplitter, 0, 1, {sym, t}, 0.4}, {squeezing, 1, 0.2, 0.1}, {displacement, 0, 0.5, 1.0},
                             {kerr, 1, 0.3}, {measure_pnr, [0, 1]}]},
    {ok, P} = pq_convert:to_piquasso(C, #{t => 0.7}),
    #{input := {gaussian, #{cov := [[V | _] | _]}}, instructions := [{'Beamsplitter', [0, 1], #{theta := 0.7}} | _]} = P,
    ?assert(abs(V - 2 * math:exp(-0.6)) < 1.0e-14),             % twice our covariance
    {ok, Back} = pq_convert:from_piquasso(P),
    {ok, L} = pq_circuit:lower(C, #{t => 0.7}),
    ?assertEqual(L#circuit.instrs, Back#circuit.instrs),
    ?assertEqual(pq_state:mean_cov(G), pq_state:mean_cov(element(2, Back#circuit.input))).
