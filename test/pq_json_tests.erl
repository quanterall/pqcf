%% FR-C09, FS Appendix B: JSON round trips and the error cases.
-module(pq_json_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

canon(C) ->
    %% the canonical form of a circuit's expressions (FR-C09 applies to circuits in canonical form)
    Can = fun(E) when is_number(E) -> E; (E) -> {ok, X} = pq_symbolic:canonicalise(E), X end,
    Instrs = [case I of
                  {unitary, Ms, {cmat, _, _, _}} = U -> _ = Ms, U;
                  {unitary, Ms, {matrix, Rows}} -> {unitary, Ms, element(2, pq_symbolic:canonicalise({matrix, Rows}))};
                  {Tag, _} -> {Tag, element(2, I)};
                  {T, M, X} -> {T, M, Can(X)};
                  {T, M, X, Y} when T =:= beamsplitter -> {T, M, X, Y};
                  {T, M, X, Y} -> {T, M, Can(X), Can(Y)};
                  {T, A, B, X, Y} -> {T, A, B, Can(X), Can(Y)}
              end || I <- C#circuit.instrs],
    C#circuit{instrs = Instrs, params = maps:map(fun(_, V) -> Can(V) end, C#circuit.params)}.

round_trip(C) ->
    {ok, J} = pq_circuit:to_json(C),
    {ok, C2} = pq_circuit:from_json(J),
    ?assertEqual(canon(C), C2),
    J.

example_b5_test() ->
    J = <<"{\"pq\": 1, \"n_modes\": 2, \"input\": {\"kind\": \"fock\", \"occupation\": [1, 0]},
           \"symbols\": {\"theta\": \"real\"}, \"params\": {\"theta\": \"pi/4\"},
           \"instructions\": [{\"op\": \"beamsplitter\", \"m1\": 0, \"m2\": 1, \"theta\": \"theta\", \"phi\": 0},
                              {\"op\": \"measure_pnr\", \"modes\": [0, 1]}]}">>,
    {ok, C} = pq_circuit:from_json(J),
    ?assertMatch(#circuit{n_modes = 2, input = {fock, [1, 0]}, symbols = #{theta := real},
                          instrs = [{beamsplitter, 0, 1, {mul, [1, {sym, theta}]}, 0}, {measure_pnr, [0, 1]}]}, C),
    ?assertEqual({ok, {mul, [{rat, 1, 4}, {const, pi}]}}, maps:find(theta, C#circuit.params)),
    {ok, #result{samples = S}} = pq_api:run(C, #{shots => 100, seed => 1}),
    ?assertEqual(100, length(S)).

round_trips_test() ->
    {ok, C0} = pq_circuit:new(3, {fock, [1, 0, 2]}),
    {ok, U} = pq_core:pack({cmat, [[{complex, 0.0, 1.0}, 0], [0, {complex, 0.6, -0.8}]]}),
    C = (C0#circuit{hbar = 2,
                    instrs = [{beamsplitter, 0, 1, {mul, [{rat, 1, 3}, {const, pi}]}, {sym, phi}},
                              {squeezing, 2, 0.25, {rat, -1, 2}}, {kerr, 1, {pow, {sym, k}, 2}},
                              {loss, 2, {sym, eta}}, {unitary, [0, 2], U},
                              {unitary, [1, 2], {matrix, [[{cos, {sym, a}}, 0], [0, 1]]}},
                              {measure_threshold, [0]}, {measure_pnr, [1, 2]}],
                    symbols = #{eta => nonnegative},
                    params = #{phi => {complex, {rat, 1, 2}, 3}, a => 0.5},
                    meta = #{name => <<"test">>, version => 3, ok => true}}),
    _ = round_trip(C),
    G = pq_state:squeezed_vacuum(1, [0.3], 2),
    {ok, Cg} = pq_circuit:new(1, {gaussian, G}),
    J = round_trip(Cg#circuit{instrs = [{measure_homodyne, 0, {mul, [{rat, 1, 2}, {const, pi}]}}]}),
    ?assertNotEqual(nomatch, binary:match(J, <<"\"kind\":\"gaussian\"">>)).

errors_test() ->
    {ok, C} = pq_circuit:new(1),
    ?assertEqual({error, {json_instruction, custom}}, pq_circuit:to_json(C#circuit{instrs = [{custom, x, #{}}]})),
    ?assertEqual({error, {json_input, meta}}, pq_circuit:to_json(C#circuit{meta = #{t => {a, b}}})),
    ?assertEqual({error, {json_version, 2}}, pq_circuit:from_json(<<"{\"pq\": 2, \"n_modes\": 1}">>)),
    ?assertMatch({error, {json_parse, _}}, pq_circuit:from_json(<<"{\"pq\": 1, ">>)),
    Bad = fun(Param) -> <<"{\"pq\":1,\"n_modes\":1,\"instructions\":[{\"op\":\"phaseshift\",\"mode\":0,\"phi\":", Param/binary, "}]}">> end,
    ?assertMatch({error, {json_parse, _}}, pq_circuit:from_json(Bad(<<"\"2 x\"">>))),
    ?assertEqual({error, {json_input, fock_state}}, pq_circuit:from_json(<<"{\"pq\":1,\"n_modes\":1,\"input\":{\"kind\":\"fock_state\"}}">>)),
    %% atom safety: a name the node has never seen is refused under allowed_symbols and creates no atom
    Fresh = <<"pqcfneveranatom", (integer_to_binary(erlang:unique_integer([positive])))/binary>>,
    ?assertEqual({error, {json_symbol, Fresh}}, pq_circuit:from_json(Bad(<<"\"", Fresh/binary, "\"">>), #{allowed_symbols => [x]})),
    ?assertError(badarg, binary_to_existing_atom(Fresh)),
    ?assertMatch({error, {json_symbol, {too_many, 3, 2}}},
                 pq_circuit:from_json(Bad(<<"\"a + b + c\"">>), #{max_json_symbols => 2})),
    ?assertEqual({error, {json_symbol, <<"sin">>}},
                 pq_circuit:from_json(<<"{\"pq\":1,\"n_modes\":1,\"symbols\":{\"sin\":\"real\"}}">>)),
    G = pq_state:vacuum(1, 2),
    {ok, Cg} = pq_circuit:new(1, {gaussian, G}),
    {ok, Jg} = pq_circuit:to_json(Cg),
    Jg2 = binary:replace(Jg, <<"\"pq\":1">>, <<"\"pq\":1,\"hbar\":1">>),
    ?assertEqual({error, {hbar_mismatch, 1, 2}}, pq_circuit:from_json(Jg2)).
