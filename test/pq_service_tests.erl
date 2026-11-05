%% The JSON service (docs/decisions/C1-python-client.md): requests and responses, atom safety, errors as values.
%% The Python client's own tests (python/tests) exercise both transports end to end.
-module(pq_service_tests).
-include_lib("eunit/include/eunit.hrl").

call(Req) -> json:decode(pq_service:handle(iolist_to_binary(json:encode(Req)))).

hom() ->
    #{<<"pq">> => 1, <<"n_modes">> => 2, <<"input">> => #{<<"kind">> => <<"fock">>, <<"occupation">> => [1, 1]},
      <<"instructions">> => [#{<<"op">> => <<"beamsplitter">>, <<"m1">> => 0, <<"m2">> => 1, <<"theta">> => <<"pi/4">>, <<"phi">> => 0},
                             #{<<"op">> => <<"measure_pnr">>, <<"modes">> => [0, 1]}]}.

setup_test_() ->
    {setup, fun() -> {ok, _} = application:ensure_all_started(pq) end, fun(_) -> ok end,
     [fun requests/0, fun atom_safety/0, fun errors/0, fun sizes/0]}.

requests() ->
    #{<<"id">> := 1, <<"ok">> := #{<<"service">> := 1}} = call(#{<<"id">> => 1, <<"op">> => <<"version">>}),
    #{<<"ok">> := #{<<"probabilities">> := P}} = call(#{<<"op">> => <<"probabilities">>, <<"circuit">> => hom()}),
    #{[2, 0] := A, [0, 2] := B} = maps:from_list([{K, V} || [K, V] <- P]),
    ?assert(abs(A - 0.5) < 1.0e-14 andalso abs(B - 0.5) < 1.0e-14),
    %% the circuit may also come as its JSON text
    #{<<"ok">> := #{<<"probabilities">> := P}} =
        call(#{<<"op">> => <<"probabilities">>, <<"circuit">> => iolist_to_binary(json:encode(hom()))}),
    #{<<"ok">> := #{<<"samples">> := S, <<"measured">> := [0, 1]}} =
        call(#{<<"op">> => <<"sample">>, <<"circuit">> => hom(), <<"options">> => #{<<"shots">> => 50, <<"seed">> => 3}}),
    ?assertEqual(50, length(S)),
    %% the same samples as the Erlang API
    {ok, R} = pq_api:sample(element(2, pq_circuit:from_json(iolist_to_binary(json:encode(hom())))), #{shots => 50, seed => 3}),
    ?assertEqual(element(4, R), S).

%% names in a request never become atoms: unknown options and ops, enum values, binding names, observables
atom_safety() ->
    Names = [<<"svc_atom_", (integer_to_binary(K))/binary>> || K <- lists:seq(1, 6)],
    [N1, N2, N3, N4, N5, N6] = Names,
    _ = call(#{<<"op">> => N1}),
    _ = call(#{<<"op">> => <<"probabilities">>, <<"circuit">> => hom(), <<"options">> => #{N2 => 1}}),
    _ = call(#{<<"op">> => <<"probabilities">>, <<"circuit">> => hom(), <<"options">> => #{<<"backend">> => N3}}),
    _ = call(#{<<"op">> => <<"probabilities">>, <<"circuit">> => hom(), <<"options">> => #{<<"bindings">> => #{N4 => 1}}}),
    _ = call(#{<<"op">> => <<"gradient">>, <<"circuit">> => hom(), <<"binding">> => #{N5 => 1},
               <<"observable">> => #{<<"mean_photon_number">> => #{<<"modes">> => [0]}}}),
    _ = call(#{<<"op">> => <<"gradient">>, <<"circuit">> => hom(), <<"binding">> => #{},
               <<"observable">> => #{<<"weighted_sum">> => [[N6, #{<<"mean_photon_number">> => #{<<"modes">> => [0]}}]]}}),
    [?assertError(badarg, binary_to_existing_atom(N)) || N <- Names].

errors() ->
    Cat = fun(Resp) -> maps:get(<<"category">>, maps:get(<<"error">>, Resp)) end,
    ?assertEqual(<<"bad_term">>, Cat(json:decode(pq_service:handle(<<"not json">>)))),
    ?assertEqual(<<"bad_term">>, Cat(json:decode(pq_service:handle(<<"[1, 2]">>)))),
    ?assertEqual(<<"bad_term">>, Cat(call(#{<<"op">> => <<"run">>}))),
    ?assertEqual(<<"bad_term">>, Cat(call(#{<<"op">> => <<"probabilities">>, <<"circuit">> => hom(),
                                          <<"options">> => #{<<"shots">> => -1}}))),
    ?assertEqual(<<"bad_term">>, Cat(call(#{<<"op">> => <<"probabilities">>, <<"circuit">> => hom(),
                                          <<"options">> => #{<<"strict">> => <<"yes">>}}))),
    %% PQCF's own errors keep their categories
    Bad = (hom())#{<<"instructions">> => [#{<<"op">> => <<"loss">>, <<"mode">> => 0, <<"eta">> => 2}]},
    ?assertEqual(<<"domain_violation">>, Cat(call(#{<<"op">> => <<"validate">>, <<"circuit">> => Bad}))),
    Sym = (hom())#{<<"instructions">> => [#{<<"op">> => <<"phaseshift">>, <<"mode">> => 0, <<"phi">> => <<"a">>},
                                          #{<<"op">> => <<"measure_pnr">>, <<"modes">> => [0]}]},
    ?assertEqual(<<"unbound_symbols">>, Cat(call(#{<<"op">> => <<"probabilities">>, <<"circuit">> => Sym}))),
    %% the request id comes back with errors too
    #{<<"id">> := <<"x">>} = call(#{<<"id">> => <<"x">>, <<"op">> => <<"nope">>}).

sizes() ->
    %% a request beyond a limit is an error response (here max_fock_dim, given as an option), not a crash
    Big = #{<<"pq">> => 1, <<"n_modes">> => 6, <<"input">> => #{<<"kind">> => <<"fock">>, <<"occupation">> => [3, 3, 3, 3, 3, 3]},
            <<"instructions">> => [#{<<"op">> => <<"kerr">>, <<"mode">> => 0, <<"kappa">> => 0.1}]},
    #{<<"error">> := #{<<"category">> := C}} = call(#{<<"op">> => <<"evolve">>, <<"circuit">> => Big,
                                                     <<"options">> => #{<<"max_fock_dim">> => 100}}),
    ?assertEqual(<<"dimension_limit">>, C).
