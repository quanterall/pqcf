%% pq_form, the optional bridge to FORM (docs/decisions/SY-FORM-comparison.md). Without FORM (the dev image) only
%% the refusal is tested; with it (the comparison image, comparison/form/Dockerfile) every result must equal
%% pq_symbolic's exactly.
-module(pq_form_tests).
-include_lib("eunit/include/eunit.hrl").

s(N) -> {sym, N}.

not_available_test() ->
    ?assertEqual({error, {symbolic_unsupported, #{form => not_available}}},
                 pq_form:expand({add, [s(x), 1]}, #{form_executable => "/nonexistent/form"})),
    ?assertNot(pq_form:available(#{form_executable => "/nonexistent/form"})).

with_form_test_() ->
    case pq_form:available() of
        false -> [];               % reported by not_available_test; the comparison image runs the rest
        true ->
            Exprs = [{pow, {add, [1, s(x), {mul, [{rat, 1, 2}, s(y)]}, {complex, 0, 1}]}, 5},
                     {mul, [{pow, s(x), -2}, {add, [s(x), 1]}, pq_symbolic:sub(s(y), {rat, 3, 7})]},
                     {pow, {add, [s(a_b), s(c), -1]}, 7},                % names with '_' (FORM reserves them)
                     pq_symbolic:divide({add, [s(x), s(y)]}, 4), 0, {add, [s(x), pq_symbolic:neg(s(x))]}],
            [?_assertEqual(pq_symbolic:canonicalise(E), pq_form:expand(E)) || E <- Exprs] ++
            [?_assertEqual(pq_symbolic:permanent(M), pq_form:permanent(M))
             || N <- [1, 2, 3, 5],
                M <- [{matrix, [[s(list_to_atom(lists:flatten(io_lib:format("m~b_~b", [I, J])))) || J <- lists:seq(1, N)]
                                || I <- lists:seq(1, N)]}]] ++
            [?_assertEqual({error, {symbolic_unsupported, #{form => exp}}}, pq_form:expand({exp, s(x)})),
             ?_assertEqual({error, {symbolic_unsupported, #{form => pow}}}, pq_form:expand({pow, s(x), {rat, 1, 2}})),
             ?_assertEqual({error, {bad_shape, permanent}}, pq_form:permanent({matrix, [[s(x), s(y)]]})),
             ?_assertMatch({error, {timeout, #{form := 1}}},
                           pq_form:expand({pow, {add, [1 | [s(list_to_atom("v" ++ integer_to_list(K))) || K <- lists:seq(1, 8)]]}, 30},
                                          #{form_timeout_ms => 1}))]
    end.
