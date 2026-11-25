%% PQCF results certified by Lean 4 through a lean-check node (examples/lean/run.sh starts both).
-module(pq_lean_demo).
-export([main/1]).

main(Node) ->
    ok = io:setopts([{encoding, unicode}]),
    Opts = #{checker => {node, Node}},
    P = fun(S) -> {ok, E} = pq_symbolic:parse(S), E end,
    Show = fun(Title, {ok, #{source := Src, verdict := V}}) ->
                   io:format("== ~s~n~ts~n-> ~p (~p ms), axioms ~p~n~n",
                             [Title, Src, maps:get(status, V), maps:get(elapsed_ms, V),
                              [{N, A} || #{name := N, axioms := A} <- maps:get(theorems, V)]])
           end,
    % 1. a factorisation: x^4 - y^4 = (x - y)(x + y)(x^2 + y^2)
    Show("factorisation", pq_lean:factor(P(<<"x^4 - y^4">>), Opts)),
    % 2. a gcd with its cofactors
    Show("gcd", pq_lean:gcd(P(<<"x^3 - x">>), P(<<"x^2 + 2*x + 1">>), Opts)),
    % 3. an exact solution of a linear system: 3 x + 2 y = 1, x - 4 y = 5/2 gives x = 9/14, y = -13/28
    Show("linear system", pq_lean:linear([{[3, 2], 1}, {[1, -4], {5, 2}}], [{9, 14}, {-13, 28}], Opts)),
    % 4. a wrong solution (one sign changed) is rejected
    {ok, #{verdict := Bad}} = pq_lean:linear([{[3, 2], 1}, {[1, -4], {5, 2}}], [{9, 14}, {13, 28}], Opts),
    io:format("== a wrong solution~n-> ~p, ~p~n  ~ts~n", [maps:get(status, Bad), maps:get(reason, Bad),
                                                      case maps:get(errors, Bad) of [#{message := M} | _] -> M; _ -> <<>> end]).
