%% Published PQCF results certified by Lean 4: the polynomial gcds (task T7, nontrivial) and factorisations (task T8,
%% nontrivial) of the FORM comparison (comparison/form, the problems of compare:gcd_problems/1 and
%% factor_problems/1: T. Ueda's polybench defaults, 5 variables, 22-30 terms of degree 22-30 per basic polynomial).
%% PQCF computes each result, writes it as a Lean statement and proof (pq_lean), and a lean-check node checks it.
-module(certify_published).
-export([main/1, main/2]).

main(Node) -> main(Node, all).

main(Node, Which) ->
    ok = io:setopts([{encoding, unicode}]),
    Opts = #{checker => {node, Node}, timeout => 1800000},
    Gcds = case Which of {gcd, K} -> [lists:nth(K, compare:gcd_problems("nontrivial"))];
                         factor -> [];
                         _ -> compare:gcd_problems("nontrivial") end,
    Facs = case Which of {gcd, _} -> []; _ -> compare:factor_problems("nontrivial") end,
    R1 = [report("T7 gcd", K, A, pq_lean:gcd(A, B, Opts)) || {K, {A, B}} <- lists:enumerate(Gcds)],
    R2 = [report("T8 factor", K, P, pq_lean:factor(P, Opts)) || {K, P} <- lists:enumerate(Facs)],
    All = R1 ++ R2,
    io:format("~nsummary: ~b of ~b valid~n", [length([v || valid <- All]), length(All)]).

report(Task, K, Input, {ok, #{source := Src, verdict := V}}) ->
    #{status := St, elapsed_ms := Ms} = V,
    Ax = lists:usort(lists:append([A || #{axioms := A} <- maps:get(theorems, V), is_list(A)])),
    io:format("~s ~b: ~b input terms, Lean source ~b bytes -> ~p in ~.1f s, axioms ~p~n",
              [Task, K, terms(Input), byte_size(Src), St, Ms / 1000, Ax]),
    case St of valid -> ok; _ -> io:format("  reason ~p~n  errors ~p~n", [maps:get(reason, V), lists:sublist(maps:get(errors, V), 2)]) end,
    St.

terms({add, Ts}) -> length(Ts);
terms(_) -> 1.
