%% Reducing a large sum of fractions to one (decision record R1, phase 4): PQCF's reconstruction over prime fields
%% (pq_recon:rational/3) against FORM 5's PolyRatFun, which adds the fractions with polynomial GCDs (FLINT) as it
%% goes. The sum is shaped like the coefficients of a physics calculation: K terms, each a polynomial over a product
%% of a few linear factors drawn from a small set (propagator-like denominators), in V variables. Both results are
%% checked against each other exactly: N_pqcf D_form - N_form D_pqcf = 0.
%%
%%   erl -pa _build/default/lib/*/ebin -pa <dir> -noshell -eval 'fractions:main("<dir>", [{20, 3}, {200, 4}])'
-module(fractions).
-export([main/2, sum/3, form_script/2]).

-define(LIMITS, #{sym_timeout_ms => 3600000, sym_max_heap_words => 4000000000, max_poly_terms => 50000000,
                  max_expand_power => 64, max_expr_nodes => 500000000}).

vars(V) -> [list_to_atom("x" ++ integer_to_list(I)) || I <- lists:seq(1, V)].

%% the sum: K terms a_i / (f_j1 f_j2 f_j3), the factors linear forms with small integer coefficients, a_i a random
%% polynomial of degree up to 3 with up to 6 terms; seeded, so FORM and PQCF get the same sum
sum(K, V, Seed) ->
    S0 = rand:seed_s(exsss, {Seed, K, V}),
    Xs = [{sym, X} || X <- vars(V)],
    {Factors, S1} = lists:mapfoldl(fun(_, S) -> linear(Xs, S) end, S0, lists:seq(1, 2 * V + 2)),
    {Terms, _} = lists:mapfoldl(fun(_, S) ->
                                        {NF, Sa} = rand:uniform_s(3, S),
                                        {Idx, Sb} = lists:mapfoldl(fun(_, Sx) -> rand:uniform_s(length(Factors), Sx) end, Sa, lists:seq(1, NF)),
                                        {Num, Sc} = numerator(Xs, Sb),
                                        {{Num, [lists:nth(I, Factors) || I <- Idx]}, Sc}
                                end, S1, lists:seq(1, K)),
    Terms.

linear(Xs, S0) ->
    {Cs, S1} = lists:mapfoldl(fun(_, S) -> {C, S2} = rand:uniform_s(7, S), {C - 4, S2} end, S0, Xs),
    {C0, S2} = rand:uniform_s(5, S1),
    {{add, [C0 | [{mul, [C, X]} || {C, X} <- lists:zip(Cs, Xs), C =/= 0]]}, S2}.

numerator(Xs, S0) ->
    {NT, S1} = rand:uniform_s(6, S0),
    lists:mapfoldl(fun(_, S) ->
                           {C, Sa} = rand:uniform_s(19, S),
                           {Es, Sb} = lists:mapfoldl(fun(_, Sx) -> rand:uniform_s(length(Xs) + 1, Sx) end, Sa, [1, 2, 3]),
                           Mono = [lists:nth(E, Xs) || E <- Es, E =< length(Xs)],
                           {{mul, [C - 10 | Mono]}, Sb}
                   end, S1, lists:seq(1, NT)).

expr(Terms) ->
    {add, [{mul, [{add, Num}, {pow, {mul, Den}, -1}]} || {Num, Den} <- Terms]}.

main(Dir, Cases) ->
    {ok, _} = application:ensure_all_started(pq),
    [begin
         Terms = sum(K, V, 1),
         E = expr(Terms),
         Vars = vars(V),
         {Tc, {ok, Prog}} = timer:tc(pq_symbolic, compile, [[E], Vars, maps:merge(?LIMITS, #{exact_constants => true})]),
         {Tr, {ok, {N, D}}} = timer:tc(pq_recon, rational, [Prog, Vars, #{}]),
         File = filename:join(Dir, io_lib:format("frac_~b_~b", [K, V])),
         ok = file:write_file(File ++ ".frm", form_script(Terms, Vars)),
         ok = file:write_file(File ++ ".pqcf", io_lib:format("~s~n~s~n", [pq_symbolic:pretty(pq_sym_norm:to_ast(N)),
                                                                          pq_symbolic:pretty(pq_sym_norm:to_ast(D))])),
         io:format("fractions K=~b V=~b: compiled in ~.2f s (~b operations); reconstructed in ~.2f s: numerator ~b terms, "
                   "denominator ~b terms~n", [K, V, Tc / 1.0e6, length(element(4, Prog)), Tr / 1.0e6, maps:size(N), maps:size(D)])
     end || {K, V} <- Cases],
    ok.

%% FORM: the same terms as rat(numerator, denominator) under PolyRatFun, summed; the result printed
form_script(Terms, Vars) ->
    P = fun(Ast) -> binary_to_list(pq_symbolic:pretty(Ast)) end,
    Body = lists:join(" +\n  ", [["rat(", P({add, Num}), ", ", P({mul, Den}), ")"] || {Num, Den} <- Terms]),
    lists:flatten(["#-\nSymbols ", lists:join(",", [atom_to_list(X) || X <- Vars]), ";\nCFunction rat;\nPolyRatFun rat;\n",
                   "Local F =\n  ", Body, ";\n.sort\n",
                   "#write \"terms %d\", termsin_(F)\nFormat nospaces;\nPrint +s;\n.end\n"]).
