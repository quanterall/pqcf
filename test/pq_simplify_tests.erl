%% Verified rewriting of exact results (pq_simplify; decision record R1, option 3): factors from an alphabet found by
%% exact division, coefficients collected by monomials, the result equal to the input and never larger.
-module(pq_simplify_tests).
-include_lib("eunit/include/eunit.hrl").

d() -> {sym, d}.
s() -> {sym, s}.
t() -> {sym, t}.

%% an AST expanded into pq_recon's map form over Vars
recon(Ast, Vars) ->
    {ok, C} = pq_symbolic:canonicalise(Ast),
    {ok, Prog} = pq_symbolic:compile([C], Vars, #{exact_constants => true}),
    {ok, P} = pq_recon:polynomial(Prog, Vars, #{}),
    P.

%% a reduction-like coefficient: (d - 4)^3 (2d - 7) (s + t)^2 / (3 t^2 (d - 6) (2d - 9)), written expanded
factored_test() ->
    Vars = [d, s, t],
    Num = {mul, [{pow, {add, [d(), -4]}, 3}, {add, [{mul, [2, d()]}, -7]}, {pow, {add, [s(), t()]}, 2}]},
    Den = {mul, [3, {pow, t(), 2}, {add, [d(), -6]}, {add, [{mul, [2, d()]}, -9]}]},
    N = recon(Num, Vars), D = recon(Den, Vars),
    {Ast, #{size_in := In, size_out := Out, form := factored}} = pq_simplify:rational(N, D, Vars, #{collect => []}),
    ?assert(pq_simplify:verify(Ast, N, D, Vars)),
    %% 15 terms over 3: 14 + 2 additions
    ?assertEqual(16, In),
    %% one each in d - 4, 2d - 7, s + t, d - 6, 2d - 9
    ?assertEqual(5, Out),
    ?assertEqual(Out, pq_simplify:adds(Ast)).

%% the expanded form when it has fewer additions; the factored one among equals
expanded_test() ->
    Vars = [d, s],
    %% d^3 + 1 = (d + 1)(d^2 - d + 1): one addition expanded, three factored
    N = recon({add, [{pow, d(), 3}, 1]}, Vars), D = recon({add, [d(), -6]}, Vars),
    {A1, #{size_in := 2, size_out := 2, form := expanded}} = pq_simplify:rational(N, D, Vars, #{collect => []}),
    ?assert(pq_simplify:verify(A1, N, D, Vars)),
    N2 = recon({mul, [{add, [d(), -4]}, s()]}, Vars),
    {A2, #{size_in := 2, size_out := 2, form := factored}} = pq_simplify:rational(N2, D, Vars, #{collect => []}),
    ?assert(pq_simplify:verify(A2, N2, D, Vars)).

%% collected by monomials in s, t, each coefficient factored in d and cancelled against D:
%% [(36d^2 - 240d + 400) s - (8d^2 - 56d + 96) t] / ((4 - d)(10 - 3d)) = 4 (3d - 10) s / (d - 4) + 8 (d - 3) t / (10 - 3d)
collected_test() ->
    Vars = [d, s, t],
    N = recon({add, [{mul, [{add, [{mul, [36, {pow, d(), 2}]}, {mul, [-240, d()]}, 400]}, s()]},
                     {mul, [-1, {add, [{mul, [8, {pow, d(), 2}]}, {mul, [-56, d()]}, 96]}, t()]}]}, Vars),
    D = recon({mul, [{add, [4, {mul, [-1, d()]}]}, {add, [10, {mul, [-3, d()]}]}]}, Vars),
    {Ast, #{size_in := 7, size_out := Out, form := collected}} = pq_simplify:rational(N, D, Vars, #{}),
    ?assert(pq_simplify:verify(Ast, N, D, Vars)),
    %% the sum of two terms; 3d - 10 over d - 4; d - 3 over 3d - 10
    ?assertEqual(5, Out),
    ?assert(is_list(lists:flatten(io_lib:format("~s", [pq_simplify:format(Ast)])))).

%% a factor outside the alphabet stays in the rest, expanded; a proposal (extra) finds it
proposal_test() ->
    Vars = [d, s, t],
    Odd = {add, [{mul, [d(), s()]}, {mul, [-1, {pow, t(), 2}]}, 5]},
    N = recon({mul, [{pow, Odd, 2}, s()]}, Vars), D = recon({add, [d(), -4]}, Vars),
    {A1, #{size_out := O1}} = pq_simplify:rational(N, D, Vars, #{}),
    {A2, #{size_out := O2}} = pq_simplify:rational(N, D, Vars, #{extra => [Odd]}),
    ?assert(pq_simplify:verify(A1, N, D, Vars)),
    ?assert(pq_simplify:verify(A2, N, D, Vars)),
    ?assert(O2 < O1),
    %% Odd (two additions), d - 4 (one)
    ?assertEqual(3, O2).

%% a wrong rewrite is caught; 0 stays 0
verify_catches_test() ->
    Vars = [d, s],
    N = recon({add, [d(), s()]}, Vars), D = recon(1, Vars),
    ?assertNot(pq_simplify:verify({add, [d(), s(), 1]}, N, D, Vars)),
    ?assertEqual({0, #{size_in => 0, size_out => 0, form => expanded}}, pq_simplify:rational(#{}, D, Vars, #{})).

%% partial fractions (Leinartas): exact for denominators with dependent affine factors, whichever variables they are
%% taken in; 1/(s t (s + t)) = 1/(t (s + t)^2) + 1/(s (s + t)^2), two factors per term
partial_fractions_test_() ->
    {timeout, 120, fun() ->
        Vars = [d, s, t],
        Cases = [{1, {mul, [s(), t(), {add, [s(), t()]}]}},
                 {{add, [d(), -3]}, {mul, [{pow, {add, [d(), -4]}, 2}, {add, [d(), -6]}, s()]}},
                 {{add, [{mul, [d(), s()]}, {pow, t(), 2}]}, {mul, [{pow, s(), 2}, {add, [s(), t()]}, {add, [d(), -4]}, {add, [{mul, [2, d()]}, -7]}]}},
                 {{mul, [s(), t()]}, {mul, [{add, [s(), {mul, [-1, t()]}]}, {add, [s(), t(), -1]}, {add, [s(), 2]}]}}],
        [begin
             N = recon(Nu, Vars), D = recon(De, Vars),
             [begin
                  {Ast, #{form := Form}} = pq_simplify:rational(N, D, Vars, #{apart => [X], alphabet => [], collect => []}),
                  ?assert(lists:member(Form, [factored, partial, expanded])),
                  ?assert(pq_simplify:verify(Ast, N, D, Vars))
              end || X <- [[d, s, t], [d], [s, t]]]
         end || {Nu, De} <- Cases],
        N1 = recon(1, Vars), D1 = recon({mul, [s(), t(), {add, [s(), t()]}]}, Vars),
        {A1, #{form := partial}} = pq_simplify:rational(N1, D1, Vars, #{apart => [[s, t]], form => partial}),
        ?assertEqual("(t^(-1)*(s + t)^(-2) + s^(-1)*(s + t)^(-2))", lists:flatten(io_lib:format("~s", [pq_simplify:format(A1)])))
    end}.

%% factors outside the alphabet, found by exact gcds: a repeated one, one shared by two coefficients
discover_test_() ->
    {timeout, 120, fun() ->
        Vars = [d, s, t],
        Odd = {add, [{mul, [d(), s()]}, {mul, [-1, {pow, t(), 2}]}, 5]},
        Other = {add, [{pow, s(), 2}, {mul, [3, d(), t()]}, 1]},
        Third = {add, [{mul, [d(), {pow, t(), 2}]}, {mul, [-1, s()]}, 2]},
        One = recon(1, Vars),
        Has = fun(Found) -> lists:member(canon(Odd), [canon(F) || F <- Found]) end,
        %% the repeated factor
        N1 = recon({mul, [{pow, Odd, 2}, Other]}, Vars), D1 = recon({add, [d(), -4]}, Vars),
        Found1 = pq_simplify:discover([{N1, D1}], Vars),
        ?assert(Has(Found1)),
        %% shared by two coefficients' numerators (each alone irreducible to the alphabet)
        ?assert(Has(pq_simplify:discover([{recon({mul, [Odd, Other]}, Vars), One}, {recon({mul, [Odd, Third]}, Vars), One}], Vars))),
        ?assertEqual([], pq_simplify:discover([{recon(Other, Vars), One}, {recon(Third, Vars), One}], Vars)),
        %% the proposals make the result smaller, verified
        {A, #{size_out := Out}} = pq_simplify:rational(N1, D1, Vars, #{extra => Found1}),
        ?assert(pq_simplify:verify(A, N1, D1, Vars)),
        {_, #{size_out := Out0}} = pq_simplify:rational(N1, D1, Vars, #{}),
        ?assert(Out < Out0)
    end}.

canon(A) -> {ok, C} = pq_symbolic:canonicalise(A), C.
