%% The one-loop massless box's differential equation in t, derived by PQCF (decision record R1): d M / dt = A(d, s, t) M
%% for the masters M = (box, s-bubble, t-bubble). Each master's derivative (ibp_family:derivative/3: external
%% momentum operators combined into d/dt) is a sum of integrals, each reduced to the masters (automatic seeds,
%% pq_recon:linear/4); every entry of A is then reconstructed as one reduced rational function (pq_recon:rational/3).
%% Written for the PINN experiment (pinn_box.py): one line per entry, "i j expression" in Python syntax.
%%
%%   box_de:run("<dir>/box_de.txt")
-module(box_de).
-export([run/1]).

-define(LIMITS, #{sym_timeout_ms => 3600000, sym_max_heap_words => 4000000000, max_expr_nodes => 500000000}).

run(File) ->
    {ok, _} = application:ensure_all_started(pq),
    Fam = ibp_family:family(box),
    Vars = [d, s, t],
    Masters = [ibp_family:canon(Fam, M) || M <- [{1, 1, 1, 1}, {1, 0, 1, 0}, {0, 1, 0, 1}]],
    Row = fun(M) ->
                  Terms = ibp_family:derivative(Fam, M, t),
                  Contribs = lists:append([case lists:member(I, Masters) of
                                               true -> [{I, C}];
                                               false -> [{Mj, {mul, [C, Cj]}} || {Mj, Cj} <- reduce(I, Masters)]
                                           end || {I, C} <- Terms]),
                  [{add, [0 | [C || {Mj, C} <- Contribs, Mj =:= M2]]} || M2 <- Masters]
          end,
    A = [Row(M) || M <- Masters],
    Exact = [[entry(E, Vars) || E <- R] || R <- A],
    Lines = [[integer_to_list(I - 1), " ", integer_to_list(J - 1), " ", python(Ast), "\n"]
             || {I, R} <- lists:enumerate(Exact), {J, Ast} <- lists:enumerate(R)],
    ok = file:write_file(File, ["# masters ", io_lib:format("~w", [Masters]), "\n" | Lines]),
    [io:format("dM~b/dt ~b: ~s~n", [I - 1, J - 1, pq_simplify:format(Ast)]) || {I, R} <- lists:enumerate(Exact), {J, Ast} <- lists:enumerate(R)],
    ok.

%% an integral in the masters: [{Master, CoefficientAst}]
reduce(I, Masters) ->
    Seeds = ibp:auto(box, I),
    #{a := A, bs := Bs, masters := Ms, target := TI} = ibp:system(Seeds, I, 1),
    {ok, [{TI, Per}]} = pq_recon:linear(A, Bs, [d, s, t], ?LIMITS#{components => [TI], rhs => length(Ms)}),
    Red = [{M, {mul, [pq_sym_norm:to_ast(N), {pow, pq_sym_norm:to_ast(D), -1}]}}
           || {C, {M, _}} <- lists:enumerate(0, Ms), {C1, {N, D}} <- Per, C1 =:= C, map_size(N) > 0],
    [] =:= [M || {M, _} <- Red, not lists:member(M, Masters)] orelse error({unexpected_masters, I, Red}),
    Red.

%% an entry as one reduced rational function, simplified
entry(Ast, Vars) ->
    {ok, Prog} = pq_symbolic:compile([Ast], Vars, #{exact_constants => true}),
    {ok, {N, D}} = pq_recon:rational(Prog, Vars, #{}),
    {S, _} = pq_simplify:rational(N, D, Vars, #{}),
    S.

python(Ast) ->
    re:replace(lists:flatten(io_lib:format("~s", [pq_simplify:format(Ast)])), "\\^", "**", [global, {return, list}]).
