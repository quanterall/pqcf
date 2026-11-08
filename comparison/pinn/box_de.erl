%% The one-loop massless box's differential equation in t, derived by PQCF (pq_ibp:differential_equation/4; decision
%% record R1): d M / dt = A(d, s, t) M for the masters M = (box, s-bubble, t-bubble). Written for the PINN experiment
%% (pinn_box.py): one line per entry, "i j expression" in Python syntax.
%%
%%   box_de:run("<dir>/box_de.txt")
-module(box_de).
-export([run/1]).

run(File) ->
    {ok, _} = application:ensure_all_started(pq),
    Fam = pq_ibp_family:example(box),
    Vars = pq_ibp_family:vars(Fam),
    Masters = [pq_ibp_family:canon(Fam, M) || M <- [{1, 1, 1, 1}, {1, 0, 1, 0}, {0, 1, 0, 1}]],
    {ok, Rows} = pq_ibp:differential_equation(Fam, Masters, t, #{}),
    Exact = [[element(1, pq_simplify:rational(N, D, Vars, #{})) || {N, D} <- R] || R <- Rows],
    Lines = [[integer_to_list(I - 1), " ", integer_to_list(J - 1), " ", python(Ast), "\n"]
             || {I, R} <- lists:enumerate(Exact), {J, Ast} <- lists:enumerate(R)],
    ok = file:write_file(File, ["# masters ", io_lib:format("~w", [Masters]), "\n" | Lines]),
    [io:format("dM~b/dt ~b: ~s~n", [I - 1, J - 1, pq_simplify:format(Ast)]) || {I, R} <- lists:enumerate(Exact), {J, Ast} <- lists:enumerate(R)],
    ok.

python(Ast) ->
    re:replace(lists:flatten(io_lib:format("~s", [pq_simplify:format(Ast)])), "\\^", "**", [global, {return, list}]).
