%% The last step of comparison/beta/run.sh: from the work directory Dir (FORM's logs res_<field>_<loops>.log, PQCF's
%% reductions <field>_<loops>.red, Forcer's values of PQCF's masters masters_<field>_<loops>.out, and Forcer's
%% Masters.prc), the files of results/ written into Out:
%%   renormalisation_form_forcer_3loops.txt  three loops, from FORM + Forcer's bare self-energies;
%%   renormalisation_pqcf_2loops.txt         two loops, every integral reduced by PQCF;
%%   form_forcer/bare_<field>_<loops>.txt    FORM + Forcer's bare self-energies (FORM's output);
%%   pqcf/bare_pqcf_<field>_<loops>.txt      PQCF's bare self-energies;
%%   pqcf/reduction_<field>_<loops>.txt      PQCF's reductions; pqcf/masters_<field>_<loops>.frm its masters for Forcer.
-module(beta_run).
-export([write_results/2]).

write_results(Dir, Out) ->
    {ok, _} = application:ensure_all_started(pq),
    ok = filelib:ensure_path(filename:join(Out, "form_forcer")),
    ok = filelib:ensure_path(filename:join(Out, "pqcf")),
    R3 = renorm:run(Dir, 3),
    ok = renorm:write(filename:join(Out, "renormalisation_form_forcer_3loops.txt"), R3,
                      "# Three-loop QCD, Landau gauge: Z_3 and Z~_3 from the FORM + Forcer bare self-energies (results/form_forcer),\n"
                      "# solved and the beta function extracted by PQCF (renorm.erl). Finite: " ++ atom_to_list(maps:get(finite, R3)) ++ "\n"),
    Tables = forcer_cmp:master_tables(filename:join(Dir, "Masters.prc")),
    Dm1 = epsexp:sinv(epsexp:exact(#{0 => #{[] => {3, 1}}, 1 => #{[] => {-2, 1}}}), 10),
    G = maps:from_list([{L, epsexp:smul(beta_pqcf:compose(Dir, "glu", L, Tables), Dm1)} || L <- [1, 2]]),
    H = maps:from_list([{L, beta_pqcf:compose(Dir, "gho", L, Tables)} || L <- [1, 2]]),
    R2 = renorm:beta(#{gluon => G, ghost => H, loops => 2}),
    ok = renorm:write(filename:join(Out, "renormalisation_pqcf_2loops.txt"), R2,
                      "# Two-loop QCD, Landau gauge: every integral reduced by PQCF (results/pqcf), the master integrals' values from\n"
                      "# Forcer, the renormalisation by PQCF. Finite: " ++ atom_to_list(maps:get(finite, R2)) ++ "\n"),
    [ok = file:write_file(filename:join(Out, io_lib:format("pqcf/bare_pqcf_~s_~b.txt", [F, L])),
                          ["# PQCF's bare self-energy, MS-bar (the gluon's divided by D - 1)\n", epsexp:to_string(maps:get(L, M))])
     || {F, M} <- [{"glu", G}, {"gho", H}], L <- [1, 2]],
    [begin
         {ok, Red} = file:read_file(filename:join(Dir, io_lib:format("~s_~b.red", [F, L]))),
         ok = file:write_file(filename:join(Out, io_lib:format("pqcf/reduction_~s_~b.txt", [F, L])), Red),
         {ok, Ms} = file:read_file(filename:join(Dir, io_lib:format("masters_~s_~b.frm", [F, L]))),
         ok = file:write_file(filename:join(Out, io_lib:format("pqcf/masters_~s_~b.frm", [F, L])), Ms)
     end || F <- ["glu", "gho"], L <- [1, 2]],
    % FORM's own output: the expression F, from its line "   F =" to the closing ";"
    [begin
         {ok, Log} = file:read_file(filename:join(Dir, io_lib:format("res_~s_~b.log", [F, L]))),
         {match, [Expr]} = re:run(Log, "(\\n   F =\\n.*?;)", [dotall, {capture, all_but_first, binary}]),
         <<"\n", Body/binary>> = Expr,
         ok = file:write_file(filename:join(Out, io_lib:format("form_forcer/bare_~s_~b.txt", [F, L])), [Body, "\n"])
     end || F <- ["glu", "gho"], L <- [1, 2, 3]],
    io:format("beta functions, three loops (FORM + Forcer): ~p~n", [maps:get(betas, R3)]),
    io:format("beta functions, two loops (PQCF): ~p~n", [maps:get(betas, R2)]),
    ok.
