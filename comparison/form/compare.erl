%% The standing comparison of pq_symbolic with FORM (docs/decisions/SY-FORM-comparison.md), run weekly by
%% .github/workflows/form-comparison.yml and by hand with comparison/form/run.sh, until the end of the project.
%%
%% Every reference task runs three ways:
%%   native   pq_symbolic (canonicalise, permanent, circuit_unitary)
%%   bridge   pq_form, FORM as an external process, the result read back into pq_symbolic's canonical form
%%   form     FORM alone on the equivalent script (its own wall time, started from here)
%% and is checked: native and bridge exactly equal (=:=), the term counts of all three equal. A failed check fails
%% the run. Results: RESULTS.md (the latest run) and history.csv (one row per task per run), in OutDir.
%%
%%   erl -pa _build/default/lib/*/ebin -pa <dir> -noshell -eval 'compare:main("comparison/form")'
-module(compare).
-export([main/1, tasks/0]).

-define(LIMITS, #{sym_timeout_ms => 900000, sym_max_heap_words => 4000000000, max_poly_terms => 50000000,
                  max_expand_power => 64, max_expr_nodes => 500000000}).

%% {Name, What, native fun, bridge fun | none, FORM script}
tasks() ->
    S = fun(N) -> {sym, list_to_atom(N)} end,
    Xs = [S("x" ++ integer_to_list(I)) || I <- lists:seq(1, 6)],
    Pow = fun(N) -> {pow, {add, [1 | Xs]}, N} end,
    P10 = {pow, {add, [1, S("x"), S("y"), S("z"), S("w")]}, 10},
    Mat = fun(N) -> {matrix, [[S(lists:flatten(io_lib:format("a~b_~b", [I, J]))) || J <- lists:seq(1, N)]
                              || I <- lists:seq(1, N)]} end,
    [{io_lib:format("T1 power n=~b", [N]), io_lib:format("(1 + x1 + ... + x6)^~b expanded", [N]),
      fun() -> pq_symbolic:canonicalise(Pow(N), ?LIMITS) end,
      fun() -> pq_form:expand(Pow(N), ?LIMITS) end,
      form_power(N)} || N <- [12, 16, 20]] ++
    [{"T2 product", "((1 + x + y + z + w)^10)^2: a dense product",
      fun() -> {ok, P} = pq_symbolic:canonicalise(P10, ?LIMITS), pq_symbolic:canonicalise({mul, [P, P]}, ?LIMITS) end,
      fun() -> pq_form:expand({mul, [P10, P10]}, ?LIMITS) end,
      "* final: F\nSymbols x,y,z,w;\nLocal P = (1+x+y+z+w)^10;\n.sort\nLocal F = P*P;\n.sort\nDrop P;\n.end\n"}] ++
    [{io_lib:format("T3 permanent n=~b", [N]), io_lib:format("permanent of a ~b x ~b matrix of symbols", [N, N]),
      fun() -> pq_symbolic:permanent(Mat(N)) end,
      fun() -> pq_form:permanent(Mat(N), ?LIMITS) end,
      form_permanent(N)} || N <- [6, 7, 8]] ++
    [{io_lib:format("T4 mesh k=~b", [K]), io_lib:format("symbolic unitary of a ~b-mode Clements mesh, every angle a symbol", [K]),
      fun() -> pq_symbolic:circuit_unitary(mesh(K), ?LIMITS) end,
      none,                                    % trigonometric entries: outside pq_form's polynomial fragment
      form_mesh(K)} || K <- [4, 6, 8]].

main(OutDir) ->
    [application:set_env(pq, K, V) || {K, V} <- [{max_symbolic_photons, 10} | maps:to_list(?LIMITS)]],
    pq_form:available() orelse fail(form_not_available),
    Version = pq_form:version(),
    Rows = [run_task(T) || T <- tasks()],
    Env = env(Version),
    write_results(OutDir, Env, Rows),
    append_history(OutDir, Env, Rows),
    [io:format("~-20s ~10s ~10s ~10s ~10s ~s~n", [N, fmt(Tn), fmt(Tb), fmt(Tf), integer_to_list(Terms), Ok])
     || #{name := N, native := Tn, bridge := Tb, form := Tf, terms := Terms, check := Ok} <- Rows],
    case [R || #{check := C} = R <- Rows, C =/= "ok"] of
        [] -> halt(0);
        Bad -> io:format("FAILED CHECKS: ~p~n", [[maps:get(name, R) || R <- Bad]]), halt(1)
    end.

fail(Why) -> io:format("comparison: ~p~n", [Why]), halt(2).

run_task({Name, What, Native, Bridge, Script}) ->
    {Tn, Rn} = timer:tc(Native),
    {Tb, Rb} = case Bridge of none -> {none, none}; _ -> timer:tc(Bridge) end,
    {Tf, Tterms} = run_form(Script),
    NTerms = case Rn of {ok, V} -> terms(V); _ -> -1 end,
    Check = if
                element(1, Rn) =/= ok -> io_lib:format("native failed: ~p", [Rn]);
                Rb =/= none andalso Rb =/= Rn -> "native and bridge differ";
                Tterms =/= NTerms -> io_lib:format("terms: native ~b, FORM ~b", [NTerms, Tterms]);
                true -> "ok"
            end,
    #{name => lists:flatten(Name), what => lists:flatten(What), native => Tn / 1.0e6,
      bridge => case Tb of none -> none; _ -> Tb / 1.0e6 end, form => Tf, terms => NTerms, check => lists:flatten(Check)}.

terms({matrix, Rows}) -> lists:sum([terms(E) || Row <- Rows, E <- Row]);
terms({add, L}) -> length(L);
terms(0) -> 0;
terms(_) -> 1.

%% FORM alone: its wall time (best of 3) and the total terms of the expressions it ends with
run_form(Script) ->
    Dir = "/tmp/compare_form_" ++ integer_to_list(erlang:unique_integer([positive])),
    ok = file:make_dir(Dir),
    F = filename:join(Dir, "t.frm"),
    ok = file:write_file(F, Script),
    Exe = os:find_executable("form"),
    %% in its own directory: FORM writes its scratch files to the current one
    Runs = [timer:tc(fun() -> os:cmd("cd " ++ Dir ++ " && " ++ Exe ++ " t.frm") end) || _ <- [1, 2, 3]],
    {T, Out} = hd(lists:keysort(1, Runs)),
    _ = file:del_dir_r(Dir),
    {T / 1.0e6, final_terms(Script, Out)}.

%% the "Terms in output" of the final expressions' last appearance; each script names them in its first line,
%% "* final: <name or prefix>" (FORM's log wraps long lines inside names, so its Drop statements cannot be read back)
final_terms(Script, Out) ->
    {match, [Prefix]} = re:run(iolist_to_binary(Script), "^\\* final: (\\S+)", [{capture, all_but_first, list}]),
    Last = lists:foldl(fun(L, Acc) ->
                               case re:run(L, "^\\s*(\\S+)\\s+Terms in output =\\s+(\\d+)", [{capture, all_but_first, list}]) of
                                   {match, [Name, N]} -> Acc#{Name => list_to_integer(N)};
                                   nomatch -> Acc
                               end
                       end, #{}, string:split(Out, "\n", all)),
    lists:sum([N || Name := N <- Last, Name =:= Prefix orelse (lists:last(Prefix) =:= $r andalso
                                                               lists:prefix(Prefix, Name))]).

%%% FORM scripts of the tasks

form_power(N) ->
    io_lib:format("* final: F\nSymbols x1,x2,x3,x4,x5,x6;\nLocal F = (1+x1+x2+x3+x4+x5+x6)^~b;\n.end\n", [N]).

form_permanent(N) ->
    A = fun(I, J) -> io_lib:format("a~b~b", [I, J]) end,
    X = fun(J) -> io_lib:format("x~b", [J]) end,
    Js = lists:seq(1, N),
    ["* final: F\nSymbols ", lists:join(",", [A(I, J) || I <- Js, J <- Js] ++ [X(J) || J <- Js]), ";\nLocal F = 1;\n.sort\n",
     [["Multiply (", lists:join("+", [[A(I, J), "*", X(J)] || J <- Js]), ");\n",
       [["id ", X(J), "^2 = 0;\n"] || J <- Js], ".sort\n"] || I <- Js],
     "id ", lists:join("*", [X(J) || J <- Js]), " = 1;\n.end\n"].

%% the beamsplitters of mesh/1 in FORM: cos t = (z + 1/z)/2, sin t = -i/2 (z - 1/z) with z = e^(i t), p = e^(i phi)
form_mesh(K) ->
    Bs = mesh_modes(K),
    N = length(Bs),
    U = fun(S, R, C) -> io_lib:format("U~br~bc~b", [S, R, C]) end,
    Ks = lists:seq(0, K - 1),
    [io_lib:format("* final: U~br\n", [N]), "Symbols ", lists:join(",", [io_lib:format("t~b", [J]) || J <- lists:seq(0, N - 1)] ++
                                 [io_lib:format("p~b", [J]) || J <- lists:seq(0, N - 1)]), ";\n",
     [[["Local ", U(0, R, C), " = ", case R =:= C of true -> "1"; false -> "0" end, ";\n"] || C <- Ks] || R <- Ks],
     ".sort\n",
     [begin
          Cos = io_lib:format("(t~b+t~b^-1)/2", [J, J]),
          Sin = io_lib:format("(-i_/2)*(t~b-t~b^-1)", [J, J]),
          Gab = io_lib:format("(-p~b^-1)*~s", [J, Sin]),
          Gba = io_lib:format("p~b*~s", [J, Sin]),
          [[[["Local ", U(J + 1, R, C), " = ",
              if R =:= M -> [Cos, "*", U(J, M, C), " + ", Gab, "*", U(J, M + 1, C)];
                 R =:= M + 1 -> [Gba, "*", U(J, M, C), " + ", Cos, "*", U(J, M + 1, C)];
                 true -> U(J, R, C)
              end, ";\n"] || R <- Ks] || C <- Ks],
           ".sort\nDrop ", lists:join(",", [U(J, R, C) || R <- Ks, C <- Ks]), ";\n.sort\n"]
      end || {J, M} <- lists:enumerate(0, Bs)],
     ".end\n"].

mesh_modes(K) -> lists:sublist(lists:append([lists:seq(L rem 2, K - 2, 2) || L <- lists:seq(0, K - 1)]), K * (K - 1) div 2).

mesh(K) ->
    Is = [{beamsplitter, M, M + 1, {sym, list_to_atom("t" ++ integer_to_list(J))}, {sym, list_to_atom("p" ++ integer_to_list(J))}}
          || {J, M} <- lists:enumerate(0, mesh_modes(K))],
    {ok, C} = pq_circuit:new(K),
    setelement(5, C, Is).

%%% Results

env(Version) ->
    Git = string:trim(os:cmd("git -C /src rev-parse --short HEAD 2>/dev/null")),
    %% the host's CPU model, from the host (a container may have no lscpu, and aarch64's /proc/cpuinfo no name)
    Cpu = case os:getenv("COMPARISON_CPU", "") of
              "" -> string:trim(os:cmd("lscpu 2>/dev/null | sed -n 's/^Model name: *//p' | head -1"));
              C -> C
          end,
    #{date => calendar:system_time_to_rfc3339(erlang:system_time(second), [{offset, "Z"}]),
      form => binary_to_list(hd(binary:split(Version, <<"  ">>))),
      pqcf => case Git of "" -> os:getenv("PQCF_COMMIT", "unknown"); _ -> Git end,
      cpu => case Cpu of "" -> "unknown"; _ -> Cpu end,
      arch => erlang:system_info(system_architecture), machine => os:getenv("COMPARISON_MACHINE", "local")}.

fmt(none) -> "-";
fmt(T) when is_float(T) -> io_lib:format("~.3f", [T]).

ratio(none, _) -> "-";
ratio(_, F) when F < 0.0005 -> "-";
ratio(T, F) -> io_lib:format("~.1f", [T / F]).

write_results(Dir, Env, Rows) ->
    #{date := D, form := F, pqcf := P, cpu := Cpu, arch := A, machine := M} = Env,
    Table = [["| ", N, " | ", W, " | ", integer_to_list(Terms), " | ", fmt(Tf), " | ", fmt(Tn), " | ", ratio(Tn, Tf),
              " | ", fmt(Tb), " | ", C, " |\n"]
             || #{name := N, what := W, terms := Terms, form := Tf, native := Tn, bridge := Tb, check := C} <- Rows],
    Doc = ["# pq_symbolic and FORM: the standing comparison\n\n",
           "Generated by `comparison/form/compare.erl`. Do not edit by hand: the weekly workflow\n",
           "(`.github/workflows/form-comparison.yml`) and `comparison/form/run.sh` rewrite this file and append to\n",
           "`history.csv`. The decision record is `docs/decisions/SY-FORM-comparison.md`.\n\n",
           "Latest run: ", D, " on ", M, " (", Cpu, ", ", A, "), PQCF ", P, ", ", F, ".\n\n",
           "All times are seconds, single-threaded. The columns are:\n\n",
           "- **FORM:** FORM alone on the equivalent script, best of three runs.\n",
           "- **native:** `pq_symbolic`.\n",
           "- **ratio:** native divided by FORM.\n",
           "- **bridge:** `pq_form`, FORM run as an external process with the result read back into PQCF's\n",
           "  canonical form. A dash means the task is outside the bridge's polynomial fragment.\n",
           "- **check:** `ok` means native and bridge are exactly equal and all term counts agree.\n\n",
           "| Task | What | Terms | FORM | native | ratio | bridge | check |\n",
           "| --- | --- | --- | --- | --- | --- | --- | --- |\n", Table],
    ok = file:write_file(filename:join(Dir, "RESULTS.md"), Doc).

append_history(Dir, Env, Rows) ->
    F = filename:join(Dir, "history.csv"),
    Header = case filelib:is_regular(F) of
                 true -> [];
                 false -> "date,machine,arch,pqcf,form,task,terms,form_s,native_s,bridge_s,check\n"
             end,
    #{date := D, form := Fv, pqcf := P, arch := A, machine := M} = Env,
    Lines = [[D, ",", M, ",", A, ",", P, ",\"", Fv, "\",\"", N, "\",", integer_to_list(T), ",", fmt(Tf), ",", fmt(Tn), ",",
              fmt(Tb), ",", case C of "ok" -> "ok"; _ -> "FAILED" end, "\n"]
             || #{name := N, terms := T, form := Tf, native := Tn, bridge := Tb, check := C} <- Rows],
    ok = file:write_file(F, [Header, Lines], [append]).
