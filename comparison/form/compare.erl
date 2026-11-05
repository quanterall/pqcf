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
    %% the classic multiplication benchmarks of the computer-algebra literature, used for FORM among others:
    %% Fateman's dense f (f + 1), f = (1 + x + y + z + t)^n, and Monagan and Pearce's sparse f g,
    %% f = (1 + x + y + 2 z^2 + 3 t^3 + 5 u^5)^n, g = (1 + u + t + 2 z^2 + 3 y^3 + 5 x^5)^n
    [{io_lib:format("T5 Fateman n=~b", [N]), io_lib:format("f (f + 1), f = (1 + x + y + z + t)^~b (Fateman)", [N]),
      fun() -> {ok, F} = pq_symbolic:canonicalise(fateman(N), ?LIMITS),
               pq_symbolic:canonicalise({mul, [F, {add, [F, 1]}]}, ?LIMITS) end,
      fun() -> pq_form:expand({mul, [fateman(N), {add, [fateman(N), 1]}]}, ?LIMITS) end,
      lists:flatten(io_lib:format("* final: F\nSymbols x,y,z,t;\nLocal P = (1+x+y+z+t)^~b;\n.sort\n"
                                  "Local F = P*(P+1);\n.sort\nDrop P;\n.end\n", [N]))} || N <- [10, 15]] ++
    [{io_lib:format("T6 sparse n=~b", [N]), io_lib:format("f g, f = (1 + x + y + 2z^2 + 3t^3 + 5u^5)^~b, g its mirror (Monagan and Pearce)", [N]),
      fun() -> {ok, F} = pq_symbolic:canonicalise(mp_f(N), ?LIMITS), {ok, G} = pq_symbolic:canonicalise(mp_g(N), ?LIMITS),
               pq_symbolic:canonicalise({mul, [F, G]}, ?LIMITS) end,
      fun() -> pq_form:expand({mul, [mp_f(N), mp_g(N)]}, ?LIMITS) end,
      lists:flatten(io_lib:format("* final: F\nSymbols x,y,z,t,u;\nLocal P = (1+x+y+2*z^2+3*t^3+5*u^5)^~b;\n"
                                  "Local Q = (1+u+t+2*z^2+3*y^3+5*x^5)^~b;\n.sort\nLocal F = P*Q;\n.sort\n"
                                  "Drop P,Q;\n.end\n", [N, N]))} || N <- [6, 8]] ++
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
    write_charts(OutDir, Env, Rows),
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

fateman(N) -> {pow, {add, [1, {sym, x}, {sym, y}, {sym, z}, {sym, t}]}, N}.

mp_f(N) ->
    {pow, {add, [1, {sym, x}, {sym, y}, {mul, [2, {pow, {sym, z}, 2}]}, {mul, [3, {pow, {sym, t}, 3}]},
                 {mul, [5, {pow, {sym, u}, 5}]}]}, N}.

mp_g(N) ->
    {pow, {add, [1, {sym, u}, {sym, t}, {mul, [2, {pow, {sym, z}, 2}]}, {mul, [3, {pow, {sym, y}, 3}]},
                 {mul, [5, {pow, {sym, x}, 5}]}]}, N}.

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
           "![The latest run](charts/comparison.svg)\n\n",
           "![Before and after the algorithm change of 2 Nov 2026](charts/improvement.svg)\n\n",
           "![The ratio over the weekly runs](charts/history.svg)\n\n",
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

%%% Charts (SVG, no dependencies: the weekly run redraws them). A dot plot on a log time axis; the palette's
%%% first three categorical slots (validated all-pairs, light and dark), each series its own marker shape;
%%% text in the ink tokens, never the series colour; a tooltip (<title>) on every marker; light and dark.

-define(W, 760).
-define(LEFT, 200).
-define(RIGHT, 70).
-define(ROW, 30).
-define(TOP, 96).

style() ->
    "<style>"
    ".bg{fill:#fcfcfb}.t1{fill:#0b0b0b}.t2{fill:#52514e}.mu{fill:#898781}.grid{stroke:#e1e0d9}.axis{stroke:#c3c2b7}"
    ".link{stroke:#c3c2b7}.s1{fill:#2a78d6}.s2{fill:#eb6834}.s3{fill:#1baf7a}.ring{stroke:#fcfcfb}"
    ".l1{stroke:#2a78d6}.l2{stroke:#eb6834}.l3{stroke:#1baf7a}"
    "@media (prefers-color-scheme: dark){.bg{fill:#1a1a19}.t1{fill:#ffffff}.t2{fill:#c3c2b7}.grid{stroke:#2c2c2a}"
    ".axis{stroke:#383835}.link{stroke:#383835}.s1{fill:#3987e5}.s2{fill:#d95926}.s3{fill:#199e70}.ring{stroke:#1a1a19}"
    ".l1{stroke:#3987e5}.l2{stroke:#d95926}.l3{stroke:#199e70}}"
    "text{font-family:system-ui,-apple-system,'Segoe UI',sans-serif}"
    "</style>".

%% a marker of series K (1 circle, 2 square, 3 diamond) at (X, Y) with a tooltip
marker(K, X0, Y0, Tip) ->
    {X, Y} = {float(X0), float(Y0)},
    Shape = case K of
                1 -> io_lib:format("<circle cx='~.1f' cy='~.1f' r='5'", [X, Y]);
                2 -> io_lib:format("<rect x='~.1f' y='~.1f' width='9' height='9' rx='1.5'", [X - 4.5, Y - 4.5]);
                3 -> io_lib:format("<path d='M~.1f ~.1f l6 6 l-6 6 l-6 -6 z'", [X, Y - 6])
            end,
    [Shape, io_lib:format(" class='s~b ring' stroke-width='2'><title>~s</title></", [K, esc(Tip)]),
     case K of 1 -> "circle"; 2 -> "rect"; 3 -> "path" end, ">"].

esc(S) -> lists:flatten(string:replace(string:replace(lists:flatten(S), "&", "&amp;", all), "<", "&lt;", all)).

%% x of a time on the log axis from 10^Lo to 10^Hi seconds
xpos(T, Lo, Hi) -> xpos(T, Lo, Hi, ?RIGHT).
xpos(T, Lo, Hi, Right) -> ?LEFT + (math:log10(max(T, math:pow(10, Lo))) - Lo) / (Hi - Lo) * (?W - ?LEFT - Right).

%% Rows: [{Label, [{SeriesK, SeriesName, Seconds | none}], RightNote}]
dot_plot(Title, Sub, Legend, Rows, Lo, Hi) -> dot_plot(Title, Sub, Legend, Rows, Lo, Hi, ?RIGHT).

%% Right: the margin kept for the notes, so no marker runs under them
dot_plot(Title, Sub, Legend, Rows, Lo, Hi, Right) ->
    H = ?TOP + length(Rows) * ?ROW + 44,
    Ticks = lists:seq(Lo, Hi),
    TickLabel = fun(E) when E < 0 -> io_lib:format("~s s", [float_to_list(math:pow(10, E), [{decimals, -E}, compact])]);
                   (E) -> io_lib:format("~b s", [round(math:pow(10, E))])
                end,
    Grid = [[io_lib:format("<line class='grid' x1='~.1f' x2='~.1f' y1='~b' y2='~b' stroke-width='1'/>",
                           [xpos(math:pow(10, E), Lo, Hi, Right), xpos(math:pow(10, E), Lo, Hi, Right), ?TOP - 10, H - 38]),
             io_lib:format("<text class='mu' x='~.1f' y='~b' font-size='11' text-anchor='middle'>~s</text>",
                           [xpos(math:pow(10, E), Lo, Hi, Right), H - 22, TickLabel(E)])] || E <- Ticks],
    Leg = lists:append([[marker(K, ?LEFT + 4 + (I - 1) * 180, 66, Name),
                         io_lib:format("<text class='t2' x='~b' y='70' font-size='12'>~s</text>",
                                       [?LEFT + 16 + (I - 1) * 180, esc(Name)])]
                        || {I, {K, Name}} <- lists:enumerate(Legend)]),
    Body = [begin
                Y = ?TOP + (I - 1) * ?ROW + ?ROW div 2,
                Vals = [V || {_, _, V} <- Ms, V =/= none],
                Link = case Vals of
                           [_, _ | _] -> io_lib:format("<line class='link' x1='~.1f' x2='~.1f' y1='~b' y2='~b' stroke-width='2'/>",
                                                       [xpos(lists:min(Vals), Lo, Hi, Right), xpos(lists:max(Vals), Lo, Hi, Right), Y, Y]);
                           _ -> ""
                       end,
                [io_lib:format("<text class='t1' x='~b' y='~b' font-size='12' text-anchor='end'>~s</text>",
                               [?LEFT - 12, Y + 4, esc(Label)]),
                 Link,
                 [marker(K, xpos(V, Lo, Hi, Right), Y, io_lib:format("~s, ~s: ~s s", [Label, Name, fmt(V)])) || {K, Name, V} <- Ms, V =/= none],
                 io_lib:format("<text class='t2' x='~b' y='~b' font-size='12' text-anchor='end'>~s</text>",
                               [?W - 8, Y + 4, esc(Note)])]
            end || {I, {Label, Ms, Note}} <- lists:enumerate(Rows)],
    ["<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 ", integer_to_list(?W), " ", integer_to_list(H),
     "' width='", integer_to_list(?W), "' role='img' aria-label='", esc(Title), "'>", style(),
     io_lib:format("<rect class='bg' width='~b' height='~b' rx='8'/>", [?W, H]),
     io_lib:format("<text class='t1' x='20' y='30' font-size='16' font-weight='600'>~s</text>", [esc(Title)]),
     io_lib:format("<text class='t2' x='20' y='48' font-size='12'>~s</text>", [esc(Sub)]),
     Grid, Leg, Body,
     io_lib:format("<text class='mu' x='~b' y='~b' font-size='11' text-anchor='end'>seconds, single-threaded, log scale</text>",
                   [?W - Right, H - 6]),
     "</svg>\n"].

write_charts(Dir, Env, Rows) ->
    ok = filelib:ensure_dir(filename:join([Dir, "charts", "x"])),
    #{date := D, machine := M, form := F} = Env,
    Latest = [{N, [{1, "FORM 5", Tf}, {2, "PQCF native", Tn}, {3, "PQCF bridge to FORM", Tb}],
               case ratio(Tn, Tf) of "-" -> ""; R -> [R, "x"] end}
              || #{name := N, form := Tf, native := Tn, bridge := Tb} <- Rows],
    ok = file:write_file(filename:join([Dir, "charts", "comparison.svg"]),
                         dot_plot("pq_symbolic and FORM 5 on the reference tasks",
                                  io_lib:format("~s, ~s, ~s. Right: PQCF native time over FORM's.",
                                                [lists:sublist(D, 10), M, hd(string:split(F, " ("))]),
                                  [{1, "FORM 5"}, {2, "PQCF native"}, {3, "PQCF bridge to FORM"}], Latest, -3, 1)),
    ok = file:write_file(filename:join([Dir, "charts", "improvement.svg"]), improvement()),
    ok = file:write_file(filename:join([Dir, "charts", "history.svg"]), history(Dir)).

%% the algorithm change of 2 Nov 2026, measured on dev-aarch64-shared (docs/decisions/SY-FORM-comparison.md);
%% "refused" (a safety limit) is drawn at the axis' end with its tooltip saying so
improvement() ->
    Rows = [{"(1 + x1..x6)^16", 0.03, 61.0, 0.64}, {"(1 + x1..x6)^20", 0.12, refused, 2.4},
            {"dense product", 0.14, 3.4, 1.3}, {"permanent 7 x 7", 0.02, 429.0, 0.02},
            {"permanent 8 x 8", 0.16, refused, 0.26}],
    Plot = [{L, [{1, "FORM 5", Tf}, {3, "PQCF before", case B of refused -> none; _ -> B end}, {2, "PQCF now", A}],
             case B of refused -> "before: refused"; _ -> [times(B / A), " faster"] end}
            || {L, Tf, B, A} <- Rows],
    dot_plot("pq_symbolic before and after the algorithm change",
             "2 Nov 2026: multinomial powers, plain products, the pruned permanent. Same answers, exactly.",
             [{1, "FORM 5"}, {3, "PQCF before"}, {2, "PQCF now"}], Plot, -2, 3, 130).

%% a speed-up for reading: 2.6x, 95x, 21,450x
times(X) when X < 10 -> io_lib:format("~.1fx", [X]);
times(X) -> [thousands(integer_to_list(round(X))), "x"].

thousands(S) when length(S) =< 3 -> S;
thousands(S) -> thousands(lists:sublist(S, length(S) - 3)) ++ "," ++ lists:nthtail(length(S) - 3, S).

%% the geometric mean over the tasks of native time over FORM's, per run, one line per machine
history(Dir) ->
    {ok, Bin} = file:read_file(filename:join(Dir, "history.csv")),
    [_ | Lines] = [L || L <- string:split(binary_to_list(Bin), "\n", all), L =/= ""],
    Recs = [begin
                [Date, Mach | _] = string:split(L, ",", all),
                Fields = string:split(L, ",", all),
                Tf = lists:nth(length(Fields) - 3, Fields), Tn = lists:nth(length(Fields) - 2, Fields),
                {lists:sublist(Date, 16), Mach, num(Tf), num(Tn)}
            end || L <- Lines],
    Runs = lists:usort([{D, M} || {D, M, _, _} <- Recs]),
    Points = [{D, M, math:pow(10, lists:sum(Ls) / length(Ls))}
              || {D, M} <- Runs, Ls <- [[math:log10(Tn / Tf) || {D2, M2, Tf, Tn} <- Recs, D2 =:= D, M2 =:= M,
                                                               is_float(Tf), is_float(Tn), Tf > 0.002]], Ls =/= []],
    Machines = lists:sublist(lists:usort([M || {_, M, _} <- Points]), 3),
    Dates = lists:usort([D || {D, _, _} <- Points]),
    {W, H, L, R, T, B} = {?W, 300, 70, 30, 80, 46},
    Xof = fun(D) -> case length(Dates) of
                        1 -> L + (W - L - R) / 2;
                        N -> L + (index(D, Dates) - 1) / (N - 1) * (W - L - R)
                    end end,
    Ys = [Y || {_, _, Y} <- Points],
    Lo = math:floor(math:log10(lists:min([1.0 | Ys]))), Hi = math:ceil(math:log10(lists:max([10.0 | Ys]))),
    Yof = fun(V) -> T + (Hi - math:log10(V)) / (Hi - Lo) * (H - T - B) end,
    Grid = [[io_lib:format("<line class='grid' x1='~b' x2='~b' y1='~.1f' y2='~.1f' stroke-width='1'/>",
                           [L, W - R, Yof(math:pow(10, E)), Yof(math:pow(10, E))]),
             io_lib:format("<text class='mu' x='~b' y='~.1f' font-size='11' text-anchor='end'>~sx</text>",
                           [L - 8, Yof(math:pow(10, E)) + 4, integer_to_list(round(math:pow(10, E)))])]
            || E <- lists:seq(round(Lo), round(Hi))],
    Series = [begin
                  Ps = [{Xof(D), Yof(V), D, V} || {D, M2, V} <- Points, M2 =:= M],
                  [case Ps of
                       [_, _ | _] -> io_lib:format("<polyline class='l~b' fill='none' stroke-width='2' stroke-linejoin='round' points='~s'/>",
                                                   [K, lists:join(" ", [io_lib:format("~.1f,~.1f", [X, Y]) || {X, Y, _, _} <- Ps])]);
                       _ -> ""
                   end,
                   [marker(K, X, Y, io_lib:format("~s, ~s: native ~sx FORM's time", [M, D, io_lib:format("~.1f", [V])])) || {X, Y, D, V} <- Ps],
                   marker(K, L + 4 + (K - 1) * 230, 58, M),
                   io_lib:format("<text class='t2' x='~b' y='62' font-size='12'>~s</text>", [L + 16 + (K - 1) * 230, esc(M)])]
              end || {K, M} <- lists:enumerate(Machines)],
    XL = case Dates of
             [] -> [];
             _ -> [io_lib:format("<text class='mu' x='~.1f' y='~b' font-size='11' text-anchor='~s'>~s</text>",
                                 [Xof(D), H - 24, A, lists:sublist(D, 10)])
                   || {D, A} <- lists:usort([{hd(Dates), "start"}, {lists:last(Dates), "end"}])]
         end,
    ["<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 ", integer_to_list(W), " ", integer_to_list(H),
     "' width='", integer_to_list(W), "' role='img' aria-label='The ratio over the weekly runs'>", style(),
     io_lib:format("<rect class='bg' width='~b' height='~b' rx='8'/>", [W, H]),
     "<text class='t1' x='20' y='30' font-size='16' font-weight='600'>PQCF native time over FORM's, run by run</text>",
     "<text class='t2' x='20' y='46' font-size='12'>The geometric mean over the reference tasks (lower is closer to FORM); one line per machine.</text>",
     Grid, Series, XL, "</svg>\n"].

index(X, L) -> length(lists:takewhile(fun(Y) -> Y =/= X end, L)) + 1.

num(S) -> try list_to_float(S) catch error:badarg -> none end.

