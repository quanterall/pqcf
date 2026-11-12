%% The cases of the Forcer comparison (comparison/forcer/README.md): FORM's own checks of Forcer
%% (form-dev/form, check/extra/forcer.frm: Forcer_example, Forcer_1, Forcer_2), reduced by PQCF and checked against
%% Forcer two ways. run.sh drives it:
%%   pqcf(Dir): every case reduced by PQCF (timed); its masters and its target written for Forcer
%%   check(Dir): after Forcer has run on those files -- PQCF's coefficients times Forcer's values of PQCF's masters
%%     against Forcer's value of the target (G-scheme, as epsilon expansions), and the MS-bar sums against the
%%     numbers FORM's checks assert. Dir/check.json; halts 1 if anything differs.
-module(forcer_suite).
-export([cases/0, pqcf/1, check/1, topo/1]).

-define(KITE, "vx(-Q,p2,p3)*vx(p1,-p2,p5)*vx(-p1,p4,Q)*vx(-p3,-p4,-p5)").

%% {Case, Vertices, Lines, Insertions, [{Coef, Powers, Numerators}]}
cases() ->
    [{"example_3loop", "vx(Q,p1,p5,p6)*vx(-p1,p2,p3)*vx(-p5,-p6,p4)*vx(-Q,-p2,-p3,-p4)", 6, [],
      [{1, ones(6), [{q, p3}, {q, p4}]}]},
     {"example_kite", ?KITE, 5, [1, 4], [{1, ones(5), []}]}]
    ++ [{zno(I, W), Vx, 11, [], terms(Sps, W)} || I <- [1, 2, 3], {Vx, Sps} <- [topo(I)], W <- [f1, f2]]
    ++ [{zno(I, f2n), Vx, 11, [], terms(Sps, f2n)} || I <- [1, 2], {Vx, Sps} <- [topo(I)]].

%% FORM's assertions: {Section, Expression, [{Case, Order}]} (order: loops plus insertions, for the MS-bar factor)
assertions() ->
    [{"Forcer_example", "F", [{"example_3loop", 3}, {"example_kite", 4}]},
     {"Forcer_1", "F1", [{zno(I, f1), 4} || I <- [1, 2, 3]]},
     {"Forcer_1", "F2", [{zno(I, f2), 4} || I <- [1, 2, 3]]},
     {"Forcer_2", "F1", [{zno(I, f1), 4} || I <- [1, 2]]},
     {"Forcer_2", "F2", [{zno(I, f2n), 4} || I <- [1, 2]]}].

zno(I, W) -> lists:flatten(io_lib:format("zno~b_~s", [I, W])).
ones(N) -> lists:duplicate(N, 1).

%% the three four-loop non-planar topologies of forcer.frm (Zno1..3) and their numerators n12..n14
topo(1) -> {"vx(-Q,p4,p5)*vx(p3,-p4,p10)*vx(p2,-p3,p9)*vx(p1,-p2,p11)*vx(-p5,p6,-p11)*vx(-p6,p7,-p10)*vx(-p7,p8,-p9)*vx(-p1,-p8,Q)",
            [{p2, p4}, {q, p2}, {q, p3}]};
topo(2) -> {"vx(-Q,p4,p5)*vx(p3,-p4,p11)*vx(p6,p7,p10)*vx(p2,-p3,-p10)*vx(p1,-p2,p9)*vx(-p5,-p6,-p9)*vx(-p7,p8,-p11)*vx(-p1,-p8,Q)",
            [{q, p2}, {p1, p4}, {q, p3}]};
topo(3) -> {"vx(-Q,p3,p4)*vx(p6,p8,p10)*vx(p5,-p10,p11)*vx(p1,-p3,-p5)*vx(-p4,-p8,p9)*vx(p7,-p9,-p11)*vx(p2,-p6,-p7)*vx(-p1,-p2,Q)",
            [{q, p6}, {q, p8}, {p3, p6}]}.

%% F1: every line once; F2: one line squared, or one numerator; f2n: the numerators only (Forcer_2's F2)
terms(_, f1) -> [{1, ones(11), []}];
terms(Sps, f2) -> [{1, [case J of K -> 2; _ -> 1 end || J <- lists:seq(1, 11)], []} || K <- lists:seq(1, 11)] ++ terms(Sps, f2n);
terms(Sps, f2n) -> [{1, ones(11), [Sp]} || Sp <- Sps].

pqcf(Dir) ->
    {ok, _} = application:ensure_all_started(pq),
    Rows = [begin
                P = forcer_cmp:problem(Vx, N, #{ex => Ex}),
                {T, Coeffs, Info} = forcer_cmp:reduce(P, Terms),
                Ms = maps:keys(Coeffs),
                ok = file:write_file(filename:join(Dir, Case ++ "_masters.frm"), forcer_cmp:masters_form(P, Ms)),
                ok = file:write_file(filename:join(Dir, Case ++ "_target.frm"), forcer_cmp:targets_form([{"T", target_form(Vx, N, Ex, Terms)}])),
                ok = file:write_file(filename:join(Dir, Case ++ ".coeffs"), io_lib:format("~p.~n", [{Ms, [maps:get(M, Coeffs) || M <- Ms]}])),
                %% the reduction itself, readable: the integral, then each master (Forcer's notation) and its exact
                %% coefficient in d
                ok = filelib:ensure_dir(filename:join([Dir, "symbolic", "x"])),
                ok = file:write_file(filename:join([Dir, "symbolic", Case ++ ".reduction.txt"]),
                                     [io_lib:format("# PQCF's reduction of ~s (d = 4 - 2 ep), ~b masters~n# integral:~n  ~s~n~n",
                                                    [Case, length(Ms), target_form(Vx, N, Ex, Terms)]),
                                      [io_lib:format("master M~b = ~s~ncoefficient = ~s~n~n",
                                                     [K, forcer_cmp:integral_text(P, M), pq_symbolic:pretty(maps:get(M, Coeffs))])
                                       || {K, M} <- lists:enumerate(Ms)]]),
                io:format("~-16s ~8.2f s  ~7b equations  ~6b needed  ~3b masters  seeds ~w~n",
                          [Case, T, maps:get(equations, Info, 0), maps:get(needed, Info, 0), length(Ms), maps:get(seeds, Info, none)]),
                #{<<"case">> => list_to_binary(Case), <<"pqcf_seconds">> => T, <<"equations">> => maps:get(equations, Info, 0),
                  <<"needed">> => maps:get(needed, Info, 0), <<"masters">> => length(Ms),
                  <<"seeds">> => list_to_binary(io_lib:format("~w", [maps:get(seeds, Info, none)]))}
            end || {Case, Vx, N, Ex, Terms} <- cases()],
    ok = file:write_file(filename:join(Dir, "pqcf.json"), json:encode(Rows)).

%% the same integrals in Forcer's input notation
target_form(Vx, N, Ex, Terms) ->
    Name = fun(q) -> "Q"; (A) -> atom_to_list(A) end,
    ExS = case Ex of [] -> ""; _ -> ["*ex(", lists:join(",", ["p" ++ integer_to_list(J) || J <- Ex]), ")"] end,
    lists:join("\n + ", [[integer_to_list(C), "*", Vx,
                          [io_lib:format("/p~b.p~b^~b", [J, J, X]) || {J, X} <- lists:enumerate(Ps), J =< N, X =/= 0],
                          [["*", Name(A), ".", Name(B)] || {A, B} <- Sps], ExS]
                         || {C, Ps, Sps} <- Terms]).

check(Dir) ->
    {ok, _} = application:ensure_all_started(pq),
    Tables = forcer_cmp:master_tables(filename:join(Dir, "Masters.prc")),
    Series = maps:from_list(
               [begin
                    {ok, [{Ms, Cs}]} = file:consult(filename:join(Dir, Case ++ ".coeffs")),
                    MV = forcer_cmp:form_exprs(filename:join(Dir, Case ++ "_masters.out")),
                    TV = forcer_cmp:form_exprs(filename:join(Dir, Case ++ "_target.out")),
                    {Eq, Ours, Theirs} = forcer_cmp:compare(maps:from_list(lists:zip(Ms, Cs)), Ms,
                                                            [maps:get("M" ++ integer_to_list(K), MV) || K <- lists:seq(1, length(Ms))],
                                                            maps:get("T", TV), Tables),
                    Prec = min(element(2, Ours), element(2, Theirs)),
                    io:format("~-16s Forcer's result, G-scheme: ~s (compared to O(ep^~b))~n",
                              [Case, case Eq of true -> "equal"; false -> "DIFFERENT" end, Prec]),
                    ok = file:write_file(filename:join([Dir, "symbolic", Case ++ ".expansion.txt"]),
                                         [io_lib:format("# ~s: PQCF's coefficients times the masters' values, G-scheme (Forcer's normalisation),~n"
                                                        "# known below O(ep^~b); equal to Forcer's own result: ~p~n", [Case, Prec, Eq]),
                                          epsexp:to_string(Ours)]),
                    {Case, {Eq, Ours}}
                end || {Case, _, _, _, _} <- cases()]),
    Asserts = [begin
                   Total = forcer_cmp:msbar([{element(2, maps:get(C, Series)), O} || {C, O} <- Parts], 4),
                   {_, E} = epsexp:series_of_form(forcer_cmp:asserted(filename:join(Dir, "forcer.frm"), Section, Expr), #{}),
                   Eq = epsexp:equal(Total, E(6)),
                   io:format("~s ~s, MS-bar, against FORM's assertion: ~s~n~s",
                             [Section, Expr, case Eq of true -> "equal"; false -> "DIFFERENT" end, epsexp:to_string(Total)]),
                   ok = file:write_file(filename:join([Dir, "symbolic", Section ++ "_" ++ Expr ++ ".msbar.txt"]),
                                        [io_lib:format("# ~s ~s: PQCF's result, MS-bar (Forcer's msbarexpand=4); equal to FORM's assertion: ~p~n",
                                                       [Section, Expr, Eq]), epsexp:to_string(Total)]),
                   #{<<"check">> => list_to_binary(Section ++ " " ++ Expr), <<"equal">> => Eq,
                     <<"result">> => iolist_to_binary(epsexp:to_string(Total))}
               end || {Section, Expr, Parts} <- assertions()],
    Components = [#{<<"case">> => list_to_binary(C), <<"equal">> => Eq} || C := {Eq, _} <- Series],
    ok = file:write_file(filename:join(Dir, "check.json"), json:encode(#{<<"components">> => Components, <<"assertions">> => Asserts})),
    case lists:all(fun(#{<<"equal">> := E}) -> E end, Components ++ Asserts) of
        true -> ok;
        false -> halt(1)
    end.
