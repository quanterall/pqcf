%% A real reduction problem for the linear-system benchmark (decision record R1, phase 5): integration-by-parts (IBP)
%% identities, Laporta-style, reduced to master integrals -- by PQCF's reconstruction from eliminations modulo primes
%% (pq_recon:linear/4), and by FORM 5 eliminating the same system symbolically under PolyRatFun. The families, their
%% identities, zero sectors and symmetries are ibp_family.erl's:
%%   - the one-loop massless box, seeds {Pos, Neg} (or {box, Pos, Neg}), target I(2,1,1,1);
%%   - the two-loop massless propagator "kite", seeds {kite, Pos, Neg}, target I(1,1,1,1,1);
%%   - FIRE's two-loop double box, seeds {doublebox, R, S} (Kira's r, s), target I(1,1,1,1,1,1,1,-1,-1).
%%
%% The system (system/3): one elimination modulo a prime at a random point learns the independent equations, the
%% masters and, by mark and sweep, the equations the target needs; those are the square system both programs solve.
%%
%%   erl -pa _build/default/lib/*/ebin -pa <dir> -noshell -eval 'ibp:main("<dir>", [{4, 0}, {5, 1}, {kite, 6, 1}])'
-module(ibp).
-export([main/2, system/3, identities/1, zero/1]).

-define(LIMITS, #{sym_timeout_ms => 3600000, sym_max_heap_words => 4000000000, max_expr_nodes => 500000000}).

-define(P, 1000000007).

%%% The identities: ibp_family.erl's, for the family's seeds. The box and the kite keep their seeds of the first
%%% benchmarks (indices in -Neg .. 2, at most Pos positive powers in total, in the non-zero sectors); the double box
%%% has Kira's (r, s).

family(4) -> ibp_family:family(box);
family(5) -> ibp_family:family(kite);
family(9) -> ibp_family:family(doublebox).

zero(A) -> ibp_family:zero(family(tuple_size(A)), A).

%% a list of equations, each a map Integral => CoefficientAst
identities({Pos, Neg}) -> identities({box, Pos, Neg});
identities({doublebox, R, S}) ->
    Fam = ibp_family:family(doublebox),
    ibp_family:identities(Fam, ibp_family:seeds(Fam, R, S));
identities({Family, Pos, Neg}) ->
    Fam = ibp_family:family(Family),
    Size = length(maps:get(props, Fam)),
    Seeds = [list_to_tuple(L) || L <- tuples(Size, lists:seq(-Neg, 2)),
                                 lists:sum([max(0, X) || X <- L]) =< Pos, not zero(list_to_tuple(L))],
    ibp_family:identities(Fam, Seeds).

tuples(0, _) -> [[]];
tuples(N, Range) -> [[X | T] || X <- Range, T <- tuples(N - 1, Range)].

%%% Laporta order and the square system

%% more complex first: positive propagators, then their total power, then numerator powers, then the tuple
weight(A) ->
    L = tuple_to_list(A),
    {length([x || X <- L, X > 0]), lists:sum([X || X <- L, X > 0]), -lists:sum([X || X <- L, X < 0]), A}.

%% The square system for Target, as a map:
%%   a, bs: A (rows as sparse maps Column => entry, 0-based) and the right-hand sides, one per master (sparse maps
%%     Master => entry): A u = Bs, the answer linear in the masters;
%%   b, vars: the same right-hand side with the masters as symbols m1, .. (for FORM);
%%   masters ([{Integral, Symbol}]), target (its column), columns (the unknowns), rows (the equations, Integral =>
%%   coefficient), point (the random point of the learning elimination), size ({equations, integrals} generated).
%%
%% One elimination mod ?P at a random point, columns in Laporta order (most complex first), learns the system
%% (pq_nif:ff_pivots/4): the independent equations, the masters (the columns without a pivot) and, by mark and
%% sweep, the equations the target needs. Those, in the pivot order of their columns, are the square system; its
%% unknowns stay in Laporta order, so that eliminating column k on row k repeats the learned elimination.
system(Seeds, Target, Seed) ->
    Eqs = [E || E <- identities(Seeds), map_size(E) > 0],
    Ints = lists:usort(lists:append([maps:keys(E) || E <- Eqs])),
    lists:member(Target, Ints) orelse error({target_not_in_system, Target}),
    Cols = [I || {_, I} <- lists:reverse(lists:sort([{weight(I), I} || I <- Ints]))],
    Index = maps:from_list(lists:zip(Cols, lists:seq(0, length(Cols) - 1))),
    rand:seed(exsss, {Seed, 1, 1}),
    Point = #{d => rand:uniform(?P - 1), s => rand:uniform(?P - 1), t => rand:uniform(?P - 1)},
    %% the coefficients evaluated mod ?P as they are; only the square system's are canonicalised (below)
    Num = [lists:sort([{maps:get(I, Index), modval(C, Point)} || I := C <- E]) || E <- Eqs],
    TCol = maps:get(Target, Index),
    {PR, PC, Needed} = pq_nif:ff_pivots(?P, Num, length(Cols), [TCol]),
    lists:member(TCol, PC) orelse error(target_is_a_master),
    PivotCol = maps:from_list(lists:zip(PR, PC)),
    EqArr = array:from_list(Eqs),
    ColArr = array:from_list(Cols),
    Rows = [maps:filter(fun(_, C) -> C =/= 0 end,
                        maps:map(fun(_, C) -> {ok, Cc} = pq_symbolic:canonicalise(C), Cc end, array:get(R, EqArr)))
            || R <- Needed],
    UCols = [array:get(maps:get(R, PivotCol), ColArr) || R <- Needed],
    UIndex = maps:from_list(lists:zip(UCols, lists:seq(0, length(UCols) - 1))),
    Masters = [I || {_, I} <- lists:usort([{maps:get(I, Index), I} || Row <- Rows, I <- maps:keys(Row), not maps:is_key(I, UIndex)])],
    MSyms = [{sym, list_to_atom("m" ++ integer_to_list(K))} || K <- lists:seq(1, length(Masters))],
    MIndex = maps:from_list(lists:zip(Masters, MSyms)),
    MPos = maps:from_list(lists:zip(Masters, lists:seq(0, length(Masters) - 1))),
    A = [maps:from_list([{maps:get(I, UIndex), C} || I := C <- Row, maps:is_key(I, UIndex)]) || Row <- Rows],
    B = [{add, [0 | [{mul, [-1, C, maps:get(I, MIndex)]} || I := C <- Row, maps:is_key(I, MIndex)]]} || Row <- Rows],
    Bs = [maps:from_list([{maps:get(I, MPos), {mul, [-1, C]}} || I := C <- Row, maps:is_key(I, MPos)]) || Row <- Rows],
    #{a => A, b => B, bs => Bs, vars => invariants(Target) ++ [X || {sym, X} <- MSyms], masters => lists:zip(Masters, MSyms),
      target => maps:get(Target, UIndex), columns => UCols, rows => Rows, point => Point,
      size => {length(Eqs), length(Cols)}}.

%% For FORM (linsys.erl's elimination): the target the last column, the rows in the pivot order of a forward
%% elimination (row k pivots on column k), so that the target is eliminated last; dense rows and B
form_system(#{rows := Rows, columns := UCols, b := B, point := Point}, Target) ->
    FCols = [I || I <- UCols, I =/= Target] ++ [Target],
    FIndex = maps:from_list(lists:zip(FCols, lists:seq(0, length(FCols) - 1))),
    Num = [lists:sort([{maps:get(I, FIndex), modval(C, Point)} || I := C <- Row, maps:is_key(I, FIndex)]) || Row <- Rows],
    {PR, PC, _} = pq_nif:ff_pivots(?P, Num, length(FCols), []),
    PC =:= lists:seq(0, length(FCols) - 1) orelse error(square_system_singular),
    RowArr = array:from_list(Rows),
    BArr = array:from_list(B),
    {[[maps:get(I, array:get(R, RowArr), 0) || I <- FCols] || R <- PR], [array:get(R, BArr) || R <- PR]}.

invariants(T) -> maps:get(vars, family(tuple_size(T))).

%% a coefficient's value mod ?P at the point (the ASTs of the identities: integers, rationals, symbols, sums,
%% products)
modval(N, _) when is_integer(N) -> (N rem ?P + ?P) rem ?P;
modval({rat, N, D}, _) -> modval(N, #{}) * inv(modval(D, #{}), ?P) rem ?P;
modval({sym, X}, Point) -> maps:get(X, Point);
modval({add, Es}, Point) -> lists:foldl(fun(E, Acc) -> (Acc + modval(E, Point)) rem ?P end, 0, Es);
modval({mul, Es}, Point) -> lists:foldl(fun(E, Acc) -> Acc * modval(E, Point) rem ?P end, 1, Es);
modval(C, Point) ->
    {ok, V} = pq_symbolic:evaluate(C, Point, #{exact => true}),
    modval(V, Point).

inv(A, M) -> pq_ff:inv(A, M).

%%% The run

main(Dir, Cases) ->
    {ok, _} = application:ensure_all_started(pq),
    [begin
         Target0 = case Seeds of
                      {kite, _, _} -> {1, 1, 1, 1, 1};
                      {doublebox, _, _} -> {1, 1, 1, 1, 1, 1, 1, -1, -1};
                      _ -> {2, 1, 1, 1}
                  end,
         %% the target as the identities name it: its canonical image under the sector symmetries
         Target = ibp_family:canon(family(tuple_size(Target0)), Target0),
         {Ts, #{a := A, bs := Bs, vars := Vars, masters := Masters, target := TI, size := {NEqs, NInts}} = Sys} =
             timer:tc(?MODULE, system, [Seeds, Target, 1]),
         N = length(A),
         %% the answer is linear in the masters: one right-hand side per master, each coefficient a rational function
         %% of d, s, t alone; one elimination per point serves them all, and the zero ones are found at once
         {T, {ok, [{TI, PerMaster}]}} = timer:tc(pq_recon, linear, [A, Bs, invariants(Target), ?LIMITS#{components => [TI], rhs => length(Masters)}]),
         PerMap = maps:from_list(PerMaster),
         Coeffs = [{M, Sym, maps:get(C, PerMap)} || {C, {M, Sym}} <- lists:enumerate(0, Masters)],
         NonZero = [{M, Sym, NumDen} || {M, Sym, {Nu, _} = NumDen} <- Coeffs, map_size(Nu) > 0],
         P = fun(E) -> binary_to_list(pq_symbolic:pretty(E)) end,
         Name = case Seeds of
                    {F, {dots, Po}, Ne} -> io_lib:format("ibp_~s_d~b_~b", [F, Po, Ne]);
                    {F, Po, Ne} -> io_lib:format("ibp_~s_~b_~b", [F, Po, Ne]);
                    {Po, Ne} -> io_lib:format("ibp_~b_~b", [Po, Ne])
                end,
         Answer = [{P(pq_sym_norm:to_ast(Nu)), P(pq_sym_norm:to_ast(De)), atom_to_list(X)} || {_, {sym, X}, {Nu, De}} <- NonZero],
         ok = file:write_file(filename:join(Dir, Name ++ ".answer"),
                              [[io_lib:format("~w", [M]), " ", Nu, " / ", De, "\n"] || {{M, _}, {Nu, De, _}} <- lists:zip([{M, S} || {M, S, _} <- NonZero], Answer)]),
         %% the FORM script (FORM eliminates from the last column: the columns and rows reversed); only for systems
         %% FORM can be expected to solve
         case N =< 4000 of
             true ->
                 {FA, FB} = form_system(Sys, Target),
                 Rev = fun(L) -> lists:reverse(L) end,
                 ok = file:write_file(filename:join(Dir, Name ++ ".frm"), form_script(Rev([Rev(Row) || Row <- FA]), Rev(FB), Vars, Answer));
             false -> ok
         end,
         io:format("IBP, seeds ~p: ~b equations in ~b integrals, ~b needed (built in ~.1f s), ~b masters; I~w reduced in ~.2f s "
                   "to ~b of them: ~s~n",
                   [Seeds, NEqs, NInts, N, Ts / 1.0e6, length(Masters), Target, T / 1.0e6, length(NonZero),
                    lists:join(", ", [io_lib:format("~w (~b/~b terms)", [M, map_size(Nu), map_size(De)]) || {M, _, {Nu, De}} <- NonZero])])
     end || Seeds <- Cases],
    ok.

%% FORM: the elimination of linsys.erl (unknown 1 the target), the masters symbols; at the end U minus PQCF's answer,
%% sum_j rat(N_j m_j, D_j), must be 0
form_script(A, B, Vars, Answer) ->
    Base = linsys:form_script(A, B, Vars, "0", "1"),
    Sum = lists:join(" + ", [["rat((", Nu, ")*", M, ",", De, ")"] || {Nu, De, M} <- Answer]),
    Check = ["Local Z = U - (", case Sum of [] -> "0"; _ -> Sum end, ");\n"],
    re:replace(Base, "Local Z = U - rat\\(\n0,\n1\\);\n", lists:flatten(Check), [{return, list}]).

