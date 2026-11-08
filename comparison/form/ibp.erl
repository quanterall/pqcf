%% The IBP benchmarks (decision record R1): reductions by PQCF (pq_ibp, pq_ibp_family) of the families' targets, each
%% coefficient simplified (pq_simplify) and, for systems FORM can be expected to solve, a FORM 5 script eliminating
%% the same square system symbolically under PolyRatFun and checking U - PQCF's answer = 0.
%%   - the one-loop massless box, seeds {Pos, Neg} (or {box, Pos, Neg}: indices in -Neg .. 2, at most Pos positive
%%     powers in total), target I(2,1,1,1);
%%   - the two-loop massless propagator "kite", seeds {kite, Pos, Neg}, target I(1,1,1,1,1);
%%   - FIRE's two-loop double box, seeds {doublebox, R, S} (Kira's r, s), target I(1,1,1,1,1,1,1,-1,-1);
%%   - {Family, {dots, D}, S}: at most D dots per sector; {auto, Family}, {auto, Family, Target} and
%%     {auto, Family, Target, SystemOpts}: the smallest complete seeds (pq_ibp:auto_seeds/3).
%%
%%   erl -pa _build/default/lib/*/ebin -pa <dir> -noshell -eval 'ibp:main("<dir>", [{5, 1}, {kite, 6, 1}, {auto, doublebox}])'
-module(ibp).
-export([main/2]).

-define(LIMITS, #{sym_timeout_ms => 3600000, sym_max_heap_words => 4000000000, max_expr_nodes => 500000000}).
-define(P, 1000000007).

family(Size) when is_integer(Size) -> pq_ibp_family:example(case Size of 4 -> box; 5 -> kite; 9 -> doublebox end);
family(Name) -> pq_ibp_family:example(Name).

target(kite) -> {1, 1, 1, 1, 1};
target(doublebox) -> {1, 1, 1, 1, 1, 1, 1, -1, -1};
target(box) -> {2, 1, 1, 1}.

%% the seed integrals of the first benchmarks: indices in -Neg .. 2, at most Pos positive powers, non-zero sectors
old_seeds(Fam, Pos, Neg) ->
    [list_to_tuple(L) || L <- tuples(pq_ibp_family:size(Fam), lists:seq(-Neg, 2)),
                         lists:sum([max(0, X) || X <- L]) =< Pos, not pq_ibp_family:zero(Fam, list_to_tuple(L))].

tuples(0, _) -> [[]];
tuples(N, Range) -> [[X | T] || X <- Range, T <- tuples(N - 1, Range)].

%% {Family, Seeds | auto, Target, SystemOpts, Name}
case_of({Pos, Neg}) when is_integer(Pos) -> case_of({box, Pos, Neg});
case_of({auto, F}) -> case_of({auto, F, target(F)});
case_of({auto, F, T}) -> case_of({auto, F, T, #{}});
case_of({auto, F, T, O}) ->
    {F, auto, T, O, io_lib:format("ibp_~s_auto_~s~s", [F, lists:join("_", [integer_to_list(X) || X <- tuple_to_list(T)]),
                                                       case O of #{order := Ord} -> ["_", atom_to_list(Ord)]; _ -> "" end])};
case_of({F, {dots, D}, S}) -> {F, {dots, D, S}, target(F), #{}, io_lib:format("ibp_~s_d~b_~b", [F, D, S])};
case_of({doublebox, R, S}) -> {doublebox, {R, S}, target(doublebox), #{}, io_lib:format("ibp_doublebox_~b_~b", [R, S])};
case_of({F, Pos, Neg}) ->
    {F, {old, Pos, Neg}, target(F), #{}, case F of box -> io_lib:format("ibp_~b_~b", [Pos, Neg]); _ -> io_lib:format("ibp_~s_~b_~b", [F, Pos, Neg]) end}.

main(Dir, Cases) ->
    {ok, _} = application:ensure_all_started(pq),
    [begin
         {F, SeedSpec, Target0, SysOpts, Name} = case_of(Case),
         Fam = family(F),
         Vars0 = pq_ibp_family:vars(Fam),
         %% the target as the identities name it: its canonical image under the sector symmetries
         Target = pq_ibp_family:canon(Fam, Target0),
         {Ts, {Spec, Sys}} =
             timer:tc(fun() ->
                              case SeedSpec of
                                  auto -> pq_ibp:auto_seeds(Fam, [Target], SysOpts);
                                  {old, Pos, Neg} -> {Case, pq_ibp:system(Fam, old_seeds(Fam, Pos, Neg), [Target], SysOpts)};
                                  _ -> {SeedSpec, pq_ibp:system(Fam, pq_ibp:seeds(Fam, SeedSpec, [Target]), [Target], SysOpts)}
                              end
                      end),
         #{a := A, bs := Bs, masters := Masters, targets := [{Target, TI}], size := {NEqs, NInts}} = Sys,
         N = length(A),
         %% the answer is linear in the masters: one right-hand side per master, each coefficient a rational function
         %% of the invariants alone; one elimination per point serves them all, and the zero ones are found at once
         {T, {ok, [{TI, PerMaster}]}} = timer:tc(pq_recon, linear, [A, Bs, Vars0, maps:merge(?LIMITS, #{components => [TI], rhs => length(Masters)})]),
         PerMap = maps:from_list(PerMaster),
         NonZero = [{M, Sym, ND} || {C, {M, Sym}} <- lists:enumerate(0, Masters), {Nu, _} = ND <- [maps:get(C, PerMap)], map_size(Nu) > 0],
         P = fun(E) -> binary_to_list(pq_symbolic:pretty(E)) end,
         Answer = [{P(pq_sym_norm:to_ast(Nu)), P(pq_sym_norm:to_ast(De)), atom_to_list(X)} || {_, {sym, X}, {Nu, De}} <- NonZero],
         ok = file:write_file(filename:join(Dir, Name ++ ".answer"),
                              [[io_lib:format("~w", [M]), " ", Nu, " / ", De, "\n"] || {{M, _, _}, {Nu, De, _}} <- lists:zip(NonZero, Answer)]),
         %% the coefficients simplified (pq_simplify), each verified against the reconstruction
         Simple = [begin
                       {Ast, St} = pq_simplify:rational(Nu, De, Vars0, #{}),
                       true = pq_simplify:verify(Ast, Nu, De, Vars0),
                       {M, Ast, St}
                   end || {M, _, {Nu, De}} <- NonZero],
         ok = file:write_file(filename:join(Dir, Name ++ ".simplified"),
                              [[io_lib:format("~w", [M]), " ", pq_simplify:format(Ast), "\n"] || {M, Ast, _} <- Simple]),
         %% the FORM script (FORM eliminates from the last column: the columns and rows reversed); only for systems
         %% FORM can be expected to solve
         case N =< 4000 of
             true ->
                 {FA, FB} = form_system(Sys, Target),
                 Rev = fun(L) -> lists:reverse(L) end,
                 FVars = Vars0 ++ [X || {_, {sym, X}} <- Masters],
                 ok = file:write_file(filename:join(Dir, Name ++ ".frm"), form_script(Rev([Rev(Row) || Row <- FA]), Rev(FB), FVars, Answer));
             false -> ok
         end,
         io:format("IBP ~s, seeds ~w: ~b equations in ~b integrals, ~b needed (built in ~.1f s), ~b masters; I~w reduced in ~.2f s "
                   "to ~b of them: ~s~n  simplified: ~b additions -> ~b~n",
                   [F, Spec, NEqs, NInts, N, Ts / 1.0e6, length(Masters), Target, T / 1.0e6, length(NonZero),
                    lists:join(", ", [io_lib:format("~w (~b/~b terms)", [M, map_size(Nu), map_size(De)]) || {M, _, {Nu, De}} <- NonZero]),
                    lists:sum([I || {_, _, #{size_in := I}} <- Simple]), lists:sum([O || {_, _, #{size_out := O}} <- Simple])])
     end || Case <- Cases],
    ok.

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

modval(C, Point) ->
    {ok, V} = pq_symbolic:evaluate(C, Point, #{exact => true}),
    case V of
        N when is_integer(N) -> (N rem ?P + ?P) rem ?P;
        {rat, N, D} -> (N rem ?P + ?P) rem ?P * pq_ff:inv((D rem ?P + ?P) rem ?P, ?P) rem ?P
    end.

%% FORM: the elimination of linsys.erl (unknown 1 the target), the masters symbols; at the end U minus PQCF's answer,
%% sum_j rat(N_j m_j, D_j), must be 0
form_script(A, B, Vars, Answer) ->
    Base = linsys:form_script(A, B, Vars, "0", "1"),
    Sum = lists:join(" + ", [["rat((", Nu, ")*", M, ",", De, ")"] || {Nu, De, M} <- Answer]),
    Check = ["Local Z = U - (", case Sum of [] -> "0"; _ -> Sum end, ");\n"],
    re:replace(Base, "Local Z = U - rat\\(\n0,\n1\\);\n", lists:flatten(Check), [{return, list}]).
