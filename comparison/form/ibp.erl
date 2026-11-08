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
-export([main/2, system/3, system/4, identities/1, zero/1, auto/2]).

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
identities({Family, R, S}) when Family =:= doublebox; is_tuple(R) ->
    Fam = ibp_family:family(Family),
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

%% more complex first: positive propagators, then their total power, then numerator powers, then the tuple (order
%% numerators: within a sector a numerator is simpler than a dot, so masters carry numerators); order dots: the
%% total of the powers and numerators first, then the numerators (masters carry dots, as FIRE's do)
weight(A, Order) ->
    L = tuple_to_list(A),
    Pos = lists:sum([X || X <- L, X > 0]),
    Neg = -lists:sum([X || X <- L, X < 0]),
    case Order of
        numerators -> {length([x || X <- L, X > 0]), Pos, Neg, A};
        dots -> {length([x || X <- L, X > 0]), Pos + Neg, Neg, A}
    end.

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
system(Seeds, Target, Seed) -> system(Seeds, Target, Seed, #{}).

%% Opts: order => numerators (default) | dots, the Laporta order that picks the masters (weight/2)
system(Seeds, Target, Seed, Opts) ->
    Eqs = [E || E <- identities(Seeds), map_size(E) > 0],
    Ints = lists:usort(lists:append([maps:keys(E) || E <- Eqs])),
    lists:member(Target, Ints) orelse error({target_not_in_system, Target}),
    Order = maps:get(order, Opts, numerators),
    Cols = [I || {_, I} <- lists:reverse(lists:sort([{weight(I, Order), I} || I <- Ints]))],
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

auto(Family, Target) -> auto(Family, Target, #{}).

auto(Family, Target, Opts) ->
    {Ta, {Chosen, Tried}} = timer:tc(fun() -> auto_seeds(Family, Target, Opts) end),
    io:format("auto seeds for ~s ~w: ~w (tried ~w, ~.1f s)~n", [Family, Target, Chosen, Tried, Ta / 1.0e6]),
    Chosen.

target(kite) -> {1, 1, 1, 1, 1};
target(doublebox) -> {1, 1, 1, 1, 1, 1, 1, -1, -1};
target(box) -> {2, 1, 1, 1}.

%% The smallest seeds that reduce the target completely, without learning: {dots, D} and numerators up to S, from
%% the target's own (its dots, its numerator powers) upward, the cheaper first. Complete: every master with a
%% non-zero coefficient (one numeric solve at a random point mod ?P) is itself a seed, its identities generated --
%% a master at the edge of the seed range, whose identities were never written, means the seeds are too small.
%% {Seeds, [Tried]}
auto_seeds(Family, Target, Opts) ->
    Fam = ibp_family:family(Family),
    L = tuple_to_list(Target),
    Dt = lists:sum([X - 1 || X <- L, X > 1]),
    St = lists:sum([-X || X <- L, X < 0]),
    Candidates = lists:sort(fun({D1, S1}, {D2, S2}) -> {D1 + S1, D1} =< {D2 + S2, D2} end,
                            [{D, S} || D <- lists:seq(Dt, Dt + 2), S <- lists:seq(St, St + 2)]),
    auto_seeds(Fam, Family, Target, Opts, Candidates, []).

auto_seeds(_, Family, _, _, [], Tried) -> error({no_seeds_reduce_the_target, Family, lists:reverse(Tried)});
auto_seeds(Fam, Family, Target, Opts, [{D, S} | More], Tried) ->
    Seeds = {Family, {dots, D}, S},
    Complete = try system(Seeds, Target, 1, Opts) of
                   Sys ->
                       SeedSet = sets:from_list([ibp_family:canon(Fam, X) || X <- ibp_family:seeds(Fam, {dots, D}, S)], [{version, 2}]),
                       lists:all(fun(M) -> sets:is_element(M, SeedSet) end, nonzero_masters(Sys))
               catch
                   error:{target_not_in_system, _} -> false;
                   error:target_is_a_master -> false
               end,
    case Complete of
        true -> {Seeds, lists:reverse([{D, S} | Tried])};
        false -> auto_seeds(Fam, Family, Target, Opts, More, [{D, S} | Tried])
    end.

%% the masters with a non-zero coefficient in the target, from one numeric solve at a random point mod ?P
nonzero_masters(#{a := A, bs := Bs, masters := Masters, target := TI, vars := Vars0}) ->
    Vars = [V || V <- Vars0, lists:member(V, [d, s, t])],
    Sp = fun(Row) -> lists:sort(maps:to_list(Row)) end,
    SA = [Sp(R) || R <- A], SB = [Sp(R) || R <- Bs],
    {ok, Prog} = pq_symbolic:compile([E || Row <- SA ++ SB, {_, E} <- Row], Vars, ?LIMITS#{exact_constants => true}),
    {ok, {ff, Ref, _, _}} = pq_ff:prepare(Prog),
    {ok, Lin} = pq_nif:ff_linsys(Ref, length(A), length(Masters), TI, 0, [[C || {C, _} <- R] || R <- SA], [[C || {C, _} <- R] || R <- SB], []),
    {ok, [Vals]} = pq_nif:ff_linsys_all(Lin, ?P, [rand:uniform(?P - 1) || _ <- Vars]),
    [M || {{M, _}, V} <- lists:zip(Masters, Vals), V =/= 0].

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
         SysOpts = case Cases0 of {auto, _, _, O} -> O; _ -> #{} end,
         Target0 = case Cases0 of
                       {auto, _, T0, _} -> T0;
                       {auto, _, T0} -> T0;
                       {auto, Fm0} -> target(Fm0);
                       {P0, _} when is_integer(P0) -> target(box);
                       _ -> target(element(1, Cases0))
                   end,
         %% the target as the identities name it: its canonical image under the sector symmetries
         Target = ibp_family:canon(family(tuple_size(Target0)), Target0),
         %% {auto, Family} or {auto, Family, Target}: the smallest seeds that reduce the target completely
         %% (auto_seeds/2)
         Seeds = case Cases0 of
                     {auto, Fm} -> auto(Fm, Target);
                     {auto, Fm, _} -> auto(Fm, Target);
                     {auto, Fm, _, _} -> auto(Fm, Target, SysOpts);
                     _ -> Cases0
                 end,
         {Ts, #{a := A, bs := Bs, vars := Vars, masters := Masters, target := TI, size := {NEqs, NInts}} = Sys} =
             timer:tc(?MODULE, system, [Seeds, Target, 1, SysOpts]),
         N = length(A),
         %% the answer is linear in the masters: one right-hand side per master, each coefficient a rational function
         %% of d, s, t alone; one elimination per point serves them all, and the zero ones are found at once
         {T, {ok, [{TI, PerMaster}]}} = timer:tc(pq_recon, linear, [A, Bs, invariants(Target), ?LIMITS#{components => [TI], rhs => length(Masters)}]),
         PerMap = maps:from_list(PerMaster),
         Coeffs = [{M, Sym, maps:get(C, PerMap)} || {C, {M, Sym}} <- lists:enumerate(0, Masters)],
         NonZero = [{M, Sym, NumDen} || {M, Sym, {Nu, _} = NumDen} <- Coeffs, map_size(Nu) > 0],
         P = fun(E) -> binary_to_list(pq_symbolic:pretty(E)) end,
         Name = case Seeds of
                    {F, {dots, Po}, Ne} when Cases0 =:= Seeds -> io_lib:format("ibp_~s_d~b_~b", [F, Po, Ne]);
                    {F, {dots, Po}, Ne} when SysOpts =/= #{} ->
                        io_lib:format("ibp_~s_d~b_~b_~s_~s", [F, Po, Ne, maps:get(order, SysOpts, numerators),
                                                              lists:join("_", [integer_to_list(X) || X <- tuple_to_list(Target)])]);
                    {F, {dots, Po}, Ne} -> io_lib:format("ibp_~s_d~b_~b_~s", [F, Po, Ne, lists:join("_", [integer_to_list(X) || X <- tuple_to_list(Target)])]);
                    {F, Po, Ne} -> io_lib:format("ibp_~s_~b_~b", [F, Po, Ne]);
                    {Po, Ne} -> io_lib:format("ibp_~b_~b", [Po, Ne])
                end,
         Answer = [{P(pq_sym_norm:to_ast(Nu)), P(pq_sym_norm:to_ast(De)), atom_to_list(X)} || {_, {sym, X}, {Nu, De}} <- NonZero],
         ok = file:write_file(filename:join(Dir, Name ++ ".answer"),
                              [[io_lib:format("~w", [M]), " ", Nu, " / ", De, "\n"] || {{M, _}, {Nu, De, _}} <- lists:zip([{M, S} || {M, S, _} <- NonZero], Answer)]),
         %% the coefficients simplified (pq_simplify: factors from an alphabet, verified against the reconstruction)
         Simple = [begin
                       {Ast, St} = pq_simplify:rational(Nu, De, invariants(Target), #{}),
                       true = pq_simplify:verify(Ast, Nu, De, invariants(Target)),
                       {M, Ast, St}
                   end || {M, _, {Nu, De}} <- NonZero],
         ok = file:write_file(filename:join(Dir, Name ++ ".simplified"),
                              [[io_lib:format("~w", [M]), " ", pq_simplify:format(Ast), "\n"] || {M, Ast, _} <- Simple]),
         io:format("  simplified: ~b additions -> ~b (forms ~w)~n",
                   [lists:sum([I || {_, _, #{size_in := I}} <- Simple]), lists:sum([O || {_, _, #{size_out := O}} <- Simple]),
                    [F || {_, _, #{form := F}} <- Simple]]),
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
     end || Cases0 <- Cases],
    ok.

%% FORM: the elimination of linsys.erl (unknown 1 the target), the masters symbols; at the end U minus PQCF's answer,
%% sum_j rat(N_j m_j, D_j), must be 0
form_script(A, B, Vars, Answer) ->
    Base = linsys:form_script(A, B, Vars, "0", "1"),
    Sum = lists:join(" + ", [["rat((", Nu, ")*", M, ",", De, ")"] || {Nu, De, M} <- Answer]),
    Check = ["Local Z = U - (", case Sum of [] -> "0"; _ -> Sum end, ");\n"],
    re:replace(Base, "Local Z = U - rat\\(\n0,\n1\\);\n", lists:flatten(Check), [{return, list}]).

