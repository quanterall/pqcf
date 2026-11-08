%% Reducing Feynman integrals to master integrals by integration-by-parts identities (decision record R1): the
%% identities of a family (pq_ibp_family) over seeds chosen automatically, one learning elimination modulo a prime
%% that keeps only the equations the targets need (mark and sweep), and the targets' coefficients reconstructed from
%% eliminations at sample points (pq_recon:linear/4: denominators guessed, the elimination learned and replayed).
%% Results are exact; the work spreads over pq_rexec's executors like any reconstruction.
%%
%%   Fam = pq_ibp_family:example(doublebox),
%%   {ok, #{T := Coefficients}, Info} = pq_ibp:reduce(Fam, [T = {1,1,1,1,1,1,1,-1,-1}], #{})
%%
%% Coefficients: [{Master, {N, D}}], N / D the reduced rational function in the family's variables (pq_recon's
%% polynomial maps; pq_simplify:rational/4 writes them in factored form). Masters and targets are integrals as
%% tuples, named by their canonical image under the family's symmetries (pq_ibp_family:canon/2).
-module(pq_ibp).
-export([reduce/3, system/4, seeds/3, auto_seeds/3, differential_equation/4]).

-define(P, 1000000007).
-define(LIMITS, #{sym_timeout_ms => 3600000, sym_max_heap_words => 4000000000, max_expr_nodes => 500000000}).
%% compiling a system's entries: exact constants, the limits raised for large systems
-define(COMPILE, #{sym_timeout_ms => 3600000, sym_max_heap_words => 4000000000, max_expr_nodes => 500000000,
                   exact_constants => true}).

-type integral() :: pq_ibp_family:integral().
-type ratfun() :: {map(), map()}.

%%% Reduction

%% Opts: seeds => auto (default) | {dots, D, S} | {R, S} (Kira's r and s) | [Seed]; order => numerators (default) |
%% dots (which integrals become masters, weight/2); and pq_recon's options (executor, workers, guess, ...).
%% {ok, #{Target => [{Master, {N, D}}]}, Info}: a zero target maps to [], a target that is a master to itself;
%% Info: #{seeds, equations, integrals, needed, masters}.
-spec reduce(pq_ibp_family:family(), [integral()], map()) -> {ok, #{integral() => [{integral(), ratfun()}]}, map()} | {error, term()}.
reduce(Fam, Targets, Opts) ->
    N = pq_ibp_family:size(Fam),
    case [T || T <- Targets, not (is_tuple(T) andalso tuple_size(T) =:= N andalso lists:all(fun is_integer/1, tuple_to_list(T)))] of
        [_ | _] = Bad -> {error, {bad_targets, Bad}};
        [] ->
            Canon = maps:from_list([{T, pq_ibp_family:canon(Fam, T)} || T <- Targets]),
            Live = lists:usort([C || _ := C <- Canon, not pq_ibp_family:zero(Fam, C)]),
            try reduce_live(Fam, Live, Opts) of
                {Red, Info} ->
                    {ok, maps:from_list([{T, case pq_ibp_family:zero(Fam, C) of
                                                 true -> [];
                                                 false -> maps:get(C, Red)
                                             end} || T := C <- Canon]), Info}
            catch
                throw:{pq_ibp, Why} -> {error, Why}
            end
    end.

reduce_live(_Fam, [], _Opts) -> {#{}, #{}};
reduce_live(Fam, Live, Opts) ->
    {SeedSpec, Sys} = case maps:get(seeds, Opts, auto) of
                          auto -> auto_seeds(Fam, Live, Opts);
                          Spec -> {Spec, system(Fam, seeds(Fam, Spec, Live), Live, Opts)}
                      end,
    #{a := A, bs := Bs, masters := Masters, targets := TIs, size := {NEqs, NInts}} = Sys,
    Vars = pq_ibp_family:vars(Fam),
    One = {#{{[], []} => {q, {1, 1}, {0, 1}}}, #{{[], []} => {q, {1, 1}, {0, 1}}}},
    %% targets that are masters themselves, and the others' columns
    Reduced = case [TI || {_, TI} <- TIs, is_integer(TI)] of
                  [] -> #{};
                  Cols ->
                      case pq_recon:linear(A, Bs, Vars, maps:merge(?LIMITS, Opts#{components => Cols, rhs => length(Masters)})) of
                          {ok, Per} -> maps:from_list(Per);
                          {error, Why} -> throw({pq_ibp, Why})
                      end
              end,
    Red = maps:from_list(
            [{T, case TI of
                     {master, _} -> [{T, One}];
                     _ -> [{M, ND} || {C, ND = {Nu, _}} <- maps:get(TI, Reduced), map_size(Nu) > 0, {M, _} <- [lists:nth(C + 1, Masters)]]
                 end} || {T, TI} <- TIs]),
    {Red, #{seeds => SeedSpec, equations => NEqs, integrals => NInts, needed => length(A), masters => length(Masters)}}.

%%% Seeds

%% the seed integrals: {dots, D, S}: in every non-zero sector, at most D powers above 1 and S numerator powers in
%% total; {R, S}: Kira's (the positive indices summing to at most R); a list: as given
-spec seeds(pq_ibp_family:family(), term(), [integral()]) -> [integral()].
seeds(Fam, {dots, D, S}, _Targets) -> pq_ibp_family:seeds(Fam, {dots, D}, S);
seeds(Fam, {R, S}, _Targets) when is_integer(R) -> pq_ibp_family:seeds(Fam, R, S);
seeds(_Fam, List, _Targets) when is_list(List) -> List.

%% The smallest seeds that reduce every target completely, without learning: {dots, D, S} from the targets' own
%% dots and numerator powers upward, the cheaper first. Complete: every master with a non-zero coefficient in a
%% target (one numeric solve at a random point mod ?P) is itself a seed, its identities generated; a master at the
%% edge of the seed range, whose identities were never written, means the seeds are too small. {Spec, System}
%% (the system as system/4 gives it)
-spec auto_seeds(pq_ibp_family:family(), [integral()], map()) -> {{dots, non_neg_integer(), non_neg_integer()}, map()}.
auto_seeds(Fam, Targets, Opts) ->
    Dt = lists:max([lists:sum([X - 1 || X <- tuple_to_list(T), X > 1]) || T <- Targets]),
    St = lists:max([lists:sum([-X || X <- tuple_to_list(T), X < 0]) || T <- Targets]),
    Candidates = lists:sort(fun({D1, S1}, {D2, S2}) -> {D1 + S1, D1} =< {D2 + S2, D2} end,
                            [{D, S} || D <- lists:seq(Dt, Dt + 2), S <- lists:seq(St, St + 2)]),
    auto_seeds(Fam, Targets, Opts, Candidates).

auto_seeds(_Fam, _Targets, _Opts, []) -> throw({pq_ibp, no_seeds_reduce_the_targets});
auto_seeds(Fam, Targets, Opts, [{D, S} | More]) ->
    Seeds = pq_ibp_family:seeds(Fam, {dots, D}, S),
    Result = try system(Fam, Seeds, Targets, Opts) of
                 Sys ->
                     SeedSet = sets:from_list([pq_ibp_family:canon(Fam, X) || X <- Seeds], [{version, 2}]),
                     case lists:all(fun(M) -> sets:is_element(M, SeedSet) end, nonzero_masters(Fam, Sys)) of
                         true -> {ok, Sys};
                         false -> incomplete
                     end
             catch
                 throw:{pq_ibp, {targets_not_in_system, _}} -> incomplete
             end,
    case Result of
        {ok, Sys2} -> {{dots, D, S}, Sys2};
        incomplete -> auto_seeds(Fam, Targets, Opts, More)
    end.

%% the masters with a non-zero coefficient in any target, from one numeric solve at a random point mod ?P
nonzero_masters(Fam, #{a := A, bs := Bs, masters := Masters, targets := TIs}) ->
    Vars = pq_ibp_family:vars(Fam),
    Sp = fun(Row) -> lists:sort(maps:to_list(Row)) end,
    SA = [Sp(R) || R <- A], SB = [Sp(R) || R <- Bs],
    {ok, Prog} = pq_symbolic:compile([E || Row <- SA ++ SB, {_, E} <- Row], Vars, ?COMPILE),
    {ok, {ff, Ref, _, _}} = pq_ff:prepare(Prog),
    PointMap = point(Vars, {3, 1, 4}),
    Point = [maps:get(V, PointMap) || V <- Vars],
    lists:usort(
      [M || {_, TI} <- TIs, is_integer(TI),
            {ok, Lin} <- [pq_nif:ff_linsys(Ref, length(A), length(Masters), TI, 0, [[C || {C, _} <- R] || R <- SA],
                                           [[C || {C, _} <- R] || R <- SB], [])],
            {ok, [Vals]} <- [pq_nif:ff_linsys_all(Lin, ?P, Point)],
            {{M, _}, V} <- lists:zip(Masters, Vals), V =/= 0]).

%%% The system

%% The square system for the targets, as a map:
%%   a, bs: A (rows as sparse maps Column => entry, 0-based) and the right-hand sides, one per master (sparse maps
%%     Master => entry): A u = Bs, the answer linear in the masters;
%%   b: the same right-hand side with the masters as symbols m1, ..;
%%   masters ([{Integral, Symbol}]), targets ([{Target, Column | {master, Integral}}]), columns (the unknowns), rows
%%   (the equations, Integral => coefficient), point (the learning elimination's), size ({equations, integrals}).
%%
%% One elimination mod ?P at a random point, columns in Laporta order (most complex first), learns the system
%% (pq_nif:ff_pivots/4): the independent equations, the masters (the columns without a pivot) and, by mark and sweep,
%% the equations the targets need. Those, in the pivot order of their columns, are the square system; its unknowns
%% stay in Laporta order, so that eliminating column k on row k repeats the learned elimination.
-spec system(pq_ibp_family:family(), [integral()], [integral()], map()) -> map().
system(Fam, Seeds, Targets, Opts) ->
    Eqs = [E || E <- pq_ibp_family:identities(Fam, Seeds), map_size(E) > 0],
    Ints = lists:usort(lists:append([maps:keys(E) || E <- Eqs])),
    case [T || T <- Targets, not lists:member(T, Ints)] of
        [] -> ok;
        Missing -> throw({pq_ibp, {targets_not_in_system, Missing}})
    end,
    Order = maps:get(order, Opts, numerators),
    Cols = [I || {_, I} <- lists:reverse(lists:sort([{weight(I, Order), I} || I <- Ints]))],
    Index = maps:from_list(lists:zip(Cols, lists:seq(0, length(Cols) - 1))),
    Point = point(pq_ibp_family:vars(Fam), {1, 1, 1}),
    %% the coefficients evaluated mod ?P as they are; only the square system's are canonicalised (below)
    Num = [lists:sort([{maps:get(I, Index), modval(C, Point)} || I := C <- E]) || E <- Eqs],
    TCols = [maps:get(T, Index) || T <- Targets],
    {PR, PC, Needed} = pq_nif:ff_pivots(?P, Num, length(Cols), TCols),
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
    TIs = [{T, case maps:find(T, UIndex) of {ok, Ix} -> Ix; error -> {master, T} end} || T <- Targets],
    #{a => A, b => B, bs => Bs, masters => lists:zip(Masters, MSyms), targets => TIs, columns => UCols, rows => Rows,
      point => Point, size => {length(Eqs), length(Cols)}}.

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

%% a coefficient's value mod ?P at the point (the identities' ASTs: integers, rationals, symbols, sums, products)
modval(N, _) when is_integer(N) -> (N rem ?P + ?P) rem ?P;
modval({rat, N, D}, _) -> modval(N, #{}) * pq_ff:inv(modval(D, #{}), ?P) rem ?P;
modval({sym, X}, Point) -> maps:get(X, Point);
modval({add, Es}, Point) -> lists:foldl(fun(E, Acc) -> (Acc + modval(E, Point)) rem ?P end, 0, Es);
modval({mul, Es}, Point) -> lists:foldl(fun(E, Acc) -> Acc * modval(E, Point) rem ?P end, 1, Es);
modval(C, Point) ->
    {ok, V} = pq_symbolic:evaluate(C, Point, #{exact => true}),
    modval(V, Point).

point(Vars, Seed) ->
    {L, _} = lists:mapfoldl(fun(V, S) -> {X, S2} = rand:uniform_s(?P - 1, S), {{V, X}, S2} end, rand:seed_s(exsss, Seed), Vars),
    maps:from_list(L).

%%% Differential equations

%% d M_i / d Var = sum_j E_ij M_j for the masters M (canonical integrals): each master's derivative
%% (pq_ibp_family:derivative/3) reduced to the masters (reduce/3, all its integrals at once), each entry reconstructed
%% as one reduced rational function. {ok, [[{N, D}]]} (rows in the masters' order) | {error, _}; an integral reducing
%% to masters outside the list is an error ({masters_incomplete, Extra}).
-spec differential_equation(pq_ibp_family:family(), [integral()], atom(), map()) -> {ok, [[ratfun()]]} | {error, term()}.
differential_equation(Fam, Masters, Var, Opts) ->
    Vars = pq_ibp_family:vars(Fam),
    Terms = [pq_ibp_family:derivative(Fam, M, Var) || M <- Masters],
    Others = lists:usort([I || Ts <- Terms, {I, _} <- Ts, not lists:member(I, Masters)]),
    case case Others of [] -> {ok, #{}, #{}}; _ -> reduce(Fam, Others, Opts) end of
        {error, _} = E -> E;
        {ok, Red, _} ->
            In = fun(I) -> case lists:member(I, Masters) of true -> [{I, 1}]; false -> [{M, ratast(ND)} || {M, ND} <- maps:get(I, Red)] end end,
            Extra = lists:usort([M || Ts <- Terms, {I, _} <- Ts, {M, _} <- In(I), not lists:member(M, Masters)]),
            case Extra of
                [_ | _] -> {error, {masters_incomplete, Extra}};
                [] ->
                    Entry = fun(Ts, Mj) -> {add, [0 | [{mul, [C, X]} || {I, C} <- Ts, {M, X} <- In(I), M =:= Mj]]} end,
                    {ok, [[rational(Entry(Ts, Mj), Vars) || Mj <- Masters] || Ts <- Terms]}
            end
    end.

ratast({N, D}) -> {mul, [pq_sym_norm:to_ast(N), {pow, pq_sym_norm:to_ast(D), -1}]}.

rational(Ast, Vars) ->
    {ok, Prog} = pq_symbolic:compile([Ast], Vars, ?COMPILE),
    {ok, ND} = pq_recon:rational(Prog, Vars, #{}),
    ND.
