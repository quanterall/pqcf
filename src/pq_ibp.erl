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
-export([reduce/3, reduce_sum/3, reduce_sums/3, system/4, seeds/3, auto_seeds/3, differential_equation/4, choose_masters/3, size_of/2]).

-define(P, 1000000007).
-define(LIMITS, #{sym_timeout_ms => 3600000, sym_max_heap_words => 4000000000, max_expr_nodes => 500000000}).
%% compiling a system's entries: exact constants, the limits raised for large systems
-define(COMPILE, #{sym_timeout_ms => 3600000, sym_max_heap_words => 4000000000, max_expr_nodes => 500000000,
                   exact_constants => true}).

-type integral() :: pq_ibp_family:integral().
-type ratfun() :: {map(), map()}.

%%% Reduction

%% Opts: seeds => auto (default) | {dots, D, S} | {R, S} (Kira's r and s) | [Seed]; order => numerators (default) |
%% dots (which integrals become masters, weight/2); sector_maps => true (default) | false (with given seeds: whether
%% the sectors that map onto another one are expanded into it instead of contributing identities); masters =>
%% [Integral] (preferred masters: each the simplest of its sector, so a master wherever the identities allow; see
%% choose_masters/3); and pq_recon's options (executor, workers, guess).
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
    {Red, #{seeds => SeedSpec, sector_maps => maps:get(sector_maps, Sys, maps:get(sector_maps, Opts, true)),
            equations => NEqs, integrals => NInts, needed => length(A), masters => length(Masters)}}.

%% The reduction of a combination sum c_T T (c_T expressions in the family's variables) as one result: {ok,
%% [{Master, {N, D}}], Info}. One extra unknown x and one equation x - sum c_T T = sum c_T M (the targets that are
%% masters themselves, on the right-hand side) are added to the targets' system, and only x is reconstructed --
%% instead of every target's coefficient of every master (reduce/3), when only their sum is wanted.
-spec reduce_sum(pq_ibp_family:family(), [{integral(), term()}], map()) -> {ok, [{integral(), ratfun()}], map()} | {error, term()}.
reduce_sum(Fam, Comb, Opts) ->
    case reduce_sums(Fam, [Comb], Opts) of
        {ok, [R], Info} -> {ok, R, Info};
        Error -> Error
    end.

%% Several combinations at once, on one system: one extra unknown and equation each, {ok, [[{Master, {N, D}}]],
%% Info} in the combinations' order (a sum split by its other factors -- colour, say -- needs no second system)
-spec reduce_sums(pq_ibp_family:family(), [[{integral(), term()}]], map()) ->
          {ok, [[{integral(), ratfun()}]], map()} | {error, term()}.
reduce_sums(Fam, Combs, Opts) ->
    %% each integral under its canonical name, zero ones dropped, equal ones merged
    Merge = fun(Comb) -> lists:foldl(fun({T, C}, Acc) ->
                                             I = pq_ibp_family:canon(Fam, T),
                                             case pq_ibp_family:zero(Fam, I) of
                                                 true -> Acc;
                                                 false -> maps:update_with(I, fun(C0) -> {add, [C0, C]} end, C, Acc)
                                             end
                                     end, #{}, Comb) end,
    Merged = [Merge(C) || C <- Combs],
    case lists:usort(lists:append([maps:keys(M) || M <- Merged])) of
        [] -> {ok, [[] || _ <- Combs], #{}};
        Live ->
            try
                {SeedSpec, Sys} = case maps:get(seeds, Opts, auto) of
                                      auto -> auto_seeds(Fam, Live, Opts);
                                      Spec -> {Spec, system(Fam, seeds(Fam, Spec, Live), Live, Opts)}
                                  end,
                #{a := A, bs := Bs0, masters := Masters0, targets := TIs, size := {NEqs, NInts}} = Sys,
                N = length(A),
                %% only the masters that can have a coefficient (known from the seed search's check), renumbered
                MasterTargets = [T || {T, {master, _}} <- TIs],
                {Masters, Bs} = case Sys of
                                    #{nonzero := NZ} ->
                                        Keep = [K || {K, {M, _}} <- lists:enumerate(0, Masters0),
                                                     lists:member(M, NZ) orelse lists:member(M, MasterTargets)],
                                        New = maps:from_list(lists:zip(Keep, lists:seq(0, length(Keep) - 1))),
                                        {[lists:nth(K + 1, Masters0) || K <- Keep],
                                         [maps:from_list([{maps:get(C, New), E} || C := E <- R, is_map_key(C, New)]) || R <- Bs0]};
                                    _ -> {Masters0, Bs0}
                                end,
                MPos0 = maps:from_list(lists:zip([M || {M, _} <- Masters], lists:seq(0, length(Masters) - 1))),
                %% a target that is a master but no master of the system (its system has none): the masters extended
                Extra = [T || T <- MasterTargets, not is_map_key(T, MPos0)],
                Ms = [M || {M, _} <- Masters] ++ Extra,
                MPos = maps:from_list(lists:zip(Ms, lists:seq(0, length(Ms) - 1))),
                K = length(Combs),
                Rows = [maps:from_list([{N + J - 1, 1} | [{TI, {mul, [-1, maps:get(T, Mg)]}} || {T, TI} <- TIs, is_integer(TI), is_map_key(T, Mg)]])
                        || {J, Mg} <- lists:enumerate(Merged)],
                Rhss = [maps:from_list([{maps:get(T, MPos), maps:get(T, Mg)} || T <- MasterTargets, is_map_key(T, Mg)]) || Mg <- Merged],
                Vars = pq_ibp_family:vars(Fam),
                Wanted = lists:seq(N, N + K - 1),
                case pq_recon:linear([maps:map(fun(_, E) -> E end, R) || R <- A] ++ Rows, Bs ++ Rhss, Vars,
                                     maps:merge(?LIMITS, Opts#{components => Wanted, rhs => length(Ms)})) of
                    {ok, Per} ->
                        Msa = list_to_tuple(Ms),
                        PerMap = maps:from_list(Per),
                        {ok, [[{element(C + 1, Msa), ND} || {C, ND = {Nu, _}} <- maps:get(W, PerMap), map_size(Nu) > 0] || W <- Wanted],
                         #{seeds => SeedSpec, sector_maps => maps:get(sector_maps, Sys, maps:get(sector_maps, Opts, true)),
                           equations => NEqs, integrals => NInts, needed => N, masters => length(Ms)}};
                    {error, Why} -> {error, Why}
                end
            catch
                throw:{pq_ibp, Why2} -> {error, Why2}
            end
    end.

%%% Seeds

%% the seed integrals: {dots, D, S}: in every non-zero sector, at most D powers above 1 and S numerator powers in
%% total; {sectors, D, S, [{Sector, Ds, Ss}]}: {dots, D, S} with the budget of some sectors raised (auto_seeds/3);
%% {R, S}: Kira's (the positive indices summing to at most R); a list: as given
-spec seeds(pq_ibp_family:family(), term(), [integral()]) -> [integral()].
seeds(Fam, {dots, D, S}, _Targets) -> pq_ibp_family:seeds(Fam, {dots, D}, S);
seeds(Fam, {sectors, D, S, Raised}, _Targets) ->
    lists:usort(pq_ibp_family:seeds(Fam, {dots, D}, S) ++ lists:append([pq_ibp_family:sector_seeds(Fam, Sec, Ds, Ss) || {Sec, Ds, Ss} <- Raised]));
seeds(Fam, {R, S}, _Targets) when is_integer(R) -> pq_ibp_family:seeds(Fam, R, S);
seeds(_Fam, List, _Targets) when is_list(List) -> List.

%% The smallest seeds that reduce every target completely, without learning: {dots, D, S} from the targets' own
%% dots and numerator powers upward, the cheaper first. Complete: every master with a non-zero coefficient in a
%% target (one numeric solve at a random point mod ?P) is itself a seed, its identities generated; a master at the
%% edge of the seed range, whose identities were never written, means the seeds are too small. {Spec, System}
%% (the system as system/4 gives it)
%%
%% First, sector by sector: from the targets' own budget, each master that is not a seed (its identities never
%% written) raises the budget of its own sector just enough to make it one, and the system is built again, until
%% every master is a seed. Most sectors keep the small budget. If that stalls (a master no raise makes a seed), the
%% global budgets below, as before.
-spec auto_seeds(pq_ibp_family:family(), [integral()], map()) -> {term(), map()}.
auto_seeds(Fam, Targets, Opts) ->
    Dt = lists:max([lists:sum([X - 1 || X <- tuple_to_list(T), X > 1]) || T <- Targets]),
    St = lists:max([lists:sum([-X || X <- tuple_to_list(T), X < 0]) || T <- Targets]),
    case maps:get(sector_seeds, Opts, true) of
        false -> global_seeds(Fam, Targets, Opts, Dt, St);
        true ->
            %% the global budgets first, while their seeds stay within three times the smallest budget's; a larger
            %% one hands over to the sector-by-sector search, which raises only the sectors that need it (four
            %% loops: R1), and if that stalls, the global budgets without the limit
            Cap = 3 * length(pq_ibp_family:seeds(Fam, {dots, Dt}, St)),
            case catch global_seeds(Fam, Targets, Opts#{seed_cap => Cap}, Dt, St) of
                {{dots, _, _}, #{} = Sys} = Found when is_map_key(a, Sys) -> Found;
                _ ->
                    case sector_seeds(Fam, Targets, Opts, Dt, St) of
                        {ok, Spec, Sys} -> {Spec, Sys};
                        _ -> global_seeds(Fam, Targets, Opts, Dt, St)
                    end
            end
    end.

%% the sector-by-sector search, with the sector maps (the smallest systems), then without
sector_seeds(Fam, Targets, Opts, Dt, St) ->
    lists:foldl(fun(_, {ok, _, _} = Done) -> Done;
                   (Mapped, _) -> raise(Fam, Targets, Opts#{sector_maps => Mapped}, Dt, St, [], 0)
                end, none, case maps:find(sector_maps, Opts) of {ok, M} -> [M]; error -> [true, false] end).

raise(_, _, _, _, _, _, 12) -> stalled;
raise(Fam, Targets, Opts, Dt, St, Raised, Round) ->
    Spec = {sectors, Dt, St, Raised},
    Seeds = seeds(Fam, Spec, Targets),
    try system(Fam, Seeds, Targets, Opts) of
        Sys ->
            SeedSet = sets:from_list([pq_ibp_family:canon(Fam, X) || X <- Seeds], [{version, 2}]),
            NZ = nonzero_masters(Fam, Sys),
            case [M || M <- NZ, not sets:is_element(M, SeedSet)] of
                [] -> {ok, Spec, Sys#{sector_maps => maps:get(sector_maps, Opts), nonzero => NZ}};
                Missing ->
                    %% each missing master's sector, its budget raised to the master's own dots and numerators
                    Need = lists:foldl(fun(M, Acc) ->
                                               L = tuple_to_list(M),
                                               Sec = [J || {J, X} <- lists:enumerate(L), X > 0],
                                               Dm = lists:sum([X - 1 || X <- L, X > 1]),
                                               Sm = lists:sum([-X || X <- L, X < 0]),
                                               {D0, S0} = maps:get(Sec, Acc, {Dt, St}),
                                               Acc#{Sec => {max(D0, Dm), max(S0, Sm)}}
                                       end, maps:from_list([{Sec, {D, S}} || {Sec, D, S} <- Raised]), Missing),
                    Raised2 = lists:sort([{Sec, D, S} || Sec := {D, S} <- Need]),
                    case Raised2 =:= lists:sort(Raised) of
                        true -> stalled;
                        false -> raise(Fam, Targets, Opts, Dt, St, Raised2, Round + 1)
                    end
            end
    catch
        throw:{pq_ibp, {targets_not_in_system, _}} -> stalled
    end.

global_seeds(Fam, Targets, Opts, Dt, St) ->
    Candidates = lists:sort(fun({D1, S1}, {D2, S2}) -> {D1 + S1, D1} =< {D2 + S2, D2} end,
                            [{D, S} || D <- lists:seq(Dt, Dt + 2), S <- lists:seq(St, St + 2)]),
    %% for each budget: with the sector maps first (the smallest system), then without (a mapped sector's integrals
    %% with dots expand into dotted integrals that the representative's identities may not reduce)
    auto_seeds(Fam, Targets, Opts, [{D, S, M} || {D, S} <- Candidates, M <- [true, false]]).

auto_seeds(_Fam, _Targets, _Opts, []) -> throw({pq_ibp, no_seeds_reduce_the_targets});
auto_seeds(Fam, Targets, Opts, [{D, S, Mapped} | More]) ->
    Seeds = pq_ibp_family:seeds(Fam, {dots, D}, S),
    %% seed_cap (from auto_seeds/3): the budgets come in growing order, so the first too large ends the search
    case length(Seeds) > maps:get(seed_cap, Opts, infinity) of
        true -> throw({pq_ibp, seed_cap});
        false -> ok
    end,
    Result = try system(Fam, Seeds, Targets, maps:remove(seed_cap, Opts#{sector_maps => Mapped})) of
                 Sys ->
                     SeedSet = sets:from_list([pq_ibp_family:canon(Fam, X) || X <- Seeds], [{version, 2}]),
                     NZ = nonzero_masters(Fam, Sys),
                     case lists:all(fun(M) -> sets:is_element(M, SeedSet) end, NZ) of
                         true -> {ok, Sys#{nonzero => NZ}};
                         false -> incomplete
                     end
             catch
                 throw:{pq_ibp, {targets_not_in_system, _}} -> incomplete
             end,
    case Result of
        {ok, Sys2} -> {{dots, D, S}, Sys2#{sector_maps => Mapped}};
        incomplete -> auto_seeds(Fam, Targets, Opts, More)
    end.

%% the masters with a non-zero coefficient in any target, from one numeric solve at a random point mod ?P: of a random
%% combination x = sum r_T T of the targets (one extra unknown and equation), whose masters are those of every
%% target but with probability about (number of masters) / ?P -- one solve instead of one per target
nonzero_masters(Fam, #{a := A, bs := Bs, masters := Masters, targets := TIs}) ->
    Vars = pq_ibp_family:vars(Fam),
    case [TI || {_, TI} <- TIs, is_integer(TI)] of
        [] -> [];
        Cols ->
            N = length(A),
            {Rs, _} = lists:mapfoldl(fun(_, S) -> rand:uniform_s(?P - 1, S) end, rand:seed_s(exsss, {2, 7, 1}), Cols),
            Row = maps:from_list([{N, 1} | [{C, -R} || {C, R} <- lists:zip(Cols, Rs)]]),
            Sp = fun(R) -> lists:sort(maps:to_list(R)) end,
            SA = [Sp(R) || R <- A] ++ [Sp(Row)], SB = [Sp(R) || R <- Bs] ++ [[]],
            {ok, Prog} = pq_symbolic:compile([E || R <- SA ++ SB, {_, E} <- R], Vars, ?COMPILE),
            {ok, {ff, Ref, _, _}} = pq_ff:prepare(Prog),
            PointMap = point(Vars, {3, 1, 4}),
            Point = [maps:get(V, PointMap) || V <- Vars],
            {ok, Lin} = pq_nif:ff_linsys(Ref, N + 1, length(Masters), N, 0, [[C || {C, _} <- R] || R <- SA],
                                         [[C || {C, _} <- R] || R <- SB], []),
            {ok, [Vals]} = pq_nif:ff_linsys_all(Lin, ?P, Point),
            lists:usort([M || {{M, _}, V} <- lists:zip(Masters, Vals), V =/= 0])
    end.

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
system(Fam0, Seeds0, Targets, Opts) ->
    Vars = pq_ibp_family:vars(Fam0),
    %% sector_maps => true (default): the sectors that map onto another one contribute no identities, their integrals
    %% expanded into the representative's; false: every sector's identities, the scalar symmetries only
    {Fam, Seeds} = case maps:get(sector_maps, Opts, true) of
                       true -> {Fam0, [A || A <- Seeds0, not pq_ibp_family:mapped(Fam0, A)]};
                       false -> {pq_ibp_family:without_sector_maps(Fam0), Seeds0}
                   end,
    %% the coefficients as exact polynomials (Exponents => {N, D}): evaluated mod ?P directly, and only the needed
    %% rows turned into ASTs
    %% in parallel chunks of the seeds, concatenated in order (the same equations in the same order as one process),
    %% the family's caches computed once here and handed to the chunks
    Warm = pq_ibp_family:warm(Fam, Seeds),
    Point = point(Vars, {1, 1, 1}),
    PointList = [maps:get(V, Point) || V <- Vars],
    %% each equation with its coefficients' values mod ?P, computed where the equation is made (one copy back)
    Valued = fun(E) -> {E, [{I, polymod(C, PointList)} || I := C <- E]} end,
    Pairs = [P || {E, _} = P <- lists:append(pmap(fun(Chunk) ->
                                                         ok = pq_ibp_family:load(Warm),
                                                         [Valued(E) || E <- pq_ibp_family:identities_poly(Fam, Chunk)]
                                                 end,
                                                 chunks(Seeds, max(64, length(Seeds) div (4 * erlang:system_info(schedulers_online)) + 1)))),
                  map_size(E) > 0]
          %% a target in a sector that maps onto another one: defined by its expansion (its own sector has no
          %% identities)
          ++ [Valued(maps:fold(fun(J, C, M) -> M#{J => pq_ibp_family:qpneg(C)} end, #{T => pq_ibp_family:qpone(Fam)},
                               maps:from_list(pq_ibp_family:expand(Fam, T))))
              || T <- Targets, pq_ibp_family:mapped(Fam, T)],
    Eqs = [E || {E, _} <- Pairs],
    Ints = lists:usort(lists:append([maps:keys(E) || E <- Eqs])),
    case [T || T <- Targets, not lists:member(T, Ints)] of
        [] -> ok;
        Missing -> throw({pq_ibp, {targets_not_in_system, Missing}})
    end,
    Order = maps:get(order, Opts, numerators),
    Prefer = sets:from_list([pq_ibp_family:canon(Fam, M) || M <- maps:get(masters, Opts, [])], [{version, 2}]),
    Weight = fun(I) ->
                     case sets:is_element(I, Prefer) of
                         true -> {length([x || X <- tuple_to_list(I), X > 0]), -1, 0, I};
                         false -> weight(I, Order)
                     end
             end,
    Cols = [I || {_, I} <- lists:reverse(lists:sort([{Weight(I), I} || I <- Ints]))],
    Index = maps:from_list(lists:zip(Cols, lists:seq(0, length(Cols) - 1))),
    Num = [lists:sort([{maps:get(I, Index), V} || {I, V} <- Vals]) || {_, Vals} <- Pairs],
    TCols = [maps:get(T, Index) || T <- Targets],
    {PR, PC, Needed} = pq_nif:ff_pivots(?P, Num, length(Cols), TCols),
    PivotCol = maps:from_list(lists:zip(PR, PC)),
    EqArr = array:from_list(Eqs),
    ColArr = array:from_list(Cols),
    Rows = [maps:map(fun(_, C) -> poly_ast(C, Vars) end, array:get(R, EqArr)) || R <- Needed],
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

%% a polynomial coefficient's value mod ?P at the point (its coordinates in the variables' order)
polymod(Poly, Xs) ->
    maps:fold(fun(Exps, {N, D}, Acc) ->
                      Mono = lists:foldl(fun({X, E}, M) -> M * pow_mod(X, E) rem ?P end, 1, lists:zip(Xs, Exps)),
                      C = (N rem ?P + ?P) rem ?P * pq_ff:inv(D rem ?P, ?P) rem ?P,
                      (Acc + C * Mono) rem ?P
              end, 0, Poly).

pow_mod(_, 0) -> 1;
pow_mod(X, 1) -> X;
pow_mod(X, E) -> X * pow_mod(X, E - 1) rem ?P.

%% a polynomial coefficient as an AST: a sum of rational multiples of monomials
poly_ast(Poly, Vars) ->
    Terms = [{mul, [case D of 1 -> N; _ -> {rat, N, D} end | [{pow, {sym, V}, E} || {V, E} <- lists:zip(Vars, Exps), E > 0]]}
             || {Exps, {N, D}} <- lists:sort(maps:to_list(Poly))],
    case Terms of [T] -> T; _ -> {add, Terms} end.

%%% Choosing the masters

%% A greedy search for a basis that makes the targets' coefficients small: each master in turn replaced by a
%% neighbour in its sector (one more dot on a line, or one more numerator power), the targets reduced again, and the
%% change kept when the coefficients' total size (size_of/2: additions after pq_simplify) falls. Every candidate is
%% an exact reduction, so the basis found is as valid as the default one; only its size is searched.
%% Opts as for reduce/3, plus rounds => R (default 2: passes over the masters). {ok, Masters, #{size_in, size_out,
%% tried, reductions}} | {error, _}.
-spec choose_masters(pq_ibp_family:family(), [integral()], map()) -> {ok, [integral()], map()} | {error, term()}.
choose_masters(Fam, Targets, Opts) ->
    Vars = pq_ibp_family:vars(Fam),
    case reduce(Fam, Targets, Opts) of
        {error, _} = E -> E;
        {ok, R0, Info} ->
            %% every trial on the default reduction's system (its seeds and sector maps): a candidate the system cannot
            %% make a master just leaves the basis as it is, and each trial costs one reduction
            Opts1 = Opts#{seeds => maps:get(seeds, Info, maps:get(seeds, Opts, auto)),
                          sector_maps => maps:get(sector_maps, Info, maps:get(sector_maps, Opts, true))},
            B0 = basis(R0),
            S0 = size_of(R0, Vars),
            %% a target is never a candidate (a basis of the targets themselves is trivially small)
            Live = [pq_ibp_family:canon(Fam, T) || T <- Targets],
            {B, S, Tried} = improve(Fam, Targets, Live, Opts1, Vars, B0, S0, maps:get(rounds, Opts, 2), 0),
            {ok, B, #{size_in => S0, size_out => S, tried => Tried, seeds => maps:get(seeds, Opts1),
                      sector_maps => maps:get(sector_maps, Opts1)}}
    end.

basis(Red) -> lists:usort([M || _ := Cs <- Red, {M, _} <- Cs]).

%% the coefficients' total size: additions after pq_simplify:rational/4
-spec size_of(#{integral() => [{integral(), ratfun()}]}, [atom()]) -> non_neg_integer().
size_of(Red, Vars) ->
    lists:sum([pq_simplify:adds(element(1, pq_simplify:rational(N, D, Vars, #{}))) || _ := Cs <- Red, {_, {N, D}} <- Cs]).

improve(_, _, _, _, _, B, S, 0, Tried) -> {B, S, Tried};
improve(Fam, Targets, Live, Opts, Vars, B0, S0, Rounds, Tried0) ->
    {B, S, Tried} =
        lists:foldl(fun(M, {B1, S1, T1}) ->
                            case lists:member(M, B1) of
                                false -> {B1, S1, T1};
                                true ->
                                    lists:foldl(fun(C, {B2, S2, T2}) ->
                                                        Try = lists:usort([C | B2 -- [M]]),
                                                        case catch reduce(Fam, Targets, Opts#{masters => Try}) of
                                                            {ok, R, _} ->
                                                                Got = basis(R),
                                                                %% a target equal to a master (by a symmetry) is no basis
                                                                Trivial = lists:any(fun(T) -> case maps:get(T, R) of
                                                                                                  [{_, {Nu, De}}] -> Nu =:= De;
                                                                                                  _ -> false
                                                                                              end end, Targets),
                                                                case lists:member(C, Got) andalso not Trivial andalso size_of(R, Vars) of
                                                                    S3 when is_integer(S3), S3 < S2 -> {Got, S3, T2 + 1};
                                                                    _ -> {B2, S2, T2 + 1}
                                                                end;
                                                            _ -> {B2, S2, T2 + 1}
                                                        end
                                                end, {B1, S1, T1}, neighbours(Fam, M) -- (B1 ++ Live))
                            end
                    end, {B0, S0, Tried0}, B0),
    case S < S0 of
        true -> improve(Fam, Targets, Live, Opts, Vars, B, S, Rounds - 1, Tried);
        false -> {B, S, Tried}
    end.

%% the integrals of M's sector one step away: a dot on one of its lines, or a numerator power on one of the others
neighbours(Fam, M) ->
    L = tuple_to_list(M),
    lists:usort([C || J <- lists:seq(1, length(L)),
                      C <- [pq_ibp_family:canon(Fam, list_to_tuple(setnth(J, L, case lists:nth(J, L) > 0 of
                                                                                       true -> lists:nth(J, L) + 1;
                                                                                       false -> lists:nth(J, L) - 1
                                                                                   end)))],
                      not pq_ibp_family:zero(Fam, C)]).

setnth(1, [_ | T], V) -> [V | T];
setnth(N, [H | T], V) -> [H | setnth(N - 1, T, V)].

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

%% F over the list, each element in a process of its own, the results in the list's order (an exception in one is
%% raised in the caller, with its class); each worker sends its outcome, the monitor catches a worker that dies
pmap(F, L) ->
    Self = self(),
    Jobs = [begin
                Tag = make_ref(),
                {Pid, Mon} = spawn_monitor(fun() -> Self ! {Tag, try {ok, F(X)} catch C:R:St -> {raise, C, R, St} end} end),
                {Tag, Pid, Mon}
            end || X <- L],
    [receive
         {Tag, Outcome} ->
             erlang:demonitor(Mon, [flush]),
             case Outcome of
                 {ok, V} -> V;
                 {raise, C, R, St} -> erlang:raise(C, R, St)
             end;
         {'DOWN', Mon, process, Pid, Why} -> error({pmap, Why})
     end || {Tag, Pid, Mon} <- Jobs].

chunks([], _) -> [];
chunks(L, N) when length(L) =< N -> [L];
chunks(L, N) -> {A, B} = lists:split(N, L), [A | chunks(B, N)].

