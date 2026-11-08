%% The reduction benchmark suite (decision record R1, "A benchmark suite"): fixed Feynman-integral reductions, the
%% ones FIRE ships as examples run by both programs, timed and checked equal. bench/ibp.sh runs it; this module is
%% PQCF's side.
%%
%%   pqcf(Dir, Label): every case reduced (pq_ibp:reduce/3, automatic seeds), timed in this node, its answer and
%%     masters written to Dir (<case>.answer, and for FIRE's cases <case>.masters.m in FIRE's format); one JSON
%%     object per case in Dir/pqcf-<Label>.json
%%   check(Dir): FIRE's cases compared with FIRE's tables in Dir (<case>.fire.tables, <case>.masters.tables), in
%%     FIRE's basis; Dir/check.json
%%   fire_cases(): the FIRE examples, as lines "Case Example Problem TargetFile" for bench/ibp.sh
-module(ibp_suite).
-export([cases/0, pqcf/2, check/1, fire_cases/1]).


%% {Case, Family, Target, Fire}: Fire = {Example, Problem} (FIRE's examples/<Example>.start, problem number) or none
cases() ->
    [{fire_box, fire_box, {3, 3, 3, 3}, {"box", 1}},
     {fire_v2, v2, {-1, 1, 1, 1, 1, 1, 1}, {"v2", 2}},
     {doublebox_1, doublebox, {1, 1, 1, 1, 1, 1, 1, -1, -1}, {"doublebox", 1}},
     {doublebox_2, doublebox, {1, 1, 1, 1, 1, 1, 1, -2, -1}, {"doublebox", 1}},
     {doublebox_3, doublebox, {1, 1, 1, 1, 1, 1, 2, -2, -2}, {"doublebox", 1}},
     {kite, kite, {2, 1, 1, 1, 1}, none},
     {massive_tadpole, tadpole, {4}, none},
     {banana3, banana3, {1, 1, 1, 2, -1, 0, 0, 0, 0}, none}].

%% FIRE's examples as pq_ibp_family families (from FIRE's notebooks: FIRE's propagators are -q^2, sign -1); the
%% others are PQCF's own
family(fire_box) ->
    %% one-loop massless box, S = T = 1: p1.p2 = -S/2, p2.p4 = -T/2, p1.p4 = (S + T)/2
    pq_ibp_family:new(#{loops => [k], ext => [p1, p2, p4], vars => [d],
                        props => [{-1, [{1, k}]}, {-1, [{1, k}, {1, p1}]}, {-1, [{1, k}, {1, p1}, {1, p2}]},
                                  {-1, [{1, k}, {1, p1}, {1, p2}, {1, p4}]}],
                        kin => #{{p1, p2} => {rat, -1, 2}, {p2, p4} => {rat, -1, 2}, {p1, p4} => 1}});
family(v2) ->
    %% two-loop massless vertex: p^2 = q^2 = 0, p.q = -1/2
    Q = fun(Ms) -> {-1, Ms} end,
    pq_ibp_family:new(#{loops => [l, r], ext => [p, q], vars => [d],
                        props => [Q([{1, l}, {-1, r}]), Q([{1, l}]), Q([{1, r}]), Q([{1, p}, {-1, l}]), Q([{1, q}, {-1, r}]),
                                  Q([{1, p}, {-1, l}, {1, r}]), Q([{1, q}, {-1, r}, {1, l}])],
                        kin => #{{p, q} => {rat, -1, 2}}});
family(tadpole) ->
    pq_ibp_family:new(#{loops => [k], ext => [], vars => [d, msq], props => [{1, [{1, k}], {sym, msq}}]});
family(banana3) ->
    pq_ibp_family:new(#{loops => [k1, k2, k3], ext => [p], vars => [d, s], kin => #{{p, p} => {sym, s}},
                        top => [1, 2, 3, 4],
                        props => [{1, [{1, k1}]}, {1, [{1, k2}]}, {1, [{1, k3}]},
                                  {1, [{1, p}, {-1, k1}, {-1, k2}, {-1, k3}]},
                                  {1, [{1, k1}, {-1, k2}]}, {1, [{1, k1}, {-1, k3}]}, {1, [{1, k2}, {-1, k3}]},
                                  {1, [{1, k1}, {-1, p}]}, {1, [{1, k2}, {-1, p}]}]});
family(Name) -> pq_ibp_family:example(Name).

fire_cases(Dir) ->
    Lines = [io_lib:format("~s ~s ~b ~s~n", [C, Ex, Pr, string:join([integer_to_list(X) || X <- tuple_to_list(T)], ",")])
             || {C, _, T, {Ex, Pr}} <- cases()],
    ok = file:write_file(filename:join(Dir, "fire_cases.txt"), Lines).

pqcf(Dir, Label) ->
    {ok, _} = application:ensure_all_started(pq),
    Rows = [run_case(Dir, C) || C <- cases()],
    ok = file:write_file(filename:join(Dir, "pqcf-" ++ Label ++ ".json"), json:encode(Rows)),
    [io:format("~-16s ~8.3f s  ~6b equations  ~5b needed  ~2b masters  ~4b additions~n",
               [maps:get(<<"case">>, R), maps:get(<<"seconds">>, R), maps:get(<<"equations">>, R),
                maps:get(<<"needed">>, R), maps:get(<<"masters">>, R), maps:get(<<"additions">>, R)]) || R <- Rows],
    ok.

run_case(Dir, {Case, F, T, Fire}) ->
    Fam = family(F),
    Vars = pq_ibp_family:vars(Fam),
    %% each case from cold caches: a fresh process, the family's caches its own
    Self = self(),
    Pid = spawn_link(fun() -> Self ! {self(), timer:tc(pq_ibp, reduce, [Fam, [T], #{}])} end),
    {Us, {ok, #{T := Coeffs}, Info}} = receive {Pid, Res} -> Res end,
    P = fun(E) -> binary_to_list(pq_symbolic:pretty(E)) end,
    ok = file:write_file(filename:join(Dir, atom_to_list(Case) ++ ".answer"),
                         [[io_lib:format("~w", [M]), " ", P(pq_sym_norm:to_ast(N)), " / ", P(pq_sym_norm:to_ast(D)), "\n"]
                          || {M, {N, D}} <- Coeffs]),
    case Fire of
        {_, Pr} ->
            ok = file:write_file(filename:join(Dir, atom_to_list(Case) ++ ".masters.m"),
                                 ["{", lists:join(",\n ", [io_lib:format("{~b, {~s}}", [Pr, string:join([integer_to_list(X) || X <- tuple_to_list(M)], ", ")])
                                                          || {M, _} <- Coeffs]), "}\n"]);
        none -> ok
    end,
    Adds = pq_ibp:size_of(#{T => Coeffs}, Vars),
    #{<<"case">> => atom_to_binary(Case), <<"family">> => atom_to_binary(F),
      <<"target">> => list_to_binary(io_lib:format("~w", [T])), <<"seconds">> => Us / 1.0e6,
      <<"seeds">> => list_to_binary(io_lib:format("~w", [maps:get(seeds, Info, none)])),
      <<"equations">> => maps:get(equations, Info, 0), <<"needed">> => maps:get(needed, Info, 0),
      <<"masters">> => length(Coeffs), <<"additions">> => Adds}.

check(Dir) ->
    {ok, _} = application:ensure_all_started(pq),
    Rows = [begin
                Fam = family(F),
                Ours = answer(filename:join(Dir, atom_to_list(Case) ++ ".answer")),
                FireT = fire_check:tables(filename:join(Dir, atom_to_list(Case) ++ ".fire.tables")),
                MT = case filelib:is_regular(filename:join(Dir, atom_to_list(Case) ++ ".masters.tables")) of
                         true -> fire_check:tables(filename:join(Dir, atom_to_list(Case) ++ ".masters.tables"));
                         false -> #{}
                     end,
                {Equal, NF, NC} = fire_check:compare(Fam, Ours, T, FireT, MT),
                io:format("~-16s FIRE ~b masters (~b up to symmetries), PQCF ~b: ~s~n",
                          [Case, NF, NC, length(Ours), case Equal of true -> "equal"; false -> "DIFFERENT" end]),
                #{<<"case">> => atom_to_binary(Case), <<"equal">> => Equal, <<"fire_masters">> => NF,
                  <<"fire_classes">> => NC}
            end || {Case, F, T, {_, _}} <- cases()],
    ok = file:write_file(filename:join(Dir, "check.json"), json:encode(Rows)),
    case lists:all(fun(#{<<"equal">> := E}) -> E end, Rows) of
        true -> ok;
        false -> halt(1)
    end.

answer(File) ->
    {ok, Bin} = file:read_file(File),
    [begin
         [Int, Rest] = string:split(Line, " "),
         [Num, Den] = string:split(Rest, " / "),
         {ok, Toks, _} = erl_scan:string(binary_to_list(Int) ++ "."),
         {ok, I} = erl_parse:parse_term(Toks),
         {ok, N} = pq_symbolic:parse(Num),
         {ok, D} = pq_symbolic:parse(Den),
         {I, {mul, [N, {pow, D, -1}]}}
     end || Line <- binary:split(Bin, <<"\n">>, [global, trim_all])].
