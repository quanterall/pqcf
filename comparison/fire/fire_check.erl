%% PQCF's reduction of the double box against FIRE's (decision record R1): both give target = sum c_j M_j, but the
%% programs may pick different masters, and FIRE's example doesn't apply every sector symmetry (12 masters, PQCF 8).
%% The comparison is in FIRE's basis: each PQCF master is one of FIRE's (up to a symmetry), or FIRE reduces it
%% (a tables file of its own); FIRE masters equal by a symmetry (pq_ibp_family:canon/2) are merged. At a random point
%% mod p, sum_j c_j [M_j in FIRE's masters] must equal FIRE's coefficients, master by master.
%%
%%   fire_check:run("<dir>/ibp_doublebox_7_2.answer", "<FIRE's doublebox1.tables>", [{Master, "<FIRE's tables>"}])
-module(fire_check).
-export([run/3]).

-define(P, 1000000007).

run(AnswerFile, TablesFile, Reduced) ->
    {ok, _} = application:ensure_all_started(pq),
    Fam = pq_ibp_family:example(doublebox),
    Canon = fun(F) -> pq_ibp_family:canon(Fam, F) end,
    Point = #{d => 1234577, s => 7654337, t => 3456791},
    Ours = [{M, value(Ast, Point)} || {M, Ast} <- answer(AnswerFile)],
    Fire = [{F, value(Ast, Point)} || {F, Ast} <- fire_target(TablesFile)],
    FireClasses = lists:usort([Canon(F) || {F, _} <- Fire]),
    %% PQCF's answer in FIRE's masters (merged by the symmetries)
    InFire = fun(M) ->
                     case lists:member(Canon(M), FireClasses) of
                         true -> [{Canon(M), 1}];
                         false ->
                             {M, File} = lists:keyfind(M, 1, Reduced),
                             [{Canon(F), value(Ast, Point)} || {F, Ast} <- fire_target(File)]
                     end
             end,
    Add = fun(K, V, Acc) -> Acc#{K => (maps:get(K, Acc, 0) + V) rem ?P} end,
    Lhs = lists:foldl(fun({M, C}, Acc) -> lists:foldl(fun({F, H}, A) -> Add(F, C * H rem ?P, A) end, Acc, InFire(M)) end, #{}, Ours),
    Rhs = lists:foldl(fun({F, V}, Acc) -> Add(Canon(F), V, Acc) end, #{}, Fire),
    Keys = lists:usort(maps:keys(Lhs) ++ maps:keys(Rhs)),
    Bad = [{K, maps:get(K, Lhs, 0), maps:get(K, Rhs, 0)} || K <- Keys, maps:get(K, Lhs, 0) =/= maps:get(K, Rhs, 0)],
    io:format("FIRE: ~b masters (~b up to symmetries); PQCF: ~b; compared at a random point mod p: ~s~n",
              [length(Fire), length(FireClasses), length(Ours),
               case Bad of [] -> "all equal"; _ -> io_lib:format("~b differ: ~p", [length(Bad), Bad]) end]),
    Bad =:= [].

%% an exact rational function evaluated at the integer point, mod p
value(Ast, Point) ->
    {ok, V} = pq_symbolic:evaluate(Ast, Point, #{exact => true}),
    case V of
        N when is_integer(N) -> (N rem ?P + ?P) rem ?P;
        {rat, N, D} -> (N rem ?P + ?P) rem ?P * pq_ff:inv((D rem ?P + ?P) rem ?P, ?P) rem ?P
    end.

%% ibp.erl's answer file: lines "Integral Numerator / Denominator"
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

%% FIRE's tables: the first entry is the target's reduction, a list of {Code, "Coefficient"}; a code is 2 followed by
%% the indices, each as 127 + index in three digits
fire_target(File) ->
    {ok, Bin} = file:read_file(File),
    [First | _] = binary:split(Bin, <<"}\n        },">>),
    {match, Ms} = re:run(First, <<"\\{(2[0-9]+),\"([^\"]+)\"\\}">>, [global, {capture, all_but_first, binary}]),
    [{decode(Code), begin {ok, Ast} = pq_symbolic:parse(Coef), Ast end} || [Code, Coef] <- Ms].

decode(<<"2", Rest/binary>>) ->
    list_to_tuple([binary_to_integer(X) - 127 || <<X:3/binary>> <= Rest]).
