%% PQCF against FORM's Forcer on FORM's own Forcer checks (comparison/forcer/README.md). An integral is given in
%% Forcer's notation: vertices vx(...) (momenta into the vertex sum to zero, Q the external one), the lines' powers,
%% numerator scalar products and ex(...) insertions (a line's power 1 + epsilon). PQCF builds the family from the
%% graph (a loop basis from momentum conservation, auxiliary invariants to complete the scalar products), reduces
%% with pq_ibp, and writes its masters back in Forcer's notation, for Forcer to give their values. The answers are
%% compared as epsilon expansions (epsexp): PQCF's coefficients times the masters' values against Forcer's own
%% result, and in MS-bar against the numbers FORM's check asserts.
-module(forcer_cmp).
-export([problem/3, reduce/2, masters_form/2, targets_form/1, compare/5, form_exprs/1, msbar/2, asserted/3, master_tables/1, integral_text/2, targets/2]).

%% Problem: #{vx => "vx(Q,p1,..)*vx(..)", lines => N, ex => [Line], terms => [{Coef, Powers, Sps}]}:
%% the integral is sum Coef * prod_i 1/(p_i^2)^Powers_i * prod {A, B} A.B, A and B momenta names (q for Q).
%% Returns #{family, map (the momenta), targets ([{Integral, CoefAst}] summed), ...}.
problem(Vx, N, Opts) ->
    Verts = parse_vx(Vx),
    {Loops, Mom} = solve(Verts, N),
    LineMoms = [maps:get(list_to_atom("p" ++ integer_to_list(I)), Mom) || I <- lists:seq(1, N)],
    %% lines with the same momentum (up to its sign: a self-energy insertion on a line) are one propagator, their
    %% powers added: Rep maps each line to its propagator
    Canon = fun(M) -> Ms = lists:sort([X || X = {C, _} <- M, C =/= 0]),
                      Neg = lists:sort([{-C, A} || {C, A} <- Ms]),
                      min(Ms, Neg) end,
    {Reps, Rep} = lists:foldl(fun({I, M}, {Rs, Map}) ->
                                      K = Canon(M),
                                      case lists:keyfind(K, 1, Rs) of
                                          {K, J} -> {Rs, Map#{I => J}};
                                          false -> J = length(Rs) + 1, {Rs ++ [{K, J}], Map#{I => J}}
                                      end
                              end, {[], #{}}, lists:enumerate(LineMoms)),
    Props = [hd([M || {I, M} <- lists:enumerate(LineMoms), maps:get(I, Rep) =:= J]) || {_, J} <- Reps],
    Aux = complete(Loops, Props),
    Ex = maps:get(ex, Opts, []),
    Fam = pq_ibp_family:new(#{loops => Loops, ext => [q], vars => [d], kin => #{{q, q} => 1},
                              props => [{1, M} || M <- Props] ++ [{1, A} || A <- Aux], top => lists:seq(1, length(Props)),
                              eps => maps:from_list([{maps:get(L, Rep), 1} || L <- Ex])}),
    #{family => Fam, mom => Mom, lines => N, rep => Rep, props => length(Props), aux => Aux, vx => Vx, ex => Ex,
      loops => length(Loops)}.

%% a term's line powers as the family's: lines sharing a propagator add their powers
fold_powers(#{rep := Rep, props := NP}, Powers) ->
    [lists:sum([X || {I, X} <- lists:enumerate(Powers), maps:get(I, Rep) =:= J]) || J <- lists:seq(1, NP)].

%% the targets of a problem: [{Coef, Powers, Sps}] as a combination of the family's integrals
targets(P = #{family := Fam, mom := Mom, aux := Aux}, Terms) ->
    M = fun(Name) -> maps:get(Name, Mom) end,
    %% numerators with rational multipliers (the kinematics are numbers), each integral's multiplier of each distinct
    %% term coefficient summed: Integral => #{Coef => Rational}; then one simplification per integral over its few
    %% distinct coefficients (symbolic merging term by term was quadratic: minutes for a graph of 1,571 terms)
    Comb = lists:foldl(
             fun({Coef, Powers, Sps}, Acc) ->
                     Base = list_to_tuple(fold_powers(P, Powers) ++ [0 || _ <- Aux]),
                     Sp = [{M(A), M(B)} || {A, B} <- Sps],
                     Exp = case pq_ibp_family:numerator_q(Fam, Base, Sp) of
                               error -> [{I, {ast, C}} || {I, C} <- pq_ibp_family:numerator(Fam, Base, Sp)];
                               Qs -> Qs
                           end,
                     lists:foldl(fun({I, R}, A2) ->
                                         maps:update_with(I, fun(Cs) -> maps:update_with(Coef, fun(R0) -> addq(R0, R) end, R, Cs) end,
                                                          #{Coef => R}, A2)
                                 end, Acc, Exp)
             end, #{}, Terms),
    [{I, C} || I := Cs <- Comb,
               {ok, C} <- [pq_symbolic:canonicalise({add, [{mul, [Coef, q_ast(R)]} || Coef := R <- Cs]})], C =/= 0].

addq({ast, A}, {ast, B}) -> {ast, {add, [A, B]}};
addq({ast, A}, B) -> {ast, {add, [A, q_ast(B)]}};
addq(A, {ast, B}) -> {ast, {add, [q_ast(A), B]}};
addq({A, B}, {C, D}) -> qnorm({A * D + C * B, B * D}).
q_ast({ast, A}) -> A;
q_ast({N, 1}) -> N;
q_ast({N, D}) -> {rat, N, D}.

%% PQCF's reduction of Terms: {Seconds, #{Master => CoefAst}, Info} (the combination's coefficients, in d; one
%% reconstruction for the whole sum, pq_ibp:reduce_sum/3)
reduce(P = #{family := Fam}, Terms) ->
    Ts = targets(P, Terms),
    {Us, {ok, Red, Info}} = timer:tc(pq_ibp, reduce_sum, [Fam, Ts, #{}]),
    {Us / 1.0e6, maps:from_list([{M, {mul, [pq_sym_norm:to_ast(Nu), {pow, pq_sym_norm:to_ast(De), -1}]}} || {M, {Nu, De}} <- Red]), Info}.

%% PQCF's masters as a FORM program for Forcer: one expression per master, M1, M2, ...
masters_form(P, Masters) ->
    ["#-\n#include- forcer.h\n",
     [io_lib:format("L M~b = ~s;~n", [K, integral_form(P, I)]) || {K, I} <- lists:enumerate(Masters)],
     "#call Forcer\nB ep;\nP +s;\n.end\n"].

targets_form(Exprs) ->
    ["#-\n#include- forcer.h\n", [io_lib:format("L ~s = ~s;~n", [N, E]) || {N, E} <- Exprs], "#call Forcer\nB ep;\nP +s;\n.end\n"].

%% an integral of the problem in Forcer's notation, as text
integral_text(P, I) -> iolist_to_binary(integral_form(P, I)).

integral_form(#{vx := Vx, lines := N, rep := Rep, props := NP, aux := Aux, ex := Ex}, I) ->
    L = tuple_to_list(I),
    {Ps, As} = lists:split(NP, L),
    %% a propagator's power on the first of its lines, 0 on the others (Forcer contracts those)
    First = maps:from_list([{J, hd([K || K <- lists:seq(1, N), maps:get(K, Rep) =:= J])} || J <- lists:seq(1, NP)]),
    Ls = [case maps:get(maps:get(K, Rep), First) of K -> lists:nth(maps:get(K, Rep), Ps); _ -> 0 end || K <- lists:seq(1, N)],
    %% a line's power: its index (ex(p) adds epsilon to the power written, whatever it is)
    Lines = [io_lib:format("*p~b.p~b^~b", [J, J, -X]) || {J, X} <- lists:enumerate(Ls), X =/= 0],
    Nums = [case X > 0 of
                true -> error({aux_propagator, I});
                false -> io_lib:format("*(~s)^~b", [square(A), -X])
            end || {A, X} <- lists:zip(Aux, As), X =/= 0],
    ExS = case Ex of [] -> ""; _ -> ["*ex(", lists:join(",", ["p" ++ integer_to_list(J) || J <- Ex]), ")"] end,
    [Vx, Lines, Nums, ExS].

%% a momentum's square in FORM: sum c_a c_b a.b
square(Mom) ->
    Name = fun(q) -> "Q"; (A) -> atom_to_list(A) end,
    lists:join("+", [io_lib:format("~b*~s.~s", [Ca * Cb, Name(A), Name(B)]) || {Ca, A} <- Mom, {Cb, B} <- Mom]).

%% FORM's printed expressions in a log: #{Name => Text} ("   Name =" up to ";")
form_exprs(File) ->
    {ok, Bin} = file:read_file(File),
    {match, Ms} = re:run(Bin, "\\n   ([A-Za-z][A-Za-z0-9]*) =\\n(.*?);", [global, dotall, {capture, all_but_first, list}]),
    maps:from_list([{N, E} || [N, E] <- Ms]).

%%% The graph

%% "vx(Q,p1,-p5)*vx(...)" -> [[{1, q}, {1, p1}, {-1, p5}], ...]
parse_vx(S) ->
    {match, Vs} = re:run(S, "vx\\(([^)]*)\\)", [global, {capture, all_but_first, list}]),
    [[case string:trim(T) of
          "-" ++ X -> {-1, name(X)};
          X -> {1, name(X)}
      end || T <- string:split(V, ",", all)] || [V] <- Vs].

name("Q") -> q;
name(X) -> list_to_atom(X).

%% momentum conservation solved: {Loops, #{p_i => [{C, Atom}]}}: the pivot momenta in terms of the free ones
%% (the loop momenta) and q
solve(Verts, N) ->
    Cols = [list_to_atom("p" ++ integer_to_list(I)) || I <- lists:seq(1, N)] ++ [q],
    Rows = [[{lists:sum([C || {C, X} <- V, X =:= Col]), 1} || Col <- Cols] || V <- Verts],
    {Red, Pivots} = rref(Rows, 1, N, []),
    Free = [lists:nth(I, Cols) || I <- lists:seq(1, N), not lists:member(I, Pivots)],
    Dep = maps:from_list(
            [{lists:nth(Pc, Cols),
              [{-int(lists:nth(Fc, Row)), lists:nth(Fc, Cols)} || Fc <- lists:seq(1, N + 1), not lists:member(Fc, Pivots),
                                                                   lists:nth(Fc, Row) =/= {0, 1}]}
             || {Pc, Row} <- lists:zip(Pivots, Red)]),
    Mom = maps:merge(maps:from_list([{F, [{1, F}]} || F <- Free]), Dep),
    {Free, Mom#{q => [{1, q}]}}.

int({X, 1}) -> X.

rref(Rows, C, N, Piv) when C > N -> {[R || R <- Rows, lists:any(fun(X) -> X =/= {0, 1} end, R)], lists:reverse(Piv)};
rref(Rows, C, N, Piv) ->
    Done = length(Piv),
    {Top, Rest} = lists:split(Done, Rows),
    case lists:splitwith(fun(R) -> lists:nth(C, R) =:= {0, 1} end, Rest) of
        {_, []} -> rref(Rows, C + 1, N, Piv);
        {Before, [P | After]} ->
            Inv = qinv(lists:nth(C, P)),
            Pn = [qmul(Inv, X) || X <- P],
            E = fun(R) -> F = lists:nth(C, R), [qsub(X, qmul(F, Y)) || {X, Y} <- lists:zip(R, Pn)] end,
            rref([E(R) || R <- Top] ++ [Pn] ++ [E(R) || R <- Before ++ After], C + 1, N, [C | Piv])
    end.

%% auxiliary invariants completing the lines to a basis of the scalar products (l_i l_j, l_i q): greedily from
%% (l_i + q)^2, (l_i + l_j)^2, (l_i - l_j)^2
complete(Loops, LineMoms) ->
    Keys = [{A, B} || {I, A} <- lists:enumerate(Loops), {J, B} <- lists:enumerate(Loops), I =< J] ++ [{A, q} || A <- Loops],
    Vec = fun(Mom) -> [lists:sum([Ca * Cb || {Ca, A} <- Mom, {Cb, B} <- Mom, {A, B} =:= K orelse {B, A} =:= K]) || K <- Keys] end,
    Cands = [[{1, A}, {1, q}] || A <- Loops]
            ++ [[{1, A}, {1, B}] || {I, A} <- lists:enumerate(Loops), {J, B} <- lists:enumerate(Loops), I < J]
            ++ [[{1, A}, {-1, B}] || {I, A} <- lists:enumerate(Loops), {J, B} <- lists:enumerate(Loops), I < J]
            ++ [[{1, A}] || A <- Loops],
    {_, Aux} = lists:foldl(fun(C, {Basis, Acc}) ->
                                   V = Vec(C),
                                   case rank(Basis ++ [V]) > rank(Basis) of
                                       true -> {Basis ++ [V], Acc ++ [C]};
                                       false -> {Basis, Acc}
                                   end
                           end, {[Vec(M) || M <- LineMoms], []}, Cands),
    Aux.

rank([]) -> 0;
rank(Vs) ->
    W = length(hd(Vs)),
    {R, _} = rref([[{X, 1} || X <- V] || V <- Vs], 1, W, []),
    length(R).

%%% The comparison

%% PQCF's combination (#{Master => CoefAst}) with the masters' values (Forcer's exact output for M1, M2, ..., as
%% text), against Forcer's exact result for the target (text); N orders. {Equal, PQCF series, Forcer series}
compare(Coeffs, Masters, MasterValues, TargetValue, Tables) ->
    N = 10,
    Val = maps:from_list([{M, element(2, epsexp:series_of_form(V, Tables))} || {M, V} <- lists:zip(Masters, MasterValues)]),
    Ours = lists:foldl(fun({M, C}, Acc) ->
                               S = epsexp:eval(C, Tables, N),
                               add(Acc, mul(S, (maps:get(M, Val))(N)))
                       end, none, maps:to_list(Coeffs)),
    {_, F} = epsexp:series_of_form(TargetValue, Tables),
    Theirs = F(N),
    {epsexp:equal(Ours, Theirs), Ours, Theirs}.

%% Forcer's msbarexpand=Depth: each part (a G-scheme series and its order: loops plus insertions) times the MS-bar
%% factor to that power, the terms with order + (power of epsilon) > Depth dropped
msbar(Parts, Depth) ->
    lists:foldl(fun({S, Order}, Acc) ->
                        F = epsexp:spow(epsexp:msbar(12), Order, 12),
                        add(Acc, epsexp:drop(mul(S, F), Depth - Order + 1))
                end, none, Parts).

%% Forcer's master integrals (forcer/Masters.prc: R.N. Lee, A.V. Smirnov, V.A. Smirnov,
%% https://arxiv.org/abs/1108.0732, and P.A. Baikov, K.G. Chetyrkin, https://arxiv.org/abs/1004.1153, in Forcer's
%% basis): #{Name => Series}, each known below the order its ALARM term marks
master_tables(File) ->
    {ok, Bin} = file:read_file(File),
    {match, Ms} = re:run(Bin, "id\\s+Master\\(([a-z0-9]+)\\)\\s*=(.*?);", [global, dotall, {capture, all_but_first, list}]),
    maps:from_list(
      [{Name, begin
                  {Text, Prec} = case re:run(Body, "\\+ep\\^([0-9]+)\\*\\(ALARM[a-z0-9]+\\)", [{capture, [0, 1], list}]) of
                                     {match, [Whole, K]} -> {string:replace(Body, Whole, "", all), list_to_integer(K)};
                                     nomatch -> {Body, 1000000}
                                 end,
                  {_, F} = epsexp:series_of_form(lists:flatten(Text), #{}),
                  epsexp:truncate(F(Prec + 2), Prec)
              end} || [Name, Body] <- Ms]).

%% the expected result FORM's check file asserts for an expression in a section: assert result("Name") =~ expr("...")
asserted(File, Section, Name) ->
    {ok, Bin0} = file:read_file(File),
    {match, [Bin]} = re:run(Bin0, "\\*--#\\[ " ++ Section ++ " :(.*?)\\*--#\\] " ++ Section ++ " :", [dotall, {capture, all_but_first, binary}]),
    {match, [E]} = re:run(Bin, "assert result\\(\"" ++ Name ++ "\"\\) =~ expr\\(\"(.*?)\"\\)", [dotall, {capture, all_but_first, list}]),
    E.

add(none, B) -> B;
add(A, B) -> epsexp:sadd(A, B).
mul(A, B) -> epsexp:smul(A, B).

qinv({A, B}) -> qnorm({B, A}).
qmul({A, B}, {C, D}) -> qnorm({A * C, B * D}).
qsub({A, B}, {C, D}) -> qnorm({A * D - C * B, B * D}).
qnorm({N, D}) when D < 0 -> qnorm({-N, -D});
qnorm({0, _}) -> {0, 1};
qnorm({N, D}) -> G = gcd(abs(N), D), {N div G, D div G}.
gcd(A, 0) -> A;
gcd(A, B) -> gcd(B, A rem B).
