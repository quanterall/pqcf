%% A real reduction problem for the linear-system benchmark (decision record R1, phase 5): integration-by-parts (IBP)
%% identities of the one-loop massless box, Laporta-style, reduced to master integrals -- by PQCF's reconstruction
%% from eliminations modulo primes (pq_recon:linear/4), and by FORM 5 eliminating the same system symbolically under
%% PolyRatFun.
%%
%% The family: propagators D1 = k^2, D2 = (k + p1)^2, D3 = (k + p1 + p2)^2, D4 = (k - p4)^2; all external momenta
%% massless, s = (p1 + p2)^2, t = (p2 + p3)^2; the integral I(a1, a2, a3, a4) = int d^d k / prod D_j^a_j (a_j <= 0:
%% numerators). A sector is scaleless, and its integrals zero, unless its positive indices include D1 and D3 (an
%% s-channel cut) or D2 and D4 (a t-channel cut). The identities: int d/dk^mu (v^mu / prod D_j^a_j) = 0 for
%% v = k, p1, p2, p4, with the scalar products written in the D's:
%%   k.k = D1, 2 k.p1 = D2 - D1, 2 k.p2 = D3 - D2 - s, 2 k.p4 = D1 - D4,
%%   p1.p2 = s/2, p1.p4 = t/2, p2.p4 = -(s + t)/2, p_i.p_i = 0.
%%
%% The system: identities for every seed with indices in {-Neg .. 2}, at most Pos positive powers in total, in the
%% nonzero sectors. One elimination modulo a prime at a random (d, s, t), columns in Laporta order (most complex
%% first, the target last), picks the independent equations and the masters (the columns without a pivot); the
%% square system of the pivot rows in the pivot columns, the masters moved to the right-hand side as symbols M1, ..,
%% is what both systems solve for the target.
%%
%% A second family, two loops: the massless propagator "kite" (below), seeds {kite, Pos, Neg}, target I(1,1,1,1,1).
%%
%%   erl -pa _build/default/lib/*/ebin -pa <dir> -noshell -eval 'ibp:main("<dir>", [{4, 0}, {5, 1}, {kite, 6, 1}])'
-module(ibp).
-export([main/2, system/3, identities/1, subst_master/2]).

-define(LIMITS, #{sym_timeout_ms => 3600000, sym_max_heap_words => 4000000000, max_expr_nodes => 500000000}).

-define(P, 1000000007).

%%% The identities

%% the scalar product of v with k + q, q = sum of external momenta (as a map momentum => coefficient), as a linear
%% combination of D's, s, t and constants: {DCoeffs (j => coefficient), ConstAst}
%% coefficients as ASTs over d, s, t
d() -> {sym, d}.
s() -> {sym, s}.
t() -> {sym, t}.

%% v . k in D's (a list of {Coefficient, j | const}); v . p in s, t
vk(k) -> [{1, 1}];
vk(p1) -> [{{rat, 1, 2}, 2}, {{rat, -1, 2}, 1}];
vk(p2) -> [{{rat, 1, 2}, 3}, {{rat, -1, 2}, 2}, {{mul, [{rat, -1, 2}, s()]}, const}];
vk(p4) -> [{{rat, 1, 2}, 1}, {{rat, -1, 2}, 4}].

pp(A, B) when A > B -> pp(B, A);
pp(p1, p2) -> {mul, [{rat, 1, 2}, s()]};
pp(p1, p4) -> {mul, [{rat, 1, 2}, t()]};
pp(p2, p4) -> {mul, [{rat, -1, 2}, {add, [s(), t()]}]};
pp(_, _) -> 0.

%% the shift q_j of propagator j: D_j = (k + q_j)^2
shift(1) -> [];
shift(2) -> [{1, p1}];
shift(3) -> [{1, p1}, {1, p2}];
shift(4) -> [{-1, p4}].

%% v . q_j: for v = k, k.q_j in D's; for v = p_i, constants
vq(k, J) -> lists:append([[{{mul, [C, Ci]}, X} || {Ci, X} <- vk(P)] || {C, P} <- shift(J)]);
vq(V, J) -> [{{mul, [C, pp(V, P)]}, const} || {C, P} <- shift(J)].

v_dot_k(k) -> vk(k);
v_dot_k(V) -> vk(V).

%% one identity: the seed's indices A (a tuple), the vector V; a list of {Integral, CoefficientAst}
identity(A, V) ->
    Div = case V of k -> [{A, d()}]; _ -> [] end,
    Deriv = lists:append(
              [begin
                   Aj = element(J, A),
                   %% -a_j D_j^(-a_j - 1) * 2 v.(k + q_j): raise index j, times the combination
                   Raised = setelement(J, A, Aj + 1),
                   [term(Raised, {mul, [-2 * Aj, C]}, X) || {C, X} <- v_dot_k(V) ++ vq(V, J)]
               end || J <- lists:seq(1, 4), element(J, A) =/= 0]),
    Div ++ Deriv.

%% a D_m (or a constant) times the integral: D_m lowers index m
term(A, C, const) -> {A, C};
term(A, C, M) -> {setelement(M, A, element(M, A) - 1), C}.

zero(A) when tuple_size(A) =:= 5 -> kite_zero(A);
zero(A) ->
    Pos = [J || J <- [1, 2, 3, 4], element(J, A) > 0],
    not ((lists:member(1, Pos) andalso lists:member(3, Pos)) orelse (lists:member(2, Pos) andalso lists:member(4, Pos))).

%% every identity for the seeds with indices in -Neg .. 2, at most Pos positive powers, in nonzero sectors: a list of
%% equations, each a map Integral => CoefficientAst (zero integrals dropped, equal integrals merged)
identities({Pos, Neg}) -> identities({box, Pos, Neg});
identities({Family, Pos, Neg}) ->
    Size = case Family of box -> 4; kite -> 5 end,
    Seeds = [list_to_tuple(L) || L <- tuples(Size, lists:seq(-Neg, 2)),
                                 lists:sum([max(0, X) || X <- L]) =< Pos, not zero(list_to_tuple(L))],
    Ids = case Family of
              box -> [identity(A, V) || A <- Seeds, V <- [k, p1, p2, p4]];
              kite -> [kite_identity(A, Dk, V) || A <- Seeds, Dk <- [k1, k2], V <- [k1, k2, p]]
          end,
    [maps:filter(fun(I, _) -> not zero(I) end,
                 lists:foldl(fun({I, C}, M) -> maps:update_with(I, fun(C0) -> {add, [C0, C]} end, C, M) end, #{}, Id))
     || Id <- Ids].

tuples(0, _) -> [[]];
tuples(N, Range) -> [[X | T] || X <- Range, T <- tuples(N - 1, Range)].

%%% The two-loop massless propagator ("kite") family
%%
%% D1 = k1^2, D2 = k2^2, D3 = (k1 - p)^2, D4 = (k2 - p)^2, D5 = (k1 - k2)^2, p^2 = s: a complete basis of the scalar
%% products (no irreducible numerators). Its graph: external vertices a (p in) and b (p out), inner X and Y; lines
%% 1: a-X, 2: X-b, 3: Y-a, 4: b-Y, 5: X-Y. A missing propagator (index =< 0) contracts its line; the integral is
%% scaleless -- zero -- if a and b merge (no momentum flows) or a remaining line becomes a self-loop (a massless
%% tadpole). The identities: d/dk_i . v for i = 1, 2 and v = k1, k2, p (six per seed).

kite_zero(A) ->
    Lines = [{1, a, x}, {2, x, b}, {3, y, a}, {4, b, y}, {5, x, y}],
    Root = lists:foldl(fun({J, U, V}, R) ->
                               case element(J, A) =< 0 of
                                   true -> union(U, V, R);
                                   false -> R
                               end
                       end, #{}, Lines),
    F = fun(V) -> find(V, Root) end,
    F(a) =:= F(b) orelse lists:any(fun({J, U, V}) -> element(J, A) > 0 andalso F(U) =:= F(V) end, Lines).

find(V, R) -> case R of #{V := P} when P =/= V -> find(P, R); _ -> V end.
union(U, V, R) -> case {find(U, R), find(V, R)} of {X, X} -> R; {X, Y} -> R#{X => Y} end.

%% scalar products in the D's: a list of {Coefficient, j | const}
kdot(k1, k1) -> [{1, 1}];
kdot(k2, k2) -> [{1, 2}];
kdot(p, p) -> [{s(), const}];
kdot(k1, k2) -> [{{rat, 1, 2}, 1}, {{rat, 1, 2}, 2}, {{rat, -1, 2}, 5}];
kdot(k1, p) -> [{{rat, 1, 2}, 1}, {{rat, -1, 2}, 3}, {{mul, [{rat, 1, 2}, s()]}, const}];
kdot(k2, p) -> [{{rat, 1, 2}, 2}, {{rat, -1, 2}, 4}, {{mul, [{rat, 1, 2}, s()]}, const}];
kdot(A, B) -> kdot(B, A).

%% d D_j / d k_i as twice a combination of momenta [{Coefficient, momentum}]
dprop(1, k1) -> [{1, k1}];
dprop(2, k2) -> [{1, k2}];
dprop(3, k1) -> [{1, k1}, {-1, p}];
dprop(4, k2) -> [{1, k2}, {-1, p}];
dprop(5, k1) -> [{1, k1}, {-1, k2}];
dprop(5, k2) -> [{-1, k1}, {1, k2}];
dprop(_, _) -> [].

kite_identity(A, Dk, V) ->
    Div = case Dk =:= V of true -> [{A, d()}]; false -> [] end,
    Deriv = lists:append(
              [begin
                   Aj = element(J, A),
                   Raised = setelement(J, A, Aj + 1),
                   [term(Raised, {mul, [-2 * Aj, C0, C]}, X) || {C0, Mom} <- dprop(J, Dk), {C, X} <- kdot(V, Mom)]
               end || J <- lists:seq(1, 5), element(J, A) =/= 0]),
    Div ++ Deriv.

%%% Laporta order and the square system

%% more complex first: positive propagators, then their total power, then numerator powers, then the tuple
weight(A) ->
    L = tuple_to_list(A),
    {length([x || X <- L, X > 0]), lists:sum([X || X <- L, X > 0]), -lists:sum([X || X <- L, X < 0]), A}.

%% {A, B, Vars, Masters, TargetIndex}: the square system for Target, the masters as symbols M1, ..
system(Seeds, Target, Seed) ->
    Eqs0 = identities(Seeds),
    Eqs = [E || E <- Eqs0, map_size(E) > 0],
    Ints = lists:usort(lists:append([maps:keys(E) || E <- Eqs])),
    lists:member(Target, Ints) orelse error({target_not_in_system, Target}),
    %% columns in Laporta order, most complex first: the elimination finds the independent equations and the
    %% masters (the columns left without a pivot, the simplest)
    Cols = [I || {_, I} <- lists:reverse(lists:sort([{weight(I), I} || I <- Ints]))],
    Index = maps:from_list(lists:zip(Cols, lists:seq(1, length(Cols)))),
    rand:seed(exsss, {Seed, 1, 1}),
    Point = #{d => rand:uniform(?P - 1), s => rand:uniform(?P - 1), t => rand:uniform(?P - 1)},
    Canon = [maps:map(fun(_, C) -> {ok, Cc} = pq_symbolic:canonicalise(C), Cc end, E) || E <- Eqs],
    Num = [maps:from_list([{maps:get(I, Index), modval(C, Point)} || I := C <- E]) || E <- Canon],
    {PivotRows, PivotCols} = eliminate(lists:enumerate(Num), length(Cols)),
    lists:member(maps:get(Target, Index), PivotCols) orelse error(target_is_a_master),
    Masters = [lists:nth(C, Cols) || C <- lists:seq(1, length(Cols)), not lists:member(C, PivotCols)],
    MSyms = [{sym, list_to_atom("m" ++ integer_to_list(K))} || K <- lists:seq(1, length(Masters))],
    MIndex = maps:from_list(lists:zip(Masters, MSyms)),
    %% the square system, its columns re-ordered with the target last, and eliminated again (it is nonsingular, so
    %% every column finds a pivot row): the rows in that pivot order, so that eliminating unknown k on row k from
    %% the first column on is valid, and the target is eliminated last
    UCols = [lists:nth(C, Cols) || C <- lists:sort(PivotCols), lists:nth(C, Cols) =/= Target] ++ [Target],
    UIndex = maps:from_list(lists:zip(UCols, lists:seq(1, length(UCols)))),
    SqRows = [lists:nth(R, Canon) || R <- PivotRows],
    SqNum = [maps:from_list([{maps:get(I, UIndex), modval(C, Point)} || I := C <- Row, maps:is_key(I, UIndex)]) || Row <- SqRows],
    {Order, Cols2} = eliminate(lists:enumerate(SqNum), length(UCols)),
    Cols2 =:= lists:seq(1, length(UCols)) orelse error(square_system_singular),
    Rows = [lists:nth(R, SqRows) || R <- Order],
    A = [[maps:get(I, Row, 0) || I <- UCols] || Row <- Rows],
    B = [{add, [0 | [{mul, [-1, C, maps:get(I, MIndex)]} || I := C <- Row, maps:is_key(I, MIndex)]]} || Row <- Rows],
    Vars = invariants(Target) ++ [X || {sym, X} <- MSyms],
    {A, B, Vars, lists:zip(Masters, MSyms), length(UCols) - 1, UCols}.

invariants(T) when tuple_size(T) =:= 5 -> [d, s];
invariants(_) -> [d, s, t].

%% Gaussian elimination mod ?P, rows sparse maps col => value, columns in order 1 .. N: {pivot rows, pivot columns}
eliminate(Rows, N) -> elim(Rows, 1, N, [], []).

elim(_Rows, C, N, PR, PC) when C > N -> {lists:reverse(PR), lists:reverse(PC)};
elim(Rows, C, N, PR, PC) ->
    case [{R, Row} || {R, Row} <- Rows, maps:get(C, Row, 0) =/= 0] of
        [] -> elim(Rows, C + 1, N, PR, PC);
        [{R, Row} | _] ->
            Inv = inv(maps:get(C, Row), ?P),
            Rows2 = [{R2, case maps:get(C, Row2, 0) of
                              0 -> Row2;
                              F -> sub(Row2, Row, F * Inv rem ?P)
                          end} || {R2, Row2} <- Rows, R2 =/= R],
            elim(Rows2, C + 1, N, [R | PR], [C | PC])
    end.

sub(A, B, F) ->
    maps:filter(fun(_, V) -> V =/= 0 end,
                maps:fold(fun(K, V, Acc) -> maps:put(K, ((maps:get(K, Acc, 0) - F * V) rem ?P + ?P) rem ?P, Acc) end, A, B)).

modval(C, Point) ->
    Bind = maps:map(fun(_, V) -> V end, Point),
    {ok, V} = pq_symbolic:evaluate(C, Bind, #{exact => true}),
    case V of
        N when is_integer(N) -> (N rem ?P + ?P) rem ?P;
        {rat, N, D} -> ((N rem ?P + ?P) rem ?P) * inv(D, ?P) rem ?P
    end.

inv(A, M) -> pq_ff:inv(A, M).

%%% The run

main(Dir, Cases) ->
    {ok, _} = application:ensure_all_started(pq),
    [begin
         Target = case Seeds of {kite, _, _} -> {1, 1, 1, 1, 1}; _ -> {2, 1, 1, 1} end,
         {Ts, {A, B, Vars, Masters, TI, _Cols}} = timer:tc(?MODULE, system, [Seeds, Target, 1]),
         N = length(B),
         %% the answer is linear in the masters: one right-hand side per master, each coefficient a rational function
         %% of d, s, t alone; one elimination per point serves them all, and the zero ones are found at once
         MSyms = [X || {_, {sym, X}} <- Masters],
         Bs = [[subst_master(Bi, X) || X <- MSyms] || Bi <- B],
         {T, {ok, [{TI, PerMaster}]}} = timer:tc(pq_recon, linear, [A, Bs, invariants(Target), ?LIMITS#{components => [TI]}]),
         Coeffs = [{M, Sym, lists:nth(C + 1, [ND || {_, ND} <- lists:sort(PerMaster)])} || {C, {M, Sym}} <- lists:enumerate(0, Masters)],
         NonZero = [{M, Sym, NumDen} || {M, Sym, {Nu, _} = NumDen} <- Coeffs, map_size(Nu) > 0],
         P = fun(E) -> binary_to_list(pq_symbolic:pretty(E)) end,
         File = filename:join(Dir, case Seeds of
                                       {F, Po, Ne} -> io_lib:format("ibp_~s_~b_~b.frm", [F, Po, Ne]);
                                       {Po, Ne} -> io_lib:format("ibp_~b_~b.frm", [Po, Ne])
                                   end),
         Rev = fun(L) -> lists:reverse(L) end,
         Answer = [{P(pq_sym_norm:to_ast(Nu)), P(pq_sym_norm:to_ast(De)), atom_to_list(X)} || {_, {sym, X}, {Nu, De}} <- NonZero],
         ok = file:write_file(File, form_script(Rev([Rev(Row) || Row <- A]), Rev(B), Vars, Answer)),
         io:format("IBP, seeds ~p: ~b x ~b system (built in ~.1f s), ~b masters; I~w reduced in ~.2f s to ~b of them: ~s~n",
                   [Seeds, N, N, Ts / 1.0e6, length(Masters), Target, T / 1.0e6, length(NonZero),
                    lists:join(", ", [io_lib:format("~w (~b/~b terms)", [M, map_size(Nu), map_size(De)]) || {M, _, {Nu, De}} <- NonZero])])
     end || Seeds <- Cases],
    ok.

%% b_i with master X set to 1 and every other master to 0
subst_master(E, X) ->
    {ok, C} = pq_symbolic:canonicalise(pq_symbolic:substitute(E, #{X => 1})),
    {ok, C0} = pq_symbolic:canonicalise(zero_masters(C)),
    C0.

zero_masters({sym, M} = S) ->
    case atom_to_list(M) of [$m | _] -> 0; _ -> S end;
zero_masters(T) when is_tuple(T) -> list_to_tuple([zero_masters(E) || E <- tuple_to_list(T)]);
zero_masters(L) when is_list(L) -> [zero_masters(E) || E <- L];
zero_masters(X) -> X.

%% FORM: the elimination of linsys.erl (unknown 1 the target), the masters symbols; at the end U minus PQCF's answer,
%% sum_j rat(N_j m_j, D_j), must be 0
form_script(A, B, Vars, Answer) ->
    Base = linsys:form_script(A, B, Vars, "0", "1"),
    Sum = lists:join(" + ", [["rat((", Nu, ")*", M, ",", De, ")"] || {Nu, De, M} <- Answer]),
    Check = ["Local Z = U - (", case Sum of [] -> "0"; _ -> Sum end, ");\n"],
    re:replace(Base, "Local Z = U - rat\\(\n0,\n1\\);\n", lists:flatten(Check), [{return, list}]).

