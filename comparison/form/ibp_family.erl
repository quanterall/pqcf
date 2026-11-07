%% Integration-by-parts identities for any family of propagators (decision record R1, phase 5): the scalar products
%% solved from the propagators' definitions, the scaleless (zero) sectors by the Lee-Pak criterion, the seeds in
%% Kira's form (r: the sum of the positive indices, s: of the numerator powers). ibp.erl's box and kite are written
%% out by hand; this module gives the same identities for them (checked) and any other family, such as FIRE's
%% double box.
%%
%% A family: #{props => [{Sign, Momentum}]} with Momentum = [{Coefficient, q}] over the loop momenta (loops) and
%% the external ones (ext); D_a = Sign_a * Momentum_a^2; kin(Pi, Pj) the external scalar products as ASTs and
%% kin_at(Pi, Pj, Point) their values mod p at a point of the invariants; top: the indices that can be positive.
%% The propagators must span the scalar products of the loop momenta (with each other and with the externals).
-module(ibp_family).
-export([family/1, identities/2, zero/2, seeds/3]).

-define(P, 1000000007).

%%% The families

d() -> {sym, d}.

family(doublebox) ->
    %% FIRE's examples/doublebox: Internal {k1, k2}, External {p1, p2, p3}, all massless, p1.p2 = s/2,
    %% p1.p3 = t/2, p2.p3 = -(s + t)/2; propagators -k1^2, -(k1+p1+p2)^2, -k2^2, -(k2+p1+p2)^2, -(k1+p1)^2,
    %% -(k1-k2)^2, -(k2-p3)^2; numerators -(k2+p1)^2, -(k1-p3)^2
    S = {sym, s}, T = {sym, t},
    Kin = fun(A, B) when A > B -> kin3(B, A, S, T); (A, B) -> kin3(A, B, S, T) end,
    #{name => doublebox, loops => [k1, k2], ext => [p1, p2, p3], vars => [d, s, t], top => [1, 2, 3, 4, 5, 6, 7],
      props => [{-1, [{1, k1}]}, {-1, [{1, k1}, {1, p1}, {1, p2}]}, {-1, [{1, k2}]}, {-1, [{1, k2}, {1, p1}, {1, p2}]},
                {-1, [{1, k1}, {1, p1}]}, {-1, [{1, k1}, {-1, k2}]}, {-1, [{1, k2}, {-1, p3}]},
                {-1, [{1, k2}, {1, p1}]}, {-1, [{1, k1}, {-1, p3}]}],
      kin => Kin};
family(box) ->
    %% ibp.erl's box: k^2, (k+p1)^2, (k+p1+p2)^2, (k-p4)^2, p1.p2 = s/2, p1.p4 = t/2, p2.p4 = -(s+t)/2
    S = {sym, s}, T = {sym, t},
    Kin = fun(A, B) when A > B -> kin3(ren(B), ren(A), S, T); (A, B) -> kin3(ren(A), ren(B), S, T) end,
    #{name => box, loops => [k], ext => [p1, p2, p4], vars => [d, s, t], top => [1, 2, 3, 4],
      props => [{1, [{1, k}]}, {1, [{1, k}, {1, p1}]}, {1, [{1, k}, {1, p1}, {1, p2}]}, {1, [{1, k}, {-1, p4}]}],
      kin => Kin};
family(kite) ->
    S = {sym, s},
    Kin = fun(p, p) -> S end,
    #{name => kite, loops => [k1, k2], ext => [p], vars => [d, s], top => [1, 2, 3, 4, 5],
      props => [{1, [{1, k1}]}, {1, [{1, k2}]}, {1, [{1, k1}, {-1, p}]}, {1, [{1, k2}, {-1, p}]}, {1, [{1, k1}, {-1, k2}]}],
      kin => Kin}.

ren(p4) -> p3;
ren(P) -> P.

kin3(p1, p2, S, _) -> {mul, [{rat, 1, 2}, S]};
kin3(p1, p3, _, T) -> {mul, [{rat, 1, 2}, T]};
kin3(p2, p3, S, T) -> {mul, [{rat, -1, 2}, {add, [S, T]}]};
kin3(_, _, _, _) -> 0.

%%% Scalar products in the D's

%% the scalar products with a loop momentum, in a fixed order: {li, lj} (i =< j) and {li, e}
sp_keys(#{loops := L, ext := E}) ->
    [{A, B} || {I, A} <- lists:enumerate(L), {J, B} <- lists:enumerate(L), I =< J] ++ [{A, B} || A <- L, B <- E].

%% q.q' for two momenta: a scalar product key, or {const, Ast} for two external ones
pair(Fam = #{loops := L}, Q1, Q2) ->
    case {lists:member(Q1, L), lists:member(Q2, L)} of
        {false, false} -> {const, (maps:get(kin, Fam))(Q1, Q2)};
        {true, false} -> {sp, {Q1, Q2}};
        {false, true} -> {sp, {Q2, Q1}};
        {true, true} -> {sp, case index(Q1, L) =< index(Q2, L) of true -> {Q1, Q2}; false -> {Q2, Q1} end}
    end.

index(X, L) -> length(lists:takewhile(fun(Y) -> Y =/= X end, L)) + 1.

%% every scalar product with a loop momentum as [{Rational, j}] (D_j) plus a constant Ast: #{Key => {Ds, Const}};
%% memoised per family
sps(Fam = #{name := Name, props := Props}) ->
    case get({sps, Name}) of
        undefined ->
            Keys = sp_keys(Fam),
            length(Keys) =:= length(Props) orelse error({not_a_basis, length(Keys), length(Props)}),
            %% D_a = sum_k M_ak sp_k + c_a
            Rows = [begin
                        Terms = [{Sign * C1 * C2, pair(Fam, Q1, Q2)} || {C1, Q1} <- Mom, {C2, Q2} <- Mom],
                        {[lists:sum([C || {C, {sp, K1}} <- Terms, K1 =:= K]) || K <- Keys],
                         {add, [0 | [{mul, [C, X]} || {C, {const, X}} <- Terms]]}}
                    end || {Sign, Mom} <- Props],
            Inv = inverse([[{C, 1} || C <- R] || {R, _} <- Rows]),
            Consts = [C || {_, C} <- Rows],
            %% sp_k = sum_a Inv_ka (D_a - c_a)
            Map = maps:from_list(
                    [{K, {[{Q, J} || {J, Q} <- lists:enumerate(InvRow), Q =/= {0, 1}],
                          {add, [0 | [{mul, [-1, rat(Q), C]} || {Q, C} <- lists:zip(InvRow, Consts), Q =/= {0, 1}]]}}}
                     || {K, InvRow} <- lists:zip(Keys, Inv)]),
            put({sps, Name}, Map),
            Map;
        Map -> Map
    end.

rat({N, 1}) -> N;
rat({N, D}) -> {rat, N, D}.

%% v.q as [{CoefficientAst, j | const}]
dot(Fam, V, Q) ->
    case pair(Fam, V, Q) of
        {const, X} -> [{X, const}];
        {sp, K} -> {Ds, C} = maps:get(K, sps(Fam)), [{rat(R), J} || {R, J} <- Ds] ++ [{C, const}]
    end.

%% the inverse of a square matrix of rationals {N, D} (Gauss-Jordan)
inverse(M) ->
    N = length(M),
    Aug = [Row ++ [case I =:= J of true -> {1, 1}; false -> {0, 1} end || J <- lists:seq(1, N)] || {I, Row} <- lists:enumerate(M)],
    [lists:nthtail(N, Row) || Row <- gj(Aug, 1, N)].

gj(Rows, C, N) when C > N -> Rows;
gj(Rows, C, N) ->
    {Before, [Piv | After]} = lists:splitwith(fun(R) -> lists:nth(C, R) =:= {0, 1} end, lists:nthtail(C - 1, Rows)),
    Done = lists:sublist(Rows, C - 1),
    Inv = qinv(lists:nth(C, Piv)),
    PivN = [qmul(Inv, X) || X <- Piv],
    Elim = fun(R) -> F = lists:nth(C, R), [qsub(X, qmul(F, Y)) || {X, Y} <- lists:zip(R, PivN)] end,
    gj([Elim(R) || R <- Done] ++ [PivN] ++ [Elim(R) || R <- Before ++ After], C + 1, N).

qnorm({N, D}) when D < 0 -> qnorm({-N, -D});
qnorm({N, D}) -> G = max(1, gcd(abs(N), D)), {N div G, D div G}.
gcd(A, 0) -> A;
gcd(A, B) -> gcd(B, A rem B).
qmul({A, B}, {C, D}) -> qnorm({A * C, B * D}).
qsub({A, B}, {C, D}) -> qnorm({A * D - C * B, B * D}).
qinv({A, B}) -> qnorm({B, A}).

%%% The identities

%% d/dl . v of the seed A (a tuple): [{Integral, CoefficientAst}], before merging
identity(Fam = #{props := Props}, A, L, V) ->
    Div = case L =:= V of true -> [{A, d()}]; false -> [] end,
    Deriv = lists:append(
              [begin
                   Aj = element(J, A),
                   Raised = setelement(J, A, Aj + 1),
                   CL = lists:sum([C || {C, Q} <- Mom, Q =:= L]),
                   %% -a_j D_j^(-a_j-1) dD_j/dl.v, dD_j/dl = 2 Sign_j CL Mom_j
                   [term(Raised, {mul, [-2 * Aj * Sign * CL * Cq, X]}, Y)
                    || CL =/= 0, {Cq, Q} <- Mom, {X, Y} <- dot(Fam, V, Q)]
               end || {J, {Sign, Mom}} <- lists:enumerate(Props), element(J, A) =/= 0]),
    Div ++ Deriv.

term(A, C, const) -> {A, C};
term(A, C, M) -> {setelement(M, A, element(M, A) - 1), C}.

%% every identity of every seed: maps Integral => CoefficientAst, zero integrals dropped, equal ones merged
identities(Fam = #{loops := Ls, ext := Es}, Seeds) ->
    [maps:filter(fun(I, _) -> not zero(Fam, I) end,
                 lists:foldl(fun({I, C}, M) -> maps:update_with(I, fun(C0) -> {add, [C0, C]} end, C, M) end, #{},
                             identity(Fam, A, L, V)))
     || A <- Seeds, L <- Ls, V <- Ls ++ Es].

%%% The seeds: every non-zero sector within top, the positive indices summing to at most R, the others (numerators)
%%% to at least -S

seeds(Fam = #{top := Top, props := Props}, R, S) ->
    N = length(Props),
    Sectors = [Sec || Sec <- subsets(Top), Sec =/= [], length(Sec) =< R],
    lists:append(
      [[list_to_tuple(fill(N, Sec, Pos, Others, Neg))
        || Pos <- compositions(Sec, R),
           Neg <- numerators(length(Others), S)]
       || Sec <- Sectors, Others <- [[J || J <- lists:seq(1, N), not lists:member(J, Sec)]],
          not zero(Fam, sector_point(N, Sec))]).

subsets([]) -> [[]];
subsets([X | Xs]) -> Rest = subsets(Xs), [[X | Ys] || Ys <- Rest] ++ Rest.

sector_point(N, Sec) -> list_to_tuple([case lists:member(J, Sec) of true -> 1; false -> 0 end || J <- lists:seq(1, N)]).

%% the powers of the sector's propagators: each >= 1, summing to at most R
compositions(Sec, R) -> parts(length(Sec), R - length(Sec), fun(E) -> E + 1 end).
%% the numerator powers on the other indices: each =< 0, summing to at least -S
numerators(K, S) -> parts(K, S, fun(E) -> -E end).

%% K non-negative integers summing to at most Budget, mapped by F
parts(0, _, _) -> [[]];
parts(K, Budget, F) -> [[F(E) | Rest] || E <- lists:seq(0, Budget), Rest <- parts(K - 1, Budget - E, F)].

fill(N, Sec, Pos, Others, Neg) ->
    M = maps:from_list(lists:zip(Sec, Pos) ++ lists:zip(Others, Neg)),
    [maps:get(J, M) || J <- lists:seq(1, N)].

%%% Scaleless sectors: the Lee-Pak criterion. With U and F the Symanzik polynomials of the sector (the x_j of the
%%% other propagators set to 0), the integral is 0 if and only if sum_j k_j x_j dG/dx_j = G has a solution in
%%% constants k_j, G = U + F: every monomial x^alpha of G has k.alpha = 1, a linear system for k. U and F are built
%%% mod p at a random point of the invariants (their supports are those of the generic ones but for a cancellation
%%% of probability about 1/p).

zero(Fam = #{name := Name}, A) ->
    Sec = [J || {J, X} <- lists:enumerate(tuple_to_list(A)), X > 0],
    case get({zero, Name, Sec}) of
        undefined -> Z = scaleless(Fam, Sec), put({zero, Name, Sec}, Z), Z;
        Z -> Z
    end.

scaleless(Fam, Sec) ->
    {U, F} = symanzik(Fam),
    In = fun(Mono) -> lists:all(fun({J, E}) -> E =:= 0 orelse lists:member(J, Sec) end, lists:enumerate(Mono)) end,
    Support = [Mono || Mono := C <- maps:merge(U, F), C =/= 0, In(Mono)],
    case [M || M := C <- U, C =/= 0, In(M)] of
        [] -> true;
        _ ->
            %% k . alpha = 1 for every monomial: rows [alpha_j (j in Sec) | 1]
            Rows = [[lists:nth(J, Mono) || J <- Sec] ++ [1] || Mono <- Support],
            consistent(Rows, length(Sec))
    end.

%% U and F mod p as maps Exponents => coefficient (memoised)
symanzik(Fam = #{name := Name, loops := Ls, props := Props, kin := Kin}) ->
    case get({symanzik, Name}) of
        undefined ->
            N = length(Props),
            rand:seed(exsss, {7, 7, 7}),
            Point = maps:from_list([{V, rand:uniform(?P - 1)} || V <- maps:get(vars, Fam)]),
            KinAt = fun(A, B) -> {ok, V} = pq_symbolic:evaluate(Kin(A, B), Point, #{exact => true}), modq(V) end,
            X = fun(J, C) -> #{[case I =:= J of true -> 1; false -> 0 end || I <- lists:seq(1, N)] => (C rem ?P + ?P) rem ?P} end,
            Coef = fun(Mom, Q) -> lists:sum([C || {C, Q1} <- Mom, Q1 =:= Q]) end,
            Ext = fun(Mom) -> [{C, Q} || {C, Q} <- Mom, not lists:member(Q, Ls)] end,
            %% M_ij = sum_a Sign_a c_ai c_aj x_a; J_i = sum_a Sign_a c_ai x_a P_a (P_a the external part, a list
            %% of {Poly, e}); J0 = sum_a Sign_a x_a P_a^2
            M = [[padd([X(J, Sign * Coef(Mom, Li) * Coef(Mom, Lj)) || {J, {Sign, Mom}} <- lists:enumerate(Props)])
                  || Lj <- Ls] || Li <- Ls],
            Js = [[{X(J, Sign * Coef(Mom, Li) * C), Q} || {J, {Sign, Mom}} <- lists:enumerate(Props), {C, Q} <- Ext(Mom)]
                  || Li <- Ls],
            J0 = padd([X(J, Sign * C1 * C2 * KinAt(Q1, Q2)) || {J, {Sign, Mom}} <- lists:enumerate(Props),
                                                              {C1, Q1} <- Ext(Mom), {C2, Q2} <- Ext(Mom)]),
            JJ = fun(Ji, Jj) -> padd([pscale(pmul(P1, P2), KinAt(Q1, Q2)) || {P1, Q1} <- Ji, {P2, Q2} <- Jj]) end,
            {U, Adj} = case M of
                           [[M11]] -> {M11, [[#{lists:duplicate(N, 0) => 1}]]};
                           [[M11, M12], [_, M22]] ->
                               {padd([pmul(M11, M22), pscale(pmul(M12, M12), ?P - 1)]),
                                [[M22, pscale(M12, ?P - 1)], [pscale(M12, ?P - 1), M11]]};
                           _ -> error(more_than_two_loops)
                       end,
            F = padd([pmul(Aij, JJ(Ji, Jj)) || {Ji, AdjRow} <- lists:zip(Js, Adj), {Jj, Aij} <- lists:zip(Js, AdjRow)]
                     ++ [pscale(pmul(J0, U), ?P - 1)]),
            put({symanzik, Name}, {U, F}),
            {U, F};
        UF -> UF
    end.

modq(N) when is_integer(N) -> (N rem ?P + ?P) rem ?P;
modq({rat, N, D}) -> modq(N) * pq_ff:inv(modq(D), ?P) rem ?P.

padd(Ps) ->
    maps:filter(fun(_, C) -> C =/= 0 end,
                lists:foldl(fun(P, Acc) -> maps:fold(fun(K, C, A) -> A#{K => (maps:get(K, A, 0) + C) rem ?P} end, Acc, P) end, #{}, Ps)).
pmul(A, B) -> padd([#{[E1 + E2 || {E1, E2} <- lists:zip(K1, K2)] => C1 * C2 rem ?P} || K1 := C1 <- A, K2 := C2 <- B]).
pscale(A, C) -> padd([maps:map(fun(_, V) -> V * C rem ?P end, A)]).

%% is the augmented system [A | b] (rows of integers, K unknowns) consistent mod p?
consistent(Rows, K) ->
    Red = reduce([[(X rem ?P + ?P) rem ?P || X <- R] || R <- Rows], 1, K),
    not lists:any(fun(R) -> lists:all(fun(X) -> X =:= 0 end, lists:sublist(R, K)) andalso lists:last(R) =/= 0 end, Red).

reduce(Rows, C, K) when C > K -> Rows;
reduce(Rows, C, K) ->
    case lists:partition(fun(R) -> lists:nth(C, R) =/= 0 end, Rows) of
        {[], _} -> reduce(Rows, C + 1, K);
        {[Piv | More], Rest} ->
            Inv = pq_ff:inv(lists:nth(C, Piv), ?P),
            Elim = fun(R) -> F = lists:nth(C, R) * Inv rem ?P, [(X - F * Y rem ?P + ?P) rem ?P || {X, Y} <- lists:zip(R, Piv)] end,
            [Piv | reduce([Elim(R) || R <- More] ++ Rest, C + 1, K)]
    end.
