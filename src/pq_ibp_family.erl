%% Integration-by-parts (IBP) identities for a family of Feynman integrals (decision record R1): the scalar products
%% solved from the propagators' definitions, the scaleless (zero) sectors by the Lee-Pak criterion, the sector
%% symmetries by Pak's criterion, the seeds in Kira's form (or with the dots per sector bounded), and derivatives
%% with respect to the invariants (for the masters' differential equations). pq_ibp builds and solves the systems.
%%
%% A family (new/1): loops and ext, the loop and the independent external momenta (atoms); props, the propagators
%% [{Sign, Momentum}], Momentum = [{IntegerCoefficient, q}] over them, D_a = Sign_a * Momentum_a^2 (numerators
%% included: they must span every scalar product with a loop momentum, as many as there are); kin, the external
%% scalar products #{{Pa, Pb} => Ast} in the invariants (a missing pair is 0); vars, the variables of the results (the
%% dimension d first, then the invariants); top, the indices that can be positive (default: all); symmetries
%% (default true). The integral I(A) = int prod D_a^-A_a for a tuple A.
%%
%% Results are cached per process by the family's id (its hash): the scalar products, the Symanzik polynomials, the
%% zero sectors, the symmetry classes and the derivative operators. Random points are drawn from local states,
%% never the caller's.
-module(pq_ibp_family).
-export([new/1, example/1, identities/2, zero/2, seeds/3, canon/2, derivative/3, vars/1, size/1, kin/3]).
%% for pq_ibp: the identities with their coefficients as polynomials (Exponents => {N, D})
-export([identities_poly/2]).
-export_type([family/0, integral/0]).

-type family() :: #{id := integer(), loops := [atom()], ext := [atom()], props := [{integer(), [{integer(), atom()}]}],
                    kin := #{{atom(), atom()} => term()}, vars := [atom()], top := [pos_integer()], symmetries := boolean()}.
-type integral() :: tuple().

-define(P, 1000000007).

%%% Families

-spec new(map()) -> family().
new(Spec) ->
    Loops = maps:get(loops, Spec),
    Ext = maps:get(ext, Spec, []),
    Props = maps:get(props, Spec),
    Kin0 = maps:get(kin, Spec, #{}),
    Vars = maps:get(vars, Spec),
    N = length(Props),
    Moms = Loops ++ Ext,
    ok = check(is_list(Loops) andalso Loops =/= [] andalso lists:all(fun is_atom/1, Moms), loops),
    ok = check(lists:all(fun({S, Mom}) -> (S =:= 1 orelse S =:= -1) andalso Mom =/= []
                                              andalso lists:all(fun({C, Q}) -> is_integer(C) andalso lists:member(Q, Moms) end, Mom);
                            (_) -> false end, Props), props),
    ok = check(lists:all(fun is_atom/1, Vars) andalso lists:member(d, Vars), vars),
    Kin = maps:from_list([{pair_key(A, B), X} || {A, B} := X <- Kin0]),
    ok = check(lists:all(fun({A, B}) -> lists:member(A, Ext) andalso lists:member(B, Ext) end, maps:keys(Kin)), kin),
    Fam0 = #{loops => Loops, ext => Ext, props => Props, kin => Kin, vars => Vars,
             top => maps:get(top, Spec, lists:seq(1, N)), symmetries => maps:get(symmetries, Spec, true)},
    Fam = Fam0#{id => erlang:phash2(Fam0)},
    ok = check(length(sp_keys(Fam)) =:= N, {not_a_basis, length(sp_keys(Fam)), N}),
    Fam.

check(true, _) -> ok;
check(false, What) -> error({bad_family, What}).

pair_key(A, B) when A =< B -> {A, B};
pair_key(A, B) -> {B, A}.

%% the external scalar product pa . pb as an AST
-spec kin(family(), atom(), atom()) -> term().
kin(#{kin := Kin}, A, B) -> maps:get(pair_key(A, B), Kin, 0).

-spec vars(family()) -> [atom()].
vars(#{vars := V}) -> V.

%% the number of indices of the family's integrals
-spec size(family()) -> pos_integer().
size(#{props := P}) -> length(P).

%% The families of the benchmarks: the one-loop massless box; the two-loop massless propagator "kite"; FIRE's
%% two-loop planar double box (examples/doublebox of FIRE 6.5: propagators -k1^2, -(k1+p1+p2)^2, -k2^2,
%% -(k2+p1+p2)^2, -(k1+p1)^2, -(k1-k2)^2, -(k2-p3)^2, numerators -(k2+p1)^2, -(k1-p3)^2)
-spec example(box | kite | doublebox) -> family().
example(box) ->
    S = {sym, s}, T = {sym, t},
    new(#{loops => [k], ext => [p1, p2, p4], vars => [d, s, t],
          props => [{1, [{1, k}]}, {1, [{1, k}, {1, p1}]}, {1, [{1, k}, {1, p1}, {1, p2}]}, {1, [{1, k}, {-1, p4}]}],
          kin => #{{p1, p2} => {mul, [{rat, 1, 2}, S]}, {p1, p4} => {mul, [{rat, 1, 2}, T]},
                   {p2, p4} => {mul, [{rat, -1, 2}, {add, [S, T]}]}}});
example(kite) ->
    new(#{loops => [k1, k2], ext => [p], vars => [d, s],
          props => [{1, [{1, k1}]}, {1, [{1, k2}]}, {1, [{1, k1}, {-1, p}]}, {1, [{1, k2}, {-1, p}]}, {1, [{1, k1}, {-1, k2}]}],
          kin => #{{p, p} => {sym, s}}});
example(doublebox) ->
    S = {sym, s}, T = {sym, t},
    new(#{loops => [k1, k2], ext => [p1, p2, p3], vars => [d, s, t], top => [1, 2, 3, 4, 5, 6, 7],
          props => [{-1, [{1, k1}]}, {-1, [{1, k1}, {1, p1}, {1, p2}]}, {-1, [{1, k2}]}, {-1, [{1, k2}, {1, p1}, {1, p2}]},
                    {-1, [{1, k1}, {1, p1}]}, {-1, [{1, k1}, {-1, k2}]}, {-1, [{1, k2}, {-1, p3}]},
                    {-1, [{1, k2}, {1, p1}]}, {-1, [{1, k1}, {-1, p3}]}],
          kin => #{{p1, p2} => {mul, [{rat, 1, 2}, S]}, {p1, p3} => {mul, [{rat, 1, 2}, T]},
                   {p2, p3} => {mul, [{rat, -1, 2}, {add, [S, T]}]}}}).

%% the per-process cache, keyed by the family's id
memo(Tag, #{id := Id}, Key, Fun) ->
    K = {?MODULE, Tag, Id, Key},
    case get(K) of
        undefined -> V = Fun(), put(K, V), V;
        V -> V
    end.

%% a random point of the variables mod ?P from a local state (the caller's random state is left alone)
point(Vars, Seed) ->
    {L, _} = lists:mapfoldl(fun(V, S) -> {X, S2} = rand:uniform_s(?P - 1, S), {{V, X}, S2} end, rand:seed_s(exsss, Seed), Vars),
    maps:from_list(L).

d() -> {sym, d}.

%%% Scalar products in the D's

%% the scalar products with a loop momentum, in a fixed order: {li, lj} (i =< j) and {li, e}
sp_keys(#{loops := L, ext := E}) ->
    [{A, B} || {I, A} <- lists:enumerate(L), {J, B} <- lists:enumerate(L), I =< J] ++ [{A, B} || A <- L, B <- E].

%% q.q' for two momenta: a scalar product key, or {const, Ast} for two external ones
pair(Fam = #{loops := L}, Q1, Q2) ->
    case {lists:member(Q1, L), lists:member(Q2, L)} of
        {false, false} -> {const, kin(Fam, Q1, Q2)};
        {true, false} -> {sp, {Q1, Q2}};
        {false, true} -> {sp, {Q2, Q1}};
        {true, true} -> {sp, case index(Q1, L) =< index(Q2, L) of true -> {Q1, Q2}; false -> {Q2, Q1} end}
    end.

index(X, L) -> length(lists:takewhile(fun(Y) -> Y =/= X end, L)) + 1.

%% every scalar product with a loop momentum as [{Rational, j}] (D_j) plus a constant Ast: #{Key => {Ds, Const}};
%% memoised per family
sps(Fam = #{props := Props}) ->
    memo(sps, Fam, all, fun() ->
            Keys = sp_keys(Fam),
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
            Map
    end).

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

%% every identity of every seed: maps Integral => CoefficientAst, zero integrals dropped, equal ones merged (equal
%% by the sector symmetries too, canon/2, unless the family has symmetries => false)
-spec identities(family(), [integral()]) -> [#{integral() => term()}].
identities(Fam = #{loops := Ls, ext := Es}, Seeds) ->
    Canon = case maps:get(symmetries, Fam) of
                true -> fun(I) -> canon(Fam, I) end;
                false -> fun(I) -> I end
            end,
    [lists:foldl(fun({I, C}, M) ->
                         case zero(Fam, I) of
                             true -> M;
                             false -> maps:update_with(Canon(I), fun(C0) -> {add, [C0, C]} end, C, M)
                         end
                 end, #{}, identity(Fam, A, L, V))
     || A <- Seeds, L <- Ls, V <- Ls ++ Es].

%%% The identities with polynomial coefficients (pq_ibp's fast path): every coefficient an exact polynomial in the
%%% family's variables, a map Exponents => {N, D} -- no general symbolic arithmetic -- and each integral's zero test
%%% and canonical name computed once per call.

-spec identities_poly(family(), [integral()]) -> [#{integral() => #{[non_neg_integer()] => {integer(), pos_integer()}}}].
identities_poly(Fam = #{loops := Ls, ext := Es}, Seeds) ->
    Sym = maps:get(symmetries, Fam),
    Name = fun(I) ->
                   case get({?MODULE, name, I}) of
                       undefined ->
                           V = case zero(Fam, I) of
                                   true -> zero;
                                   false when Sym -> canon(Fam, I);
                                   false -> I
                               end,
                           put({?MODULE, name, I}, V),
                           V;
                       V -> V
                   end
           end,
    Eqs = [lists:foldl(fun({I, C}, M) ->
                               case Name(I) of
                                   zero -> M;
                                   J -> case M of
                                            #{J := C0} -> M#{J => qpadd(C0, C)};
                                            _ -> M#{J => C}
                                        end
                               end
                       end, #{}, identity_poly(Fam, A, L, V))
           || A <- Seeds, L <- Ls, V <- Ls ++ Es],
    [maps:filter(fun(_, C) -> map_size(C) > 0 end, E) || E <- Eqs].

identity_poly(Fam = #{props := Props, vars := Vars}, A, L, V) ->
    Div = case L =:= V of
              true -> [{A, #{[case X of d -> 1; _ -> 0 end || X <- Vars] => {1, 1}}}];
              false -> []
          end,
    Deriv = lists:append(
              [begin
                   Aj = element(J, A),
                   Raised = setelement(J, A, Aj + 1),
                   CL = lists:sum([C || {C, Q} <- Mom, Q =:= L]),
                   [term(Raised, qpscale(X, {-2 * Aj * Sign * CL * Cq, 1}), Y)
                    || CL =/= 0, {Cq, Q} <- Mom, {X, Y} <- dot_poly(Fam, V, Q)]
               end || {J, {Sign, Mom}} <- lists:enumerate(Props), element(J, A) =/= 0]),
    Div ++ Deriv.

%% dot/3 with polynomial coefficients (memoised per pair)
dot_poly(Fam = #{vars := Vars}, V, Q) ->
    memo(dot_poly, Fam, {V, Q}, fun() ->
            [{P, Y} || {X, Y} <- dot(Fam, V, Q), P <- [qpoly(X, Vars)], map_size(P) > 0]
    end).

%% an AST of numbers, symbols, sums, products and non-negative integer powers as a polynomial
qpoly(N, Vars) when is_integer(N) -> qpconst({N, 1}, Vars);
qpoly({rat, N, D}, Vars) -> qpconst(qnorm({N, D}), Vars);
qpoly({sym, X}, Vars) -> #{[case V of X -> 1; _ -> 0 end || V <- Vars] => {1, 1}};
qpoly({add, Es}, Vars) -> lists:foldl(fun(E, Acc) -> qpadd(Acc, qpoly(E, Vars)) end, #{}, Es);
qpoly({mul, Es}, Vars) -> lists:foldl(fun(E, Acc) -> qpmul(Acc, qpoly(E, Vars)) end, qpconst({1, 1}, Vars), Es);
qpoly({pow, B, K}, Vars) when is_integer(K), K >= 0 ->
    lists:foldl(fun(_, Acc) -> qpmul(Acc, qpoly(B, Vars)) end, qpconst({1, 1}, Vars), lists:seq(1, K)).

qpconst({0, _}, _) -> #{};
qpconst(Q, Vars) -> #{[0 || _ <- Vars] => Q}.

qpadd(A, B) ->
    maps:fold(fun(K, V, Acc) ->
                      case qadd(maps:get(K, Acc, {0, 1}), V) of
                          {0, _} -> maps:remove(K, Acc);
                          S -> Acc#{K => S}
                      end
              end, A, B).
qpmul(A, B) ->
    lists:foldl(fun(T, Acc) -> qpadd(Acc, T) end, #{},
                [#{[X + Y || {X, Y} <- lists:zip(Ka, Kb)] => qmul(Va, Vb)} || Ka := Va <- A, Kb := Vb <- B]).
qpscale(_, {0, _}) -> #{};
qpscale(A, Q) -> maps:map(fun(_, V) -> qmul(V, Q) end, A).

qadd({A, B}, {C, D}) -> qnorm({A * D + C * B, B * D}).

%%% The seeds: every non-zero sector within top, the positive indices summing to at most R, the others (numerators)
%%% to at least -S

%%% With R = {dots, D} instead: every sector's powers at most D above 1 in total (Kira's r lets a small sector carry
%%% many dots, which a reduction of the top sector rarely needs).

-spec seeds(family(), pos_integer() | {dots, non_neg_integer()}, non_neg_integer()) -> [integral()].
seeds(Fam = #{top := Top, props := Props}, R, S) ->
    N = length(Props),
    Budget = fun(Sec) -> case R of {dots, D} -> length(Sec) + D; _ -> R end end,
    Sectors = [Sec || Sec <- subsets(Top), Sec =/= [], length(Sec) =< Budget(Sec)],
    lists:append(
      [[list_to_tuple(fill(N, Sec, Pos, Others, Neg))
        || Pos <- compositions(Sec, Budget(Sec)),
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

-spec zero(family(), integral()) -> boolean().
zero(Fam, A) ->
    Sec = [J || {J, X} <- lists:enumerate(tuple_to_list(A)), X > 0],
    memo(zero, Fam, Sec, fun() -> scaleless(Fam, Sec) end).

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
symanzik(Fam = #{loops := Ls, props := Props, vars := Vars}) ->
    memo(symanzik, Fam, all, fun() ->
            N = length(Props),
            Point = point(Vars, {7, 7, 7}),
            KinAt = fun(A, B) -> {ok, V} = pq_symbolic:evaluate(kin(Fam, A, B), Point, #{exact => true}), modq(V) end,
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
                           _ -> error({not_supported, more_than_two_loops})
                       end,
            F = padd([pmul(Aij, JJ(Ji, Jj)) || {Ji, AdjRow} <- lists:zip(Js, Adj), {Jj, Aij} <- lists:zip(Js, AdjRow)]
                     ++ [pscale(pmul(J0, U), ?P - 1)]),
            {U, F}
    end).

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

%%% Sector symmetries (Pak's criterion): two sectors are equivalent if U + F of one is U + F of the other with the
%%% Feynman parameters permuted (compared mod p at a random point of the invariants, as for the zero sectors). The
%%% Feynman-parameter integral of an integral without numerators depends on U, F and the powers alone, so then
%%% I_S(a) = I_R(a permuted). Every such integral maps to the smallest of its images (over the equivalent sectors
%%% and each sector's own automorphisms); an integral with numerators stays as it is (its identities reduce it to
%%% ones without).

-spec canon(family(), integral()) -> integral().
canon(Fam, A) ->
    L = tuple_to_list(A),
    Sec = [J || {J, X} <- lists:enumerate(L), X > 0],
    case lists:any(fun(X) -> X < 0 end, L) orelse maps:get(Sec, classes(Fam), none) of
        true -> A;
        none -> A;
        {_Rep, Perms} -> lists:min([image(A, Sec, Perm, length(L)) || Perm <- Perms])
    end.

%% the integral with index j of the sector moved to Perm's image of j
image(A, Sec, Perm, N) ->
    M = maps:from_list([{R, element(J, A)} || {J, R} <- lists:zip(Sec, Perm)]),
    list_to_tuple([maps:get(K, M, 0) || K <- lists:seq(1, N)]).

%% Sector => {Representative, [Perm]}: Perm the images in the representative of the sector's indices, in order, for
%% every map that carries U + F of one onto the other (memoised)
classes(Fam = #{top := Top, props := Props}) ->
    memo(classes, Fam, all, fun() ->
            N = length(Props),
            {U, F} = symanzik(Fam),
            G = maps:merge(U, F),
            Sectors = lists:sort(fun(X, Y) -> {length(X), X} =< {length(Y), Y} end,
                                 [S || S <- subsets(Top), S =/= [], not zero(Fam, sector_point(N, S))]),
            {Classes, _} =
                lists:foldl(
                  fun(S, {Acc, Reps}) ->
                          GS = restrict(G, S),
                          Sig = signature(S, GS),
                          case first_match(S, GS, [Rp || {_, Sig1, _} = Rp <- Reps, Sig1 =:= Sig]) of
                              {Rep, Perms} -> {Acc#{S => {Rep, Perms}}, Reps};
                              none -> {Acc#{S => {S, matches(S, GS, S, GS)}}, Reps ++ [{S, Sig, GS}]}
                          end
                  end, {#{}, []}, Sectors),
            Classes
    end).

first_match(_, _, []) -> none;
first_match(S, GS, [{R, _, GR} | Rest]) ->
    case matches(S, GS, R, GR) of
        [] -> first_match(S, GS, Rest);
        Perms -> {R, Perms}
    end.

%% U + F with the parameters outside the sector set to 0, as exponents over the sector's indices
restrict(G, S) ->
    maps:from_list([{[lists:nth(J, E) || J <- S], C} || E := C <- G, C =/= 0,
                    lists:all(fun({J, X}) -> X =:= 0 orelse lists:member(J, S) end, lists:enumerate(E))]).

signature(S, GS) -> {length(S), map_size(GS), lists:sort(maps:values(GS))}.

%% the bijections S -> R carrying GS onto GR
matches(S, GS, R, GR) ->
    [Perm || Perm <- perms(R),
             begin
                 %% the exponent of R's k-th index is that of the S index mapped to it
                 Pos = maps:from_list(lists:zip(Perm, lists:seq(1, length(S)))),
                 maps:from_list([{[lists:nth(maps:get(Rk, Pos), E) || Rk <- R], C} || E := C <- GS]) =:= GR
             end].

perms([]) -> [[]];
perms(L) -> [[X | P] || X <- L, P <- perms(L -- [X])].

%%% Derivatives with respect to an invariant (for the masters' differential equations). An operator
%%% O_ij = p_i . d/dp_j (p_i, p_j external) acts on the integrand through the propagators, like an IBP identity; on
%%% a function of the scalar products q_ab = p_a . p_b it moves them by O_ij q_ab = [j = a] q_ib + [j = b] q_ai. A
%%% combination sum c_ij O_ij that moves every q_ab exactly as Var does (the others fixed, the legs on shell) is
%%% d/dVar: the c_ij solve a small linear system (exact, pq_recon:linear/4), the operators left out set to 0.

%% d I(A) / d Var as [{Integral, CoefficientAst}], zero integrals dropped, equal ones merged (the symmetries too)
-spec derivative(family(), integral(), atom()) -> [{integral(), term()}].
derivative(Fam, A, Var) ->
    Cs = operator(Fam, Var),
    Terms = lists:append([[{I, {mul, [C, X]}} || {I, X} <- ext_op(Fam, A, Pi, Pj)] || {{Pi, Pj}, C} <- Cs]),
    maps:to_list(lists:foldl(fun({I, C}, M) ->
                                     case zero(Fam, I) of
                                         true -> M;
                                         false -> maps:update_with(canon(Fam, I), fun(C0) -> {add, [C0, C]} end, C, M)
                                     end
                             end, #{}, Terms)).

%% p_i . d/dp_j of 1 / prod D^a: -a_k D_k^(-a_k-1) p_i . dD_k/dp_j, dD_k/dp_j = 2 Sign_k c_(k, p_j) Mom_k
ext_op(Fam = #{props := Props}, A, Pi, Pj) ->
    lists:append(
      [begin
           Ak = element(K, A),
           Raised = setelement(K, A, Ak + 1),
           Cj = lists:sum([C || {C, Q} <- Mom, Q =:= Pj]),
           [term(Raised, {mul, [-2 * Ak * Sign * Cj * Cq, X]}, Y) || Cj =/= 0, {Cq, Q} <- Mom, {X, Y} <- dot(Fam, Pi, Q)]
       end || {K, {Sign, Mom}} <- lists:enumerate(Props), element(K, A) =/= 0]).

%% [{{Pi, Pj}, CoefficientAst}]: sum c_ij O_ij = d/dVar on the scalar products of the externals (memoised)
operator(Fam = #{ext := Es, vars := Vars}, Var) ->
    memo(operator, Fam, Var, fun() ->
            Pairs = [{A, B} || {I, A} <- lists:enumerate(Es), {J, B} <- lists:enumerate(Es), I =< J],
            Ops = [{Pi, Pj} || Pi <- Es, Pj <- Es],
            Q = fun(A, B) -> kin(Fam, A, B) end,
            Act = fun({Pi, Pj}, {A, B}) ->
                          {add, [0] ++ [Q(Pi, B) || Pj =:= A] ++ [Q(A, Pi) || Pj =:= B]}
                  end,
            Rows = [[Act(O, P) || O <- Ops] || P <- Pairs],
            Rhs = [diff(Q(A, B), Var) || {A, B} <- Pairs],
            %% the operators kept: a set of as many as there are pairs, the matrix nonsingular at a random point
            Keep = independent(Rows, length(Pairs), Vars),
            A = [[lists:nth(K, Row) || K <- Keep] || Row <- Rows],
            {ok, Sol} = pq_recon:linear(A, Rhs, Vars, #{}),
            Cs = [{lists:nth(K, Ops), {mul, [pq_sym_norm:to_ast(N), {pow, pq_sym_norm:to_ast(D), -1}]}}
                  || {C, {N, D}} <- Sol, K <- [lists:nth(C + 1, Keep)], map_size(N) > 0],
            Cs
    end).

%% the first N columns (in order) that are independent at a random point mod p
independent(Rows, N, Vars) ->
    Point = point(Vars, {5, 5, 5}),
    Num = [[begin {ok, X} = pq_symbolic:evaluate(E, Point, #{exact => true}), modq(X) end || E <- Row] || Row <- Rows],
    Cols = length(hd(Num)),
    pick(Num, lists:seq(1, Cols), [], N).

pick(_, _, Chosen, N) when length(Chosen) =:= N -> lists:reverse(Chosen);
pick(_, [], Chosen, _) -> error({operators_dependent, lists:reverse(Chosen)});
pick(Num, [C | Cs], Chosen, N) ->
    Try = lists:reverse([C | Chosen]),
    Sub = [[lists:nth(K, Row) || K <- Try] || Row <- Num],
    case rank(Sub) =:= length(Try) of
        true -> pick(Num, Cs, [C | Chosen], N);
        false -> pick(Num, Cs, Chosen, N)
    end.

%% the rank of a matrix mod p (rows of residues)
rank(Rows) ->
    K = length(hd(Rows)),
    length([R || R <- reduce([R ++ [0] || R <- Rows], 1, K), lists:any(fun(X) -> X =/= 0 end, lists:sublist(R, K))]).

%% d/dVar of an AST built from numbers, symbols, sums and products
diff(N, _) when is_integer(N) -> 0;
diff({rat, _, _}, _) -> 0;
diff({sym, X}, X) -> 1;
diff({sym, _}, _) -> 0;
diff({add, Es}, V) -> {add, [diff(E, V) || E <- Es]};
diff({mul, Es}, V) ->
    {add, [{mul, [diff(E, V) | [F || {J, F} <- lists:enumerate(Es), J =/= I]]} || {I, E} <- lists:enumerate(Es)]}.

