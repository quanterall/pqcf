%% Integration-by-parts (IBP) identities for a family of Feynman integrals (decision record R1): the scalar products
%% solved from the propagators' definitions, the scaleless (zero) sectors by the Lee-Pak criterion, the sector
%% symmetries by Pak's criterion, the seeds in Kira's form (or with the dots per sector bounded), and derivatives
%% with respect to the invariants (for the masters' differential equations). pq_ibp builds and solves the systems.
%%
%% A family (new/1): loops and ext, the loop and the independent external momenta (atoms); props, the propagators
%% [{Sign, Momentum} | {Sign, Momentum, MassSquared}], Momentum = [{IntegerCoefficient, q}] over them, D_a =
%% Sign_a * (Momentum_a^2 - MassSquared_a), the mass squared an AST in the variables, 0 if absent (numerators
%% included: they must span every scalar product with a loop momentum, as many as there are); kin, the external
%% scalar products #{{Pa, Pb} => Ast} in the invariants (a missing pair is 0); vars, the variables of the results (the
%% dimension d first, then the invariants); top, the indices that can be positive (default: all); symmetries
%% (default true). The integral I(A) = int prod D_a^-A_a for a tuple A.
%%
%% Results are cached per process by the family's id (its hash): the scalar products, the Symanzik polynomials, the
%% zero sectors, the symmetry classes and the derivative operators. Random points are drawn from local states,
%% never the caller's.
-module(pq_ibp_family).
-export([new/1, example/1, identities/2, zero/2, seeds/3, canon/2, derivative/3, vars/1, size/1, kin/3, numerator/3]).
%% for pq_ibp: the identities with their coefficients as polynomials (Exponents => {N, D})
-export([identities_poly/2, warm/2, load/1, expand/2, mapped/2, without_sector_maps/1, qpone/1, qpneg/1]).
-export_type([family/0, integral/0]).

-type family() :: #{id := integer(), loops := [atom()], ext := [atom()], props := [{integer(), [{integer(), atom()}]}],
                    masses := [term()],
                    kin := #{{atom(), atom()} => term()}, vars := [atom()], top := [pos_integer()], symmetries := boolean(),
                    sector_maps => boolean(), eps => #{pos_integer() => pos_integer()}}.
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
    Prop = fun({S, Mom}) -> {S, Mom, 0}; ({S, Mom, M2}) -> {S, Mom, M2}; (X) -> X end,
    Props3 = [Prop(P) || P <- Props],
    ok = check(lists:all(fun({S, Mom, _}) -> (S =:= 1 orelse S =:= -1) andalso Mom =/= []
                                                 andalso lists:all(fun({C, Q}) -> is_integer(C) andalso lists:member(Q, Moms) end, Mom);
                            (_) -> false end, Props3), props),
    ok = check(lists:all(fun is_atom/1, Vars) andalso lists:member(d, Vars), vars),
    Kin = maps:from_list([{pair_key(A, B), X} || {A, B} := X <- Kin0]),
    ok = check(lists:all(fun({A, B}) -> lists:member(A, Ext) andalso lists:member(B, Ext) end, maps:keys(Kin)), kin),
    %% eps => #{Line => M}: the line's power is its index plus M epsilon (d = 4 - 2 epsilon), as a one-loop insertion
    %% leaves it (Forcer's notation, https://arxiv.org/abs/1704.06650); such a line is in every sector, and the
    %% symmetries, which would have to respect the insertions, are off
    Eps = maps:get(eps, Spec, #{}),
    Top = maps:get(top, Spec, lists:seq(1, N)),
    ok = check(lists:all(fun({J, M}) -> lists:member(J, Top) andalso is_integer(M) andalso M > 0 end, maps:to_list(Eps)), eps),
    Fam0 = #{loops => Loops, ext => Ext, props => [{S, Mom} || {S, Mom, _} <- Props3], masses => [M2 || {_, _, M2} <- Props3],
             kin => Kin, vars => Vars, eps => Eps,
             top => Top, symmetries => maps:get(symmetries, Spec, true) andalso map_size(Eps) =:= 0},
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

%% The family's caches as computed in this process, after computing what identities over these seeds need (the
%% scalar products, the Symanzik polynomials, every sector's zero test and symmetry class): for load/1 in the
%% processes that generate identities in parallel, so that none computes them again
-spec warm(family(), [integral()]) -> [{term(), term()}].
warm(Fam = #{id := Id, top := Top, props := Props}, _Seeds) ->
    N = length(Props),
    _ = sps(Fam),
    _ = classes(Fam),
    _ = sector_maps(Fam),
    _ = [zero(Fam, sector_point(N, S)) || S <- subsets(lists:seq(1, N))],
    _ = [{zero(Fam, sector_point(N, S))} || S <- subsets(Top)],
    [{K, V} || {{?MODULE, _Tag, Id2, _} = K, V} <- get(), Id2 =:= Id].

%% the caches warm/2 returned, in this process
-spec load([{term(), term()}]) -> ok.
load(Entries) ->
    lists:foreach(fun({K, V}) -> put(K, V) end, Entries).

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
sps(Fam = #{props := Props, masses := Masses}) ->
    memo(sps, Fam, all, fun() ->
            Keys = sp_keys(Fam),
            %% D_a = sum_k M_ak sp_k + c_a, the mass in the constant: c_a = Sign_a (P_a^2 - m_a^2)
            Rows = [begin
                        Terms = [{Sign * C1 * C2, pair(Fam, Q1, Q2)} || {C1, Q1} <- Mom, {C2, Q2} <- Mom],
                        {[lists:sum([C || {C, {sp, K1}} <- Terms, K1 =:= K]) || K <- Keys],
                         {add, [0, {mul, [-Sign, M2]} | [{mul, [C, X]} || {C, {const, X}} <- Terms]]}}
                    end || {{Sign, Mom}, M2} <- lists:zip(Props, Masses)],
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

%% The integral A times a product of scalar products u.v (u, v momenta as [{Integer, Momentum}]), as integrals of the
%% family: [{Integral, CoefficientAst}], each scalar product written in the D's (a numerator D_j lowers index j)
-spec numerator(family(), integral(), [{[{integer(), atom()}], [{integer(), atom()}]}]) -> [{integral(), term()}].
numerator(Fam, A, Sps) ->
    Expand = fun({U, V}) -> [{{mul, [Cu * Cv, X]}, Y} || {Cu, Qu} <- U, {Cv, Qv} <- V, {X, Y} <- dot(Fam, Qu, Qv)] end,
    Terms = lists:foldl(fun(Sp, Acc) ->
                                [{{mul, [C, X]}, case Y of const -> I; _ -> setelement(Y, I, element(Y, I) - 1) end}
                                 || {C, I} <- Acc, {X, Y} <- Expand(Sp)]
                        end, [{1, A}], Sps),
    Merged = lists:foldl(fun({C, I}, M) -> maps:update_with(I, fun(C0) -> {add, [C0, C]} end, C, M) end, #{}, Terms),
    [{I, C} || I := C0 <- Merged, {ok, C} <- [pq_symbolic:canonicalise(C0)], C =/= 0].

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
    Eps = maps:get(eps, Fam, #{}),
    Div = case L =:= V of true -> [{A, d()}]; false -> [] end,
    Deriv = lists:append(
              [begin
                   Aj = element(J, A),
                   Raised = setelement(J, A, Aj + 1),
                   CL = lists:sum([C || {C, Q} <- Mom, Q =:= L]),
                   %% the power a_j: the index, plus M epsilon on a line with an insertion
                   Pow = case Eps of
                             #{J := M} -> {add, [Aj, {mul, [M, {rat, 1, 2}, {add, [4, {mul, [-1, d()]}]}]}]};
                             _ -> Aj
                         end,
                   %% -a_j D_j^(-a_j-1) dD_j/dl.v, dD_j/dl = 2 Sign_j CL Mom_j
                   [term(Raised, {mul, [-2 * Sign * CL * Cq, Pow, X]}, Y)
                    || CL =/= 0, {Cq, Q} <- Mom, {X, Y} <- dot(Fam, V, Q)]
               end || {J, {Sign, Mom}} <- lists:enumerate(Props), element(J, A) =/= 0 orelse is_map_key(J, Eps)]),
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
identities_poly(Fam = #{loops := Ls, ext := Es, vars := Vars}, Seeds0) ->
    One = qpconst({1, 1}, Vars),
    %% an integral's name in the equations: [] if it is zero, its expansion if its sector maps onto another one
    %% (numerators included), else itself under the scalar symmetries (canon/2)
    Name = fun(I) -> memo(name, Fam, I, fun() -> case zero(Fam, I) of
                                                      true -> [];
                                                      false -> expand_or_canon(Fam, I, One)
                                                  end end) end,
    Seeds = Seeds0,
    Eqs = [lists:foldl(fun({I, C}, M) ->
                               lists:foldl(fun({J, Cj}, M2) ->
                                                   Cc = qpmul(C, Cj),
                                                   case M2 of
                                                       #{J := C0} -> M2#{J => qpadd(C0, Cc)};
                                                       _ -> M2#{J => Cc}
                                                   end
                                           end, M, Name(I))
                       end, #{}, identity_poly(Fam, A, L, V))
           || A <- Seeds, L <- Ls, V <- Ls ++ Es],
    [maps:filter(fun(_, C) -> map_size(C) > 0 end, E) || E <- Eqs].

expand_or_canon(Fam = #{symmetries := true}, I, One) ->
    case mapped(Fam, I) of
        true -> expand(Fam, I);
        false -> [{canon(Fam, I), One}]
    end;
expand_or_canon(_Fam, I, One) -> [{I, One}].

identity_poly(Fam = #{props := Props, vars := Vars}, A, L, V) ->
    Eps = maps:get(eps, Fam, #{}),
    Div = case L =:= V of
              true -> [{A, #{[case X of d -> 1; _ -> 0 end || X <- Vars] => {1, 1}}}];
              false -> []
          end,
    Deriv = lists:append(
              [begin
                   Aj = element(J, A),
                   Raised = setelement(J, A, Aj + 1),
                   CL = lists:sum([C || {C, Q} <- Mom, Q =:= L]),
                   %% the power a_j: the index, plus M epsilon = M (4 - d) / 2 on a line with an insertion
                   Pow = case Eps of
                             #{J := M} -> qpadd(qpconst({Aj + 2 * M, 1}, Vars),
                                                #{[case X of d -> 1; _ -> 0 end || X <- Vars] => {-M, 2}});
                             _ -> qpconst({Aj, 1}, Vars)
                         end,
                   [term(Raised, qpscale(qpmul(X, Pow), {-2 * Sign * CL * Cq, 1}), Y)
                    || CL =/= 0, map_size(Pow) > 0, {Cq, Q} <- Mom, {X, Y} <- dot_poly(Fam, V, Q)]
               end || {J, {Sign, Mom}} <- lists:enumerate(Props), element(J, A) =/= 0 orelse is_map_key(J, Eps)]),
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

%% the polynomial 1 in the family's variables, and a polynomial's negative (for pq_ibp)
qpone(#{vars := Vars}) -> qpconst({1, 1}, Vars).
qpneg(P) -> qpscale(P, {-1, 1}).

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
    %% a line with an insertion is in every sector
    Ins = maps:keys(maps:get(eps, Fam, #{})),
    Sectors = [Sec || Sec <- subsets(Top), Sec =/= [], length(Sec) =< Budget(Sec), Ins -- Sec =:= []],
    Plain = lists:append(
              [[list_to_tuple(fill(N, Sec, Pos, Others, Neg))
                || Pos <- compositions(Sec, Budget(Sec)),
                   Neg <- numerators(length(Others), S)]
               || Sec <- Sectors, Others <- [[J || J <- lists:seq(1, N), not lists:member(J, Sec)]],
                  not zero(Fam, sector_point(N, Sec))]),
    %% a line with an insertion stays in the sector at any integer part: also the seeds with it lowered, down to -S
    lists:usort(lists:foldl(fun(J, Acc) ->
                                    [setelement(J, T, X) || T <- Acc, X <- lists:seq(-S, element(J, T))]
                            end, Plain, Ins)).

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
    Eps = maps:get(eps, Fam, #{}),
    Sec = [J || {J, X} <- lists:enumerate(tuple_to_list(A)), X > 0 orelse is_map_key(J, Eps)],
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
symanzik(Fam = #{loops := Ls, props := Props, masses := Masses, vars := Vars}) ->
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
            MassAt = fun(M2) -> {ok, V} = pq_symbolic:evaluate(M2, Point, #{exact => true}), modq(V) end,
            J0 = padd([X(J, Sign * C1 * C2 * KinAt(Q1, Q2)) || {J, {Sign, Mom}} <- lists:enumerate(Props),
                                                              {C1, Q1} <- Ext(Mom), {C2, Q2} <- Ext(Mom)]
                      ++ [X(J, -Sign * MassAt(M2)) || {J, {{Sign, _}, M2}} <- lists:enumerate(lists:zip(Props, Masses)), M2 =/= 0]),
            JJ = fun(Ji, Jj) -> padd([pscale(pmul(P1, P2), KinAt(Q1, Q2)) || {P1, Q1} <- Ji, {P2, Q2} <- Jj]) end,
            One = #{lists:duplicate(N, 0) => 1},
            U = det(M, One),
            Adj = adjugate(M, One),
            F = padd([pmul(Aij, JJ(Ji, Jj)) || {Ji, AdjRow} <- lists:zip(Js, Adj), {Jj, Aij} <- lists:zip(Js, AdjRow)]
                     ++ [pscale(pmul(J0, U), ?P - 1)]),
            {U, F}
    end).

%% the determinant and the adjugate of a square matrix of polynomials mod p (Laplace expansion: the matrices are
%% L x L, L the number of loops)
det([[A]], _One) -> A;
det([Row | _] = M, One) ->
    padd([pscale(pmul(A, det(minor(M, 1, J), One)), case J rem 2 of 1 -> 1; 0 -> ?P - 1 end)
          || {J, A} <- lists:enumerate(Row), map_size(A) > 0]).

adjugate([[_]], One) -> [[One]];
adjugate(M, One) ->
    L = length(M),
    [[pscale(det(minor(M, J, I), One), case (I + J) rem 2 of 0 -> 1; 1 -> ?P - 1 end) || J <- lists:seq(1, L)]
     || I <- lists:seq(1, L)].

%% M without row I and column J
minor(M, I, J) ->
    [[X || {Jx, X} <- lists:enumerate(Row), Jx =/= J] || {Ix, Row} <- lists:enumerate(M), Ix =/= I].

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
canon(#{symmetries := false}, A) -> A;
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

%% the bijections S -> R carrying GS onto GR, by backtracking: index i of S can go only to an index k of R with the
%% same invariant (the coefficients and exponents of the terms containing it), and each assignment must keep the
%% pairwise invariants with those made before; a complete assignment is then checked on every term. The same maps
%% as trying every permutation (as this did, 11! for a four-loop top sector), found without enumerating them.
matches(S, GS, R, GR) ->
    N = length(S),
    Ts = [{list_to_tuple(E), C} || E := C <- GS],
    Tr = [{list_to_tuple(E), C} || E := C <- GR],
    One = fun(Terms, I) -> lists:sort([{C, element(I, E)} || {E, C} <- Terms, element(I, E) > 0]) end,
    Two = fun(Terms, I, J) -> lists:sort([{C, element(I, E), element(J, E)} || {E, C} <- Terms, element(I, E) + element(J, E) > 0]) end,
    InvR = maps:from_list([{K, One(Tr, K)} || K <- lists:seq(1, N)]),
    Cands = [{I, [K || K <- lists:seq(1, N), maps:get(K, InvR) =:= One(Ts, I)]} || I <- lists:seq(1, N)],
    case lists:any(fun({_, Ks}) -> Ks =:= [] end, Cands) of
        true -> [];
        false ->
            Order = [I || {_, I} <- lists:sort([{length(Ks), I} || {I, Ks} <- Cands])],
            CandMap = maps:from_list(Cands),
            PairS = fun(I, J) -> memo_pair({s, I, J}, fun() -> Two(Ts, I, J) end) end,
            PairR = fun(K, L) -> memo_pair({r, K, L}, fun() -> Two(Tr, K, L) end) end,
            erase_pairs(),
            Found = assign(Order, CandMap, [], PairS, PairR),
            erase_pairs(),
            Rt = list_to_tuple(R),
            [Perm || Asg <- Found,
                     Perm <- [[element(maps:get(I, maps:from_list(Asg)), Rt) || I <- lists:seq(1, N)]],
                     begin
                         Pos = maps:from_list(lists:zip(Perm, lists:seq(1, N))),
                         Idx = [maps:get(Rk, Pos) || Rk <- R],
                         same_terms([{E, C} || {E, C} <- Ts], Idx, GR)
                     end]
    end.

assign([], _, Asg, _, _) -> [Asg];
assign([I | Rest], Cands, Asg, PairS, PairR) ->
    Used = [K || {_, K} <- Asg],
    lists:append([assign(Rest, Cands, [{I, K} | Asg], PairS, PairR)
                  || K <- maps:get(I, Cands), not lists:member(K, Used),
                     lists:all(fun({J, L}) -> PairS(I, J) =:= PairR(K, L) end, Asg)]).

%% the pair invariants of one matches/4 call, cached in the process dictionary for its duration
memo_pair(Key, Fun) ->
    case get({pq_ibp_pair, Key}) of
        undefined -> V = Fun(), put({pq_ibp_pair, Key}, V), V;
        V -> V
    end.

erase_pairs() -> _ = [erase(K) || {{pq_ibp_pair, _} = K, _} <- get()], ok.

%% every term, its exponents permuted, a term of GR with the same coefficient (the sizes are equal: the signature)
same_terms([], _, _) -> true;
same_terms([{Et, C} | Rest], Idx, GR) ->
    case maps:find([element(I, Et) || I <- Idx], GR) of
        {ok, C} -> same_terms(Rest, Idx, GR);
        _ -> false
    end.


%%% Derivatives with respect to an invariant (for the masters' differential equations). An operator
%%% O_ij = p_i . d/dp_j (p_i, p_j external) acts on the integrand through the propagators, like an IBP identity; on
%%% a function of the scalar products q_ab = p_a . p_b it moves them by O_ij q_ab = [j = a] q_ib + [j = b] q_ai. A
%%% combination sum c_ij O_ij that moves every q_ab exactly as Var does (the others fixed, the legs on shell) is
%%% d/dVar: the c_ij solve a small linear system (exact, pq_recon:linear/4), the operators left out set to 0.

%% d I(A) / d Var as [{Integral, CoefficientAst}], zero integrals dropped, equal ones merged (the symmetries too)
-spec derivative(family(), integral(), atom()) -> [{integral(), term()}].
derivative(Fam = #{props := Props, masses := Masses}, A, Var) ->
    Cs = operator(Fam, Var),
    Terms0 = lists:append([[{I, {mul, [C, X]}} || {I, X} <- ext_op(Fam, A, Pi, Pj)] || {{Pi, Pj}, C} <- Cs]),
    %% the masses: d D_k^-a / dVar = a Sign_k (d m_k^2 / dVar) D_k^(-a-1)
    MassTerms = [{setelement(K, A, Ak + 1), {mul, [Ak * Sign, Dm]}}
                 || {K, {{Sign, _}, M2}} <- lists:enumerate(lists:zip(Props, Masses)), M2 =/= 0,
                    Ak <- [element(K, A)], Ak =/= 0, Dm <- [diff(M2, Var)], Dm =/= 0],
    Terms = Terms0 ++ MassTerms,
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
operator(#{ext := []}, _Var) -> [];
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
diff({pow, B, K}, V) when is_integer(K) -> {mul, [K, {pow, B, K - 1}, diff(B, V)]};
diff({mul, Es}, V) ->
    {add, [{mul, [diff(E, V) | [F || {J, F} <- lists:enumerate(Es), J =/= I]]} || {I, E} <- lists:enumerate(Es)]}.

%%% Sector maps: the changes of loop momenta behind the symmetries (as in LiteRed, Kira and FIRE). For a sector S
%%% equivalent to its representative R by a permutation of the propagators (Pak's criterion, classes/1), a change
%%% k -> A k + B p of the loop momenta (A integer, det A = +-1, so the measure is unchanged) with
%%% Momentum_i(A k + B p) = +-Momentum_sigma(i)(k) for every propagator i of S carries S's integrand onto R's. Under
%%% it a numerator D_j becomes a linear combination of the family's D's, so every integral of S, numerators
%%% included, equals a combination of integrals of R and its subsectors (expand/2). Such a sector needs no identities
%%% of its own. A symmetry with no such change (it would need to move the external momenta) keeps the scalar
%%% mapping only.

%% the family with the scalar symmetries only (no sector expanded into another), as a family of its own (its own
%% caches)
-spec without_sector_maps(family()) -> family().
without_sector_maps(Fam = #{id := Id}) -> Fam#{sector_maps => false, id => erlang:phash2({without_sector_maps, Id})}.

%% does the integral's sector map onto another one?
-spec mapped(family(), integral()) -> boolean().
mapped(#{sector_maps := false}, _) -> false;
mapped(Fam = #{symmetries := true}, A) ->
    is_map_key([J || {J, X} <- lists:enumerate(tuple_to_list(A)), X > 0], sector_maps(Fam));
mapped(_, _) -> false.

%% Sector => {Representative, Perm, {A, B}} for the sectors that map onto their representative (memoised)
sector_maps(Fam) ->
    memo(sector_maps, Fam, all, fun() ->
            maps:from_list([{S, {R, Perm, T}} || S := {R, Perms} <- classes(Fam), R =/= S,
                                                 {Perm, T} <- [first_transform(Fam, S, Perms)], T =/= none])
    end).

first_transform(_Fam, _S, []) -> {none, none};
first_transform(Fam, S, [Perm | Perms]) ->
    case first_ext(Fam, S, Perm, ext_maps(Fam)) of
        none -> first_transform(Fam, S, Perms);
        T -> {Perm, T}
    end.

first_ext(_Fam, _S, _Perm, []) -> none;
first_ext(Fam, S, Perm, [C | Cs]) ->
    case transform(Fam, S, Perm, C) of
        none -> first_ext(Fam, S, Perm, Cs);
        {A, B} -> {A, B, C}
    end.

%% The maps of the external momenta that keep every external scalar product (squares included) the same polynomial
%% in the invariants: under one, every integral of the family keeps its value. Each p_e goes to +-p_f or to
%% +-(the dependent momentum -sum p), as a row of coefficients over the independent ones; the identity first.
%% (memoised)
ext_maps(Fam = #{ext := Es, vars := Vars}) ->
    memo(ext_maps, Fam, all, fun() ->
            E = length(Es),
            Unit = fun(F) -> [case G of F -> 1; _ -> 0 end || G <- lists:seq(1, E)] end,
            Images = lists:append([[Unit(F), [-X || X <- Unit(F)]] || F <- lists:seq(1, E)])
                     ++ case E > 1 of true -> [lists:duplicate(E, -1), lists:duplicate(E, 1)]; false -> [] end,
            Identity = [Unit(F) || F <- lists:seq(1, E)],
            Kin = fun(Ra, Rb) ->
                          qpoly({add, [0 | [{mul, [Ca * Cb, kin(Fam, Pa, Pb)]} || {Ca, Pa} <- lists:zip(Ra, Es), Ca =/= 0,
                                                                                  {Cb, Pb} <- lists:zip(Rb, Es), Cb =/= 0]]}, Vars)
                  end,
            Keeps = fun(C) ->
                            lists:all(fun({I, J}) -> Kin(lists:nth(I, C), lists:nth(J, C)) =:= Kin(lists:nth(I, Identity), lists:nth(J, Identity)) end,
                                      [{I, J} || I <- lists:seq(1, E), J <- lists:seq(I, E)])
                    end,
            Cands = [C || C <- choices(E, Images), C =/= Identity, Keeps(C)],
            [Identity | Cands]
    end).

choices(0, _) -> [[]];
choices(K, Images) -> [[I | Rest] || I <- Images, Rest <- choices(K - 1, Images)].

%% {A, B} with Momentum_i(A k + B p, C p) = eps_i Momentum_sigma(i)(k, p) for i in S, some signs eps, or none (C
%% the externals' map, rows of coefficients)
transform(#{loops := Ls, ext := Es, props := Props}, S, Perm, Cmap) ->
    Coef = fun(Mom, Q) -> lists:sum([C || {C, Q1} <- Mom, Q1 =:= Q]) end,
    Mom = fun(J) -> element(2, lists:nth(J, Props)) end,
    L = length(Ls), E = length(Es),
    %% unknowns: a_lj (l, j loops) then b_le (l loop, e external), row-major
    Ua = fun(Li, Lj) -> (Li - 1) * L + Lj end,
    Ub = fun(Li, Ee) -> L * L + (Li - 1) * E + Ee end,
    NU = L * L + L * E,
    Pairs = lists:zip(S, Perm),
    try_signs(Pairs, [], fun(Signs) ->
            Rows = lists:append(
                     [[row(NU, [{Ua(Li, Lj), Coef(Mom(I), lists:nth(Li, Ls))} || Li <- lists:seq(1, L)],
                           Eps * Coef(Mom(R), lists:nth(Lj, Ls)))
                       || Lj <- lists:seq(1, L)]
                      ++ [row(NU, [{Ub(Li, Ee), Coef(Mom(I), lists:nth(Li, Ls))} || Li <- lists:seq(1, L)],
                              Eps * Coef(Mom(R), lists:nth(Ee, Es))
                              - lists:sum([Coef(Mom(I), lists:nth(Ef, Es)) * lists:nth(Ee, lists:nth(Ef, Cmap)) || Ef <- lists:seq(1, E)]))
                          || Ee <- lists:seq(1, E)]
                      || {{I, R}, Eps} <- lists:zip(Pairs, Signs)]),
            %% where the conditions leave freedom: the identity (a = 1 on the diagonal, the rest 0)
            Default = [case U =< L * L andalso (U - 1) div L =:= (U - 1) rem L of true -> 1; false -> 0 end || U <- lists:seq(1, NU)],
            case solve_q(Rows, NU, Default) of
                none -> none;
                Sol ->
                    case lists:all(fun({_, 1}) -> true; (_) -> false end, Sol) of
                        false -> none;
                        true ->
                            Ints = [N || {N, 1} <- Sol],
                            A = [lists:sublist(Ints, (Li - 1) * L + 1, L) || Li <- lists:seq(1, L)],
                            B = [lists:sublist(Ints, L * L + (Li - 1) * E + 1, E) || Li <- lists:seq(1, L)],
                            case abs(idet(A)) of 1 -> {A, B}; _ -> none end
                    end
            end
    end).

%% every choice of signs, the first with a result (the first sign +1: k -> -k covers the other half)
try_signs([], Signs, F) -> F(lists:reverse(Signs));
try_signs([_ | Rest], [], F) -> try_signs(Rest, [1], F);
try_signs([_ | Rest], Signs, F) ->
    case try_signs(Rest, [1 | Signs], F) of
        none -> try_signs(Rest, [-1 | Signs], F);
        T -> T
    end.

row(NU, Coefs, Rhs) ->
    M = lists:foldl(fun({U, C}, Acc) -> Acc#{U => maps:get(U, Acc, 0) + C} end, #{}, Coefs),
    [{maps:get(U, M, 0), 1} || U <- lists:seq(1, NU)] ++ [{Rhs, 1}].

%% the solution of the rational system (rows [coefficients | rhs]), free unknowns at their defaults; none if
%% inconsistent
solve_q(Rows, NU, Default) ->
    {Red, Pivots} = rref([[qnorm(X) || X <- R] || R <- Rows], 1, NU, [], []),
    case lists:any(fun(R) -> lists:all(fun(X) -> X =:= {0, 1} end, lists:sublist(R, NU)) andalso lists:last(R) =/= {0, 1} end, Red) of
        true -> none;
        false ->
            Free = [U || U <- lists:seq(1, NU), not lists:member(U, Pivots)],
            Vals0 = maps:from_list([{U, {lists:nth(U, Default), 1}} || U <- Free]),
            Vals = lists:foldl(fun({R, U}, V) ->
                                       %% x_U = rhs - sum over free unknowns
                                       X = lists:foldl(fun(F, Acc) -> qsub(Acc, qmul(lists:nth(F, R), maps:get(F, V))) end,
                                                       lists:last(R), Free),
                                       V#{U => X}
                               end, Vals0, lists:zip(lists:sublist(Red, length(Pivots)), Pivots)),
            [maps:get(U, Vals) || U <- lists:seq(1, NU)]
    end.

%% reduced row echelon form over the rationals: {Rows (pivot rows first, in pivot order), Pivot columns}
rref(Rows, C, NU, Done, Pivots) when C > NU -> {lists:reverse(Done) ++ Rows, lists:reverse(Pivots)};
rref(Rows, C, NU, Done, Pivots) ->
    case lists:splitwith(fun(R) -> lists:nth(C, R) =:= {0, 1} end, Rows) of
        {_, []} -> rref(Rows, C + 1, NU, Done, Pivots);
        {Before, [Piv | After]} ->
            Inv = qinv(lists:nth(C, Piv)),
            PivN = [qmul(Inv, X) || X <- Piv],
            Elim = fun(R) -> F = lists:nth(C, R), case F of {0, 1} -> R; _ -> [qsub(X, qmul(F, Y)) || {X, Y} <- lists:zip(R, PivN)] end end,
            rref([Elim(R) || R <- Before ++ After], C + 1, NU, [PivN | [Elim(R) || R <- Done]], [C | Pivots])
    end.

idet([[A]]) -> A;
idet([Row | _] = M) ->
    lists:sum([case J rem 2 of 1 -> 1; 0 -> -1 end * X * idet(minor(M, 1, J)) || {J, X} <- lists:enumerate(Row), X =/= 0]).

%% An integral of a mapped sector as a combination of integrals that need no mapping: [{Integral, Poly}], zero ones
%% dropped, scalar ones canonical (memoised)
-spec expand(family(), integral()) -> [{integral(), #{[non_neg_integer()] => {integer(), pos_integer()}}}].
expand(Fam = #{props := Props, masses := Masses, vars := Vars}, A) ->
    memo(expand, Fam, A, fun() ->
            One = qpconst({1, 1}, Vars),
            N = length(Props),
            Sec = [J || {J, X} <- lists:enumerate(tuple_to_list(A)), X > 0],
            case maps:find(Sec, sector_maps(Fam)) of
                error -> [{canon(Fam, A), One}];
                {ok, {_R, Perm, T}} ->
                    Base = maps:from_list(lists:zip(Perm, [element(J, A) || J <- Sec])),
                    %% each numerator D_j(T k) as a linear form in the D's: #{j | const => Poly}
                    Factors = lists:append([lists:duplicate(-element(J, A), moved(Fam, J, T, lists:nth(J, Props), lists:nth(J, Masses)))
                                            || J <- lists:seq(1, N), element(J, A) < 0]),
                    Terms = lists:foldl(fun(F, Acc) ->
                                                [{Lower#{K => maps:get(K, Lower, 0) + 1}, qpmul(C, Ck)} || {Lower, C} <- Acc, {K, Ck} <- maps:to_list(F), K =/= const]
                                                ++ [{Lower, qpmul(C, Cc)} || {Lower, C} <- Acc, #{const := Cc} <- [F]]
                                        end, [{#{}, One}], Factors),
                    Ints = [{list_to_tuple([maps:get(K, Base, 0) - maps:get(K, Lower, 0) || K <- lists:seq(1, N)]), C} || {Lower, C} <- Terms],
                    %% each integral named in turn (zero, mapped again, or canonical), equal ones merged
                    Merged = lists:foldl(fun({I, C}, M) ->
                                                 Named = case zero(Fam, I) of
                                                             true -> [];
                                                             false -> case mapped(Fam, I) of
                                                                          true -> expand(Fam, I);
                                                                          false -> [{canon(Fam, I), One}]
                                                                      end
                                                         end,
                                                 lists:foldl(fun({J, Cj}, M2) ->
                                                                     Cc = qpmul(C, Cj),
                                                                     case M2 of #{J := C0} -> M2#{J => qpadd(C0, Cc)}; _ -> M2#{J => Cc} end
                                                             end, M, Named)
                                         end, #{}, Ints),
                    [{I, C} || I := C <- Merged, map_size(C) > 0]
            end
    end).

%% D_j(A k + B p) as #{j | const => Poly}: the momentum moved, squared through the scalar products, the mass
moved(Fam = #{loops := Ls, ext := Es, vars := Vars}, _J, {A, B, C}, {Sign, Mom}, M2) ->
    %% the moved momentum as [{Coefficient, q}]: each loop l -> sum_j a_lj k_j + sum_e b_le p_e
    Moved = lists:append([case lists:member(Q, Ls) of
                              true ->
                                  Li = index(Q, Ls),
                                  [{Co * X, Lj} || {Lj, X} <- lists:zip(Ls, lists:nth(Li, A)), X =/= 0]
                                  ++ [{Co * X, Pe} || {Pe, X} <- lists:zip(Es, lists:nth(Li, B)), X =/= 0];
                              false ->
                                  %% an external momentum goes to its image under the externals' map
                                  [{Co * X, Pe} || {Pe, X} <- lists:zip(Es, lists:nth(index(Q, Es), C)), X =/= 0]
                          end || {Co, Q} <- Mom]),
    Terms = [{qpscale(P, {Sign * C1 * C2, 1}), K} || {C1, Q1} <- Moved, {C2, Q2} <- Moved, {P, K} <- dot_poly(Fam, Q1, Q2)]
            ++ [{qpscale(qpoly(M2, Vars), {-Sign, 1}), const} || M2 =/= 0],
    lists:foldl(fun({P, K}, Acc) -> Acc#{K => qpadd(maps:get(K, Acc, #{}), P)} end, #{}, Terms).

