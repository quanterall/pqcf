%% The symbolic engine (FS 4.2, 5.2; SY). Pure: no process state, callable without the application started.
%% Every operation returns the canonical form of SY 3-4, so no simplify call is needed after another one.
%% Phase 0 scope: construction, canonicalise/simplify, subs, evaluate, shape, the matrix operations (expand
%% mode), printing and parsing, parameter resolution, bounded execution. Differentiation, equiv, compile and
%% the circuit-level functions come in later phases.
-module(pq_symbolic).
-compile({no_auto_import, [size/1]}).
-export([sym/1, num/1, i/0, pi/0, add/1, mul/1, sub/2, neg/1, divide/2, pow/2, sqrt/1, sin/1, cos/1, exp/1,
         conj/1, matrix/1,
         canonicalise/1, canonicalise/2, simplify/1, simplify/2, subs/2, subs/3, evaluate/2, evaluate/3,
         differentiate/2, differentiate/3, gradient/2, gradient/3, jacobian/2, jacobian/3,
         matmul/2, kron/2, smul/2, madd/1, dagger/1, transpose/1,
         shape/1, free_symbols/1, size/1, hash/1, domain_of/2,
         pretty/1, pretty/2, to_latex/1, parse/1, parse/2,
         resolve_params/2, resolve_params/3, circuit_unitary/1, circuit_unitary/2, gaussian_evolve/2, is_unitary/1, is_unitary/2, permanent/1, bs_probability/3, compile/2, compile/3, glynn_reference/2, eval_slp/2, eval_batch/2, equiv/3,
         substitute/2, gcd/2, gcd/3, factor/1, factor/2, cancel/1, cancel/2, numer_denom/1, numer_denom/2,
         series/3, series/4, limit/3, limit/4]).

-export_type([sym_expr/0, binding/0, opts/0]).

-type sym_expr() :: term().
-type binding() :: #{atom() => integer() | float() | {rat, integer(), pos_integer()} | {complex, term(), term()}}.
-type opts() :: map().
-type result(T) :: {ok, T} | {error, pq_core:reason()}.

%%% Construction (plain term builders)

sym(X) when is_atom(X) -> {sym, X}.
num(N) when is_number(N) -> N;
num({rat, N, D} = R) when is_integer(N), is_integer(D), D > 0 -> R.
i() -> {complex, 0, 1}.
pi() -> {const, pi}.
add(L) -> {add, L}.
mul(L) -> {mul, L}.
sub(A, B) -> {add, [A, {mul, [-1, B]}]}.
neg(A) -> {mul, [-1, A]}.
divide(A, B) -> {mul, [A, {pow, B, -1}]}.
pow(A, E) -> {pow, A, E}.
sqrt(A) -> {pow, A, {rat, 1, 2}}.
sin(A) -> {sin, A}.
cos(A) -> {cos, A}.
exp(A) -> {exp, A}.
conj(A) -> {conj, A}.
matrix(Rows) -> {matrix, Rows}.

%%% Normalisation

-spec canonicalise(sym_expr()) -> result(sym_expr()).
canonicalise(E) -> canonicalise(E, #{}).

-spec canonicalise(sym_expr(), opts()) -> result(sym_expr()).
canonicalise(E, Opts) ->
    %% rpoly: large plain products and powers stay in Rust until to_ast (K2)
    finish(run(E, Opts, fun(Ctx) -> output_deferred(pq_sym_norm:norm(E, Ctx#{rpoly => maps:get(poly_kernel, Ctx)}), Ctx) end)).

simplify(E) -> canonicalise(E, #{}).
simplify(E, Opts) -> canonicalise(E, Opts).

%% Eager, simultaneous substitution (SY 5.1); a replacement must fit the symbol's domain.
-spec subs(sym_expr(), #{atom() => sym_expr()}) -> result(sym_expr()).
subs(E, Map) -> subs(E, Map, #{}).

subs(E, Map, Opts) ->
    Domains = maps:get(domains, Opts, #{}),
    Bad = [{X, D} || X := R <- Map, D <- [maps:get(X, Domains, real)], not within(domain_of(R, Domains), D)],
    case lists:sort(Bad) of
        [] -> canonicalise(substitute(E, Map), Opts);
        [{X, D} | _] -> {error, {domain_violation, {X, D}}}
    end.

%% the raw replacement, without normalising (also used for {subs, E, Map} nodes)
substitute({sym, X} = S, Map) -> maps:get(X, Map, S);
substitute({subs, E, Inner}, Map) -> substitute(substitute(E, Inner), Map);
substitute({matrix, Rows}, Map) -> {matrix, [[substitute(E, Map) || E <- R] || R <- Rows]};
substitute({Op, L}, Map) when Op =:= add; Op =:= mul; Op =:= madd -> {Op, [substitute(X, Map) || X <- L]};
substitute({pow, B, E}, Map) -> {pow, substitute(B, Map), E};
substitute({Op, A}, Map) when Op =:= sin; Op =:= cos; Op =:= exp; Op =:= conj; Op =:= transpose;
                              Op =:= dagger -> {Op, substitute(A, Map)};
substitute({Op, A, B}, Map) when Op =:= smul; Op =:= matmul; Op =:= kron ->
    {Op, substitute(A, Map), substitute(B, Map)};
substitute(X, _Map) -> X.

%%% Differentiation (SY 5.3, FR-S05): canonical results; option order => N for higher derivatives

-spec differentiate(sym_expr(), atom()) -> result(sym_expr()).
differentiate(E, X) -> differentiate(E, X, #{}).

differentiate(E, X, Opts) when is_atom(X) ->
    N = maps:get(order, Opts, 1),
    run(E, Opts, fun(Ctx) ->
                         V = pq_sym_norm:norm(E, Ctx),
                         output(diff_value(V, X, N, Ctx), Ctx)
                 end).

diff_value({scalar, P}, X, N, Ctx) ->
    {scalar, lists:foldl(fun(_, Q) -> pq_sym_diff:diff(Q, X, Ctx) end, P, lists:seq(1, N))};
diff_value({identity, K}, _X, _N, _Ctx) -> {mat, K, K, lists:duplicate(K, lists:duplicate(K, #{}))};
diff_value({mat, R, C, Rows}, X, N, Ctx) ->
    {mat, R, C, [[element(2, diff_value({scalar, E}, X, N, Ctx)) || E <- Row] || Row <- Rows]}.

-spec gradient(sym_expr(), [atom()]) -> result([sym_expr()]).
gradient(E, Xs) -> gradient(E, Xs, #{}).

gradient(E, Xs, Opts) ->
    run(E, Opts, fun(Ctx) ->
                         {scalar, P} = scalar_value(pq_sym_norm:norm(E, Ctx)),
                         [output({scalar, pq_sym_diff:diff(P, X, Ctx)}, Ctx) || X <- Xs]
                 end).

%% the Jacobian of a scalar (one row) or of a column vector (one row per entry), one column per symbol
-spec jacobian(sym_expr(), [atom()]) -> result(sym_expr()).
jacobian(E, Xs) -> jacobian(E, Xs, #{}).

jacobian(E, Xs, Opts) ->
    run(E, Opts, fun(Ctx) ->
                         Entries = case pq_sym_norm:norm(E, Ctx) of
                                       {scalar, P} -> [P];
                                       {mat, _, 1, Rows} -> [P || [P] <- Rows];
                                       _ -> throw({bad_shape, {jacobian, E}})
                                   end,
                         Xs =/= [] orelse throw({bad_shape, {jacobian, no_symbols}}),
                         output({mat, length(Entries), length(Xs),
                                 [[pq_sym_diff:diff(P, X, Ctx) || X <- Xs] || P <- Entries]}, Ctx)
                 end).

scalar_value({scalar, _} = V) -> V;
scalar_value(_) -> throw({bad_shape, scalar_expected}).

%%% Evaluation

-spec evaluate(sym_expr(), binding()) -> result(term()).
evaluate(E, B) -> evaluate(E, B, #{}).

evaluate(E, B, Opts) ->
    run(E, Opts, fun(Ctx) -> pq_sym_eval:eval(pq_sym_norm:norm(E, Ctx), B, Opts) end).

%%% Matrices (expand mode, SY 4.6)

matmul(A, B) -> canonicalise({matmul, A, B}).
kron(A, B) -> canonicalise({kron, A, B}).
smul(S, M) -> canonicalise({smul, S, M}).
madd(Ms) -> canonicalise({madd, Ms}).
dagger(M) -> canonicalise({dagger, M}).
transpose(M) -> canonicalise({transpose, M}).

%%% Circuit-level (SY 6.1, FR-S06)

%% The symbolic single-photon unitary of a passive circuit: rows of the running matrix updated per instruction
%% (beamsplitter, phaseshift, unitary), each entry normalised at once; symbols take the circuit's domains.
-spec circuit_unitary(tuple()) -> result(sym_expr()).
circuit_unitary(C) -> circuit_unitary(C, #{}).

circuit_unitary(C, Opts) ->
    N = element(2, C),              % #circuit.n_modes
    Instrs = element(5, C),         % #circuit.instrs
    Domains = element(7, C),        % #circuit.symbols
    Max = pq_core:get(max_symbolic_modes, Opts),
    try
        N =< Max orelse throw({symbolic_limit, #{n_modes => N, max_symbolic_modes => Max}}),
        Ctx0 = context(Opts#{domains => maps:merge(Domains, maps:get(domains, Opts, #{}))}),
        %% the entries go straight to value_to_ast: large products may stay in Rust (K2)
        Ctx = Ctx0#{rpoly => maps:get(poly_kernel, Ctx0)},
        Id = list_to_tuple([list_to_tuple([pq_sym_poly:const(pq_sym_exact:c(case I =:= J of true -> 1; false -> 0 end))
                                           || J <- lists:seq(1, N)]) || I <- lists:seq(1, N)]),
        U = lists:foldl(fun(I, Acc) -> unitary_step(I, Acc, Ctx) end, Id, Instrs),
        {ok, pq_sym_norm:value_to_ast({mat, N, N, [tuple_to_list(R) || R <- tuple_to_list(U)]})}
    catch
        throw:{Cat, _} = R when is_atom(Cat) -> {error, R}
    end.

unitary_step({beamsplitter, A, B, T, P}, U, Ctx) ->
    E = fun(X) -> scalar_poly(X, Ctx) end,
    I = {complex, 0, 1},
    Cos = E({cos, T}),
    Gab = E({mul, [-1, {exp, {mul, [{mul, [-1, I]}, P]}}, {sin, T}]}),
    Gba = E({mul, [{exp, {mul, [I, P]}}, {sin, T}]}),
    rows(U, [A, B], [[Cos, Gab], [Gba, Cos]], Ctx);
unitary_step({phaseshift, M, P}, U, Ctx) ->
    rows(U, [M], [[scalar_poly({exp, {mul, [{complex, 0, 1}, P]}}, Ctx)]], Ctx);
unitary_step({unitary, Ms, {cmat, _, _, _} = Packed}, U, Ctx) ->
    {ok, {cmat, G}} = pq_core:unpack(Packed),
    rows(U, Ms, [[scalar_poly(Z, Ctx) || Z <- Row] || Row <- G], Ctx);
unitary_step({unitary, Ms, {matrix, G}}, U, Ctx) ->
    rows(U, Ms, [[scalar_poly(Z, Ctx) || Z <- Row] || Row <- G], Ctx);
unitary_step(I, _U, _Ctx) -> throw({unsupported_instruction, {element(1, I), symbolic_unitary}}).

scalar_poly(E, Ctx) ->
    case pq_sym_norm:norm(E, Ctx) of
        {scalar, P} -> P;
        _ -> throw({bad_shape, {scalar_expected, E}})
    end.

%% U <- G U on the rows Ms
rows(U, Ms, G, Ctx) ->
    Old = [element(M + 1, U) || M <- Ms],
    New = [list_to_tuple([pq_sym_poly:sum([pq_sym_poly:mul(Gij, element(Col, Row), Ctx) || {Gij, Row} <- lists:zip(GRow, Old)])
                          || Col <- lists:seq(1, tuple_size(element(1, U)))]) || GRow <- G],
    lists:foldl(fun({M, R}, Acc) -> setelement(M + 1, Acc, R) end, U, lists:zip(Ms, New)).

%% SY 6.2, FR-G04: the symbolic mean (a 2m x 1 matrix) and covariance (2m x 2m), xxpp, of a Gaussian circuit from
%% Mean0 and Cov0 (canonical, or numbers), with Hbar a number. The maps of A.2-A.6 are applied locally (rows of the
%% mean, rows and columns of the covariance on the instruction's modes): displacement, squeezing (cosh and sinh as
%% exponentials), phase shift, beamsplitter, passive unitary (the real symplectic [[Re U, -Im U], [Im U, Re U]],
%% Re and Im through conj, exact for real symbols) and loss (sqrt(eta) a power 1/2; declare eta nonnegative).
-spec gaussian_evolve(tuple(), map()) -> {ok, {sym_expr(), sym_expr()}} | {error, pq_core:reason()}.
gaussian_evolve(C, Opts) ->
    N = element(2, C),              % #circuit.n_modes
    Instrs = element(5, C),         % #circuit.instrs
    Domains = element(7, C),        % #circuit.symbols
    #{mean := Mean0, cov := Cov0, hbar := H} = Opts,
    Max = pq_core:get(max_symbolic_modes, Opts),
    try
        N =< Max orelse throw({symbolic_limit, #{n_modes => N, max_symbolic_modes => Max}}),
        Ctx = context(Opts#{domains => maps:merge(Domains, maps:get(domains, Opts, #{}))}),
        Mu = list_to_tuple([scalar_poly(X, Ctx) || X <- Mean0]),
        Sig = list_to_tuple([list_to_tuple([scalar_poly(X, Ctx) || X <- Row]) || Row <- Cov0]),
        {Mu2, Sig2} = lists:foldl(fun(I, Acc) -> gaussian_step(I, Acc, N, H, Ctx) end, {Mu, Sig}, Instrs),
        {ok, {pq_sym_norm:value_to_ast({mat, 2 * N, 1, [[X] || X <- tuple_to_list(Mu2)]}),
              pq_sym_norm:value_to_ast({mat, 2 * N, 2 * N, [tuple_to_list(R) || R <- tuple_to_list(Sig2)]})}}
    catch
        throw:{Cat, _} = R when is_atom(Cat) -> {error, R}
    end.

gaussian_step({beamsplitter, A, B, T, P}, St, N, _H, Ctx) ->
    I = {complex, 0, 1},
    U = [[{cos, T}, {mul, [-1, {exp, {mul, [{mul, [-1, I]}, P]}}, {sin, T}]}],
         [{mul, [{exp, {mul, [I, P]}}, {sin, T}]}, {cos, T}]],
    passive_step(U, [A, B], St, N, Ctx);
gaussian_step({phaseshift, M, P}, St, N, _H, Ctx) ->
    passive_step([[{exp, {mul, [{complex, 0, 1}, P]}}]], [M], St, N, Ctx);
gaussian_step({unitary, Ms, {cmat, _, _, _} = Packed}, St, N, _H, Ctx) ->
    {ok, {cmat, G}} = pq_core:unpack(Packed),
    passive_step(G, Ms, St, N, Ctx);
gaussian_step({unitary, Ms, {matrix, G}}, St, N, _H, Ctx) ->
    passive_step(G, Ms, St, N, Ctx);
gaussian_step({displacement, M, R, P}, {Mu, Sig}, N, H, Ctx) ->
    K = {mul, [{pow, 2 * H, {rat, 1, 2}}, R]},
    X = M + 1, Pp = N + M + 1,
    Mu2 = setelement(Pp, setelement(X, Mu, pq_sym_poly:sum([element(X, Mu), scalar_poly({mul, [K, {cos, P}]}, Ctx)])),
                     pq_sym_poly:sum([element(Pp, Mu), scalar_poly({mul, [K, {sin, P}]}, Ctx)])),
    {Mu2, Sig};
gaussian_step({squeezing, M, R, P}, St, N, _H, Ctx) ->
    Ch = {mul, [{rat, 1, 2}, {add, [{exp, R}, {exp, {mul, [-1, R]}}]}]},
    Sh = {mul, [{rat, 1, 2}, {add, [{exp, R}, {mul, [-1, {exp, {mul, [-1, R]}}]}]}]},
    E = fun(X) -> scalar_poly(X, Ctx) end,
    S = [[E({add, [Ch, {mul, [-1, {cos, P}, Sh]}]}), E({mul, [-1, {sin, P}, Sh]})],
         [E({mul, [-1, {sin, P}, Sh]}), E({add, [Ch, {mul, [{cos, P}, Sh]}]})]],
    symplectic_step(S, [M], St, N, Ctx);
gaussian_step({loss, M, Eta}, {Mu, Sig}, N, H, Ctx) ->
    R = scalar_poly({pow, Eta, {rat, 1, 2}}, Ctx),
    Add = scalar_poly({mul, [{add, [1, {mul, [-1, Eta]}]}, H, {rat, 1, 2}]}, Ctx),
    Idx = [M + 1, N + M + 1],
    Mu2 = lists:foldl(fun(K, A) -> setelement(K, A, pq_sym_poly:mul(R, element(K, A), Ctx)) end, Mu, Idx),
    %% rows and columns of the mode times sqrt(eta) (its own block by eta), then (1 - eta) hbar/2 on its diagonal
    Scaled = lists:foldl(fun(K, A) -> scale_row_col(A, K, R, Ctx) end, Sig, Idx),
    Sig2 = lists:foldl(fun(K, A) -> set2(A, K, K, pq_sym_poly:sum([get2(A, K, K), Add])) end, Scaled, Idx),
    {Mu2, Sig2};
gaussian_step(I, _St, _N, _H, _Ctx) -> throw({unsupported_instruction, {element(1, I), symbolic_gaussian}}).

passive_step(U0, Ms, St, N, Ctx) ->
    U = [[scalar_poly(Z, Ctx) || Z <- Row] || Row <- U0],
    Half = pq_sym_exact:c(pq_sym_exact:rat(1, 2)),
    MinusHalfI = pq_sym_exact:cq(pq_sym_exact:rat(0, 1), pq_sym_exact:rat(-1, 2)),
    Re = fun(Z) -> pq_sym_poly:scale(pq_sym_poly:sum([Z, pq_sym_poly:conj(Z, Ctx)]), Half) end,
    %% Im z = (z - conj z) / (2i) = -(i/2)(z - conj z)
    Im = fun(Z) -> pq_sym_poly:scale(pq_sym_poly:sum([Z, pq_sym_poly:scale(pq_sym_poly:conj(Z, Ctx), pq_sym_exact:c(-1))]), MinusHalfI) end,
    Neg = fun(Z) -> pq_sym_poly:scale(Z, pq_sym_exact:c(-1)) end,
    S = [[Re(Z) || Z <- Row] ++ [Neg(Im(Z)) || Z <- Row] || Row <- U] ++
        [[Im(Z) || Z <- Row] ++ [Re(Z) || Z <- Row] || Row <- U],
    symplectic_step(S, Ms, St, N, Ctx).

%% mu <- S mu on the indices of Ms (x then p), Sigma <- S Sigma S^T on those rows and columns
symplectic_step(S, Ms, {Mu, Sig}, N, Ctx) ->
    Idx = [M + 1 || M <- Ms] ++ [N + M + 1 || M <- Ms],
    Dot = fun(Row, Vals) -> pq_sym_poly:sum([pq_sym_poly:mul(A, B, Ctx) || {A, B} <- lists:zip(Row, Vals)]) end,
    OldMu = [element(K, Mu) || K <- Idx],
    Mu2 = lists:foldl(fun({K, Row}, A) -> setelement(K, A, Dot(Row, OldMu)) end, Mu, lists:zip(Idx, S)),
    Dim = tuple_size(Sig),
    %% rows: Sigma[Idx, :] <- S Sigma[Idx, :]
    OldRows = [element(K, Sig) || K <- Idx],
    NewRows = [list_to_tuple([Dot(Row, [element(J, R) || R <- OldRows]) || J <- lists:seq(1, Dim)]) || Row <- S],
    Sig1 = lists:foldl(fun({K, R}, A) -> setelement(K, A, R) end, Sig, lists:zip(Idx, NewRows)),
    %% columns: Sigma[:, Idx] <- Sigma[:, Idx] S^T
    Sig2 = list_to_tuple([begin
                              Old = [element(K, R) || K <- Idx],
                              lists:foldl(fun({K, SRow}, A) -> setelement(K, A, Dot(SRow, Old)) end, R, lists:zip(Idx, S))
                          end || R <- tuple_to_list(Sig1)]),
    {Mu2, Sig2}.

scale_row_col(Sig, K, F, Ctx) ->
    Row = list_to_tuple([pq_sym_poly:mul(F, X, Ctx) || X <- tuple_to_list(element(K, Sig))]),
    Sig1 = setelement(K, Sig, Row),
    list_to_tuple([setelement(K, R, pq_sym_poly:mul(F, element(K, R), Ctx)) || R <- tuple_to_list(Sig1)]).

get2(T, I, J) -> element(J, element(I, T)).
set2(T, I, J, V) -> setelement(I, T, setelement(J, element(I, T), V)).

%% SY 6.3: the permanent of a symbolic k x k matrix by Glynn's formula in canonical arithmetic,
%% 2^(k-1) terms, k at most max_symbolic_photons
-spec permanent(sym_expr()) -> result(sym_expr()).
permanent(M) ->
    finish(run(M, #{}, fun(Ctx) ->
                        {mat, K, K2, Rows} = case pq_sym_norm:norm(M, Ctx) of
                                                 {mat, _, _, _} = X -> X;
                                                 _ -> throw({bad_shape, permanent})
                                             end,
                        K =:= K2 orelse throw({bad_shape, {permanent, K, K2}}),
                        %% in Rust when the entries allow (pq_sym_poly:kernel_permanent/3), made in the caller
                        check_photons(K),
                        P = case K > 0 andalso pq_sym_poly:kernel_permanent(Rows, K, Ctx) of
                                {ok, H} -> H;
                                _ -> glynn(Rows, K, Ctx)
                            end,
                        output_deferred({scalar, P}, Ctx)
                end)).

%% The permanent row by row over the subsets of columns already used (the expansion FORM does with its pruning,
%% comparison of 2 Nov 2026): for each set S of i columns, the sum over every assignment of rows 1..i to S; a term
%% that would use a column twice is never made, so the intermediates stay the size of partial results -- where
%% Glynn's formula expanded 2^(k-1) products of k linear forms and relied on their cancelling (7 x 7 symbolic:
%% 429 s with Glynn). Each entry is used as it is, repeated rows or columns included (bs_probability).
check_photons(K) ->
    Max = pq_core:get(max_symbolic_photons),
    K =< Max orelse throw({symbolic_limit, #{photons => K, max_symbolic_photons => Max}}).

glynn(Rows, K, Ctx) ->
    check_photons(K),
    case K of
        0 -> pq_sym_poly:one();
        _ ->
            Entries = list_to_tuple([list_to_tuple(R) || R <- Rows]),
            Final = lists:foldl(
                      fun(I, Layer) ->
                              Row = element(I, Entries),
                              maps:fold(
                                fun(Used, P, Acc) ->
                                        lists:foldl(
                                          fun(J, A2) ->
                                                  Bit = 1 bsl (J - 1),
                                                  E = element(J, Row),
                                                  case Used band Bit =:= 0 andalso map_size(E) > 0 of
                                                      true ->
                                                          T = pq_sym_poly:mul(P, E, Ctx),
                                                          maps:update_with(Used bor Bit, fun(Q) -> pq_sym_poly:add(Q, T) end, T, A2);
                                                      false -> A2
                                                  end
                                          end, Acc, lists:seq(1, K))
                                end, #{}, Layer)
                      end, #{0 => pq_sym_poly:one()}, lists:seq(1, K)),
            maps:get((1 bsl K) - 1, Final, #{})
    end.

%% Glynn's formula, kept as the test oracle of the subset expansion (pq_symbolic_tests)
glynn_reference(Rows, K) ->
    run({add, [0]}, #{}, fun(Ctx) -> output({scalar, glynn_formula([[pq_sym_norm_entry(E, Ctx) || E <- R] || R <- Rows], K, Ctx)}, Ctx) end).

pq_sym_norm_entry(E, Ctx) -> scalar_poly(E, Ctx).

glynn_formula(Rows, K, Ctx) ->
    case K of
        0 -> pq_sym_poly:one();
        _ ->
            Cols = transpose_lists(Rows),
            Terms = [begin
                         Deltas = [1 | [case (Bits bsr (I - 1)) band 1 of 1 -> -1; 0 -> 1 end || I <- lists:seq(1, K - 1)]],
                         Sign = lists:foldl(fun erlang:'*'/2, 1, Deltas),
                         Prod = lists:foldl(fun(Col, Acc) ->
                                                    S = pq_sym_poly:sum([pq_sym_poly:scale(E, pq_sym_exact:c(D)) || {D, E} <- lists:zip(Deltas, Col)]),
                                                    pq_sym_poly:mul(Acc, S, Ctx)
                                            end, pq_sym_poly:one(), Cols),
                         pq_sym_poly:scale(Prod, pq_sym_exact:c(Sign))
                     end || Bits <- lists:seq(0, (1 bsl (K - 1)) - 1)],
            pq_sym_poly:scale(pq_sym_poly:sum(Terms), pq_sym_exact:c(pq_sym_exact:rat(1, 1 bsl (K - 1))))
    end.

transpose_lists([[] | _]) -> [];
transpose_lists(Rows) -> [[hd(R) || R <- Rows] | transpose_lists([tl(R) || R <- Rows])].

%% SY 6.3: |perm(U_(T,S))|^2 / (prod S! prod T!) for a symbolic unitary (rows repeated by T, columns by S)
-spec bs_probability(sym_expr(), [non_neg_integer()], [non_neg_integer()]) -> result(sym_expr()).
bs_probability(U, S, T) ->
    run(U, #{}, fun(Ctx) ->
                        {mat, M, M, Rows} = case pq_sym_norm:norm(U, Ctx) of
                                                {mat, R, R, _} = X when length(S) =:= R, length(T) =:= R -> X;
                                                _ -> throw({bad_shape, bs_probability})
                                            end,
                        _ = M,
                        lists:sum(S) =:= lists:sum(T) orelse throw({bad_shape, {photon_numbers, S, T}}),
                        RowIdx = lists:append([lists:duplicate(Tj, J) || {J, Tj} <- lists:enumerate(1, T)]),
                        ColIdx = lists:append([lists:duplicate(Si, I) || {I, Si} <- lists:enumerate(1, S)]),
                        Sub = [[lists:nth(C, lists:nth(R, Rows)) || C <- ColIdx] || R <- RowIdx],
                        P = glynn(Sub, length(RowIdx), Ctx),
                        Fact = lists:foldl(fun(X, A) -> A * fact(X) end, 1, S ++ T),
                        Prob = pq_sym_poly:mul(P, pq_sym_poly:conj(P, Ctx), Ctx),
                        output({scalar, pq_sym_poly:scale(Prob, pq_sym_exact:c(pq_sym_exact:rat(1, Fact)))}, Ctx)
                end).

fact(0) -> 1;
fact(N) -> N * fact(N - 1).

%% SY 6.1: U U^dagger against the identity by canonical form: true or false for entries in E (SY 3.5; circuits of
%% beamsplitters and phase shifts), unknown when an entry holds an opaque base or a float
-spec is_unitary(sym_expr()) -> {ok, true | false | unknown} | {error, pq_core:reason()}.
is_unitary(M) -> is_unitary(M, #{}).

is_unitary(M, Opts) ->
    case canonicalise({matmul, M, {dagger, M}}, Opts) of
        {ok, {matrix, Rows}} ->
            N = length(Rows),
            Id = [[case I =:= J of true -> 1; false -> 0 end || J <- lists:seq(1, N)] || I <- lists:seq(1, N)],
            case Rows =:= Id of
                true -> {ok, true};
                false ->
                    case decisive(Rows) of
                        true -> {ok, false};
                        false -> {ok, unknown}
                    end
            end;
        {ok, {identity, _}} -> {ok, true};
        {ok, _} -> {error, {bad_shape, is_unitary}};
        {error, _} = E -> E
    end.

%% no opaque base and no float anywhere: the canonical form is faithful (SY 3.5)
decisive(T) when is_float(T) -> false;
decisive({pow, B, E}) -> not is_tuple(B) orelse element(1, B) =/= add andalso decisive(B) andalso decisive(E);
decisive(T) when is_tuple(T) -> lists:all(fun decisive/1, tuple_to_list(T));
decisive(L) when is_list(L) -> lists:all(fun decisive/1, L);
decisive(_) -> true.

%%% Straight-line programs (SY 5.5, FR-S10) and equivalence (SY 5.4, FR-S11)

%% Exprs (scalars, or one matrix) compiled for repeated evaluation; Symbols are the inputs
-spec compile([sym_expr()], [atom()]) -> result(tuple()).
compile(Exprs, Symbols) -> compile(Exprs, Symbols, #{}).

%% Opts as for evaluate/3 (domains above all): evaluate(E, B, Opts) is eval_slp(compile([E], free_symbols(E), Opts), B)
-spec compile([sym_expr()], [atom()], opts()) -> result(tuple()).
compile(Exprs, Symbols, Opts) ->
    run({add, [0 | [E || E <- Exprs, sh_ok(E)]]}, Opts, fun(Ctx) ->
                                                       Vs = [pq_sym_norm:norm(E, Ctx) || E <- Exprs],
                                                       pq_sym_slp:compile(Vs, Symbols, maps:with([optimize, order_search, exact_constants], Opts))
                                               end).

sh_ok(E) -> shape(E) =:= {ok, scalar}.

-spec eval_slp(tuple(), binding()) -> result(term()).
eval_slp(P, B) -> pq_sym_slp:eval(P, B).

-spec eval_batch(tuple(), [binding()]) -> result(list()).
eval_batch(P, Bs) -> pq_sym_slp:eval_batch(P, Bs).

%% The greatest common divisor of two polynomials with rational coefficients (decision record SY2): primitive, with
%% integer coefficients, the leading term positive (terms ordered lexicographically by the sorted symbols); the gcd
%% of 0 and 0 is 0. By the modular method (pq_sym_gcd). Complex coefficients, pi, exponentials and sums in
%% denominators are refused: {error, {symbolic_unsupported, gcd}}.
-spec gcd(sym_expr(), sym_expr()) -> result(sym_expr()).
gcd(A, B) -> gcd(A, B, #{}).

-spec gcd(sym_expr(), sym_expr(), opts()) -> result(sym_expr()).
gcd(A, B, Opts) ->
    %% a sum of monomials with exact coefficients is read as it is; anything else is canonicalised first
    Vars = lists:usort(raw_symbols(A) ++ raw_symbols(B)),
    Direct = run({add, [A, B]}, Opts, fun(Ctx) ->
                                              Pa = pq_sym_gcd:from_ast(A, Vars),
                                              Pb = pq_sym_gcd:from_ast(B, Vars),
                                              gcd_out(pq_sym_gcd:gcd(Pa, Pb), Vars, Ctx)
                                      end),
    case Direct of
        {error, {symbolic_unsupported, gcd}} ->
            case {canonicalise(A, Opts), canonicalise(B, Opts)} of
                {{ok, Ca}, {ok, Cb}} ->
                    finish(run({add, [Ca, Cb]}, Opts, fun(Ctx) ->
                                                              G = pq_sym_gcd:gcd(pq_sym_gcd:from_ast(Ca, Vars), pq_sym_gcd:from_ast(Cb, Vars)),
                                                              gcd_out(G, Vars, Ctx)
                                                      end));
                {{error, _} = E, _} -> E;
                {_, E} -> E
            end;
        _ -> finish(Direct)
    end.

%% The factorisation of a polynomial with rational coefficients into irreducibles over the rationals (decision record
%% SY2): {ok, {Constant, [{Factor, Multiplicity}]}}, each factor canonical, primitive with integer coefficients and
%% its leading term positive, sorted by multiplicity, then by factor; the product is the input. By Zassenhaus's and
%% Wang's methods (pq_sym_factor). Refused as gcd/2,3 refuses: {error, {symbolic_unsupported, factor}}.
-spec factor(sym_expr()) -> result({sym_expr(), [{sym_expr(), pos_integer()}]}).
factor(E) -> factor(E, #{}).

-spec factor(sym_expr(), opts()) -> result({sym_expr(), [{sym_expr(), pos_integer()}]}).
factor(E, Opts) ->
    Vars = lists:usort(raw_symbols(E)),
    F = fun(X) -> run(X, Opts#{bounded => maps:get(bounded, Opts, true)},
                      fun(_Ctx) ->
                              P = try pq_sym_gcd:from_ast(X, Vars)
                                  catch throw:{symbolic_unsupported, gcd} -> throw({symbolic_unsupported, factor})
                                  end,
                              case map_size(P) of
                                  0 -> {0, []};
                                  _ -> {C, Fs} = pq_sym_factor:factor(P), pq_sym_factor:result(C, Fs, Vars)
                              end
                      end) end,
    case F(E) of
        {error, {symbolic_unsupported, factor}} ->
            case canonicalise(E, Opts) of
                {ok, Ce} -> F(Ce);
                Err -> Err
            end;
        R -> R
    end.

%% A rational function of symbols with rational coefficients in normal form (decision record SY2): N * D^-1 with N
%% and D coprime polynomials, D primitive with integer coefficients and a positive leading coefficient, so equal
%% rational functions give equal results (the canonical form alone keeps a sum in a denominator as an opaque atom).
%% Anything else: {error, {symbolic_unsupported, rational}}.
-spec cancel(sym_expr()) -> result(sym_expr()).
cancel(E) -> cancel(E, #{}).

-spec cancel(sym_expr(), opts()) -> result(sym_expr()).
cancel(E, Opts) ->
    with_ratfun(E, Opts, fun(R, Vars, Ctx) -> output_deferred(pq_sym_norm:norm(pq_sym_ratfun:to_ast(R, Vars), Ctx), Ctx) end).

%% {N, D} of cancel/1,2's normal form, each canonical
-spec numer_denom(sym_expr()) -> result({sym_expr(), sym_expr()}).
numer_denom(E) -> numer_denom(E, #{}).

-spec numer_denom(sym_expr(), opts()) -> result({sym_expr(), sym_expr()}).
numer_denom(E, Opts) ->
    case with_ratfun(E, Opts, fun(R, Vars, _) -> {R, Vars} end) of
        {ok, {{N, D}, Vars}} ->
            {ok, Na} = canonicalise(pq_sym_gcd:to_ast(N, Vars), Opts),
            {ok, Da} = canonicalise(pq_sym_gcd:to_ast(D, Vars), Opts),
            {ok, {Na, Da}};
        Err -> Err
    end.

with_ratfun(E, Opts, Out) ->
    case canonicalise(E, Opts) of
        {ok, C} ->
            Vars = free_symbols(C),
            finish(run(C, Opts, fun(Ctx) -> Out(pq_sym_ratfun:from_ast(C, Vars), Vars, Ctx) end));
        Err -> Err
    end.

%% equal as rational functions (true | false), or none when either is not a rational function of symbols
ratfun_equiv(X, Y) ->
    Vars = lists:usort(free_symbols(X) ++ free_symbols(Y)),
    try pq_sym_ratfun:equal(pq_sym_ratfun:from_ast(X, Vars), pq_sym_ratfun:from_ast(Y, Vars))
    catch throw:{symbolic_unsupported, rational} -> none
    end.

%% The series of E in the symbol X up to (not including) X^N (decision record SY2, phase 2): {ok, Sum}, Sum a sum of
%% terms C X^e log(X)^m, e rational, sorted by e then m, each C exact (rational functions of the other symbols, with
%% opaque atoms such as {log, 2}, {const, pi}, {const, euler_gamma}, {zeta, 3}, {gamma, {rat, 1, 3}}). Opts: at => 0 |
%% A | infinity, and the expression limits. The functions: those of pq_series. Refused:
%% {error, {symbolic_unsupported, #{series => Why}}}.
-spec series(sym_expr(), atom(), integer()) -> result(sym_expr()).
series(E, X, N) -> series(E, X, N, #{}).

-spec series(sym_expr(), atom(), integer(), opts()) -> result(sym_expr()).
series(E, X, N, Opts) when is_atom(X), is_integer(N) ->
    case run(E, Opts, fun(_) -> pq_series:series(E, X, N, Opts) end) of
        {ok, {ok, S}} -> {ok, S};
        {ok, {error, _} = Err} -> Err;
        Err -> Err
    end.

%% The limit of E as X tends to A (0, an expression free of X, or infinity), from its series' leading term: {ok, L}
%% for a finite limit, {ok, {diverges, Term}} (Term the leading term C X^e log(X)^m) otherwise.
-spec limit(sym_expr(), atom(), term()) -> result(term()).
limit(E, X, A) -> limit(E, X, A, #{}).

-spec limit(sym_expr(), atom(), term(), opts()) -> result(term()).
limit(E, X, A, Opts) ->
    limit_order(E, X, A, Opts, 1).

limit_order(_, _, _, _, N) when N > 64 -> {error, {symbolic_unsupported, #{limit => no_leading_term}}};
limit_order(E, X, A, Opts, N) ->
    case run(E, Opts, fun(_) -> pq_series:terms(E, X, N, Opts#{at => A}) end) of
        {ok, {ok, {[], _}}} -> limit_order(E, X, A, Opts, 2 * N);
        {ok, {ok, {Ts, _}}} ->
            Key = fun({Ex, M, _}) -> {-ex_value(Ex, A), M} end,
            % the dominant term: at infinity the highest power, at a point the lowest; then the highest log power
            {Ex, M, C} = hd(lists:sort(fun(T1, T2) -> Key(T1) >= Key(T2) end, Ts)),
            case {ex_sign(Ex, A), M} of
                {pos, _} -> {ok, 0};
                {zero, 0} -> {ok, C};
                _ -> {ok, {diverges, {mul, [C] ++ [{pow, {sym, X}, Ex} || Ex =/= 0] ++ [{pow, {log, {sym, X}}, M} || M > 0]}}}
            end;
        {ok, {error, _} = Err} -> Err;
        Err -> Err
    end.

ex_value(K, A) when is_integer(K) -> ex_value({rat, K, 1}, A);
ex_value({rat, P, Q}, infinity) -> -P / Q;
ex_value({rat, P, Q}, _) -> P / Q.

ex_sign(Ex, A) ->
    V = ex_value(Ex, A),
    if V > 0 -> pos; V == 0 -> zero; true -> neg end.

gcd_out(G, Vars, Ctx) -> output_deferred(pq_sym_norm:norm(pq_sym_gcd:to_ast(G, Vars), Ctx), Ctx).

%% true: identical canonical forms (sound for all inputs); false: different forms of two expressions in E (SY 3.5),
%% where the form is faithful -- no opaque base, no float, no exponential of pi (finding: exp(i pi/5) equals a
%% radical expression that E also contains), and not a constant-only difference with e or pi (transcendence);
%% unknown otherwise, then with numeric => N, N random real bindings may show a difference: {false, Binding}.
-spec equiv(sym_expr(), sym_expr(), opts()) -> {ok, true | false | unknown | {false, binding()}} | {error, pq_core:reason()}.
equiv(A, B, Opts) ->
    case {canonicalise(A, Opts), canonicalise(B, Opts)} of
        {{ok, X}, {ok, X}} -> {ok, true};
        {{ok, X}, {ok, Y}} ->
            {ok, D} = canonicalise({add, [X, {mul, [-1, Y]}]}, Opts),
            Faithful = faithful(X) andalso faithful(Y) andalso not (free_symbols(D) =:= [] andalso has_transcendental(D)),
            case Faithful of
                true -> {ok, false};
                false ->
                    %% rational functions of symbols: decided by their normal forms (SY2)
                    case ratfun_equiv(X, Y) of
                        none -> numeric_equiv(D, maps:get(numeric, Opts, 0), Opts);
                        Eq -> {ok, Eq}
                    end
            end;
        {{error, _} = E, _} -> E;
        {_, E} -> E
    end.

faithful(T) when is_float(T) -> false;
faithful({pow, B, _}) when is_tuple(B), element(1, B) =:= add -> false;
faithful({pow, {mul, _} , _}) -> false;
faithful({conj, {pow, _, _}}) -> false;
faithful({exp, Arg}) -> not lists:member({const, pi}, flatten(Arg)) andalso faithful(Arg);
faithful(T) when is_tuple(T) -> lists:all(fun faithful/1, tuple_to_list(T));
faithful(L) when is_list(L) -> lists:all(fun faithful/1, L);
faithful(_) -> true.

has_transcendental(D) -> lists:any(fun(X) -> X =:= {const, pi} orelse (is_tuple(X) andalso element(1, X) =:= exp) end, flatten(D)).

flatten(T) when is_tuple(T) -> [T | lists:append([flatten(X) || X <- tuple_to_list(T)])];
flatten(L) when is_list(L) -> lists:append([flatten(X) || X <- L]);
flatten(X) -> [X].

numeric_equiv(_D, 0, _Opts) -> {ok, unknown};
numeric_equiv(D, N, Opts) ->
    Syms = free_symbols(D),
    Domains = maps:get(domains, Opts, #{}),
    Draw = fun(K) ->
                   maps:from_list([{S, case maps:get(S, Domains, real) of
                                           real -> 2 * frac(K, I) - 1;
                                           complex -> {complex, 2 * frac(K, I) - 1, 2 * frac(K + 7, I) - 1};
                                           _ -> 0.1 + frac(K, I)
                                       end} || {I, S} <- lists:enumerate(Syms)])
           end,
    Diff = [B || K <- lists:seq(1, N), B <- [Draw(K)],
                 case evaluate(D, B, Opts) of
                     {ok, V} -> magnitude(V) > 1.0e-9;
                     _ -> false
                 end],
    case Diff of
        [B | _] -> {ok, {false, B}};
        [] -> {ok, unknown}
    end.

%% deterministic points in [0, 1) (a fixed low-discrepancy sequence: no rand)
frac(K, I) -> X = K * 0.6180339887498949 + I * 0.4142135623730951, X - math:floor(X).

magnitude({complex, R, I}) -> math:sqrt(R * R + I * I);
magnitude(X) -> abs(X).

%%% Inspection

%% FR-S08: scalar | {Rows, Cols}, without normalising
-spec shape(sym_expr()) -> {ok, scalar | {pos_integer(), pos_integer()}} | {error, {bad_shape, term()}}.
shape(E) ->
    try {ok, sh(E)}
    catch throw:{bad_shape, _} = R -> {error, R}
    end.

sh(N) when is_number(N) -> scalar;
sh({rat, _, _}) -> scalar;
sh({complex, _, _}) -> scalar;
sh({const, _}) -> scalar;
sh({sym, _}) -> scalar;
sh({Op, L} = N) when Op =:= add; Op =:= mul -> [scalar_or(sh(X), N) || X <- L], scalar;
sh({pow, B, _} = N) -> scalar_or(sh(B), N);
sh({Op, A} = N) when Op =:= sin; Op =:= cos; Op =:= exp; Op =:= conj -> scalar_or(sh(A), N);
sh({matrix, [R | _] = Rows} = N) when is_list(R), R =/= [] ->
    lists:all(fun(X) -> is_list(X) andalso length(X) =:= length(R) end, Rows) orelse throw({bad_shape, N}),
    [scalar_or(sh(X), N) || Row <- Rows, X <- Row],
    {length(Rows), length(R)};
sh({identity, K}) when is_integer(K), K > 0 -> {K, K};
sh({madd, [M | Ms]} = N) ->
    S = matrix_or(sh(M), N),
    [S =:= matrix_or(sh(X), N) orelse throw({bad_shape, N}) || X <- Ms],
    S;
sh({smul, S, M} = N) -> scalar_or(sh(S), N), matrix_or(sh(M), N);
sh({matmul, A, B} = N) ->
    {R, K} = matrix_or(sh(A), N),
    case matrix_or(sh(B), N) of
        {K, C} -> {R, C};
        _ -> throw({bad_shape, N})
    end;
sh({kron, A, B} = N) -> {R1, C1} = matrix_or(sh(A), N), {R2, C2} = matrix_or(sh(B), N), {R1 * R2, C1 * C2};
sh({Op, M} = N) when Op =:= transpose; Op =:= dagger -> {R, C} = matrix_or(sh(M), N), {C, R};
sh({subs, E, _}) -> sh(E);
sh(N) -> throw({bad_shape, N}).

scalar_or(scalar, _) -> scalar;
scalar_or(_, N) -> throw({bad_shape, N}).

matrix_or({_, _} = S, _) -> S;
matrix_or(_, N) -> throw({bad_shape, N}).

%% symbols of the canonical form (so x - x has none)
-spec free_symbols(sym_expr()) -> [atom()].
free_symbols(E) ->
    case run(E, #{}, fun(Ctx) -> pq_sym_eval:free_symbols(pq_sym_norm:norm(E, Ctx)) end) of
        {ok, L} -> L;
        {error, _} -> lists:usort(raw_symbols(E))
    end.

raw_symbols({sym, X}) -> [X];
raw_symbols(T) when is_tuple(T) -> raw_symbols(tuple_to_list(T));
raw_symbols(L) when is_list(L) -> lists:append([raw_symbols(X) || X <- L]);
raw_symbols(M) when is_map(M) -> raw_symbols(maps:values(M));
raw_symbols(_) -> [].

%% node count of the unfolded tree (SY 7.3: sharing is not counted)
-spec size(sym_expr()) -> non_neg_integer().
size(E) -> count(E, infinity).

count(E, Cap) ->
    try nodes(E, 0, Cap)
    catch throw:{over, N} -> N
    end.

nodes(_, N, Cap) when is_integer(Cap), N > Cap -> throw({over, N});
nodes(T, N, Cap) when is_tuple(T) -> nodes_list(tl(tuple_to_list(T)), N + 1, Cap);
nodes(L, N, Cap) when is_list(L) -> nodes_list(L, N, Cap);
nodes(M, N, Cap) when is_map(M) -> nodes_list(maps:values(M), N, Cap);
nodes(_, N, _Cap) -> N + 1.

nodes_list([], N, _Cap) -> N;
nodes_list([H | T], N, Cap) -> nodes_list(T, nodes(H, N, Cap), Cap).

%% SHA-256 of the canonical form (SY 7.4)
-spec hash(sym_expr()) -> binary().
hash(E) ->
    {ok, V} = run(E, #{}, fun(Ctx) -> pq_sym_norm:norm(E, Ctx) end),
    crypto:hash(sha256, term_to_binary(V, [deterministic])).

%% Conservative domain inference (SY 5.1), positive < nonnegative < real < complex.
-spec domain_of(sym_expr(), #{atom() => pq_core:domain()}) -> pq_core:domain().
domain_of(N, _) when is_number(N) -> sign_domain(N);
domain_of({rat, N, _}, _) -> sign_domain(N);
domain_of({complex, R, I}, D) ->
    case I == 0 orelse I =:= {rat, 0, 1} of
        true -> domain_of(R, D);
        false -> complex
    end;
domain_of({const, i}, _) -> complex;
domain_of({const, _}, _) -> positive;
domain_of({sym, X}, D) -> maps:get(X, D, real);
domain_of({add, L}, D) -> join([domain_of(X, D) || X <- L]);
domain_of({mul, L}, D) ->
    case join([domain_of(X, D) || X <- L]) of
        nonnegative -> nonnegative;
        Other -> Other
    end;
domain_of({pow, B, E}, D) ->
    case {domain_of(B, D), E} of
        {positive, _} -> positive;
        {nonnegative, _} -> nonnegative;
        {real, N} when is_integer(N), N rem 2 =:= 0 -> nonnegative;
        {real, N} when is_integer(N) -> real;
        _ -> complex
    end;
domain_of({exp, A}, D) -> case within(domain_of(A, D), real) of true -> positive; false -> complex end;
domain_of({F, A}, D) when F =:= sin; F =:= cos -> case within(domain_of(A, D), real) of true -> real; false -> complex end;
domain_of({conj, A}, D) -> domain_of(A, D);
domain_of({subs, E, M}, D) -> domain_of(substitute(E, M), D);
domain_of(_, _) -> complex.

sign_domain(N) when N > 0 -> positive;
sign_domain(N) when N == 0 -> nonnegative;
sign_domain(_) -> real.

rank(positive) -> 0;
rank(nonnegative) -> 1;
rank(real) -> 2;
rank(complex) -> 3.

join([]) -> positive;
join(Ds) -> lists:foldl(fun(A, B) -> case rank(A) >= rank(B) of true -> A; false -> B end end, positive, Ds).

within(A, B) -> rank(A) =< rank(B).

%%% Printing and parsing (SY 5.6)

-spec pretty(sym_expr()) -> binary().
pretty(E) -> pretty(E, #{}).

pretty(E, Opts) ->
    case canonicalise(E, Opts) of
        {ok, C} -> pq_sym_print:ascii(C);
        {error, _} -> pq_sym_print:ascii(E)
    end.

-spec to_latex(sym_expr()) -> binary().
to_latex(E) ->
    case canonicalise(E) of
        {ok, C} -> pq_sym_print:latex(C);
        {error, _} -> pq_sym_print:latex(E)
    end.

-spec parse(binary()) -> result(sym_expr()).
parse(Bin) -> parse(Bin, #{}).
parse(Bin, Opts) -> pq_sym_parse:parse(Bin, Opts).

%%% Parameter resolution (SY 6.5, FR-S15): the binding overrides the defaults; the remaining defaults are resolved
%%% in dependency order (exactly, until the final conversion); a cycle is param_cycle, a symbol left without a
%%% value unbound_symbols. Returns every value as a float or {complex, Re, Im}.

-spec resolve_params(#{atom() => sym_expr()}, binding()) -> result(#{atom() => float() | {complex, float(), float()}}).
resolve_params(Defaults, Binding) -> resolve_params(Defaults, Binding, #{}).

resolve_params(Defaults, Binding, Opts) ->
    try
        maps:foreach(fun(_, V) -> number_value(V) orelse throw({bad_expr, {binding, V}}) end, Binding),
        Pending = maps:without(maps:keys(Binding), Defaults),
        Deps = maps:map(fun(_, E) -> [S || S <- free_symbols(E), maps:is_key(S, Pending)] end, Pending),
        Order = topo(Deps),
        Env = lists:foldl(
                fun(P, Env0) ->
                        {ok, V} = ok_or_throw(canonicalise(substitute(maps:get(P, Pending), Env0), Opts)),
                        Env0#{P => V}
                end, Binding, Order),
        Unbound = lists:usort(lists:append([free_symbols(V) || V <- maps:values(Env)])),
        Unbound =:= [] orelse throw({unbound_symbols, Unbound}),
        {ok, maps:map(fun(_, V) -> {ok, X} = ok_or_throw(evaluate(V, #{}, Opts)), X end, Env)}
    catch
        throw:{_, _} = R -> {error, R}
    end.

number_value(N) when is_number(N) -> true;
number_value({rat, N, D}) -> is_integer(N) andalso is_integer(D) andalso D > 0;
number_value({complex, _, _}) -> true;
number_value(_) -> false.

ok_or_throw({ok, _} = R) -> R;
ok_or_throw({error, R}) -> throw(R).

%% Kahn's algorithm, smallest name first, so the order is deterministic
topo(Deps) -> topo(Deps, []).

topo(Deps, Acc) when map_size(Deps) =:= 0 -> lists:reverse(Acc);
topo(Deps, Acc) ->
    case lists:sort([P || P := [] <- Deps]) of
        [] -> throw({param_cycle, lists:sort(maps:keys(Deps))});
        [P | _] ->
            Rest = maps:map(fun(_, Ds) -> lists:delete(P, Ds) end, maps:remove(P, Deps)),
            topo(Rest, [P | Acc])
    end.

%%% Bounded execution (SY 7.2, FR-S13): a call on a large input runs in a monitored process with a heap limit and
%%% a wall-clock limit; small inputs (fewer than 100 nodes) and calls with bounded => false run inline.

run(E, Opts, F) ->
    Ctx = context(Opts),
    MaxNodes = maps:get(max_expr_nodes, Ctx),
    case count(E, MaxNodes) of
        N when N > MaxNodes -> {error, {symbolic_limit, #{reason => max_expr_nodes, max => MaxNodes}}};
        N when N < 100 -> guarded(F, Ctx);
        _ ->
            case maps:get(bounded, Opts, true) of
                false -> guarded(F, Ctx);
                true -> bounded(F, Ctx)
            end
    end.

context(Opts) ->
    Keys = [max_poly_terms, max_expand_power, max_expr_nodes, sym_timeout_ms, sym_max_heap_words],
    maps:merge(maps:from_list([{K, pq_core:get(K, Opts)} || K <- Keys]),
               #{domains => maps:get(domains, Opts, #{}), poly_kernel => maps:get(poly_kernel, Opts, true)}).

guarded(F, Ctx) ->
    try F(Ctx) of
        V -> {ok, V}
    catch
        throw:{Cat, _} = R when is_atom(Cat) -> {error, R};
        error:badarith -> {error, {non_finite, arith}}
    end.

bounded(F, Ctx) ->
    Parent = self(),
    Tag = make_ref(),
    {Pid, MRef} = spawn_opt(fun() -> Parent ! {Tag, guarded(F, Ctx#{bounded => false})} end,
                            [monitor, {max_heap_size, #{size => maps:get(sym_max_heap_words, Ctx), kill => true,
                                                        error_logger => false}}]),
    Timeout = maps:get(sym_timeout_ms, Ctx),
    receive
        {Tag, Result} ->
            erlang:demonitor(MRef, [flush]),
            Result;
        {'DOWN', MRef, process, Pid, killed} -> {error, {symbolic_limit, #{reason => heap}}};
        {'DOWN', MRef, process, Pid, Reason} -> {error, {symbolic_limit, #{reason => {crashed, Reason}}}}
    after Timeout ->
        exit(Pid, kill),
        receive {'DOWN', MRef, process, Pid, _} -> ok end,
        receive {Tag, _} -> ok after 0 -> ok end,
        {error, {symbolic_limit, #{reason => timeout}}}
    end.

%% the canonical AST, checked against max_expr_nodes
output(V, Ctx) ->
    Ast = pq_sym_norm:value_to_ast(V),
    Max = maps:get(max_expr_nodes, Ctx),
    count(Ast, Max) > Max andalso throw({symbolic_limit, #{reason => max_expr_nodes, max => Max}}),
    Ast.

%% output/2 for a result that is one polynomial kept in Rust (K2), made in the caller (finish/1): the bounded worker
%% returns the handle, not millions of terms copied out of its heap. finish/1 counts the terms' nodes in Rust and
%% refuses a result over max_expr_nodes before building it, as output/2 does after. The handle's own size and the
%% caller's conversion are bounded by max_poly_terms, which the worker has already checked.
output_deferred({scalar, {rpoly, Vars, R}} = V, Ctx) ->
    case lists:all(fun({sym, _}) -> true; ({expi, _}) -> true; (_) -> false end, Vars) of
        true -> {deferred_rpoly, Vars, R, maps:get(max_expr_nodes, Ctx)};
        false -> output(V, Ctx)
    end;
output_deferred(V, Ctx) -> output(V, Ctx).

finish({ok, {deferred_rpoly, Vars, R, Max}}) ->
    Weights = [case V of {expi, X} -> {true, size({sym, X})}; _ -> {false, size(V)} end || V <- Vars],
    case pq_nif:rpoly_nodes(R, Weights) of
        N when N > Max -> {error, {symbolic_limit, #{reason => max_expr_nodes, max => Max}}};
        _ -> {ok, pq_nif:rpoly_to_ast(R, Vars)}
    end;
finish(Result) -> Result.
