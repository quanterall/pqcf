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
         matmul/2, kron/2, smul/2, madd/1, dagger/1, transpose/1,
         shape/1, free_symbols/1, size/1, hash/1, domain_of/2,
         pretty/1, pretty/2, to_latex/1, parse/1, parse/2,
         resolve_params/2, resolve_params/3,
         substitute/2]).

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
    run(E, Opts, fun(Ctx) -> output(pq_sym_norm:norm(E, Ctx), Ctx) end).

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
               #{domains => maps:get(domains, Opts, #{})}).

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
