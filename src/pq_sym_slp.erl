%% Straight-line programs (SY 5.5, FR-S10): canonical expressions compiled into ops on slots, with common
%% subexpressions shared through a map from structural key to slot, so every atom power, exponential factor and
%% monomial is computed once for all outputs. The program is plain data (it can be sent to other nodes); eval runs
%% the ops in order, so the operation order is fixed (D1).
%%
%%   {slp, 1, Inputs, Ops, Outputs}
%%   op() :: {const, Slot, {Re, Im}} | {input, Slot, Sym} | {mul, Slot, A, B} | {add, Slot, A, B}
%%         | {pow, Slot, A, {N, D}} | {exp, Slot, A} | {conj, Slot, A} | {sub, Slot, A, B} (A - B; Horner only)
-module(pq_sym_slp).
-export([compile/2, compile/3, eval/2, eval_batch/2]).

-record(st, {ops = [], keys = #{}, next = 0, ctx, opt = none, order = []}).

-spec compile([pq_sym_norm:value()], [atom()]) -> {slp, 1, [atom()], [tuple()], [non_neg_integer()] | {matrix, pos_integer(), pos_integer(), [non_neg_integer()]}}.
compile(Values, Inputs) -> compile(Values, Inputs, #{}).

%% Opts: optimize => none (the default: the program of compile/2, operation for operation) | horner (1.3;
%% docs/decisions/H1-horner.md): every polynomial as a multivariate Horner scheme -- P = R + x^k Q on a variable that
%% several terms share, recursively, in one occurrence order fixed for the polynomial (most terms first, ties by
%% canonical order, so the program is deterministic; it gave 4-6 % fewer operations than choosing greedily at every
%% step on the resultant benchmarks), identical sub-polynomials computed once. Fewer operations, a different
%% operation order: results equal within rounding, not bit for bit, to the default's. cse | horner_cse (1.3;
%% docs/decisions/H2-cse.md): the default program or the Horner scheme, then common subexpressions across it
%% (pq_sym_cse: sums and products as n-ary nodes, the operand pair most of them share made one operation, until
%% none is shared) -- FORM's O2 step; also within rounding. order_search => #{trials => N, seed, batch, nodes} with
%% horner | horner_cse (1.3; docs/decisions/H3-order-search.md): N candidate variable orders, one for all the
%% polynomials, searched in parallel batches (pq_sym_search), the order of the fewest operations kept -- the same
%% program for any number of cores or nodes.
-spec compile([pq_sym_norm:value()], [atom()], map()) -> tuple().
compile(Values0, Inputs, Opts) ->
    Values = [expand_identity(V) || V <- Values0],
    Opt = maps:get(optimize, Opts, none),
    Search = maps:get(order_search, Opts, #{}),
    Order = case (Opt =:= horner orelse Opt =:= horner_cse) andalso maps:get(trials, Search, 0) > 0 of
                true ->
                    %% one order for every polynomial, its tail each polynomial's own occurrence order
                    Score = fun(O) -> ops(build(Values, Inputs, Opt, O)) end,
                    {Best, Cost} = pq_sym_search:run(Score, global_order(Values), Search),
                    %% never worse than each polynomial's own order
                    case Cost < Score([]) of true -> Best; false -> [] end;
                false -> []
            end,
    build(Values, Inputs, Opt, Order).

%% the operations of a program (inputs and constants are not operations)
ops({slp, 1, _, Ops, _}) -> length([x || Op <- Ops, element(1, Op) =/= input, element(1, Op) =/= const]).

%% the atoms Horner splits on, over all the outputs, by the number of terms they occur in
global_order(Values) ->
    Ps = lists:append([case V of {scalar, P} -> [P]; {mat, _, _, Rows} -> lists:append(Rows) end || V <- Values]),
    Counts = lists:foldl(fun(P, Acc) -> maps:fold(fun(K, N, A) -> maps:update_with(K, fun(M) -> M + N end, N, A) end,
                                                 Acc, counts(P)) end, #{}, Ps),
    [K || {_, K} <- lists:sort([{-N, K} || K := N <- Counts])].

build(Values, Inputs, Opt, Order) ->
    St0 = lists:foldl(fun(X, S) -> element(2, slot({input, X}, fun(Sl) -> {input, Sl, X} end, S)) end,
                      #st{opt = Opt, order = Order}, Inputs),
    %% horner_cse: the Horner scheme, then the pair search over it
    {Outs, St} = lists:mapfoldl(fun out/2, St0, Values),
    Outputs = case {Values, Outs} of
                  {[{mat, R, C, _}], [Slots]} -> {matrix, R, C, Slots};
                  _ -> lists:append(Outs)
              end,
    Prog = {slp, 1, Inputs, lists:reverse(St#st.ops), Outputs},
    case St#st.opt of
        cse -> pq_sym_cse:run(Prog);
        horner_cse -> pq_sym_cse:run(Prog);
        _ -> Prog
    end.

expand_identity({identity, N}) ->
    {mat, N, N, [[case I =:= J of true -> pq_sym_poly:one(); false -> #{} end || J <- lists:seq(1, N)] || I <- lists:seq(1, N)]};
expand_identity(V) -> V.

out({scalar, P}, S) -> {Sl, S2} = poly(P, S), {[Sl], S2};
out({mat, _, _, Rows}, S) -> lists:mapfoldl(fun(P, Acc) -> poly(P, Acc) end, S, lists:append(Rows)).

%% a slot for a structural key, made once
slot(Key, Make, #st{keys = K} = S) ->
    case K of
        #{Key := Sl} -> {Sl, S};
        _ ->
            Sl = S#st.next,
            {Sl, S#st{ops = [Make(Sl) | S#st.ops], keys = K#{Key => Sl}, next = Sl + 1}}
    end.

const(Z, S) -> slot({const, Z}, fun(Sl) -> {const, Sl, Z} end, S).

%% a polynomial: the sum of coefficient times monomial, in canonical key order (or a Horner scheme)
poly(P, S) when map_size(P) =:= 0 -> const({0.0, 0.0}, S);
poly(P, #st{opt = Opt} = S) when map_size(P) > 1, Opt =:= horner orelse Opt =:= horner_cse -> horner(P, S#st.order ++ order(P), S);
poly(P, S) ->
    case maps:find({poly, P}, S#st.keys) of
        {ok, Sl} -> {Sl, S};
        error ->
            {Terms, S2} = lists:mapfoldl(fun({M, C}, Acc) -> term(M, C, Acc) end, S, pq_sym_poly:sorted(P)),
            {Sum, S3} = lists:foldl(fun(T, {none, A}) -> {T, A};
                                       (T, {Prev, A}) -> slot({add, Prev, T}, fun(Sl) -> {add, Sl, Prev, T} end, A)
                                    end, {none, S2}, Terms),
            {Sum, S3#st{keys = (S3#st.keys)#{{poly, P} => Sum}}}
    end.

%%% Horner schemes

%% Every node of the scheme is a primitive polynomial: P = G Prim, G its content (pq_sym_poly:content/2, the sign
%% included, so the first term of Prim is positive), and the slot computes Prim; -(a - b), 2 (a - b) and a - b share
%% one slot (FORM extracts the same content at every level of its Horner trees). A parent folds the signs into add /
%% sub and multiplies by a content only where it is not +-1.
horner(P, Order, S) ->
    {Sl, G, S2} = hc(P, Order, S),
    scaled(G, Sl, S2).

hc(P, _Order, S) when map_size(P) =:= 0 -> {Sl, S2} = const({0.0, 0.0}, S), {Sl, {1, 1}, S2};
hc(P, Order, S) ->
    {G, Prim} = pq_sym_poly:content(P, true),
    case maps:find({hpoly, Prim}, S#st.keys) of
        {ok, Sl} -> {Sl, G, S};
        error ->
            {Sl, S2} = hprim(Prim, Order, S),
            {Sl, G, S2#st{keys = (S2#st.keys)#{{hpoly, Prim} => Sl}}}
    end.

hprim(P, _Order, S) when map_size(P) =:= 1 -> [{M, C}] = maps:to_list(P), term(M, C, S);
hprim(P, Order, S) ->
    case split_var(P, Order) of
        none -> sum_terms(P, S);
        X ->
            {With, Without} = maps:fold(fun({A, _} = M, C, {W, O}) ->
                                                case int_exp(X, A) of
                                                    0 -> {W, O#{M => C}};
                                                    _ -> {W#{M => C}, O}
                                                end
                                        end, {#{}, #{}}, P),
            K = lists:min([int_exp(X, A) || {A, _} <- maps:keys(With)]),
            Q = maps:fold(fun({A, E}, C, Acc) -> Acc#{{lower(X, K, A), E} => C} end, #{}, With),
            {Sx, S3} = atom_pow({X, {K, 1}}, S),
            {Sq, Gq, S4} = hc(Q, Order, S3),
            {Prod, S5} = slot({mul, Sx, Sq}, fun(Sl0) -> {mul, Sl0, Sx, Sq} end, S4),
            case map_size(Without) of
                0 -> scaled(Gq, Prod, S5);
                _ ->
                    {Sr, Gr, S6} = hc(Without, Order, S5),
                    combine(Gr, Sr, Gq, Prod, S6)
            end
    end.

%% Gr Sr + Gq Sp as one slot: |G| multiplied in where it is not 1, the signs as add or sub
combine(Gr, Sr, Gq, Sp, S) ->
    {SignR, Ar, S1} = magnitude(Gr, Sr, S),
    {SignQ, Aq, S2} = magnitude(Gq, Sp, S1),
    case {SignR, SignQ} of
        {1, 1} -> slot({add, Ar, Aq}, fun(Sl) -> {add, Sl, Ar, Aq} end, S2);
        {1, -1} -> slot({sub, Ar, Aq}, fun(Sl) -> {sub, Sl, Ar, Aq} end, S2);
        {-1, 1} -> slot({sub, Aq, Ar}, fun(Sl) -> {sub, Sl, Aq, Ar} end, S2);
        {-1, -1} ->
            {T, S3} = slot({add, Ar, Aq}, fun(Sl) -> {add, Sl, Ar, Aq} end, S2),
            scaled({-1, 1}, T, S3)
    end.

magnitude({N, D}, Sl, S) ->
    Sign = case N < 0 of true -> -1; false -> 1 end,
    {Sl2, S2} = scaled({abs(N), D}, Sl, S),
    {Sign, Sl2, S2}.

%% G times the slot (G a rational)
scaled({1, 1}, Sl, S) -> {Sl, S};
scaled({N, D}, Sl, S) ->
    {Cs, S2} = const({N / D, 0.0}, S),
    slot({mul, Cs, Sl}, fun(T) -> {mul, T, Cs, Sl} end, S2).

%% the atoms with a positive integer exponent, by the number of terms they occur in (most first; ties by canonical order)
order(P) -> [K || {_, K} <- lists:sort([{-N, K} || K := N <- counts(P)])].

counts(P) ->
    maps:fold(fun({A, _}, _, Acc) ->
                      lists:foldl(fun({K, {N, 1}}, Ac) when is_integer(N), N > 0 -> maps:update_with(K, fun(V) -> V + 1 end, 1, Ac);
                                     (_, Ac) -> Ac
                                  end, Acc, A)
              end, #{}, P).

%% the first atom of the order that at least two terms of P have
split_var(P, Order) ->
    C = counts(P),
    case [K || K <- Order, maps:get(K, C, 0) >= 2] of
        [] -> none;
        [K | _] -> K
    end.

int_exp(X, A) ->
    case lists:keyfind(X, 1, A) of
        {X, {N, 1}} when is_integer(N), N > 0 -> N;
        _ -> 0
    end.

lower(X, K, A) ->
    [Y || Y <- [case Y0 of {X, {N, 1}} -> {X, {N - K, 1}}; _ -> Y0 end || Y0 <- A], element(2, Y) =/= {0, 1}].

%% the sum of the terms of a Horner leaf (no atom left to split on), in canonical key order: the terms with a
%% positive real coefficient added, then those with a negative one subtracted with its magnitude (a subtraction
%% costs what an addition does; a coefficient -1 no longer costs a product)
sum_terms(P, S) ->
    {Signed, S2} = lists:mapfoldl(fun({M, {q, {N, D}, {0, 1}}}, Acc) when N < 0 ->
                                          {T, Acc2} = term(M, {q, {-N, D}, {0, 1}}, Acc),
                                          {{-1, T}, Acc2};
                                     ({M, C}, Acc) ->
                                          {T, Acc2} = term(M, C, Acc),
                                          {{1, T}, Acc2}
                                  end, S, pq_sym_poly:sorted(P)),
    Pos = [T || {1, T} <- Signed], Neg = [T || {-1, T} <- Signed],
    Add = fun(T, {none, A}) -> {T, A};
             (T, {Prev, A}) -> slot({add, Prev, T}, fun(Sl) -> {add, Sl, Prev, T} end, A)
          end,
    case lists:foldl(Add, {none, S2}, Pos) of
        {none, S3} ->
            {Sum, S4} = lists:foldl(Add, {none, S3}, Neg),
            scaled({-1, 1}, Sum, S4);
        Acc0 ->
            lists:foldl(fun(T, {Prev, A}) -> slot({sub, Prev, T}, fun(Sl) -> {sub, Sl, Prev, T} end, A) end, Acc0, Neg)
    end.

term(M, C, S) ->
    Z = pq_sym_exact:to_complex(C),
    case M =:= pq_sym_poly:one_mono() of
        true -> const(Z, S);
        false ->
            {Ms, S2} = mono(M, S),
            case Z of
                {1.0, +0.0} -> {Ms, S2};
                _ ->
                    {Cs, S3} = const(Z, S2),
                    slot({mul, Cs, Ms}, fun(Sl) -> {mul, Sl, Cs, Ms} end, S3)
            end
    end.

mono({Alg, Exp} = M, S) ->
    case maps:find({mono, M}, S#st.keys) of
        {ok, Sl} -> {Sl, S};
        error ->
            {Fs, S2} = lists:mapfoldl(fun atom_pow/2, S, Alg),
            {Es, S3} = lists:mapfoldl(fun exp_factor/2, S2, Exp),
            {Prod, S4} = lists:foldl(fun(F, {none, A}) -> {F, A};
                                        (F, {Prev, A}) -> slot({mul, Prev, F}, fun(Sl) -> {mul, Sl, Prev, F} end, A)
                                     end, {none, S3}, Fs ++ Es),
            {Prod, S4#st{keys = (S4#st.keys)#{{mono, M} => Prod}}}
    end.

atom_pow({{cbase, B}, E}, S) ->
    {Bs, S2} = poly(B, S),
    {P, S3} = slot({pow, Bs, E}, fun(Sl) -> {pow, Sl, Bs, E} end, S2),
    slot({conj, P}, fun(Sl) -> {conj, Sl, P} end, S3);
atom_pow({K, E}, S) ->
    {B, S2} = atom(K, S),
    case E of
        {1, 1} -> {B, S2};
        _ -> slot({pow, B, E}, fun(Sl) -> {pow, Sl, B, E} end, S2)
    end.

atom({sym, X}, S) -> slot({input, X}, fun(Sl) -> {input, Sl, X} end, S);
atom({conj, X}, S) ->
    {I, S2} = atom({sym, X}, S),
    slot({conj, I}, fun(Sl) -> {conj, Sl, I} end, S2);
atom(pi, S) -> const({math:pi(), 0.0}, S);
atom({rad, P}, S) -> const({math:sqrt(P), 0.0}, S);
atom({base, B}, S) -> poly(B, S).

exp_factor({B, C}, S) ->
    {Bs, S2} = case B =:= pq_sym_poly:one_mono() of
                   true -> const({1.0, 0.0}, S);
                   false -> mono(B, S)
               end,
    {Cs, S3} = const(pq_sym_exact:to_complex(C), S2),
    {Arg, S4} = slot({mul, Cs, Bs}, fun(Sl) -> {mul, Sl, Cs, Bs} end, S3),
    slot({exp, Arg}, fun(Sl) -> {exp, Sl, Arg} end, S4).

%%% Evaluation

-spec eval(tuple(), map()) -> {ok, [float() | {complex, float(), float()}] | pq_core:packed_cmatrix()} | {error, term()}.
eval({slp, 1, Inputs, Ops, Outputs}, Binding) ->
    case [X || X <- Inputs, not maps:is_key(X, Binding)] of
        [] ->
            try
                Vals = run(Ops, Binding, #{}),
                {ok, result(Outputs, Vals)}
            catch
                error:badarith -> {error, {non_finite, arith}};
                throw:{non_finite, _} = R -> {error, R}
            end;
        Missing -> {error, {unbound_symbols, lists:sort(Missing)}}
    end.

-spec eval_batch(tuple(), [map()]) -> {ok, list()} | {error, term()}.
eval_batch(P, Bindings) ->
    Rs = [eval(P, B) || B <- Bindings],
    case [E || {error, _} = E <- Rs] of
        [] -> {ok, [V || {ok, V} <- Rs]};
        [E | _] -> E
    end.

run([], _B, V) -> V;
run([Op | T], B, V) -> run(T, B, step(Op, B, V)).

step({const, S, Z}, _B, V) -> V#{S => Z};
step({input, S, X}, B, V) -> V#{S => pq_sym_eval:value(maps:get(X, B))};
step({mul, S, A, C}, _B, V) -> {Ar, Ai} = maps:get(A, V), {Cr, Ci} = maps:get(C, V), V#{S => {Ar * Cr - Ai * Ci, Ar * Ci + Ai * Cr}};
step({add, S, A, C}, _B, V) -> {Ar, Ai} = maps:get(A, V), {Cr, Ci} = maps:get(C, V), V#{S => {Ar + Cr, Ai + Ci}};
step({sub, S, A, C}, _B, V) -> {Ar, Ai} = maps:get(A, V), {Cr, Ci} = maps:get(C, V), V#{S => {Ar - Cr, Ai - Ci}};
step({pow, S, A, E}, _B, V) -> V#{S => pq_sym_eval:rpow(maps:get(A, V), E)};
step({exp, S, A}, _B, V) -> {R, I} = maps:get(A, V), M = math:exp(R), V#{S => {M * math:cos(I), M * math:sin(I)}};
step({conj, S, A}, _B, V) -> {R, I} = maps:get(A, V), V#{S => {R, -I}}.

result({matrix, R, C, Slots}, V) ->
    Rows = chunk([begin {Re, Im} = maps:get(S, V), {complex, Re, Im} end || S <- Slots], C),
    _ = R,
    {ok, P} = pq_core:pack({cmat, Rows}),
    P;
result(Slots, V) -> [out_value(maps:get(S, V)) || S <- Slots].

out_value({Re, Im}) when Im == 0 -> Re + 0.0;
out_value({Re, Im}) -> {complex, Re + 0.0, Im + 0.0}.

chunk([], _) -> [];
chunk(L, K) -> {A, B} = lists:split(K, L), [A | chunk(B, K)].
