%% Straight-line programs (SY 5.5, FR-S10): canonical expressions compiled into ops on slots, with common
%% subexpressions shared through a map from structural key to slot, so every atom power, exponential factor and
%% monomial is computed once for all outputs. The program is plain data (it can be sent to other nodes); eval runs
%% the ops in order, so the operation order is fixed (D1).
%%
%%   {slp, 1, Inputs, Ops, Outputs}
%%   op() :: {const, Slot, {Re, Im}} | {input, Slot, Sym} | {mul, Slot, A, B} | {add, Slot, A, B}
%%         | {pow, Slot, A, {N, D}} | {exp, Slot, A} | {conj, Slot, A}
-module(pq_sym_slp).
-export([compile/2, eval/2, eval_batch/2]).

-record(st, {ops = [], keys = #{}, next = 0, ctx}).

-spec compile([pq_sym_norm:value()], [atom()]) -> {slp, 1, [atom()], [tuple()], [non_neg_integer()] | {matrix, pos_integer(), pos_integer(), [non_neg_integer()]}}.
compile(Values0, Inputs) ->
    Values = [expand_identity(V) || V <- Values0],
    St0 = lists:foldl(fun(X, S) -> element(2, slot({input, X}, fun(Sl) -> {input, Sl, X} end, S)) end, #st{}, Inputs),
    {Outs, St} = lists:mapfoldl(fun out/2, St0, Values),
    Outputs = case {Values, Outs} of
                  {[{mat, R, C, _}], [Slots]} -> {matrix, R, C, Slots};
                  _ -> lists:append(Outs)
              end,
    {slp, 1, Inputs, lists:reverse(St#st.ops), Outputs}.

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

%% a polynomial: the sum of coefficient times monomial, in canonical key order
poly(P, S) when map_size(P) =:= 0 -> const({0.0, 0.0}, S);
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
