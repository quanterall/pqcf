%% Common subexpressions across a straight-line program (decision record H2; FORM's O2 step after Horner): the
%% sums and products of pq_sym_slp's program as n-ary nodes, then, greedily, the pair of operands that most nodes
%% share made one operation and substituted into all of them -- sum(a, b, c) and sum(a, b, d) become t = a + b, then
%% t + c and t + d -- until no pair is shared. Sums and products are commutative and associative for complex numbers,
%% so the values are the program's within rounding, not bit for bit. Deterministic: the pair is the one shared by
%% the most nodes, ties broken by term order; nodes are visited in slot order (D1).
-module(pq_sym_cse).
-export([run/1]).

%% a node with more distinct operands than this takes no part in the pair search (its pairs are quadratic in
%% number); it is still emitted, and still gains from pairs found elsewhere
-define(WIDE, 64).

-spec run(tuple()) -> tuple().
run({slp, 1, Inputs, Ops, Outputs}) ->
    OutSlots = case Outputs of {matrix, _, _, Ss} -> Ss; Ss -> Ss end,
    Uses = uses(Ops, OutSlots),
    Nodes = flatten(Ops, Uses),
    {Nodes2, _Next} = extract(Nodes, length(Ops)),
    {NewOps, NewOuts} = emit(Ops, Nodes2, OutSlots),
    {slp, 1, Inputs, NewOps,
     case Outputs of {matrix, R, C, _} -> {matrix, R, C, NewOuts}; _ -> NewOuts end}.

%% how often every slot is an operand or an output
uses(Ops, Outs) ->
    Count = fun(S, M) -> maps:update_with(S, fun(N) -> N + 1 end, 1, M) end,
    lists:foldl(fun(Op, M) -> lists:foldl(Count, M, operands(Op)) end, lists:foldl(Count, #{}, Outs), Ops).

operands({mul, _, A, B}) -> [A, B];
operands({add, _, A, B}) -> [A, B];
operands({sub, _, A, B}) -> [A, B];
operands({pow, _, A, _}) -> [A];
operands({exp, _, A}) -> [A];
operands({conj, _, A}) -> [A];
operands(_) -> [].

%% Slot => {add | mul, SortedOperands} | {op, Op}: a sum (product) used once, only by a sum (product), merges into
%% it; slots ascend, so an operand is flattened before its user
flatten(Ops, Uses) ->
    lists:foldl(
      fun({K, S, A, B}, N) when K =:= add; K =:= mul ->
              %% a merged operand leaves the map: its user holds its operands now (keeping it made a chain of n
              %% sums n nodes of 1..n operands)
              Part = fun(X, {Ls, M}) -> case {M, Uses} of
                                            {#{X := {K, L}}, #{X := 1}} -> {[L | Ls], maps:remove(X, M)};
                                            _ -> {[[X] | Ls], M}
                                        end end,
              {Ls, N2} = lists:foldl(Part, {[], N}, [A, B]),
              N2#{S => {K, lists:merge(Ls)}};
         (Op, N) -> N#{element(2, Op) => {op, Op}}
      end, #{}, Ops).

%%% The pair search

%% Index: {K, A, B} (A =< B) => #{Node => true}; Heap: gb_set of {-Count, Pair} for the pairs of count >= 2, so
%% its smallest element is the pair of the highest count, the smallest in term order among those
extract(Nodes, Next) ->
    Idx = lists:foldl(fun(S, I) ->
                              lists:foldl(fun(P, I2) -> maps:update_with(P, fun(M) -> M#{S => true} end, #{S => true}, I2) end,
                                          I, pairs(maps:get(S, Nodes)))
                      end, #{}, lists:sort(maps:keys(Nodes))),
    Heap = maps:fold(fun(P, Set, H) -> heap_put(P, map_size(Set), H) end, gb_sets:new(), Idx),
    loop(Nodes, Idx, Heap, Next).

loop(Nodes, Idx, Heap, Next) ->
    case gb_sets:is_empty(Heap) of
        true -> {Nodes, Next};
        false ->
            {_, {K, A, B} = P} = gb_sets:smallest(Heap),
            New = Next,
            {Nodes2, Idx2, Heap2} =
                lists:foldl(fun(S, {Ns, I, H}) ->
                                    {K, L} = Old = maps:get(S, Ns),
                                    Upd = {K, lists:merge([New], L -- [A, B])},
                                    {I2, H2} = update(S, Old, Upd, I, H),
                                    {Ns#{S => Upd}, I2, H2}
                            end, {Nodes, Idx, Heap}, lists:sort(maps:keys(maps:get(P, Idx)))),
            loop(Nodes2#{New => {K, [A, B]}}, Idx2, Heap2, Next + 1)
    end.

%% the node's pairs that it lost or gained, and only those, in the index and the heap
update(S, Old, New, Idx, Heap) ->
    PO = lists:usort(pairs(Old)), PN = lists:usort(pairs(New)),
    Change = fun(F) -> fun(P, {I, H}) ->
                               M0 = maps:get(P, I, #{}),
                               M = F(M0),
                               {I#{P => M}, heap_put(P, map_size(M), heap_del(P, map_size(M0), H))}
                       end end,
    {I1, H1} = lists:foldl(Change(fun(M) -> maps:remove(S, M) end), {Idx, Heap}, ordsets:subtract(PO, PN)),
    lists:foldl(Change(fun(M) -> M#{S => true} end), {I1, H1}, ordsets:subtract(PN, PO)).

pairs({K, L}) when K =:= add; K =:= mul ->
    U = lists:usort(L),
    case length(U) > ?WIDE of
        true -> [];
        false -> [{K, A, B} || A <- U, B <- U, A < B] ++ [{K, A, A} || A <- U, length([x || X <- L, X =:= A]) >= 2]
    end;
pairs(_) -> [].

heap_put(P, N, H) when N >= 2 -> gb_sets:add({-N, P}, H);
heap_put(_, _, H) -> H.

heap_del(P, N, H) when N >= 2 -> gb_sets:delete({-N, P}, H);
heap_del(_, _, H) -> H.

%%% Back to a program

%% the inputs first, in their slots' order, then every output's operations depth first; an n-ary node as a chain
%% of binary operations over its operands in slot order, a node of one operand as that operand; identical
%% operations shared
emit(Ops, Nodes, Outs) ->
    St0 = lists:foldl(fun({input, S, X}, St) -> {_, St2} = put_op(S, fun(T) -> {input, T, X} end, {input, X}, St), St2;
                         (_, St) -> St
                      end, {#{}, #{}, [], 0}, Ops),
    {NewOuts, {_, _, Rev, _}} = lists:mapfoldl(fun(S, St) -> node(S, Nodes, St) end, St0, Outs),
    {lists:reverse(Rev), NewOuts}.

%% St = {Old => New, Key => New, ReversedOps, Next}
node(S, Nodes, {Done, _, _, _} = St) ->
    case Done of
        #{S := T} -> {T, St};
        _ ->
            {T, {D, K, R, N}} = build(S, maps:get(S, Nodes), Nodes, St),
            {T, {D#{S => T}, K, R, N}}
    end.

build(S, {op, Op}, Nodes, St) ->
    {New, St2} = lists:mapfoldl(fun(X, Acc) -> node(X, Nodes, Acc) end, St, operands(Op)),
    case {Op, New} of
        {{const, _, Z}, []} -> put_op(S, fun(T) -> {const, T, Z} end, {const, Z}, St2);
        {{input, _, X}, []} -> put_op(S, fun(T) -> {input, T, X} end, {input, X}, St2);
        {{pow, _, _, E}, [A]} -> put_op(S, fun(T) -> {pow, T, A, E} end, {pow, A, E}, St2);
        {{sub, _, _, _}, [A, B]} -> put_op(S, fun(T) -> {sub, T, A, B} end, {sub, A, B}, St2);
        {{exp, _, _}, [A]} -> put_op(S, fun(T) -> {exp, T, A} end, {exp, A}, St2);
        {{conj, _, _}, [A]} -> put_op(S, fun(T) -> {conj, T, A} end, {conj, A}, St2)
    end;
build(S, {K, L}, Nodes, St) ->
    {[First | Rest], St2} = lists:mapfoldl(fun(X, Acc) -> node(X, Nodes, Acc) end, St, L),
    lists:foldl(fun(B, {A, Acc}) ->
                        put_op(S, fun(T) -> {K, T, A, B} end, {K, min(A, B), max(A, B)}, Acc)
                end, {First, St2}, Rest).

put_op(_S, Make, Key, {D, Keys, R, N}) ->
    case Keys of
        #{Key := T} -> {T, {D, Keys, R, N}};
        _ -> {N, {D, Keys#{Key => N}, [Make(N) | R], N + 1}}
    end.
