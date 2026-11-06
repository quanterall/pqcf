%% The search over Horner variable orders (decision record H3; FORM's O3/O4 step): a batched local search over
%% permutations, every batch scored in parallel -- one process per candidate, on the given nodes -- and the same
%% result whatever the number of cores or nodes (D1):
%%
%%   - round R makes Batch candidates from the current order, candidate I by one to three random moves (swap two
%%     atoms, or move one) drawn from a generator seeded with {Seed, R, I} alone;
%%   - every candidate is scored by Score(Order), a deterministic integer (pq_sym_slp: the operations of the
%%     compiled program);
%%   - the candidate of the lowest {score, index} replaces the current order if its score is no worse (equal
%%     scores move along plateaus), and the best order seen is returned with its score.
%%
%% Score runs in the worker processes; on another node, its module must be loaded there (pq_sym_slp is, wherever pq
%% runs). The workers are linked to the caller, so a caller killed for its timeout takes them along.
-module(pq_sym_search).
-export([run/3]).

-type order() :: [term()].

-spec run(fun((order()) -> integer()), order(), map()) -> {order(), integer()}.
run(Score, Order0, Opts) ->
    Trials = maps:get(trials, Opts, 0),
    Batch = max(1, maps:get(batch, Opts, 16)),
    Seed = maps:get(seed, Opts, 1),
    Nodes = case maps:get(nodes, Opts, [node()]) of [] -> [node()]; Ns -> Ns end,
    [Cost0] = pmap(Score, [Order0], Nodes),
    case length(Order0) < 2 of
        true -> {Order0, Cost0};
        false -> rounds(0, (Trials + Batch - 1) div Batch, Batch, Seed, Score, Nodes, {Order0, Cost0}, {Order0, Cost0})
    end.

rounds(R, Rounds, _Batch, _Seed, _Score, _Nodes, _Cur, Best) when R >= Rounds -> Best;
rounds(R, Rounds, Batch, Seed, Score, Nodes, {Order, Cost}, Best) ->
    Cands = [neighbour(Order, {Seed, R, I}) || I <- lists:seq(1, Batch)],
    Costs = pmap(Score, Cands, Nodes),
    [{C, _, O} | _] = lists:sort(lists:zip3(Costs, lists:seq(1, Batch), Cands)),
    Cur = case C =< Cost of true -> {O, C}; false -> {Order, Cost} end,
    Best2 = case C < element(2, Best) of true -> {O, C}; false -> Best end,
    rounds(R + 1, Rounds, Batch, Seed, Score, Nodes, Cur, Best2).

%% one to three moves, each a swap of two positions or a move of one element to another position
neighbour(Order, {A, B, C}) ->
    S0 = rand:seed_s(exsss, {erlang:phash2(A), B + 1, C + 1}),
    N = length(Order),
    {K, S1} = rand:uniform_s(3, S0),
    {O, _} = lists:foldl(fun(_, {L, S}) ->
                                 {I, Sa} = rand:uniform_s(N, S),
                                 {J, Sb} = rand:uniform_s(N, Sa),
                                 {Kind, Sc} = rand:uniform_s(2, Sb),
                                 {move(Kind, I, J, L), Sc}
                         end, {Order, S1}, lists:seq(1, K)),
    O.

move(_, I, I, L) -> L;
move(1, I, J, L) ->
    T = list_to_tuple(L),
    tuple_to_list(setelement(I, setelement(J, T, element(I, T)), element(J, T)));
move(2, I, J, L) ->
    X = lists:nth(I, L),
    {H, T} = lists:split(J - 1, lists:delete(X, L)),
    H ++ [X | T].

%% Score over the items, one linked process each, the items spread over the nodes in turn; results in item order
pmap(Score, Items, Nodes) ->
    Parent = self(),
    NN = length(Nodes),
    Refs = [begin
                Ref = make_ref(),
                Node = lists:nth((I rem NN) + 1, Nodes),
                spawn_link(Node, fun() -> Parent ! {Ref, Score(X)} end),
                Ref
            end || {I, X} <- lists:enumerate(0, Items)],
    [receive {Ref, V} -> V end || Ref <- Refs].
