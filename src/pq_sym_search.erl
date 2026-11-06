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
-export_type([order/0]).

-type order() :: [term()].

-spec run(fun((order()) -> integer()), order(), map()) -> {order(), integer()}.
run(Score, Order0, Opts) ->
    case maps:get(strategy, Opts, local) of
        local -> local(Score, Order0, Opts);
        mcts -> mcts(Score, Order0, Opts)
    end.

local(Score, Order0, Opts) ->
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

%%% Monte Carlo tree search over order prefixes (FORM's O3; Kuipers, Plaat, Vermaseren, van den Herik 2013)
%%
%% A node is a prefix of the order (the atoms chosen so far, first to last); its children add one of the remaining
%% atoms. A playout descends from the root by UCT -- mean reward + C sqrt(2 ln N / n), the reward of a cost c being
%% (Max - c) / (Max - Min) over the costs seen so far -- to the first child never visited (children in the order
%% of Order0), and completes the order by a shuffle of the rest seeded with {Seed, playout}. A batch of Batch
%% playouts is chosen one after the other, each pending one counted as a visit of the worst reward so that they
%% spread, then scored in parallel, then applied in playout order: the same tree whatever the cores or nodes.

%% direction backward: the tree chooses the order from its end (the innermost Horner levels first): the search runs
%% on reversed orders, scored reversed back
mcts(Score, Order0, Opts) ->
    case maps:get(direction, Opts, forward) of
        forward -> mcts1(Score, Order0, Opts);
        backward ->
            {O, C} = mcts1(fun(Rev) -> Score(lists:reverse(Rev)) end, lists:reverse(Order0), Opts),
            {lists:reverse(O), C}
    end.

mcts1(Score, Order0, Opts) ->
    Trials = maps:get(trials, Opts, 0),
    Batch = max(1, maps:get(batch, Opts, 16)),
    Seed = maps:get(seed, Opts, 1),
    C = maps:get(exploration, Opts, 1.0),
    Nodes = case maps:get(nodes, Opts, [node()]) of [] -> [node()]; Ns -> Ns end,
    [Cost0] = pmap(Score, [Order0], Nodes),
    St0 = #{tree => #{}, min => Cost0, max => Cost0, best => {Order0, Cost0}},
    St = mcts_rounds(0, (Trials + Batch - 1) div Batch, Batch, Seed, C, Order0, Score, Nodes, St0),
    maps:get(best, St).

mcts_rounds(R, Rounds, _B, _Seed, _C, _O0, _Score, _Nodes, St) when R >= Rounds -> St;
mcts_rounds(R, Rounds, B, Seed, C, O0, Score, Nodes, St) ->
    %% choose the batch: each path's nodes get a pending visit
    {Plays, Pend} = lists:mapfoldl(fun(I, Pending) ->
                                           Path = select(O0, [], St, Pending, C),
                                           Order = playout(O0, Path, {Seed, R * B + I}),
                                           {{Path, Order}, lists:foldl(fun(P, M) -> maps:update_with(P, fun(V) -> V + 1 end, 1, M) end,
                                                                       Pending, prefixes(Path))}
                                   end, #{}, lists:seq(1, B)),
    _ = Pend,
    Costs = pmap(Score, [O || {_, O} <- Plays], Nodes),
    St2 = lists:foldl(fun({{Path, Order}, Cost}, S) -> backprop(Path, Order, Cost, S) end, St, lists:zip(Plays, Costs)),
    mcts_rounds(R + 1, Rounds, B, Seed, C, O0, Score, Nodes, St2).

%% the prefixes of a path, the root's empty one excluded
prefixes(Path) -> [lists:sublist(Path, K) || K <- lists:seq(1, length(Path))].

%% descend by UCT to the first unvisited child
select(O0, Path, St, Pending, C) ->
    case O0 -- Path of
        [] -> Path;
        Rest ->
            Tree = maps:get(tree, St),
            Visits = fun(P) -> {N, _} = maps:get(P, Tree, {0, 0}), N + maps:get(P, Pending, 0) end,
            Kids = [Path ++ [X] || X <- Rest],
            case [K || K <- Kids, Visits(K) =:= 0] of
                [K | _] -> K;
                [] ->
                    Total = lists:sum([Visits(K) || K <- Kids]),
                    {Min, Max} = {maps:get(min, St), maps:get(max, St)},
                    Value = fun(K) ->
                                    {N, Sum} = maps:get(K, Tree, {0, 0}),
                                    P = maps:get(K, Pending, 0),
                                    Mean = (Sum + P * Max) / (N + P),
                                    Reward = case Max > Min of true -> (Max - Mean) / (Max - Min); false -> 0.5 end,
                                    Reward + C * math:sqrt(2 * math:log(Total) / (N + P))
                            end,
                    {_, Best} = lists:max([{Value(K), -I} || {I, K} <- lists:enumerate(Kids)]),
                    select(O0, lists:nth(-Best, Kids), St, Pending, C)
            end
    end.

playout(O0, Path, {A, B}) ->
    S0 = rand:seed_s(exsss, {erlang:phash2(A), B + 1, 7}),
    {Ks, _} = lists:mapfoldl(fun(_, S) -> rand:uniform_s(S) end, S0, O0 -- Path),
    Path ++ [X || {_, X} <- lists:sort(lists:zip(Ks, O0 -- Path))].

backprop(Path, Order, Cost, St) ->
    Tree = lists:foldl(fun(P, T) -> maps:update_with(P, fun({N, S}) -> {N + 1, S + Cost} end, {1, Cost}, T) end,
                       maps:get(tree, St), prefixes(Path)),
    Best = case Cost < element(2, maps:get(best, St)) of true -> {Order, Cost}; false -> maps:get(best, St) end,
    St#{tree := Tree, min := min(Cost, maps:get(min, St)), max := max(Cost, maps:get(max, St)), best := Best}.
