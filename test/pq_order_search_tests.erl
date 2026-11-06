%% The search over Horner variable orders (compile/3 with order_search; pq_sym_search; docs/decisions/H3-order-search.md):
%% never more operations than the order it starts from, values within rounding, and the same program whatever the
%% nodes the candidates are scored on.
-module(pq_order_search_tests).
-include_lib("eunit/include/eunit.hrl").

ops({slp, 1, _, Ops, _}) -> length([x || Op <- Ops, element(1, Op) =/= input, element(1, Op) =/= const]).

close(A, B) -> abs(A - B) =< 1.0e-12 * max(1.0, abs(A)).

%% the search alone, on a score with a known optimum: the inversions against a sorted order
inversions(L) -> length([x || {I, A} <- lists:enumerate(L), {J, B} <- lists:enumerate(L), I < J, A > B]).

search_finds_better_orders_test() ->
    Start = lists:reverse(lists:seq(1, 8)),
    Score = fun inversions/1,
    {Best, Cost} = pq_sym_search:run(Score, Start, #{trials => 400, batch => 8, seed => 3}),
    ?assertEqual(Cost, inversions(Best)),
    ?assertEqual(lists:seq(1, 8), lists:sort(Best)),
    ?assert(Cost < inversions(Start) div 4),
    %% deterministic: the same order twice; another seed, maybe another order, never a worse start
    ?assertEqual({Best, Cost}, pq_sym_search:run(Score, Start, #{trials => 400, batch => 8, seed => 3})),
    {_, Cost2} = pq_sym_search:run(Score, Start, #{trials => 400, batch => 8, seed => 4}),
    ?assert(Cost2 < inversions(Start)),
    %% no trials, or a single atom: the start
    ?assertEqual({Start, inversions(Start)}, pq_sym_search:run(Score, Start, #{trials => 0})),
    ?assertEqual({[x], 0}, pq_sym_search:run(fun(_) -> 0 end, [x], #{trials => 50})).

resultant() ->
    {ok, E} = pq_symbolic:parse(<<"a0^2*b2^3 - a0*a1*b1*b2^2 + a1^2*b0*b2^2 - 2*a0*a2*b0*b2^2 + a0*a2*b1^2*b2 "
                                  "- a1*a2*b0*b1*b2 + a2^2*b0^2*b2 + a0^2*b1*b2*b3 + a3*a1*b0*b3^2 + a3^2*b0^3 "
                                  "- a3*a2*b0^2*b1 + a3*a0*b1^3 - a2*a0*b0*b3^2">>),
    E.

%% compile/3: never more operations than Horner (+ CSE) without the search, the same values
compile_search_test() ->
    E = resultant(),
    Syms = pq_symbolic:free_symbols(E),
    B = maps:from_list([{X, 0.1 * K + 0.3} || {K, X} <- lists:enumerate(Syms)]),
    {ok, V0} = pq_symbolic:evaluate(E, B),
    [begin
         {ok, Plain} = pq_symbolic:compile([E], Syms, #{optimize => O}),
         {ok, S} = pq_symbolic:compile([E], Syms, #{optimize => O, order_search => #{trials => 64, batch => 8}}),
         ?assert(ops(S) =< ops(Plain)),
         {ok, [V]} = pq_symbolic:eval_slp(S, B),
         ?assert(close(V0, V)),
         %% no trials: the program without the search
         ?assertEqual({ok, Plain}, pq_symbolic:compile([E], Syms, #{optimize => O, order_search => #{trials => 0}}))
     end || O <- [horner, horner_cse]],
    %% the option without Horner is ignored
    ?assertEqual(pq_symbolic:compile([E], Syms), pq_symbolic:compile([E], Syms, #{order_search => #{trials => 64}})).

%% the same program with the candidates scored on this node or spread over a peer and this node
nodes_test_() ->
    {setup,
     fun() ->
             os:cmd("epmd -daemon"),
             timer:sleep(300),
             {ok, _} = net_kernel:start(pqsearch, #{name_domain => shortnames}),
             Paths = lists:append([["-pa", P] || P <- code:get_path(), string:find(P, "_build") =/= nomatch]),
             {ok, Pid, Node} = peer:start_link(#{name => pqsearchw, args => Paths}),
             {Pid, Node}
     end,
     fun({Pid, _}) -> catch peer:stop(Pid), net_kernel:stop() end,
     fun({_, Node}) ->
             [?_test(begin
                         E = resultant(),
                         Syms = pq_symbolic:free_symbols(E),
                         Opts = fun(Ns) -> #{optimize => horner_cse,
                                             order_search => #{trials => 48, batch => 6, seed => 7, nodes => Ns}} end,
                         {ok, Local} = pq_symbolic:compile([E], Syms, Opts([node()])),
                         {ok, Spread} = pq_symbolic:compile([E], Syms, Opts([Node, node()])),
                         {ok, Remote} = pq_symbolic:compile([E], Syms, Opts([Node])),
                         ?assertEqual(Local, Spread),
                         ?assertEqual(Local, Remote),
                         %% the tree search too
                         M = fun(Ns) -> #{optimize => horner_cse,
                                          order_search => #{trials => 48, batch => 6, seed => 7, strategy => mcts, nodes => Ns}} end,
                         ?assertEqual(pq_symbolic:compile([E], Syms, M([node()])), pq_symbolic:compile([E], Syms, M([Node, node()])))
                     end)]
     end}.

%% the tree search (strategy mcts), both directions: better orders on the known optimum, deterministic
mcts_test() ->
    Start = lists:reverse(lists:seq(1, 8)),
    [begin
         Opts = #{trials => 400, batch => 8, seed => 5, strategy => mcts, exploration => 0.25, direction => D},
         {Best, Cost} = pq_sym_search:run(fun inversions/1, Start, Opts),
         ?assertEqual(Cost, inversions(Best)),
         ?assertEqual(lists:seq(1, 8), lists:sort(Best)),
         ?assert(Cost < inversions(Start) div 4),
         ?assertEqual({Best, Cost}, pq_sym_search:run(fun inversions/1, Start, Opts))
     end || D <- [forward, backward]],
    %% compile/3 with it: never worse than without, the same values
    E = resultant(),
    Syms = pq_symbolic:free_symbols(E),
    {ok, Plain} = pq_symbolic:compile([E], Syms, #{optimize => horner_cse}),
    {ok, S} = pq_symbolic:compile([E], Syms, #{optimize => horner_cse, order_search => #{trials => 64, batch => 8, strategy => mcts}}),
    ?assert(ops(S) =< ops(Plain)),
    B = maps:from_list([{X, 0.1 * K + 0.3} || {K, X} <- lists:enumerate(Syms)]),
    {ok, [V1]} = pq_symbolic:eval_slp(Plain, B),
    {ok, [V2]} = pq_symbolic:eval_slp(S, B),
    ?assert(close(V1, V2)).
