%% The print order by binary keys (pq_sym_norm:to_ast/1; decision record K1): the same order as the general path,
%% which a term with pi forces -- so the output for P must equal the output for P + pi, the pi term removed.
-module(pq_print_order_tests).
-include_lib("eunit/include/eunit.hrl").

-define(PI_TERM, {mul, [7, {const, pi}]}).

general(P) ->
    {add, Ts} = pq_sym_norm:to_ast(pq_sym_poly:add(P, #{{[{pi, {1, 1}}], []} => {q, {7, 1}, {0, 1}}})),
    case [T || T <- Ts, T =/= ?PI_TERM] of [T] -> T; L -> {add, L} end.

%% names chosen to stress the order: prefixes, digits, capitals, underscores, non-ASCII
names() -> [a, ab, b, a1, a10, a2, 'A', z, x_1, x, 'λ', aa].

random_poly(N, Lo, Hi) ->
    Ns = names(),
    maps:from_list([{{lists:usort([{{sym, lists:nth(rand:uniform(length(Ns)), Ns)}, {E, 1}}
                                   || _ <- lists:seq(1, rand:uniform(4)), E <- [Lo + rand:uniform(Hi - Lo + 1) - 1], E =/= 0])
                      , []}, {q, {rand:uniform(19) - 10, 1}, {0, 1}}} || _ <- lists:seq(1, N)]).

%% usort above keeps one exponent per name only if they're equal; drop monomials with a name twice
clean(P) -> maps:filter(fun({Alg, _}, {q, {C, _}, _}) -> C =/= 0 andalso length(lists:ukeysort(1, Alg)) =:= length(Alg) end, P).

same_order_test() ->
    rand:seed(exsss, {5, 6, 7}),
    [begin
         P = clean(random_poly(2 + K * 3, case K rem 2 of 0 -> -3; 1 -> 0 end, 5)),
         case map_size(P) of
             0 -> ok;
             _ -> ?assertEqual(general(P), pq_sym_norm:to_ast(P))
         end
     end || K <- lists:seq(1, 80)].

%% a constant term, and a polynomial of one term
edges_test() ->
    P = #{{[], []} => {q, {3, 1}, {0, 1}}, {[{{sym, x}, {2, 1}}], []} => {q, {1, 1}, {0, 1}},
          {[{{sym, x}, {-1, 1}}], []} => {q, {2, 1}, {0, 1}}},
    ?assertEqual(general(P), pq_sym_norm:to_ast(P)),
    ?assertEqual({add, [3, {mul, [2, {pow, {sym, x}, -1}]}, {mul, [1, {pow, {sym, x}, 2}]}]}, pq_sym_norm:to_ast(P)),
    ?assertEqual({mul, [1, {sym, y}]}, pq_sym_norm:to_ast(#{{[{{sym, y}, {1, 1}}], []} => {q, {1, 1}, {0, 1}}})).
