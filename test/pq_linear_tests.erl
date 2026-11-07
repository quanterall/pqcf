%% Linear systems solved by reconstruction (pq_recon:linear/4; decision record R1, phase 5): each solution component
%% a reduced rational function from eliminations mod p, never a symbolic elimination; checked exactly: A u = b.
-module(pq_linear_tests).
-include_lib("eunit/include/eunit.hrl").

x() -> {sym, x}.
y() -> {sym, y}.
z() -> {sym, z}.

%% A u = b exactly at random rational points: every entry and every N, D evaluated exactly (the canonical form can't
%% combine fractions with different denominators, so a symbolic A u - b wouldn't reduce to 0)
satisfied(A, B, Sol) ->
    Vars = lists:usort([X || {sym, X} <- lists:flatten([tuple_leaves(E) || E <- lists:append(A) ++ B])] ++ [x]),
    U = [{N, D} || {_, {N, D}} <- lists:sort(Sol)],
    [begin
         Bind = maps:from_list([{X, {rat, rand:uniform(97) - 48, rand:uniform(13)}} || X <- Vars]),
         Ev = fun(E) -> {ok, V} = pq_symbolic:evaluate(E, Bind, #{exact => true}), q(V) end,
         Us = [qdiv(Ev(pq_sym_norm:to_ast(N)), Ev(pq_sym_norm:to_ast(D))) || {N, D} <- U],
         [?assertEqual(Ev(Bi), lists:foldl(fun({Aij, Uj}, Acc) -> qadd(Acc, qmul(Ev(Aij), Uj)) end, {0, 1}, lists:zip(Row, Us)))
          || {Row, Bi} <- lists:zip(A, B)]
     end || _ <- lists:seq(1, 4)].

tuple_leaves(T) when is_tuple(T) -> [T | lists:flatten([tuple_leaves(E) || E <- tuple_to_list(T)])];
tuple_leaves(L) when is_list(L) -> lists:flatten([tuple_leaves(E) || E <- L]);
tuple_leaves(X) -> [X].

q(N) when is_integer(N) -> {N, 1};
q({rat, N, D}) -> {N, D}.
qadd(A, B) -> pq_sym_exact:radd(A, B).
qmul(A, B) -> pq_sym_exact:rmul(A, B).
qdiv(A, B) -> pq_sym_exact:rmul(A, pq_sym_exact:rinv(B)).

two_by_two_test() ->
    A = [[x(), 1], [1, y()]], B = [1, 0],
    {ok, Sol} = pq_recon:linear(A, B, [x, y], #{}),
    satisfied(A, B, Sol),
    %% u_0 = y / (x y - 1): one-term numerator, two-term denominator
    {0, {N0, D0}} = lists:keyfind(0, 1, Sol),
    ?assertEqual({1, 2}, {maps:size(N0), maps:size(D0)}).

%% a sparse 5 x 5 system with linear entries, like an integration-by-parts system, in three variables
sparse_test_() ->
    {timeout, 300, fun() ->
                           rand:seed(exsss, {2, 0, 2}),
                           Lin = fun() -> {add, [rand:uniform(7) - 4, {mul, [rand:uniform(5) - 3, x()]}, {mul, [rand:uniform(5) - 3, y()]},
                                                  {mul, [rand:uniform(5) - 3, z()]}]} end,
                           N = 5,
                           A = [[case I =:= J orelse rand:uniform(3) =:= 1 of true -> Lin(); false -> 0 end
                                 || J <- lists:seq(1, N)] || I <- lists:seq(1, N)],
                           B = [rand:uniform(5) - 3 || _ <- lists:seq(1, N)],
                           case pq_recon:linear(A, B, [x, y, z], #{}) of
                               {ok, Sol} -> satisfied(A, B, Sol), ?assertEqual(N, length(Sol));
                               Other -> ?assertEqual(ok, Other)
                           end
                   end}.

%% only the components asked for; the same on peer-free executors with other pools
components_test() ->
    A = [[x(), y(), 0], [0, x(), z()], [y(), 0, {add, [x(), 1]}]], B = [1, {rat, 1, 2}, z()],
    {ok, [{2, _}] = One} = pq_recon:linear(A, B, [x, y, z], #{components => [2]}),
    {ok, All} = pq_recon:linear(A, B, [x, y, z], #{workers => 1}),
    ?assertEqual(lists:keyfind(2, 1, All), hd(One)),
    satisfied(A, B, All).

%% a system that is singular for every x has no solution to reconstruct; a malformed one is refused
refusals_test() ->
    ?assertMatch({error, _}, pq_recon:linear([[x(), x()], [x(), x()]], [1, 2], [x], #{})),
    ?assertMatch({error, {bad_shape, _}}, pq_recon:linear([[x(), 1]], [1], [x], #{})).

%% several right-hand sides at once: one elimination per point for all of them, the zero ones found by the shared
%% test, each solution the one a single right-hand side gives
multi_rhs_test() ->
    A = [[x(), y(), 0], [0, x(), z()], [y(), 0, {add, [x(), 1]}]],
    B1 = [1, {rat, 1, 2}, z()], B2 = [0, 0, 0], B3 = [x(), 0, 1],
    Bs = [[P, Q, R] || {P, Q, R} <- lists:zip3(B1, B2, B3)],
    {ok, [{2, PerRhs}]} = pq_recon:linear(A, Bs, [x, y, z], #{components => [2]}),
    ?assertEqual([0, 1, 2], [C || {C, _} <- PerRhs]),
    {ok, [{2, S1}]} = pq_recon:linear(A, B1, [x, y, z], #{components => [2]}),
    {ok, [{2, S3}]} = pq_recon:linear(A, B3, [x, y, z], #{components => [2]}),
    ?assertEqual({S1, S3}, {element(2, lists:keyfind(0, 1, PerRhs)), element(2, lists:keyfind(2, 1, PerRhs))}),
    ?assertEqual({#{}, #{{[], []} => {q, {1, 1}, {0, 1}}}}, element(2, lists:keyfind(1, 1, PerRhs))).
