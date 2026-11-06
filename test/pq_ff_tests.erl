%% Prime fields for reconstruction (pq_ff, native/pq_nif/src/ffslp.rs; decision record R1).
-module(pq_ff_tests).
-include_lib("eunit/include/eunit.hrl").

primes_test() ->
    ?assertEqual(9223372036854775783, pq_ff:prime(1)),          % the largest prime below 2^63
    Ps = [pq_ff:prime(K) || K <- lists:seq(1, 6)],
    ?assertEqual(Ps, lists:reverse(lists:usort(Ps))),
    ?assert(lists:all(fun pq_ff:is_prime/1, Ps)),
    ?assert(pq_ff:is_prime((1 bsl 61) - 1)),
    [?assertNot(pq_ff:is_prime(C)) || C <- [1, 561, 1105, 3215031751, 3825123056546413051]],   % Carmichael, strong pseudoprimes
    ?assert(pq_ff:is_prime(2)).

%% fractions back from their residues: one prime when they are small, two when they are not
ratrecon_test() ->
    rand:seed(exsss, {1, 2, 3}),
    P1 = pq_ff:prime(1), P2 = pq_ff:prime(2),
    [begin
         %% one 63-bit prime determines |N|, D below about 2^31; two primes, below about 2^62
         N = rand:uniform(1 bsl 30) - (1 bsl 29), D = rand:uniform(1 bsl 30),
         {N1, D1} = pq_exact_frac(N, D),
         Res = fun(P) -> N1 * pq_ff:inv(D1, P) rem P end,
         ?assertEqual({N1, D1}, pq_ff:ratrecon(Res(P1), P1)),
         Big = {N1 * (1 bsl 25) + 1, D1 * 7},
         {Nb, Db} = pq_exact_frac(element(1, Big), element(2, Big)),
         RB = fun(P) -> pq_mod(Nb * pq_ff:inv(Db, P), P) end,
         {A12, M12} = pq_ff:crt({RB(P1), P1}, {RB(P2), P2}),
         ?assertEqual({Nb, Db}, pq_ff:ratrecon(A12, M12))
     end || _ <- lists:seq(1, 50)],
    %% whatever comes back is consistent with the residue and within the bound
    [case pq_ff:ratrecon(A, 1000003) of
         false -> ok;
         {N, D} -> ?assertEqual(0, pq_mod(N - A * D, 1000003)), ?assert(abs(N) * abs(N) =< 1000003 div 2), ?assert(D * D =< 1000003 div 2)
     end || A <- [123456789, 1, 999999, 500001]].

pq_exact_frac(N, D) -> G = pq_sym_exact:gcd(abs(N), D), {N div G, D div G}.
pq_mod(A, M) -> ((A rem M) + M) rem M.

%% an independent oracle: the expression tree itself evaluated mod P (throws bad on a division by zero)
oracle(N, _B, P) when is_integer(N) -> pq_mod(N, P);
oracle({rat, N, D}, _B, P) -> pq_mod(N * inverse(D, P), P);
oracle({sym, X}, B, P) -> maps:get(X, B) rem P;
oracle({add, L}, B, P) -> lists:foldl(fun(E, A) -> (A + oracle(E, B, P)) rem P end, 0, L);
oracle({mul, L}, B, P) -> lists:foldl(fun(E, A) -> A * oracle(E, B, P) rem P end, 1, L);
oracle({pow, E, K}, B, P) when K >= 0 -> powm(oracle(E, B, P), K, P);
oracle({pow, E, K}, B, P) -> powm(inverse(oracle(E, B, P), P), -K, P).

inverse(A, P) -> case pq_mod(A, P) of 0 -> throw(bad); X -> pq_ff:inv(X, P) end.
powm(_, 0, _) -> 1;
powm(B, K, P) -> B * powm(B, K - 1, P) rem P.

syms() -> [x, y, z, w].

leaf() ->
    case rand:uniform(6) of
        1 -> rand:uniform(9) - 5;
        2 -> {rat, rand:uniform(9) - 5, rand:uniform(7)};
        _ -> {sym, lists:nth(rand:uniform(4), syms())}
    end.

%% random rational functions: sums and products of leaves, powers (negative ones of sums too)
expr(0) -> leaf();
expr(D) ->
    case rand:uniform(5) of
        1 -> {add, [expr(D - 1) || _ <- lists:seq(1, 2 + rand:uniform(3))]};
        2 -> {mul, [expr(D - 1) || _ <- lists:seq(1, 2 + rand:uniform(2))]};
        3 -> {pow, {add, [expr(D - 1), expr(D - 1), 1]}, rand:uniform(3)};
        4 -> {pow, {add, [expr(D - 1), {sym, lists:nth(rand:uniform(4), syms())}]}, -rand:uniform(2)};
        _ -> leaf()
    end.

eval_test_() ->
    {timeout, 120, fun() ->
                           rand:seed(exsss, {9, 9, 9}),
                           Primes = [pq_ff:prime(1), 1000003, 101],
                           [begin
                                E = expr(1 + K rem 3),
                                Syms = syms(),
                                [begin
                                     {ok, Prog} = pq_symbolic:compile([E], Syms, #{optimize => O, exact_constants => true}),
                                     {ok, FF} = pq_ff:prepare(Prog),
                                     [begin
                                          Points = [[rand:uniform(P) - 1 || _ <- Syms] || _ <- lists:seq(1, 5)],
                                          {ok, Vals} = case pq_ff:eval(FF, P, Points) of
                                                           bad_prime -> {ok, skip};
                                                           R -> R
                                                       end,
                                          Vals =:= skip orelse
                                              [begin
                                                   B = maps:from_list(lists:zip(Syms, X)),
                                                   Want = try [oracle(E, B, P)] catch throw:bad -> bad end,
                                                   %% the canonical form may cancel a factor the tree divides by:
                                                   %% only a point good for both is compared
                                                   ?assert(Want =:= bad orelse V =:= bad orelse V =:= Want)
                                               end || {X, V} <- lists:zip(Points, Vals)]
                                      end || P <- Primes]
                                 end || O <- [none, horner_cse]]
                            end || K <- lists:seq(1, 60)]
                   end}.

%% a fraction's value: (x + 1/2)^2 / (y - 3) at a point, and the bad point y = 3
fraction_test() ->
    E = {mul, [{pow, {add, [{sym, x}, {rat, 1, 2}]}, 2}, {pow, {add, [{sym, y}, -3]}, -1}]},
    {ok, Prog} = pq_symbolic:compile([E], [x, y], #{exact_constants => true}),
    {ok, FF} = pq_ff:prepare(Prog),
    P = 1000003,
    {ok, [[V], bad]} = pq_ff:eval(FF, P, [[4, 5], [4, 3]]),
    ?assertEqual(pq_mod(81 * pq_ff:inv(4 * 2, P), P), V).

%% what has no meaning mod p is refused
refusals_test() ->
    [begin
         {ok, Prog} = pq_symbolic:compile([E], [x], #{exact_constants => true}),
         ?assertMatch({error, {symbolic_unsupported, _}}, pq_ff:prepare(Prog))
     end || E <- [{add, [{const, pi}, {sym, x}]}, {pow, {sym, x}, {rat, 1, 2}}, {exp, {sym, x}},
                  {mul, [{complex, 0, 1}, {sym, x}]}]],
    %% without exact constants, a float constant is refused too
    {ok, P2} = pq_symbolic:compile([{add, [{rat, 1, 3}, {sym, x}]}], [x]),
    ?assertMatch({error, _}, pq_ff:prepare(P2)).

%% exact constants change nothing in the float evaluation
exact_constants_test() ->
    rand:seed(exsss, {4, 4, 4}),
    [begin
         E = expr(2),
         B = maps:from_list([{X, 0.3 + rand:uniform()} || X <- syms()]),
         [begin
              {ok, A} = pq_symbolic:compile([E], syms(), #{optimize => O}),
              {ok, X} = pq_symbolic:compile([E], syms(), #{optimize => O, exact_constants => true}),
              ?assertEqual(pq_symbolic:eval_slp(A, B), pq_symbolic:eval_slp(X, B))
          end || O <- [none, horner_cse]]
     end || _ <- lists:seq(1, 30)].
