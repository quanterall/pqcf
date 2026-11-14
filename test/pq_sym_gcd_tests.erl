%% The modular multivariate gcd (pq_sym_gcd, decision record SY2): known gcds by construction, agreement with the
%% exact Erlang gcd of pq_simplify (pseudo-remainder sequences), and exact division.
-module(pq_sym_gcd_tests).
-include_lib("eunit/include/eunit.hrl").

%% x = [1,0,0], y = [0,1,0], z = [0,0,1]
p(Terms) -> maps:from_list([{E, C} || {E, C} <- Terms, C =/= 0]).

mul(A, B) ->
    lists:foldl(fun({K, V}, Acc) ->
                        case maps:get(K, Acc, 0) + V of 0 -> maps:remove(K, Acc); S -> Acc#{K => S} end
                end, #{},
                [{[X + Y || {X, Y} <- lists:zip(Ka, Kb)], Va * Vb} || Ka := Va <- A, Kb := Vb <- B]).

prim(P) -> element(2, pq_sym_gcd:primitive(P)).

small_test() ->
    % (x - 1)(x + 2), (x - 1)(x + 3): x - 1
    ?assertEqual(p([{[1], 1}, {[0], -1}]), pq_sym_gcd:gcd(p([{[2], 1}, {[1], 1}, {[0], -2}]), p([{[2], 1}, {[1], 2}, {[0], -3}]))),
    % 6 x^2 y and 4 x y^3: x y (contents are units over Q)
    ?assertEqual(p([{[1, 1], 1}]), pq_sym_gcd:gcd(p([{[2, 1], 6}]), p([{[1, 3], 4}]))),
    % coprime, and with zero
    ?assertEqual(p([{[0, 0], 1}]), pq_sym_gcd:gcd(p([{[1, 0], 1}, {[0, 1], 1}]), p([{[1, 0], 1}, {[0, 1], -1}]))),
    ?assertEqual(p([{[1, 0], 1}, {[0, 0], 2}]), pq_sym_gcd:gcd(#{}, p([{[1, 0], -3}, {[0, 0], -6}]))),
    % rational coefficients: x/2 + 1/3 and x^2 - 4/9 share 3x + 2
    ?assertEqual(p([{[1], 3}, {[0], 2}]), pq_sym_gcd:gcd(p([{[1], {1, 2}}, {[0], {1, 3}}]), p([{[2], 1}, {[0], {-4, 9}}]))).

%% a random polynomial in N variables: Terms terms, exponents up to D, coefficients in [-C, C]
rnd(N, Terms, D, C) ->
    p([{[rand:uniform(D + 1) - 1 || _ <- lists:seq(1, N)], rand:uniform(2 * C + 1) - C - 1} || _ <- lists:seq(1, Terms)]).

known_gcds_test_() ->
    {timeout, 120, fun() ->
        rand:seed(exsss, {7, 11, 13}),
        lists:foreach(
          fun(N) ->
                  lists:foreach(fun(_) ->
                                        G = rnd(N, 4, 3, 50), A = rnd(N, 5, 3, 50), B = rnd(N, 5, 3, 50),
                                        case map_size(G) > 0 andalso map_size(A) > 0 andalso map_size(B) > 0 of
                                            false -> ok;
                                            true ->
                                                Got = pq_sym_gcd:gcd(mul(A, G), mul(B, G)),
                                                % gcd(a g, b g) = g gcd(a, b)
                                                ?assertEqual(prim(mul(G, pq_sym_gcd:gcd(A, B))), Got),
                                                {ok, _} = pq_sym_gcd:divide(mul(A, G), Got),
                                                {ok, _} = pq_sym_gcd:divide(mul(B, G), Got)
                                        end
                                end, lists:seq(1, 15))
          end, [1, 2, 3, 4, 5])
    end}.

agrees_with_prs_test_() ->
    {timeout, 120, fun() ->
        rand:seed(exsss, {3, 5, 8}),
        Q = fun(P) -> maps:map(fun(_, C) -> {C, 1} end, P) end,
        lists:foreach(fun(_) ->
                              G = rnd(3, 3, 2, 9), A = mul(rnd(3, 3, 2, 9), G), B = mul(rnd(3, 3, 2, 9), G),
                              case map_size(A) > 0 andalso map_size(B) > 0 of
                                  false -> ok;
                                  true -> ?assertEqual(Q(pq_sym_gcd:gcd(A, B)), pq_simplify:gcd_prs(Q(A), Q(B)))
                              end
                      end, lists:seq(1, 25))
    end}.

repeated_factor_test() ->
    % P = (x + y)^3 (x - 2 z): gcd(P, dP/dx) = (x + y)^2
    S = p([{[1, 0, 0], 1}, {[0, 1, 0], 1}]),
    P = mul(mul(S, mul(S, S)), p([{[1, 0, 0], 1}, {[0, 0, 1], -2}])),
    Dp = maps:fold(fun([0 | _], _, Acc) -> Acc;
                      ([K | R], C, Acc) -> Acc#{[K - 1 | R] => maps:get([K - 1 | R], Acc, 0) + K * C} end,
                   #{}, P),
    ?assertEqual(mul(S, S), pq_sym_gcd:gcd(P, maps:filter(fun(_, C) -> C =/= 0 end, Dp))).

divide_test() ->
    A = p([{[2, 0], 1}, {[0, 2], -1}]),  % x^2 - y^2
    ?assertEqual({ok, p([{[1, 0], 1}, {[0, 1], -1}])}, pq_sym_gcd:divide(A, p([{[1, 0], 1}, {[0, 1], 1}]))),
    ?assertEqual({ok, p([{[1, 0], {1, 2}}, {[0, 1], {-1, 2}}])}, pq_sym_gcd:divide(A, p([{[1, 0], 2}, {[0, 1], 2}]))),
    ?assertEqual(false, pq_sym_gcd:divide(A, p([{[1, 0], 1}, {[0, 1], 2}]))).

large_test_() ->
    % a gcd of 3 variables, degree 10 each: (1 + x + y + z)^10 times coprime cofactors
    {timeout, 120, fun() ->
        L = p([{[0, 0, 0], 1}, {[1, 0, 0], 1}, {[0, 1, 0], 1}, {[0, 0, 1], 1}]),
        G = lists:foldl(fun(_, Acc) -> mul(Acc, L) end, #{[0, 0, 0] => 1}, lists:seq(1, 10)),
        A = mul(G, p([{[3, 0, 0], 1}, {[0, 2, 1], -5}, {[0, 0, 0], 7}])),
        B = mul(G, p([{[0, 3, 0], 2}, {[1, 0, 2], 3}, {[0, 0, 0], -1}])),
        ?assertEqual(G, pq_sym_gcd:gcd(A, B))
    end}.

api_test() ->
    X = {sym, x}, Y = {sym, y},
    % (x + y)(x - y) and (x + y)^2: x + y, canonical
    {ok, G} = pq_symbolic:gcd({mul, [{add, [X, Y]}, {add, [X, {mul, [-1, Y]}]}]}, {pow, {add, [X, Y]}, 2}),
    ?assertEqual(pq_symbolic:canonicalise({add, [X, Y]}), {ok, G}),
    % a plain sum of monomials, read directly: 6 x^2 y - 4 x y and 3 x y^2: x y
    {ok, G2} = pq_symbolic:gcd({add, [{mul, [6, X, X, Y]}, {mul, [-4, X, Y]}]}, {mul, [3, X, {pow, Y, 2}]}),
    ?assertEqual(pq_symbolic:canonicalise({mul, [X, Y]}), {ok, G2}),
    ?assertEqual({ok, 1}, pq_symbolic:gcd({add, [X, 1]}, {add, [X, 2]})),
    ?assertEqual({ok, 0}, pq_symbolic:gcd(0, 0)),
    ?assertEqual({error, {symbolic_unsupported, gcd}}, pq_symbolic:gcd({mul, [{complex, 0, 1}, X]}, X)),
    ?assertEqual({error, {symbolic_unsupported, gcd}}, pq_symbolic:gcd({exp, X}, X)).

squarefree_test() ->
    % 12 (x + y)^3 (x - 2 z)^2 (y z + 1) x^2: content 12, factors (y z + 1), x (x - 2 z), (x + y)
    S = p([{[1, 0, 0], 1}, {[0, 1, 0], 1}]),
    T = p([{[1, 0, 0], 1}, {[0, 0, 1], -2}]),
    U = p([{[0, 1, 1], 1}, {[0, 0, 0], 1}]),
    X = p([{[1, 0, 0], 1}]),
    P = lists:foldl(fun mul/2, #{[0, 0, 0] => 12}, [S, S, S, T, T, U, X, X]),
    ?assertEqual({{12, 1}, [{U, 1}, {mul(X, T), 2}, {S, 3}]}, pq_sym_gcd:squarefree(P)),
    % a square-free polynomial, and a constant
    ?assertEqual({{1, 1}, [{mul(S, T), 1}]}, pq_sym_gcd:squarefree(mul(S, T))),
    ?assertEqual({{-3, 2}, []}, pq_sym_gcd:squarefree(#{[0, 0, 0] => {-3, 2}})),
    % the content's own powers: (y + 1)^2 (x + y)
    Y1 = p([{[0, 1, 0], 1}, {[0, 0, 0], 1}]),
    ?assertEqual({{1, 1}, [{S, 1}, {Y1, 2}]}, pq_sym_gcd:squarefree(mul(mul(Y1, Y1), S))).

squarefree_random_test_() ->
    {timeout, 120, fun() ->
        rand:seed(exsss, {5, 5, 5}),
        lists:foreach(fun(_) ->
                              A = rnd(3, 3, 2, 9), B = rnd(3, 3, 2, 9), C = rnd(3, 2, 2, 9),
                              case lists:all(fun(M) -> map_size(M) > 0 end, [A, B, C]) of
                                  false -> ok;
                                  true ->
                                      P = lists:foldl(fun mul/2, A, [B, B, C, C, C]),
                                      {K, Fs} = pq_sym_gcd:squarefree(P),
                                      % the product back is P, and each factor is square-free
                                      Back = lists:foldl(fun({F, M}, Acc) -> lists:foldl(fun(_, X) -> mul(X, F) end, Acc, lists:seq(1, M)) end,
                                                         #{[0, 0, 0] => 1}, Fs),
                                      ?assertEqual(prim(P), prim(Back)),
                                      ?assertEqual(element(1, pq_sym_gcd:primitive(P)), K),
                                      [?assertEqual({{1, 1}, [{F, 1}]}, pq_sym_gcd:squarefree(F)) || {F, _} <- Fs]
                              end
                      end, lists:seq(1, 20))
    end}.
