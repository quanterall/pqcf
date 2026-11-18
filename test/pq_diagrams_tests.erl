%% pq_diagrams (plan "theoretical physics instruments", F4): Feynman graphs and their weights 1/S against
%% zero-dimensional field theory, where every propagator and vertex is 1 and the sums of 1/S over graphs are the
%% coefficients of generating functions (Cvitanovic, Lautrup and Pearson, Phys. Rev. D 18, 1939 (1978),
%% doi:10.1103/PhysRevD.18.1939): ln Z for connected vacuum graphs, the cumulants of Z(J) for connected n-point graphs,
%% the self-energy Sigma = 1 - 1/G and the vertex G_c^(4) / G^4 for one-particle-irreducible ones. The generating
%% functions are computed here in exact rationals, independently of the graphs.
-module(pq_diagrams_tests).
-include_lib("eunit/include/eunit.hrl").

phi4() -> #{vertices => [{v4, [phi, phi, phi, phi]}], propagators => [{phi, phi}]}.
phi3() -> #{vertices => [{v3, [phi, phi, phi]}], propagators => [{phi, phi}]}.

%%% Zero-dimensional generating functions, truncated power series in the coupling (lists of exact fractions)

%% <x^(2m)> for the Gaussian: (2m - 1)!!
moment(M) when M rem 2 =:= 1 -> 0;
moment(M) -> dfact(M - 1).
dfact(N) when N =< 1 -> 1;
dfact(N) -> N * dfact(N - 2).
fact(0) -> 1;
fact(N) -> N * fact(N - 1).

%% Z_k = sum_n g^n / (n! (K!)^n) <x^(K n + k)>, K the vertex degree, to order N
zk(K, Kext, N) -> [q(moment(K * Nv + Kext), fact(Nv) * pow(fact(K), Nv)) || Nv <- lists:seq(0, N)].

pow(_, 0) -> 1;
pow(X, K) -> X * pow(X, K - 1).

q(A, B) -> norm(A, B).
norm(0, _) -> {0, 1};
norm(A, B) -> G = gcd(abs(A), abs(B)), S = case B < 0 of true -> -1; false -> 1 end, {S * A div G, S * B div G}.
gcd(A, 0) -> A;
gcd(A, B) -> gcd(B, A rem B).
add({A, B}, {C, D}) -> norm(A * D + C * B, B * D).
sub(X, {C, D}) -> add(X, {-C, D}).
mul({A, B}, {C, D}) -> norm(A * C, B * D).
dvd({A, B}, {C, D}) -> norm(A * D, B * C).

smul(A, B) -> [lists:foldl(fun(I, S) -> add(S, mul(lists:nth(I + 1, A), lists:nth(K - I + 1, B))) end, {0, 1}, lists:seq(0, K))
               || K <- lists:seq(0, length(A) - 1)].
sdiv(A, B) ->
    % A / B for B(0) =/= 0, term by term
    lists:foldl(fun(K, Acc) ->
                        S = lists:foldl(fun(I, S0) -> add(S0, mul(lists:nth(I + 1, B), lists:nth(K - I + 1, Acc))) end,
                                        {0, 1}, lists:seq(1, K)),
                        Acc ++ [dvd(sub(lists:nth(K + 1, A), S), hd(B))]
                end, [], lists:seq(0, length(A) - 1)).
slog(Z) ->
    % ln Z with Z(0) = 1: L' = Z'/Z
    N = length(Z) - 1,
    Dz = [mul({K, 1}, lists:nth(K + 1, Z)) || K <- lists:seq(1, N)] ++ [{0, 1}],
    R = sdiv(Dz, Z),
    [{0, 1} | [dvd(lists:nth(K, R), {K, 1}) || K <- lists:seq(1, N)]].
ssub(A, B) -> [sub(X, Y) || {X, Y} <- lists:zip(A, B)].
one(N) -> [{1, 1} | lists:duplicate(N, {0, 1})].

%%% The checks

%% the sum of 1/S of the graphs with Nv vertices: phi^4 with E legs has Nv = L - 1 + E/2 vertices
sum(T, Ext, L, Opts) -> pq_diagrams:sum_weights(pq_diagrams:graphs(T, Ext, L, Opts)).

phi4_vacuum_test() ->
    Ln = slog(zk(4, 0, 5)),
    % L loops <-> L - 1 vertices
    [?assertEqual(lists:nth(L, Ln), sum(phi4(), [], L, [])) || L <- lists:seq(2, 6)],
    ?assertEqual([1, 2, 4, 10, 28], [length(pq_diagrams:graphs(phi4(), [], L, [])) || L <- lists:seq(2, 6)]).

phi4_two_point_test() ->
    N = 4,
    Z0 = zk(4, 0, N), Z2 = zk(4, 2, N),
    G = sdiv(Z2, Z0),                               % connected: Z1 = 0 by symmetry
    Sigma = ssub(one(N), sdiv(one(N), G)),          % G = 1 / (1 - Sigma)
    % L loops <-> L vertices
    [?assertEqual(lists:nth(L + 1, G), sum(phi4(), [phi, phi], L, [])) || L <- lists:seq(0, N)],
    [?assertEqual(lists:nth(L + 1, Sigma), sum(phi4(), [phi, phi], L, [one_pi])) || L <- lists:seq(1, N)].

phi4_four_point_test() ->
    N = 3,
    Z0 = zk(4, 0, N), Z2 = zk(4, 2, N), Z4 = zk(4, 4, N),
    G = sdiv(Z2, Z0),
    Gc4 = ssub(sdiv(Z4, Z0), [mul({3, 1}, X) || X <- smul(G, G)]),
    G4 = smul(smul(G, G), smul(G, G)),
    Gamma4 = sdiv(Gc4, G4),
    % L loops <-> L + 1 vertices
    [?assertEqual(lists:nth(L + 2, Gc4), sum(phi4(), [phi, phi, phi, phi], L, [])) || L <- lists:seq(0, N - 1)],
    [?assertEqual(lists:nth(L + 2, Gamma4), sum(phi4(), [phi, phi, phi, phi], L, [one_pi])) || L <- lists:seq(0, N - 1)].

phi3_vacuum_test_() -> {timeout, 120, fun phi3_vacuum/0}.
phi3_vacuum() ->
    % phi^3: 2 n vertices at n + 1 loops
    Ln = slog(zk(3, 0, 8)),
    [?assertEqual(lists:nth(2 * (L - 1) + 1, Ln), sum(phi3(), [], L, [])) || L <- lists:seq(2, 4)],
    ?assertEqual([2, 5, 17], [length(pq_diagrams:graphs(phi3(), [], L, [])) || L <- lists:seq(2, 4)]).

%% two fields: phi^2 chi (g/2!), both propagators 1; vacuum graphs equal ln of
%% Z = sum_n g^n / (n! 2^n) <x^(2n)> <y^n>, x and y independent Gaussians
two_fields_test() ->
    T = #{vertices => [{v, [phi, phi, chi]}], propagators => [{phi, phi}, {chi, chi}]},
    N = 6,
    Z = [q(moment(2 * Nv) * moment(Nv), fact(Nv) * pow(2, Nv)) || Nv <- lists:seq(0, N)],
    Ln = slog(Z),
    % 2 n vertices at n + 1 loops
    [?assertEqual(lists:nth(2 * (L - 1) + 1, Ln), sum(T, [], L, [])) || L <- lists:seq(2, 4)].

%% the classic small cases: phi^4 self-energy at two loops is the sunset (1/6) and the double tadpole (1/4)
named_graphs_test() ->
    Gs = pq_diagrams:graphs(phi4(), [phi, phi], 2, [one_pi]),
    ?assertEqual([{1, 6}, {1, 4}], lists:sort(fun({A, B}, {C, D}) -> A * D =< C * B end, [W || {_, W} <- Gs])),
    % the second has a snail (a line from a vertex to itself), not a tadpole: no_tadpoles keeps it, no_snails does not
    ?assertEqual(2, length(pq_diagrams:graphs(phi4(), [phi, phi], 2, [one_pi, no_tadpoles]))),
    ?assertEqual([{1, 6}], [W || {_, W} <- pq_diagrams:graphs(phi4(), [phi, phi], 2, [one_pi, no_snails])]),
    % a tadpole: the one-loop phi^3 one-point graph hangs on one line; in the self-energy at two loops no_tadpoles
    % removes the graphs with a tadpole on an internal line
    ?assertEqual([], pq_diagrams:graphs(phi3(), [phi, phi], 2, [no_tadpoles, one_pi]) --
                     pq_diagrams:graphs(phi3(), [phi, phi], 2, [one_pi])).
