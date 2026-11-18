%% pq_diagrams (plan "theoretical physics instruments", F4): Feynman graphs and their weights 1/S against
%% zero-dimensional field theory, where every propagator and vertex is 1 and the sums of 1/S over graphs are the
%% coefficients of generating functions (Cvitanovic, Lautrup and Pearson, Phys. Rev. D 18, 1939 (1978),
%% doi:10.1103/PhysRevD.18.1939): ln Z for connected vacuum graphs, the cumulants of Z(J) for connected n-point graphs,
%% the self-energy Sigma = 1 - 1/G and the vertex G_c^(4) / G^4 for one-particle-irreducible ones. The generating
%% functions are computed here in exact rationals, independently of the graphs. QED and QCD diagrams against FORM 5's
%% generator (test/oracle/gt_diagrams.sh, test/fixtures/gt_diagrams.json): the same number of diagrams and the same
%% multiset of weights sign/S; the orderly method against the Wick pairings, graph by graph.
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

%%% Orderly generation and fermions

%% phi^3 vacuum graphs at five loops (eight vertices), out of reach of the pairings, against ln Z
phi3_five_loops_test_() ->
    {timeout, 120, fun() ->
        Ln = slog(zk(3, 0, 8)),
        Gs = pq_diagrams:graphs(phi3(), [], 5, []),
        ?assertEqual(71, length(Gs)),
        ?assertEqual(lists:nth(9, Ln), pq_diagrams:sum_weights(Gs))
    end}.

%% the orderly method gives the same graphs and weights as Wick's pairings
orderly_equals_wick_test_() ->
    {timeout, 300, fun() ->
        T2 = #{vertices => [{v, [phi, phi, chi]}], propagators => [{phi, phi}, {chi, chi}]},
        Cases = [{phi4(), [], 5, []}, {phi4(), [phi, phi], 3, [one_pi]}, {phi4(), [phi, phi, phi, phi], 2, []},
                 {phi3(), [], 4, []}, {phi3(), [phi, phi], 3, [one_pi, no_tadpoles]}, {T2, [], 4, []}],
        [?assertEqual(lists:sort(pq_diagrams:graphs(T, E, L, O, #{method => wick})), lists:sort(pq_diagrams:graphs(T, E, L, O)))
         || {T, E, L, O} <- Cases]
    end}.

qed() -> #{vertices => [{v, ['ELE', ele, pho]}], propagators => [{ele, 'ELE'}, {pho, pho}], fermions => [ele, 'ELE']}.

%% zero-dimensional QED with N flavours: Z = <(1 - e A)^N>, so ln Z = sum over vacuum graphs of sign N^(fermion loops)/S
qed_zero_dimensional_test_() ->
    {timeout, 120, fun() ->
        M = 8,
        [begin
             Z = [mul({binom(N, K) * pow(-1, K), 1}, {moment(K), 1}) || K <- lists:seq(0, M)],
             Ln = slog(Z),
             [begin
                  Gs = pq_diagrams:graphs(qed(), [], L, []),
                  Sum = lists:foldl(fun({G, {A, B}}, Acc) -> add(Acc, {A * pow(N, pq_diagrams:fermion_loops(G, qed())), B}) end,
                                    {0, 1}, Gs),
                  ?assertEqual(lists:nth(2 * L - 1, Ln), Sum)
              end || L <- lists:seq(2, 5)]
         end || N <- [1, 2, 3]]
    end}.

binom(N, K) when K > N -> 0;
binom(N, K) -> fact(N) div (fact(K) * fact(N - K)).

%% QED and QCD against FORM 5: the same number of diagrams and the same weights
form_test_() ->
    {timeout, 300, fun() ->
        [begin
             Ext = [conj(binary_to_atom(F)) || F <- In] ++ [binary_to_atom(F) || F <- Out],
             Opts = case O of <<"1pi">> -> [one_pi]; _ -> [] end,
             Gs = pq_diagrams:graphs(model(M), Ext, L, Opts),
             ?assertEqual({Name, L, N}, {Name, L, length(Gs)}),
             ?assertEqual({Name, L, sort_q([frac(C) || C <- Cs])}, {Name, L, sort_q([W || {_, W} <- Gs])})
         end || #{<<"name">> := Name, <<"model">> := M, <<"in">> := In, <<"out">> := Out, <<"loops">> := L,
                  <<"options">> := O, <<"count">> := N, <<"coefficients">> := Cs} <- form_fixture()]
    end}.

model(<<"QED">>) -> qed();
model(<<"QCD">>) ->
    #{vertices => [{qqg, ['QUA', qua, glu]}, {ccg, ['GHO', gho, glu]}, {g3, [glu, glu, glu]}, {g4, [glu, glu, glu, glu]}],
      propagators => [{qua, 'QUA'}, {gho, 'GHO'}, {glu, glu}], fermions => [qua, 'QUA', gho, 'GHO']}.

%% FORM's incoming leg of a charged field carries the conjugate on its external node
conj(ele) -> 'ELE';
conj(qua) -> 'QUA';
conj(gho) -> 'GHO';
conj(X) -> X.

frac(B) ->
    case binary:split(B, <<"/">>) of
        [N] -> {binary_to_integer(N), 1};
        [N, D] -> {binary_to_integer(N), binary_to_integer(D)}
    end.

sort_q(L) -> lists:sort(fun({A, B}, {C, D}) -> A * D =< C * B end, L).

form_fixture() ->
    Dir = filename:dirname(code:which(?MODULE)),
    Paths = [filename:join([Dir, "..", "test", "fixtures", "gt_diagrams.json"]),
             filename:join([Dir, "..", "..", "..", "..", "test", "fixtures", "gt_diagrams.json"]),
             "test/fixtures/gt_diagrams.json"],
    [P | _] = [P || P <- Paths, filelib:is_regular(P)],
    {ok, B} = file:read_file(P),
    json:decode(B).
