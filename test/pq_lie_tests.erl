%% pq_lie (plan "theoretical physics instruments", F1): Wigner symbols against SymPy (test/oracle/gt_wigner.py,
%% test/fixtures/gt_wigner.json) exactly, and the SU(N) colour invariants against their published closed forms
%% (van Ritbergen, Schellekens, Vermaseren, Int. J. Mod. Phys. A 14, 41 (1999), arXiv:hep-ph/9802376).
-module(pq_lie_tests).
-include_lib("eunit/include/eunit.hrl").

-define(L, pq_lie).

fixture() ->
    Dir = filename:dirname(code:which(?MODULE)),
    Paths = [filename:join([Dir, "..", "test", "fixtures", "gt_wigner.json"]),
             filename:join([Dir, "..", "..", "..", "..", "test", "fixtures", "gt_wigner.json"]), "test/fixtures/gt_wigner.json"],
    [P | _] = [P || P <- Paths, filelib:is_regular(P)],
    {ok, B} = file:read_file(P),
    json:decode(B).

q(B) ->
    case binary:split(B, <<"/">>) of
        [N, D] -> {binary_to_integer(N), binary_to_integer(D)};
        [N] -> {binary_to_integer(N), 1}
    end.

%% {S, Q} against SymPy's sign and square (zero: sign 0)
same({0, _}, 0, _) -> true;
same({S, Q}, S, Q) -> true;
same(_, _, _) -> false.

%% every 3j with j =< 3 (1,684 with the random ones up to j = 10), exactly
wigner_3j_test() ->
    #{<<"3j">> := Cases} = fixture(),
    [?assert(same(?L:wigner_3j(J1, J2, J3, M1, M2, M3), S, q(Sq)))
     || [J1, J2, J3, M1, M2, M3, S, Sq] <- Cases].

clebsch_gordan_test() ->
    #{<<"cg">> := Cases} = fixture(),
    [?assert(same(?L:clebsch_gordan(J1, M1, J2, M2, J, M), S, q(Sq))) || [J1, M1, J2, M2, J, M, S, Sq] <- Cases].

%% every 6j with j =< 3 and a sample up to j = 8 (3,418), exactly
wigner_6j_test() ->
    #{<<"6j">> := Cases} = fixture(),
    [?assert(same(?L:wigner_6j(J1, J2, J3, J4, J5, J6), S, q(Sq))) || [J1, J2, J3, J4, J5, J6, S, Sq] <- Cases].

%% 150 9j symbols (sums of square roots) to 40 digits
wigner_9j_test_() ->
    {timeout, 120, fun() ->
        #{<<"9j">> := Cases} = fixture(),
        [begin
             E = ?L:wigner_9j(J1, J2, J3, J4, J5, J6, J7, J8, J9),
             {ok, V} = pq_symbolic:evaluate_digits(E, #{}, 45),
             Want = pq_mp:from_string(binary_to_list(W), 160),
             D = pq_mp:sub(pq_mp:from_string(V, 160), Want, 160),
             ?assert(element(2, D) =:= 0 orelse abs(pq_mp:to_float(D)) < 1.0e-40)
         end || [J1, J2, J3, J4, J5, J6, J7, J8, J9, W] <- Cases]
    end}.

%% the SU(N) invariants with N symbolic, by Fierz reduction of traces, equal to their closed forms exactly
colour_test_() ->
    {timeout, 120, fun() ->
        N = {sym, n},
        P = fun(K) -> {pow, N, K} end,
        Want = #{tf => {rat, 1, 2},
                 cf => {mul, [{add, [P(2), -1]}, {pow, {mul, [2, N]}, -1}]},
                 ca => N,
                 na => {add, [P(2), -1]},
                 d33 => {mul, [{add, [P(2), -4]}, {add, [P(2), -1]}, {pow, N, -1}]},
                 d44_ff => {mul, [{add, [P(4), {mul, [-6, P(2)]}, 18]}, {pow, {mul, [96, P(2)]}, -1}]},
                 d44_fa => {mul, [N, {add, [P(2), 6]}, {rat, 1, 48}]},
                 d44_aa => {mul, [P(2), {add, [P(2), 36]}, {rat, 1, 24}]}},
        [begin
             {ok, Z} = pq_symbolic:cancel({add, [?L:invariant(I), {mul, [-1, W]}]}),
             ?assertEqual({I, 0}, {I, Z})
         end || I := W <- Want]
    end}.

%% a Fierz identity as a check of the reduction itself: T^a T^b T^a = -(1/(2N)) T^b, so
%% Tr(T^a T^b T^a T^b) = -(1/(2N)) Tr(T^b T^b) = -(N^2 - 1)/(4N)
fierz_test() ->
    {ok, Z} = pq_symbolic:cancel({add, [?L:colour(?L:contract([{1, [[a, b, a, b]]}])),
                                        {mul, [{add, [{pow, {sym, n}, 2}, -1]}, {pow, {mul, [4, {sym, n}]}, -1}]}]}),
    ?assertEqual(0, Z).
