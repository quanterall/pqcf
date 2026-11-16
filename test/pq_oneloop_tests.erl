%% pq_oneloop (plan "theoretical physics instruments", F3): the two-point functions against two independent references
%% (test/oracle/gt_oneloop.sh, test/fixtures/gt_oneloop.json): mpmath's Feynman-parameter integrals at 40 digits
%% (B0, B1, B11 to 35 digits) and LoopTools 2.16 (B0, B1, B00, B11, DB0 to 1e-11 relative), in 60 configurations over
%% every kinematic region.
-module(pq_oneloop_tests).
-include_lib("eunit/include/eunit.hrl").

fixture() ->
    Dir = filename:dirname(code:which(?MODULE)),
    Paths = [filename:join([Dir, "..", "test", "fixtures", "gt_oneloop.json"]),
             filename:join([Dir, "..", "..", "..", "..", "test", "fixtures", "gt_oneloop.json"]), "test/fixtures/gt_oneloop.json"],
    [P | _] = [P || P <- Paths, filelib:is_regular(P)],
    {ok, B} = file:read_file(P),
    json:decode(B).

q(Bin) ->
    case binary:split(Bin, <<"/">>) of
        [N] -> binary_to_integer(N);
        [N, D] -> {binary_to_integer(N), binary_to_integer(D)}
    end.

%% |x - y| =< tol max(1, |y|), x a pq_mp number, y a decimal string or float
close(X, Y, Digits) ->
    P = 200,
    Yv = case Y of B when is_binary(B) -> pq_mp:from_string(binary_to_list(B), P); F -> pq_mp:from(float(F), P) end,
    Diff = abs(pq_mp:to_float(pq_mp:sub(X, Yv, P))),
    Diff =< math:pow(10, -Digits) * max(1.0, abs(pq_mp:to_float(Yv))).

%% the cases where the function is not finite or not defined here, and the error PQCF gives
expected_error(F, S, A, B) ->
    Zero = fun(V) -> V =:= 0 end,
    if
        F =/= b0, F =/= b00, S =:= 0 -> {pq_oneloop, zero_momentum};
        F =:= db0 -> db0_error(S, A, B, Zero);
        true -> none
    end.

db0_error(S, A, B, Zero) ->
    {Sv, Av, Bv} = {val(S), val(A), val(B)},
    Lam = Sv * Sv + Av * Av + Bv * Bv - 2 * (Sv * Av + Sv * Bv + Av * Bv),
    case {Zero(A), Zero(B)} of
        {true, true} -> none;
        {true, false} when Bv == Sv -> {pq_oneloop, infrared_divergent};
        {false, true} when Av == Sv -> {pq_oneloop, infrared_divergent};
        {false, false} when Lam == 0, Av + Bv - Sv =< 0 -> {pq_oneloop, threshold};
        _ -> none
    end.

val({N, D}) -> N / D;
val(N) -> float(N).

%% 60 configurations: mpmath to 35 digits, LoopTools to 1e-11
two_point_test_() ->
    {timeout, 300, fun() ->
        lists:foreach(
          fun(#{<<"s">> := Sb, <<"m1sq">> := Ab, <<"m2sq">> := Bb, <<"mp">> := Mp, <<"looptools">> := Lt}) ->
                  {S, A, B} = {q(Sb), q(Ab), q(Bb)},
                  [case expected_error(F, S, A, B) of
                       none ->
                           {R, I} = pq_oneloop:F(S, A, B, 160),
                           Name = atom_to_binary(F),
                           [Lr, Li] = maps:get(Name, Lt),
                           ?assert(close(R, Lr, 11), {F, Sb, Ab, Bb, looptools, re}),
                           ?assert(close(I, Li, 11), {F, Sb, Ab, Bb, looptools, im}),
                           case maps:find(Name, Mp) of
                               {ok, [Mr, Mi]} ->
                                   ?assert(close(R, Mr, 35), {F, Sb, Ab, Bb, mpmath, re}),
                                   ?assert(close(I, Mi, 35), {F, Sb, Ab, Bb, mpmath, im});
                               error -> ok
                           end;
                       Err ->
                           ?assertError(Err, pq_oneloop:F(S, A, B, 160))
                   end || F <- [b0, b1, b00, b11, db0]]
          end, fixture())
    end}.

%% A0, the scale mu^2, and closed forms
closed_forms_test() ->
    P = 200,
    V = fun(X) -> pq_mp:to_float(X) end,
    {A0, Z} = pq_oneloop:a0(2, P),
    ?assert(abs(V(A0) - 2 * (1 - math:log(2))) < 1.0e-15),
    ?assertEqual(0.0, V(Z)),
    % B0(m^2; m^2, m^2) = 2 - pi/sqrt(3)
    {B, _} = pq_oneloop:b0(1, 1, 1, P),
    ?assert(abs(V(B) - (2 - math:pi() / math:sqrt(3))) < 1.0e-15),
    % mu^2: B0 depends on log mu^2 as + log mu^2; A0 as + m^2 log mu^2
    {B2, _} = pq_oneloop:b0(5, 1, 2, P, #{mu2 => 3}),
    {B1, _} = pq_oneloop:b0(5, 1, 2, P),
    ?assert(abs(V(B2) - V(B1) - math:log(3)) < 1.0e-15),
    {A3, _} = pq_oneloop:a0(2, P, #{mu2 => 3}),
    ?assert(abs(V(A3) - V(A0) - 2 * math:log(3)) < 1.0e-15),
    % the massless bubble is scaleless at p^2 = 0
    ?assertEqual({0.0, 0.0}, begin {R0, I0} = pq_oneloop:b0(0, 0, 0, P), {V(R0), V(I0)} end),
    % to many digits: B0(1; 1, 1) to 60 digits
    {Bs, _} = pq_oneloop:value(b0, [1, 1, 1], 60),
    Want = pq_mp:sub(pq_mp:from(2, 300), pq_mp:divide(pq_mp:pi(300), pq_mp:sqrt(pq_mp:from(3, 300), 300), 300), 300),
    ?assertEqual(pq_mp:to_string(Want, 60), Bs).

%% atan, added to pq_mp for B0 between the thresholds
atan_test() ->
    P = 400,
    ?assertEqual(pq_mp:to_string(pq_mp:ldexp(pq_mp:pi(P), -2), 100), pq_mp:to_string(pq_mp:atan(1, P), 100)),
    % atan(sqrt 3) = pi/3, atan(-x) = -atan x, large arguments
    ?assertEqual(pq_mp:to_string(pq_mp:divide(pq_mp:pi(P), pq_mp:from(3, P), P), 100),
                 pq_mp:to_string(pq_mp:atan(pq_mp:sqrt(pq_mp:from(3, P), P), P), 100)),
    ?assertEqual(pq_mp:to_string(pq_mp:neg(pq_mp:atan({1, 7}, P)), 100), pq_mp:to_string(pq_mp:atan({-1, 7}, P), 100)),
    ?assert(abs(pq_mp:to_float(pq_mp:atan(1.0e6, P)) - math:atan(1.0e6)) < 1.0e-15).
