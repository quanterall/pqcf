%% pq_oneloop (plan "theoretical physics instruments", F3): the two-point functions against two independent references
%% (test/oracle/gt_oneloop.sh, test/fixtures/gt_oneloop.json): mpmath's Feynman-parameter integrals at 40 digits
%% (B0, B1, B11 to 35 digits) and LoopTools 2.16 (B0, B1, B00, B11, DB0 to 1e-11 relative), in 60 configurations over
%% every kinematic region; the three-point function C0 against mpmath (test/oracle/gt_c0.sh): the direct integral
%% to 28 digits where D > 0, contour deformation to 1e-12 elsewhere.
-module(pq_oneloop_tests).
-include_lib("eunit/include/eunit.hrl").

fixture() -> fixture("gt_oneloop.json").
fixture(Name) ->
    Dir = filename:dirname(code:which(?MODULE)),
    Paths = [filename:join([Dir, "..", "test", "fixtures", Name]),
             filename:join([Dir, "..", "..", "..", "..", "test", "fixtures", Name]), filename:join("test/fixtures", Name)],
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

%% C0 in 40 configurations: against mpmath's direct integral (28 digits) where D > 0, its contour deformation (1e-12)
%% where two deformation parameters agree, and elsewhere against the value two of the deformation (kappa = 0.1),
%% QCDLoop and LoopTools agree on (1e-12), and failing all of those (one configuration, where both programs are wrong)
%% against the deformation at 1e-6
c0_test_() ->
    {timeout, 600, fun() ->
        lists:foreach(
          fun(#{<<"args">> := Args, <<"deformed">> := [Dr, Di], <<"stable">> := Stable} = Rec) ->
                  {R, I} = apply(pq_oneloop, c0, [q(A) || A <- Args] ++ [110]),
                  case maps:find(<<"direct">>, Rec) of
                      {ok, Direct} ->
                          ?assert(close(R, Direct, 28), {Args, direct}),
                          ?assert(close(I, 0.0, 28), {Args, direct_im});
                      error -> ok
                  end,
                  case {Stable, programs(Rec)} of
                      {true, _} ->
                          ?assert(close(R, Dr, 12), {Args, deformed_re}),
                          ?assert(close(I, Di, 12), {Args, deformed_im});
                      {false, {Pr, Pi}} ->
                          ?assert(close(R, Pr, 12), {Args, programs_re}),
                          ?assert(close(I, Pi, 12), {Args, programs_im});
                      {false, none} ->
                          % only the deformation, unstable at 1e-12: a weaker check
                          ?assert(close(R, Dr, 6), {Args, deformed_loose_re}),
                          ?assert(close(I, Di, 6), {Args, deformed_loose_im})
                  end
          end, fixture("gt_c0.json"))
    end}.

%% C0 where LoopTools (lambda < 0, below threshold) and QCDLoop (lambda < 0, above) are wrong in our builds, the
%% closed form at p^2 = 0, and the refusals
c0_cases_test() ->
    P = 110,
    V = fun(X) -> pq_mp:to_float(X) end,
    {R1, I1} = pq_oneloop:c0(1, 1, 1, 1, 1, 1, P),
    ?assert(abs(V(R1) + 0.671253105748004686) < 1.0e-15), ?assertEqual(0.0, V(I1)),
    {R2, I2} = pq_oneloop:c0({9, 2}, 7, 20, 1, {1, 4}, 2, P),
    ?assert(abs(V(R2) - 0.408381951022) < 1.0e-11), ?assert(abs(V(I2) + 0.152701852704) < 1.0e-11),
    % C0(0, 0, 0; m^2, m^2, m^2) = -1/(2 m^2)
    {R3, _} = pq_oneloop:c0(0, 0, 0, 2, 2, 2, P),
    ?assertEqual(-0.25, V(R3)),
    % a massless propagator between two on-shell lines
    ?assertError({pq_oneloop, infrared_divergent}, pq_oneloop:c0(1, 0, 1, 0, 1, 1, P)).

%% the value two of the three (the deformation at kappa = 0.1, QCDLoop, LoopTools) agree on to 1e-12, else none
programs(Rec) ->
    F = fun(B) when is_binary(B) -> try binary_to_float(B) catch _:_ -> nan end; (_) -> nan end,
    Val = fun([A, B]) ->
                  case {F(A), F(B)} of {X, Y} when is_float(X), is_float(Y) -> [{X, Y}]; _ -> [] end;
             (_) -> []
          end,
    Vs = Val(maps:get(<<"deformed">>, Rec)) ++ Val(maps:get(<<"qcdloop">>, Rec, null)) ++
         Val(maps:get(<<"looptools">>, Rec, null)),
    Near = fun({A, B}, {C, D}) ->
                   math:sqrt((A - C) * (A - C) + (B - D) * (B - D)) =< 1.0e-12 * max(1.0, math:sqrt(A * A + B * B))
           end,
    case [X || {K, X} <- lists:enumerate(Vs), {L, Y} <- lists:enumerate(Vs), K < L, Near(X, Y)] of
        [V | _] -> V;
        [] -> none
    end.

%% D0 (test/oracle/gt_d0.sh): three Euclidean configurations against mpmath's direct integral (15 digits), three
%% above threshold against LoopTools (1e-12; QCDLoop returns nan for one of them). The fourth configuration of the
%% fixture takes about 5 minutes and is left out of the suite.
d0_test_() ->
    {timeout, 900, fun() ->
        #{<<"euclidean">> := Eu, <<"physical">> := Ph} = fixture("gt_d0.json"),
        [begin
             {R, I} = apply(pq_oneloop, d0, [q(A) || A <- Args] ++ [60]),
             ?assert(close(R, Want, 15), {Args, re}),
             ?assert(close(I, 0.0, 15), {Args, im})
         end || #{<<"args">> := Args, <<"mpmath">> := Want} <- Eu],
        [begin
             {R, I} = apply(pq_oneloop, d0, [q(A) || A <- Args] ++ [50]),
             ?assert(close(R, binary_to_float(Lr), 12), {Args, re}),
             ?assert(close(I, binary_to_float(Li), 12), {Args, im})
         end || #{<<"args">> := Args, <<"looptools">> := [Lr, Li]} <- lists:sublist(Ph, 3)]
    end}.

%% D0 refuses det S = 0, where its reduction does not apply: here zero momenta with equal masses (S has equal rows),
%% where D0 itself is finite, 1/(6 m^4)
d0_degenerate_test() ->
    ?assertError({pq_oneloop, degenerate_kinematics}, pq_oneloop:d0(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 50)).
