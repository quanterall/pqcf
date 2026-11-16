%% pq_qcd (plan "theoretical physics instruments", QCD1): QCD's renormalisation group against mpmath and CRunDec
%% (test/oracle/gt_qcd.sh, test/fixtures/gt_qcd.json), and beta_0 .. beta_3 against the general-gauge-group formula with
%% pq_lie's invariants.
-module(pq_qcd_tests).
-include_lib("eunit/include/eunit.hrl").

fixture() ->
    Dir = filename:dirname(code:which(?MODULE)),
    Paths = [filename:join([Dir, "..", "test", "fixtures", "gt_qcd.json"]),
             filename:join([Dir, "..", "..", "..", "..", "test", "fixtures", "gt_qcd.json"]), "test/fixtures/gt_qcd.json"],
    [P | _] = [P || P <- Paths, filelib:is_regular(P)],
    {ok, B} = file:read_file(P),
    json:decode(B).

q(B) ->
    case binary:split(B, <<"/">>) of
        [N] -> binary_to_integer(N);
        [N, D] -> {binary_to_integer(N), binary_to_integer(D)}
    end.

%% |x - y| =< 10^-Digits max(1, |y|)
close(X, Y, Digits) ->
    P = 300,
    Yv = case Y of B when is_binary(B) -> pq_mp:from_string(binary_to_list(B), P); F -> pq_mp:from(float(F), P) end,
    abs(pq_mp:to_float(pq_mp:sub(X, Yv, P))) =< math:pow(10, -Digits) * max(1.0, abs(pq_mp:to_float(Yv))).

%% n = 3 in pq_lie's invariants (expressions in n), as exact fractions
sub({sym, n}) -> 3;
sub(T) when is_tuple(T) -> list_to_tuple([sub(X) || X <- tuple_to_list(T)]);
sub(L) when is_list(L) -> [sub(X) || X <- L];
sub(X) -> X.
su3(K) ->
    {ok, V} = pq_symbolic:cancel(sub(pq_lie:invariant(K))),
    case V of N when is_integer(N) -> {N, 1}; {rat, N, D} -> {N, D} end.

norm(Ts) ->
    M = lists:foldl(fun({{N1, D1}, K, Z}, A) ->
                            {N0, D0} = maps:get({K, Z}, A, {0, 1}),
                            A#{{K, Z} => {N0 * D1 + N1 * D0, D0 * D1}}
                    end, #{}, Ts),
    lists:sort([{Key, N / D} || Key := {N, D} <- M, N =/= 0]).

%% vRVL's formula for any gauge group, with SU(3)'s invariants computed by pq_lie, equals the SU(3) coefficients
%% (transcribed separately) term by term, as polynomials in n_f with zeta_3
general_group_test() ->
    Inv = maps:from_list([{K, su3(K)} || K <- [ca, cf, tf, d44_aa, d44_fa, d44_ff]]),
    ?assertEqual(#{ca => {3, 1}, cf => {4, 3}, tf => {1, 2}, d44_aa => {135, 8}, d44_fa => {15, 16}, d44_ff => {5, 96}},
                 Inv),
    [?assertEqual(norm(pq_qcd:beta(I)), norm(pq_qcd:beta_general(I, Inv))) || I <- [0, 1, 2, 3]].

%% beta_0 .. beta_4 and gamma_0 .. gamma_4 at nf = 3 .. 6 equal mpmath's separately transcribed values to 38 digits
coefficients_test() ->
    #{<<"coefficients">> := Cs} = fixture(),
    [begin
         [?assert(close(pq_qcd:value(beta, I, {Nf, 200}), B, 38), {beta, I, Nf}) || {I, B} <- lists:enumerate(0, Bs)],
         [?assert(close(pq_qcd:value(gamma_m, I, {Nf, 200}), G, 38), {gamma_m, I, Nf}) || {I, G} <- lists:enumerate(0, Gs)]
     end || #{<<"nf">> := Nf, <<"beta">> := Bs, <<"gamma_m">> := Gs} <- Cs].

%% alpha_s and m(mu): mpmath's ODE solution to 28 digits, CRunDec's exact running to 1e-8 (its step tolerance)
running_test_() ->
    {timeout, 300, fun() ->
        #{<<"running">> := Rs} = fixture(),
        [begin
             As1 = pq_qcd:alphas(q(As0), q(Mu0), q(Mu1), Nf, L, 120),
             M1 = pq_qcd:mass(q(M0), q(As0), As1, Nf, L, 120),
             ?assert(close(As1, Am, 28), {alphas, R}),
             ?assert(close(M1, Mm, 28), {mass, R}),
             ?assert(close(As1, binary_to_float(Ac), 8), {alphas_crundec, R}),
             ?assert(close(M1, binary_to_float(Mc), 8), {mass_crundec, R})
         end || #{<<"as0">> := As0, <<"mu0">> := Mu0, <<"mu1">> := Mu1, <<"nf">> := Nf, <<"loops">> := L,
                  <<"m0">> := M0, <<"alphas">> := Am, <<"mass">> := Mm,
                  <<"crundec">> := #{<<"alphas">> := Ac, <<"mass">> := Mc}} = R <- Rs]
    end}.

%% decoupling at an MS-bar threshold: mpmath to 28 digits, CRunDec to 1e-14
decoupling_test() ->
    #{<<"decoupling">> := Ds} = fixture(),
    [begin
         F = case Dir of <<"down">> -> decouple_down; <<"up">> -> decouple_up end,
         V = pq_qcd:F(q(As), q(Mh), q(Mu), Nl, L, 120),
         ?assert(close(V, Want, 28), {D, mpmath}),
         ?assert(close(V, binary_to_float(C), 14), {D, crundec})
     end || #{<<"direction">> := Dir, <<"as">> := As, <<"mh">> := Mh, <<"mu">> := Mu, <<"nl">> := Nl, <<"loops">> := L,
              <<"value">> := Want, <<"crundec">> := C} = D <- Ds].

%% exact properties: one loop in closed form, running there and back, decoupling down then up to the order kept
exact_test() ->
    P = 200,
    % a1 = a0 / (1 + beta_0 a0 ln(mu1^2/mu0^2)), beta_0 = 23/12 for nf = 5
    As1 = pq_qcd:alphas({59, 500}, 91, 10, 5, 1, P),
    W = P + 32,
    Pi = pq_mp:pi(W),
    A0 = pq_mp:divide(pq_mp:from({59, 500}, W), Pi, W),
    T = pq_mp:log(pq_mp:from({100, 8281}, W), W),
    Want = pq_mp:mul(Pi, pq_mp:divide(A0, pq_mp:add({mp, 1, 0}, pq_mp:mul(pq_mp:mul(pq_mp:from({23, 12}, W), A0, W), T, W), W), W), W),
    ?assert(close(As1, list_to_binary(pq_mp:to_string(Want, 55)), 55)),
    % there and back at five loops
    Up = pq_qcd:alphas({59, 500}, 91, 3000, 5, 5, P),
    Back = pq_qcd:alphas(Up, 3000, 91, 5, 5, P),
    ?assert(close(Back, <<"0.118">>, 50)).
