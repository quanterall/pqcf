%% pq_perturb (plan "theoretical physics instruments", QM1): the anharmonic oscillator's perturbation series against
%% Bender and Wu (Phys. Rev. 184, 1231 (1969)), exactly at low order and by their large-order formula; and an exactly
%% solvable perturbation.
-module(pq_perturb_tests).
-include_lib("eunit/include/eunit.hrl").

%% H = p^2 + x^2 + g x^4. Bender and Wu's H' = p^2 + x^2/4 + lambda x^4/4 is (1/2)(p^2 + y^2 + 2 lambda y^4) after
%% x = sqrt(2) y, so their coefficients are E'_n = 2^(n-1) E_n: 3/4, -21/8, 333/16, -30885/128, 916731/256,
%% -65518401/1024, 2723294673/2048 (their Table I)
bender_wu_low_orders_test() ->
    Es = pq_perturb:ground_state(#{4 => 1}, 7),
    Want = [{3, 4}, {-21, 8}, {333, 16}, {-30885, 128}, {916731, 256}, {-65518401, 1024}, {2723294673, 2048}],
    [?assertEqual(W, norm(N * (1 bsl (K - 1)), D)) || {K, {N, D}, W} <- lists:zip3(lists:seq(1, 7), Es, Want)].

%% Bender and Wu's large-order formula, in this normalisation E_n ~ 2 (-1)^(n+1) sqrt(6)/pi^(3/2) (3/2)^n
%% Gamma(n + 1/2) (1 - 95/(72 n) + O(1/n^2)): the ratio at n = 50 and 100 within 1e-3 and 3e-4 of 1, the remainder
%% falling like 1/n^2
bender_wu_large_order_test_() ->
    {timeout, 300, fun() ->
        Es = pq_perturb:ground_state(#{4 => 1}, 100),
        R = fun(N) ->
                    {Num, Den} = lists:nth(N, Es),
                    ?assertEqual(case N rem 2 of 0 -> -1; 1 -> 1 end, sign(Num)),
                    LogA = math:log(2) + 0.5 * math:log(6) - 1.5 * math:log(math:pi()) + N * math:log(1.5) + lgamma(N + 0.5),
                    math:exp(ilog(abs(Num)) - ilog(Den) - LogA) / (1 - 95 / (72 * N))
            end,
        D50 = abs(R(50) - 1), D100 = abs(R(100) - 1),
        ?assert(D50 < 1.0e-3), ?assert(D100 < 3.0e-4),
        ?assert(D50 / D100 > 3.0 andalso D50 / D100 < 5.0)
    end}.

%% V = x^2: H = p^2 + (1 + g) x^2 has E = sqrt(1 + g), so E_n is the binomial coefficient C(1/2, n), exactly
quadratic_test() ->
    Es = pq_perturb:ground_state(#{2 => 1}, 30),
    {Want, _} = lists:mapfoldl(fun(N, C) -> C2 = qmul(C, {2 * (1 - N) + 1, 2 * N}), {C2, C2} end, {1, 1}, lists:seq(1, 30)),
    ?assertEqual(Want, Es).

qmul({A, B}, {C, D}) -> norm(A * C, B * D).
norm(N, D) when D < 0 -> norm(-N, -D);
norm(N, D) -> G = gcd(abs(N), D), {N div G, D div G}.
gcd(A, 0) -> A;
gcd(A, B) -> gcd(B, A rem B).
sign(X) when X > 0 -> 1;
sign(_) -> -1.
ilog(X) -> B = bitlen(X, 0), case B > 60 of true -> (B - 60) * math:log(2) + math:log(X bsr (B - 60)); false -> math:log(X) end.
bitlen(0, A) -> A;
bitlen(X, A) -> bitlen(X bsr 1, A + 1).
lgamma(X) -> (X - 0.5) * math:log(X) - X + 0.5 * math:log(2 * math:pi()) + 1 / (12 * X) - 1 / (360 * X * X * X) + 1 / (1260 * math:pow(X, 5)).
