%% Rayleigh-Schroedinger perturbation theory to high order, exactly (plan "theoretical physics instruments", QM1).
%%
%% For H = p^2 + x^2 + g V(x), V a polynomial with rational coefficients, the ground state's energy as a power series
%% E(g) = 1 + sum_n E_n g^n. Writing psi = e^(-x^2/2) sum_n g^n B_n(x), the Schroedinger equation -psi'' + x^2 psi +
%% g V psi = E psi becomes, order by order (Bender and Wu, Phys. Rev. 184, 1231 (1969)),
%%   L B_n = sum_(k=1..n) E_k B_(n-k) - V B_(n-1),   L = -d^2/dx^2 + 2x d/dx,   B_0 = 1,
%% with L x^m = 2m x^m - m(m-1) x^(m-2): triangular on monomials, so B_n is found from its highest degree down and E_n
%% from the constant term (L kills constants; B_n(0) = 0 for n >= 1, the intermediate normalisation). Every
%% coefficient is an exact rational; hundreds of orders take seconds.
-module(pq_perturb).
-export([ground_state/2]).

%% [E_1, ..., E_N] for the potential V (a map Power => coefficient, the coefficients integers or {Num, Den})
ground_state(V, N) ->
    Vq = maps:map(fun(_, C) -> q(C) end, V),
    {Es, _} = lists:foldl(fun(K, {Es0, Bs}) ->
                                  {E, B} = order(K, Es0, Bs, Vq),
                                  {Es0 ++ [E], Bs ++ [B]}
                          end, {[], [#{0 => {1, 1}}]}, lists:seq(1, N)),
    Es.

%% E_n and B_n from E_1 .. E_(n-1) and B_0 .. B_(n-1)
order(N, Es, Bs, V) ->
    % R = sum_(k=1..n-1) E_k B_(n-k) - V B_(n-1); then L B_n = E_n B_0 + R = E_n + R
    R0 = lists:foldl(fun(K, Acc) -> padd(Acc, pscale(lists:nth(N - K + 1, Bs), lists:nth(K, Es))) end, #{}, lists:seq(1, N - 1)),
    R = padd(R0, pscale(pmul(V, lists:nth(N, Bs)), {-1, 1})),
    % solve L B = R from the top degree down: the coefficient of x^m in L B is 2m b_m - (m+2)(m+1) b_(m+2)
    Top = lists:max([0 | maps:keys(R)]),
    {B, Rest} = lists:foldl(fun(M, {Bacc, _}) when M =:= 0 -> {Bacc, ok};
                               (M, {Bacc, Ok}) ->
                                    Rm = maps:get(M, R, {0, 1}),
                                    Up = q_mul({(M + 2) * (M + 1), 1}, maps:get(M + 2, Bacc, {0, 1})),
                                    Bm = q_mul(q_add(Rm, Up), {1, 2 * M}),
                                    {case Bm of {0, _} -> Bacc; _ -> Bacc#{M => Bm} end, Ok}
                            end, {#{}, ok}, lists:seq(Top, 0, -1)),
    ok = Rest,
    % the constant term: 0 * b_0 - 2 b_2 = E_n + R_0, so E_n = -2 b_2 - R_0
    En = q_add(q_mul({-2, 1}, maps:get(2, B, {0, 1})), q_mul({-1, 1}, maps:get(0, R, {0, 1}))),
    {En, B}.

%%% polynomials in x: maps Power => {Num, Den}

padd(A, B) -> maps:filter(fun(_, {N, _}) -> N =/= 0 end,
                          maps:fold(fun(K, V, M) -> M#{K => case M of #{K := W} -> q_add(V, W); _ -> V end} end, A, B)).
pscale(A, C) -> maps:filter(fun(_, {N, _}) -> N =/= 0 end, maps:map(fun(_, V) -> q_mul(V, C) end, A)).
pmul(A, B) -> maps:fold(fun(K1, V1, M) -> maps:fold(fun(K2, V2, M2) ->
                                                            K = K1 + K2, V = q_mul(V1, V2),
                                                            M2#{K => case M2 of #{K := W} -> q_add(V, W); _ -> V end}
                                                    end, M, B) end, #{}, A).

q({N, D}) -> norm(N, D);
q(N) when is_integer(N) -> {N, 1}.
q_add({A, B}, {C, D}) -> norm(A * D + C * B, B * D).
q_mul({A, B}, {C, D}) -> norm(A * C, B * D).
norm(N, D) when D < 0 -> norm(-N, -D);
norm(N, D) -> G = gcd(abs(N), D), {N div G, D div G}.
gcd(A, 0) -> A;
gcd(A, B) -> gcd(B, A rem B).
