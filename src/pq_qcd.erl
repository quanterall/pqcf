%% The renormalisation group of QCD (plan "theoretical physics instruments", QCD1): the beta function and the quark
%% mass anomalous dimension to five loops, exactly; the running of alpha_s and of MS-bar quark masses at arbitrary
%% precision; alpha_s across a heavy-quark threshold (MS-bar mass) to three loops.
%%
%% Conventions (those of RunDec and CRunDec, Herren and Steinhauser, arXiv:1703.03751): a = alpha_s / pi,
%%   da / d ln mu^2 = -sum_(i = 0..4) beta_i a^(i + 2),   d ln m / d ln mu^2 = -sum_(i = 0..4) gamma_i a^(i + 1).
%% Sources: beta_0 .. beta_2 (Gross, Wilczek, Politzer; Caswell, Jones; Tarasov, Vladimirov, Zharkov), beta_3 (van
%% Ritbergen, Vermaseren and Larin, hep-ph/9701390), beta_4 (Baikov, Chetyrkin and Kuehn, arXiv:1606.08659; Herzog,
%% Ruijl, Ueda, Vermaseren and Vogt, arXiv:1701.01404); gamma_0 .. gamma_3 (Chetyrkin, hep-ph/9703278; Vermaseren,
%% Larin and van Ritbergen, hep-ph/9703284), gamma_4 (Baikov, Chetyrkin and Kuehn, arXiv:1402.6611); decoupling to
%% three loops (Chetyrkin, Kniehl and Steinhauser, hep-ph/9706430).
%%
%% A coefficient is a list of terms {Q, Nf^K, Zetas}: Q a fraction {N, D}, K the power of n_f, Zetas a sorted list of
%% zeta arguments (e.g. [3, 3] for zeta_3^2). beta_general/2 gives beta_0 .. beta_3 for any gauge group in the colour
%% invariants of pq_lie (vRVL's formula), so that the SU(3) values above can be checked against an independent
%% transcription.
-module(pq_qcd).
-export([beta/1, gamma_m/1, coefficient/3, value/3, beta_general/2,
         alphas/6, mass/6, decouple_down/6, decouple_up/6]).

%%% The coefficients for SU(3), polynomials in n_f

beta(0) -> scale([{{11, 1}, 0, []}, {{-2, 3}, 1, []}], {1, 4});
beta(1) -> scale([{{102, 1}, 0, []}, {{-38, 3}, 1, []}], {1, 16});
beta(2) -> scale([{{2857, 2}, 0, []}, {{-5033, 18}, 1, []}, {{325, 54}, 2, []}], {1, 64});
beta(3) -> scale([{{149753, 6}, 0, []}, {{3564, 1}, 0, [3]},
                  {{-1078361, 162}, 1, []}, {{-6508, 27}, 1, [3]},
                  {{50065, 162}, 2, []}, {{6472, 81}, 2, [3]},
                  {{1093, 729}, 3, []}], {1, 256});
beta(4) -> scale([{{8157455, 16}, 0, []}, {{621885, 2}, 0, [3]}, {{-88209, 2}, 0, [4]}, {{-288090, 1}, 0, [5]},
                  {{-336460813, 1944}, 1, []}, {{-4811164, 81}, 1, [3]}, {{33935, 6}, 1, [4]}, {{1358995, 27}, 1, [5]},
                  {{25960913, 1944}, 2, []}, {{698531, 81}, 2, [3]}, {{-10526, 9}, 2, [4]}, {{-381760, 81}, 2, [5]},
                  {{-630559, 5832}, 3, []}, {{-48722, 243}, 3, [3]}, {{1618, 27}, 3, [4]}, {{460, 9}, 3, [5]},
                  {{1205, 2916}, 4, []}, {{-152, 81}, 4, [3]}], {1, 1024}).

gamma_m(0) -> [{{1, 1}, 0, []}];
gamma_m(1) -> scale([{{202, 3}, 0, []}, {{-20, 9}, 1, []}], {1, 16});
gamma_m(2) -> scale([{{1249, 1}, 0, []}, {{-2216, 27}, 1, []}, {{-160, 3}, 1, [3]}, {{-140, 81}, 2, []}], {1, 64});
gamma_m(3) -> scale([{{4603055, 162}, 0, []}, {{135680, 27}, 0, [3]}, {{-8800, 1}, 0, [5]},
                     {{-91723, 27}, 1, []}, {{-34192, 9}, 1, [3]}, {{880, 1}, 1, [4]}, {{18400, 9}, 1, [5]},
                     {{5242, 243}, 2, []}, {{800, 9}, 2, [3]}, {{-160, 3}, 2, [4]},
                     {{-332, 243}, 3, []}, {{64, 27}, 3, [3]}], {1, 256});
gamma_m(4) -> scale([{{99512327, 162}, 0, []}, {{46402466, 243}, 0, [3]}, {{96800, 1}, 0, [3, 3]},
                     {{-698126, 9}, 0, [4]}, {{-231757160, 243}, 0, [5]}, {{242000, 1}, 0, [6]}, {{412720, 1}, 0, [7]},
                     {{-150736283, 1458}, 1, []}, {{-12538016, 81}, 1, [3]}, {{-75680, 9}, 1, [3, 3]},
                     {{2038742, 27}, 1, [4]}, {{49876180, 243}, 1, [5]}, {{-638000, 9}, 1, [6]}, {{-1820000, 27}, 1, [7]},
                     {{1320742, 729}, 2, []}, {{2010824, 243}, 2, [3]}, {{46400, 27}, 2, [3, 3]}, {{-166300, 27}, 2, [4]},
                     {{-264040, 81}, 2, [5]}, {{92000, 27}, 2, [6]},
                     {{91865, 1458}, 3, []}, {{12848, 81}, 3, [3]}, {{448, 9}, 3, [4]}, {{-5120, 27}, 3, [5]},
                     {{-260, 243}, 4, []}, {{-320, 243}, 4, [3]}, {{64, 27}, 4, [4]}], {1, 1024}).

scale(Terms, F) -> [{qmul(Q, F), K, Z} || {Q, K, Z} <- Terms].

%% a coefficient at a given n_f: the terms {Q, Zetas} with n_f substituted, like terms collected
coefficient(Which, I, Nf) ->
    Terms = case Which of beta -> beta(I); gamma_m -> gamma_m(I) end,
    collect([{qmul(Q, qpow({Nf, 1}, K)), Z} || {Q, K, Z} <- Terms]).

collect(Terms) ->
    M = lists:foldl(fun({Q, Z}, Acc) -> Acc#{Z => qadd(maps:get(Z, Acc, {0, 1}), Q)} end, #{}, Terms),
    lists:sort([{Q, Z} || Z := Q <- M, element(1, Q) =/= 0]).

%% its value at P bits
value(Which, I, {Nf, P}) -> num(coefficient(Which, I, Nf), P);
value(Which, I, Nf) -> num(coefficient(Which, I, Nf), 200).

num(Terms, P) ->
    W = P + 32,
    S = lists:foldl(fun({Q, Z}, Acc) ->
                            pq_mp:add(Acc, lists:foldl(fun(K, X) -> pq_mp:mul(X, zeta(K, W), W) end, pq_mp:from(Q, W), Z), W)
                    end, {mp, 0, 0}, Terms),
    pq_mp:from(S, P).

zeta(K, W) ->
    Key = {?MODULE, zeta, K, W},
    case get(Key) of
        undefined -> V = pq_mp:zeta(K, W), put(Key, V), V;
        V -> V
    end.

%%% beta_0 .. beta_3 for any gauge group (van Ritbergen, Vermaseren and Larin, hep-ph/9701390, eq. (8)), in a = alpha_s/pi:
%%% vRVL write beta(a_s) with a_s = alpha_s/(4 pi), so beta_i here = vRVL's beta_i / 4^(i + 1).
%%% Inv: #{ca, cf, tf, d44_aa, d44_fa, d44_ff} (the d44 divided by N_A, as pq_lie gives them), exact fractions;
%%% returns the coefficient at n_f symbolic as {Q, Nf^K, Zetas} terms.
beta_general(I, Inv) ->
    #{ca := Ca, cf := Cf, tf := Tf} = Inv,
    T = fun(Q, K, Z, Factors) -> {qmul(Q, lists:foldl(fun qmul/2, {1, 1}, Factors)), K, Z} end,
    Raw = case I of
              0 -> [T({11, 3}, 0, [], [Ca]), T({-4, 3}, 1, [], [Tf])];
              1 -> [T({34, 3}, 0, [], [Ca, Ca]), T({-4, 1}, 1, [], [Cf, Tf]), T({-20, 3}, 1, [], [Ca, Tf])];
              2 -> [T({2857, 54}, 0, [], [Ca, Ca, Ca]), T({2, 1}, 1, [], [Cf, Cf, Tf]),
                    T({-205, 9}, 1, [], [Cf, Ca, Tf]), T({-1415, 27}, 1, [], [Ca, Ca, Tf]),
                    T({44, 9}, 2, [], [Cf, Tf, Tf]), T({158, 27}, 2, [], [Ca, Tf, Tf])];
              3 ->
                  #{d44_aa := Daa, d44_fa := Dfa, d44_ff := Dff} = Inv,
                  [T({150653, 486}, 0, [], [Ca, Ca, Ca, Ca]), T({-44, 9}, 0, [3], [Ca, Ca, Ca, Ca]),
                   T({-80, 9}, 0, [], [Daa]), T({704, 3}, 0, [3], [Daa]),
                   T({-39143, 81}, 1, [], [Ca, Ca, Ca, Tf]), T({136, 3}, 1, [3], [Ca, Ca, Ca, Tf]),
                   T({7073, 243}, 1, [], [Ca, Ca, Cf, Tf]), T({-656, 9}, 1, [3], [Ca, Ca, Cf, Tf]),
                   T({-4204, 27}, 1, [], [Ca, Cf, Cf, Tf]), T({352, 9}, 1, [3], [Ca, Cf, Cf, Tf]),
                   T({46, 1}, 1, [], [Cf, Cf, Cf, Tf]),
                   T({512, 9}, 1, [], [Dfa]), T({-1664, 3}, 1, [3], [Dfa]),
                   T({7930, 81}, 2, [], [Ca, Ca, Tf, Tf]), T({224, 9}, 2, [3], [Ca, Ca, Tf, Tf]),
                   T({1352, 27}, 2, [], [Cf, Cf, Tf, Tf]), T({-704, 9}, 2, [3], [Cf, Cf, Tf, Tf]),
                   T({17152, 243}, 2, [], [Ca, Cf, Tf, Tf]), T({448, 9}, 2, [3], [Ca, Cf, Tf, Tf]),
                   T({-704, 9}, 2, [], [Dff]), T({512, 3}, 2, [3], [Dff]),
                   T({424, 243}, 3, [], [Ca, Tf, Tf, Tf]), T({1232, 243}, 3, [], [Cf, Tf, Tf, Tf])]
          end,
    scale(Raw, {1, 1 bsl (2 * (I + 1))}).

%%% Running

%% alpha_s^(nf)(mu1) from alpha_s^(nf)(mu0), with Loops = 1..5 terms of beta: solve
%%   ln(mu1^2 / mu0^2) = - int_(a0)^(a1) da / beta(a)
%% by Newton's method (the derivative is -1/beta(a1)), the integral by pq_quad. Arguments exact (integers,
%% fractions) or pq_mp numbers; the result a pq_mp number at P bits (alpha_s, not a).
alphas(As0, Mu0, Mu1, Nf, Loops, P) ->
    W = P + 32,
    Pi = pq_mp:pi(W),
    A0 = pq_mp:divide(mp(As0, W), Pi, W),
    Bs = [value(beta, I, {Nf, W}) || I <- lists:seq(0, Loops - 1)],
    T = pq_mp:log(pq_mp:divide(sq(mp(Mu1, W), W), sq(mp(Mu0, W), W), W), W),
    % start: the one-loop solution
    Guess = pq_mp:divide(A0, pq_mp:add({mp, 1, 0}, pq_mp:mul(pq_mp:mul(hd(Bs), A0, W), T, W), W), W),
    A1 = newton(A0, T, Bs, Guess, W, 0),
    pq_mp:from(pq_mp:mul(A1, Pi, W), P).

newton(A0, T, Bs, A, W, I) ->
    I > 60 andalso error({pq_qcd, no_convergence}),
    % f(A) = int_(A0)^A da / beta(a) + T  (= 0 at the solution); f'(A) = 1 / beta(A)
    {Int, _} = pq_quad:integrate(fun(X) -> pq_mp:divide({mp, 1, 0}, betaf(Bs, X, W), W) end, A0, A, W - 32),
    F = pq_mp:add(Int, T, W),
    Step = pq_mp:mul(F, betaf(Bs, A, W), W),
    A2 = pq_mp:sub(A, Step, W),
    case small(Step, A2, W - 40) of
        true -> A2;
        false -> newton(A0, T, Bs, A2, W, I + 1)
    end.

%% beta(a) = sum beta_i a^(i + 2)
betaf(Bs, A, W) -> poly(Bs, A, 2, W).

poly(Cs, A, Shift, W) ->
    {S, _} = lists:foldl(fun(C, {Acc, Pw}) -> {pq_mp:add(Acc, pq_mp:mul(C, Pw, W), W), pq_mp:mul(Pw, A, W)} end,
                         {{mp, 0, 0}, ipow(A, Shift, W)}, Cs),
    S.

%% m^(nf)(mu1) from m^(nf)(mu0): ln(m1/m0) = int_(a0)^(a1) gamma(a) / beta(a) da, gamma(a) = sum gamma_i a^(i + 1)
%% (Loops terms of each), a0 and a1 the couplings at mu0 and mu1 (alpha_s, not a)
mass(M0, As0, As1, Nf, Loops, P) ->
    W = P + 32,
    Pi = pq_mp:pi(W),
    [A0, A1] = [pq_mp:divide(mp(X, W), Pi, W) || X <- [As0, As1]],
    Bs = [value(beta, I, {Nf, W}) || I <- lists:seq(0, Loops - 1)],
    Gs = [value(gamma_m, I, {Nf, W}) || I <- lists:seq(0, Loops - 1)],
    {Int, _} = pq_quad:integrate(fun(X) -> pq_mp:divide(poly(Gs, X, 1, W), betaf(Bs, X, W), W) end, A0, A1, P + 16),
    pq_mp:from(pq_mp:mul(mp(M0, W), pq_mp:exp(Int, W), W), P).

%%% Decoupling at an MS-bar threshold (CKS, hep-ph/9706430): with L = ln(mu^2 / m_h(mu)^2), m_h the MS-bar mass of
%%% the heavy quark in the n_l + 1 flavour theory, a = alpha_s^(nl + 1)(mu) / pi,
%%%   alpha_s^(nl) = alpha_s^(nl + 1) (1 + c1 a + c2 a^2 + c3 a^3),
%%% Loops as in RunDec: Loops = 1 no change, 2 .. 4 up to c1 .. c3. Upward the inverse series in a = alpha_s^(nl)/pi:
%%% d1 = -c1, d2 = 2 c1^2 - c2, d3 = -5 c1^3 + 5 c1 c2 - c3.

dec_coeffs(L, Nl, W) ->
    Z3 = zeta(3, W),
    Lp = fun(Cs) -> poly(Cs, L, 0, W) end,
    C1 = Lp([{mp, 0, 0}, pq_mp:from({-1, 6}, W)]),
    C2 = Lp([pq_mp:from({11, 72}, W), pq_mp:from({-11, 24}, W), pq_mp:from({1, 36}, W)]),
    Nlm = pq_mp:from(Nl, W),
    C3 = pq_mp:add(pq_mp:add(Lp([pq_mp:from({564731, 124416}, W), pq_mp:from({-955, 576}, W),
                                 pq_mp:from({53, 576}, W), pq_mp:from({-1, 216}, W)]),
                             pq_mp:mul(Nlm, Lp([pq_mp:from({-2633, 31104}, W), pq_mp:from({67, 576}, W),
                                                pq_mp:from({-1, 36}, W)]), W), W),
                   pq_mp:mul(pq_mp:from({-82043, 27648}, W), Z3, W), W),
    [C1, C2, C3].

%% alpha_s^(nl)(mu) from alpha_s^(nl + 1)(mu); Mh = m_h^(nl + 1)(mu), Mu = mu
decouple_down(As, Mh, Mu, Nl, Loops, P) ->
    W = P + 32,
    {A, _L, Cs} = dec_setup(As, Mh, Mu, Nl, W),
    series(As, A, lists:sublist(Cs, Loops - 1), P, W).

decouple_up(As, Mh, Mu, Nl, Loops, P) ->
    W = P + 32,
    {A, _L, [C1, C2, C3]} = dec_setup(As, Mh, Mu, Nl, W),
    M = fun(X, Y) -> pq_mp:mul(X, Y, W) end,
    D1 = pq_mp:neg(C1),
    D2 = pq_mp:sub(M(pq_mp:from(2, W), M(C1, C1)), C2, W),
    D3 = pq_mp:sub(pq_mp:add(M(pq_mp:from(-5, W), M(C1, M(C1, C1))), M(pq_mp:from(5, W), M(C1, C2)), W), C3, W),
    series(As, A, lists:sublist([D1, D2, D3], Loops - 1), P, W).

dec_setup(As, Mh, Mu, Nl, W) ->
    Asm = mp(As, W),
    A = pq_mp:divide(Asm, pq_mp:pi(W), W),
    L = pq_mp:log(pq_mp:divide(sq(mp(Mu, W), W), sq(mp(Mh, W), W), W), W),
    {A, L, dec_coeffs(L, Nl, W)}.

series(As, A, Cs, P, W) ->
    S = poly([{mp, 1, 0} | Cs], A, 0, W),
    pq_mp:from(pq_mp:mul(mp(As, W), S, W), P).

%%% Helpers

mp({mp, _, _} = X, W) -> pq_mp:from(X, W);
mp(X, W) when is_float(X) -> pq_mp:from(X, W);
mp(X, W) -> pq_mp:from(q(X), W).

sq(X, W) -> pq_mp:mul(X, X, W).

ipow(_, 0, _) -> {mp, 1, 0};
ipow(X, K, W) -> pq_mp:mul(X, ipow(X, K - 1, W), W).

small(Step, A, Bits) ->
    case Step of
        {mp, 0, _} -> true;
        _ -> mag(Step) < mag(A) - Bits
    end.

mag({mp, 0, _}) -> -1 bsl 30;
mag({mp, M, E}) -> E + length(integer_to_list(abs(M), 2)).

q(N) when is_integer(N) -> {N, 1};
q({N, D}) -> {N, D}.

qmul({A, B}, {C, D}) -> norm(A * C, B * D).
qadd({A, B}, {C, D}) -> norm(A * D + C * B, B * D).
qpow(_, 0) -> {1, 1};
qpow(X, K) -> qmul(X, qpow(X, K - 1)).
norm(N, D) when D < 0 -> norm(-N, -D);
norm(N, D) -> G = max(1, gcd(abs(N), D)), {N div G, D div G}.
gcd(A, 0) -> A;
gcd(A, B) -> gcd(B, A rem B).
