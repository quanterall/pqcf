%% Integrals on an interval at arbitrary precision: tanh-sinh (double exponential) quadrature, Takahasi and Mori
%% (Publ. RIMS 9, 721 (1974), doi:10.2977/prims/1195192451), on pq_mp numbers. Integrable singularities at the ends
%% (logarithms, x^-1/2) cost nothing extra, so an integrand with known singular points is integrated piecewise
%% between them (integrate/5 with breakpoints).
%%
%% x = c + h2 tanh(pi/2 sinh u), weight h2 (pi/2) cosh u / cosh^2(pi/2 sinh u); the step h is halved (the new nodes
%% are the odd multiples) until two successive sums agree to the precision. Nodes are placed as end -+ distance
%% (exact near the ends) and dropped closer than 2^-(2P + 48) (hi - lo): there an integrable singularity as strong as
%% x^-1/2 contributes below 2^-P. So the nodes carry bits(P) = 2P + 96 bits (more on a piece short against its ends),
%% and the integrand should compute at that precision where it takes differences with a singular point. It gets a
%% point x (pq_mp) and returns a pq_mp number or {Re, Im} of pq_mp numbers; the integral is {Re, Im} at P bits.
-module(pq_quad).
-export([integrate/4, integrate/5, bits/1]).

-define(MAX_LEVEL, 14).

%% the precision of the nodes, and the one the integrand should use near its singular points
bits(P) -> 2 * P + 96.

integrate(F, Lo, Hi, P) -> integrate(F, Lo, Hi, [], P).

%% the sum over the pieces between Lo, the breakpoints inside (Lo, Hi) and Hi
integrate(F, Lo0, Hi0, Breaks, P) ->
    W = bits(P),
    Lo = pq_mp:from(Lo0, W), Hi = pq_mp:from(Hi0, W),
    Inside = lists:usort(fun(A, B) -> pq_mp:cmp(A, B) =< 0 end,
                         [B || B0 <- Breaks, B <- [pq_mp:from(B0, W)], pq_mp:cmp(B, Lo) > 0, pq_mp:cmp(B, Hi) < 0]),
    Pts = [Lo | Inside] ++ [Hi],
    Parts = [piece(F, A, B, P, W + max(0, max(mag(A), mag(B)) - mag(pq_mp:sub(B, A, W))))
             || {A, B} <- lists:zip(lists:droplast(Pts), tl(Pts))],
    Sum = lists:foldl(fun(X, Acc) -> cadd(Acc, X, W) end, {zero(), zero()}, Parts),
    {R, I} = Sum,
    {pq_mp:from(R, P), pq_mp:from(I, P)}.

piece(F, Lo, Hi, P, W) ->
    C = pq_mp:ldexp(pq_mp:add(Lo, Hi, W), -1),
    H2 = pq_mp:ldexp(pq_mp:sub(Hi, Lo, W), -1),
    Pi2 = pq_mp:ldexp(pq_mp:pi(W), -1),
    Cut = pq_mp:ldexp(pq_mp:sub(Hi, Lo, W), -(2 * P + 48)),
    Node = fun(U) -> node(F, U, {Lo, Hi}, H2, Pi2, Cut, W) end,
    % level 0: u = 0, +-1, +-2, ..
    S0 = cadd(node0(F, C, H2, Pi2, W), sum_side(Node, 1, 1, W), W),
    refine(Node, S0, 1, P, W).

node0(F, C, H2, Pi2, W) -> cscale(val(F(C)), pq_mp:mul(H2, Pi2, W), W).

%% the nodes u = k h, k = Start, Start + Step, .. on both sides, until both fall below the cut
sum_side(Node, Start, Step, W) -> sum_side(Node, Start, Step, W, {zero(), zero()}).
sum_side(Node, K, Step, W, Acc) ->
    case Node(K) of
        stop -> Acc;
        V -> sum_side(Node, K + Step, Step, W, cadd(Acc, V, W))
    end.

%% level L: h = 2^-L; the sum h * sum f w, the new nodes the odd multiples of h
refine(Node, Prev, L, P, W) ->
    L > ?MAX_LEVEL andalso error({pq_quad, no_convergence}),
    H = pq_mp:ldexp({mp, 1, 0}, -L),
    New = sum_side(fun(K) -> Node(pq_mp:mul(pq_mp:from(K, W), H, W)) end, 1, 2, W),
    % Prev is the sum at step 2h without its factor: S_L = (S_(L-1) + 2 h new) where S_(L-1) = 2h * Prev
    Sum = cadd(Prev, New, W),
    PrevV = cscale(Prev, pq_mp:ldexp(H, 1), W),
    SumV = cscale(Sum, H, W),
    D = cabs1(csub(SumV, PrevV, W)),
    Scale = max(0, mag(cabs1(SumV))),
    case L >= 3 andalso pq_mp:cmp(D, pq_mp:ldexp({mp, 1, 0}, Scale - P - 4)) < 0 of
        true -> SumV;
        false -> refine(Node, Sum, L + 1, P, W)
    end.

%% the two nodes at +-u (u > 0) with their weights, without the step factor; stop when both are beyond the cut
node(F, K, Ends, H2, Pi2, Cut, W) when is_integer(K) -> node(F, pq_mp:from(K, W), Ends, H2, Pi2, Cut, W);
node(F, U, {Lo, Hi}, H2, Pi2, Cut, W) ->
    Eu = pq_mp:exp(U, W),
    Sinh = pq_mp:ldexp(pq_mp:sub(Eu, pq_mp:divide({mp, 1, 0}, Eu, W), W), -1),
    Cosh = pq_mp:ldexp(pq_mp:add(Eu, pq_mp:divide({mp, 1, 0}, Eu, W), W), -1),
    V = pq_mp:mul(Pi2, Sinh, W),
    Ev = pq_mp:exp(V, W),
    Inv = pq_mp:divide({mp, 1, 0}, Ev, W),
    Ch = pq_mp:add(Ev, Inv, W),                         % 2 cosh v
    % distance of the nodes from the ends: h2 (1 - tanh v) = h2 2 e^-v / (2 cosh v)
    Dist = pq_mp:divide(pq_mp:mul(H2, pq_mp:ldexp(Inv, 1), W), Ch, W),
    case pq_mp:cmp(Dist, Cut) < 0 of
        true -> stop;
        false ->
            Wt = pq_mp:divide(pq_mp:mul(pq_mp:mul(H2, Pi2, W), pq_mp:ldexp(Cosh, 2), W), pq_mp:mul(Ch, Ch, W), W),
            % the nodes as end -+ distance: exact near the ends, where the singularities are
            Fp = val(F(pq_mp:sub(Hi, Dist, W))),
            Fm = val(F(pq_mp:add(Lo, Dist, W))),
            cscale(cadd(Fp, Fm, W), Wt, W)
    end.

val({mp, _, _} = X) -> {X, zero()};
val({R, I}) -> {R, I}.

zero() -> {mp, 0, 0}.
cadd({A, B}, {C, D}, W) -> {pq_mp:add(A, C, W), pq_mp:add(B, D, W)}.
csub({A, B}, {C, D}, W) -> {pq_mp:sub(A, C, W), pq_mp:sub(B, D, W)}.
cscale({A, B}, S, W) -> {pq_mp:mul(A, S, W), pq_mp:mul(B, S, W)}.
%% |re| + |im|, enough for a convergence test
cabs1({A, B}) -> pq_mp:add(absmp(A), absmp(B), 64).
absmp({mp, M, E}) -> {mp, abs(M), E}.
mag({mp, 0, _}) -> -1 bsl 30;
mag({mp, M, E}) -> E + length(integer_to_list(abs(M), 2)).
