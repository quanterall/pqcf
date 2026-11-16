%% Integrals on an interval at arbitrary precision: tanh-sinh (double exponential) quadrature, Takahasi and Mori
%% (Publ. RIMS 9, 721 (1974), doi:10.2977/prims/1195192451), on pq_mp numbers. Integrable singularities at the ends
%% (logarithms, x^-1/2) cost nothing extra, so an integrand with known singular points is integrated piecewise
%% between them (integrate/5 with breakpoints).
%%
%% x = c + h2 tanh(pi/2 sinh u), weight h2 (pi/2) cosh u / cosh^2(pi/2 sinh u); the step h is halved (the new nodes
%% are the odd multiples) until two successive sums agree to the precision. Nodes are placed as end -+ distance
%% (exact near the ends) and dropped closer than 2^-(2P + 48) (hi - lo): there an integrable singularity as strong as
%% x^-1/2 contributes below 2^-P. Option ends => 3 goes to 2^-(3P + 48), for singularities up to x^-2/3 (the
%% six-dimensional box has such points where two singular curves meet). The nodes carry bits(P) = 2P + 96 bits
%% (bits(P, 3) = 3P + 96 with ends => 3; more on a piece short against its ends),
%% and the integrand should compute at that precision where it takes differences with a singular point. It gets a
%% point x (pq_mp) and returns a pq_mp number or {Re, Im} of pq_mp numbers; the integral is {Re, Im} at P bits.
-module(pq_quad).
-export([integrate/4, integrate/5, integrate/6, bits/1, bits/2]).

-define(MAX_LEVEL, 14).

%% the precision of the nodes, and the one the integrand should use near its singular points
bits(P) -> bits(P, 2).

%% the same when the nodes go to 2^-(E P + 48) from the ends (option ends => E)
bits(P, E) -> E * P + 96.

integrate(F, Lo, Hi, P) -> integrate(F, Lo, Hi, [], P).

%% the sum over the pieces between Lo, the breakpoints inside (Lo, Hi) and Hi; Lo > Hi gives minus the integral from
%% Hi to Lo
integrate(F, Lo0, Hi0, Breaks, P) -> integrate(F, Lo0, Hi0, Breaks, P, #{}).

%% Opts: #{parallel => true} evaluates each level's new nodes in parallel processes (for expensive integrands, such
%% as an integral inside an integral); the sums are the same, in the same order. For an integral inside an integral,
%% #{on_no_convergence => last, max_level => L} returns the last level's sum instead of failing after level L: it
%% fails where its singular points merge, at outer nodes so close to the outer breakpoint that their weight makes any
%% value negligible; the outer integral keeps its own convergence test.
integrate(F, Lo0, Hi0, Breaks, P, Opts) ->
    Lim = {maps:get(max_level, Opts, ?MAX_LEVEL), maps:get(on_no_convergence, Opts, error)},
    Eval = case maps:get(parallel, Opts, false) of
               true -> fun(Xs) -> pmap(fun(X) -> val(F(X)) end, Xs) end;
               false -> fun(Xs) -> [val(F(X)) || X <- Xs] end
           end,
    Ends = maps:get(ends, Opts, 2),
    W = bits(P, Ends),
    case pq_mp:cmp(pq_mp:from(Lo0, W), pq_mp:from(Hi0, W)) of
        0 -> {{mp, 0, 0}, {mp, 0, 0}};
        1 ->
            % the limits reversed: minus the integral the other way
            {R, I} = integrate1(Eval, Lim, Hi0, Lo0, Breaks, P, Ends, W),
            {pq_mp:neg(R), pq_mp:neg(I)};
        -1 -> integrate1(Eval, Lim, Lo0, Hi0, Breaks, P, Ends, W)
    end.

integrate1(Eval, Lim, Lo0, Hi0, Breaks, P, Ends, W) ->
    Lo = pq_mp:from(Lo0, W), Hi = pq_mp:from(Hi0, W),
    Inside = lists:usort(fun(A, B) -> pq_mp:cmp(A, B) =< 0 end,
                         [B || B0 <- Breaks, B <- [pq_mp:from(B0, W)], pq_mp:cmp(B, Lo) > 0, pq_mp:cmp(B, Hi) < 0]),
    Pts = [Lo | Inside] ++ [Hi],
    Parts = [piece(Eval, Lim, Ends, A, B, P, W + max(0, max(mag(A), mag(B)) - mag(pq_mp:sub(B, A, W))))
             || {A, B} <- lists:zip(lists:droplast(Pts), tl(Pts))],
    Sum = lists:foldl(fun(X, Acc) -> cadd(Acc, X, W) end, {zero(), zero()}, Parts),
    {R, I} = Sum,
    {pq_mp:from(R, P), pq_mp:from(I, P)}.

piece(Eval, Lim, Ends, Lo, Hi, P, W) ->
    C = pq_mp:ldexp(pq_mp:add(Lo, Hi, W), -1),
    H2 = pq_mp:ldexp(pq_mp:sub(Hi, Lo, W), -1),
    Pi2 = pq_mp:ldexp(pq_mp:pi(W), -1),
    Cut = pq_mp:ldexp(pq_mp:sub(Hi, Lo, W), -(Ends * P + 48)),
    Node = fun(U) -> node(U, {Lo, Hi}, H2, Pi2, Cut, W) end,
    % level 0: u = 0, +-1, +-2, ..
    [F0] = Eval([C]),
    S0 = cadd(cscale(F0, pq_mp:mul(H2, Pi2, W), W), level(Eval, Node, 1, 1, W), W),
    refine(Eval, Lim, Node, S0, 1, P, W).

%% the nodes u = k h, k = Start, Start + Step, .. on both sides until they fall below the cut: their positions first
%% (cheap), then the integrand at all of them (Eval, in parallel if asked), then the weighted sum
level(Eval, Node, Start, Step, W) ->
    Pos = positions(Node, Start, Step, []),
    Vals = Eval(lists:append([[Xp, Xm] || {Xp, Xm, _} <- Pos])),
    sum_pairs(Pos, Vals, W, {zero(), zero()}).

positions(Node, K, Step, Acc) ->
    case Node(K) of
        stop -> lists:reverse(Acc);
        Pt -> positions(Node, K + Step, Step, [Pt | Acc])
    end.

sum_pairs([], [], _, Acc) -> Acc;
sum_pairs([{_, _, Wt} | Pos], [Fp, Fm | Vals], W, Acc) ->
    sum_pairs(Pos, Vals, W, cadd(Acc, cscale(cadd(Fp, Fm, W), Wt, W), W)).

%% level L: h = 2^-L; the sum h * sum f w, the new nodes the odd multiples of h
refine(_, {Max, last}, _, Prev, L, _, W) when L > Max ->
    cscale(Prev, pq_mp:ldexp({mp, 1, 0}, -(L - 1)), W);
refine(_, {Max, _}, _, _, L, _, _) when L > Max -> error({pq_quad, no_convergence});
refine(Eval, Lim, Node, Prev, L, P, W) ->
    H = pq_mp:ldexp({mp, 1, 0}, -L),
    New = level(Eval, fun(K) -> Node(pq_mp:mul(pq_mp:from(K, W), H, W)) end, 1, 2, W),
    % Prev is the sum at step 2h without its factor: S_L = (S_(L-1) + 2 h new) where S_(L-1) = 2h * Prev
    Sum = cadd(Prev, New, W),
    PrevV = cscale(Prev, pq_mp:ldexp(H, 1), W),
    SumV = cscale(Sum, H, W),
    D = cabs1(csub(SumV, PrevV, W)),
    Scale = max(0, mag(cabs1(SumV))),
    case L >= 3 andalso pq_mp:cmp(D, pq_mp:ldexp({mp, 1, 0}, Scale - P - 4)) < 0 of
        true -> SumV;
        false -> refine(Eval, Lim, Node, Sum, L + 1, P, W)
    end.

%% the two nodes at +-u (u > 0), as end -+ distance (exact near the ends, where the singularities are), and their
%% weight without the step factor; stop when they are beyond the cut
node(K, Ends, H2, Pi2, Cut, W) when is_integer(K) -> node(pq_mp:from(K, W), Ends, H2, Pi2, Cut, W);
node(U, {Lo, Hi}, H2, Pi2, Cut, W) ->
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
            {pq_mp:sub(Hi, Dist, W), pq_mp:add(Lo, Dist, W), Wt}
    end.

%% F at each X in parallel processes, at most one per scheduler at a time; results in order; a crash is re-raised
pmap(F, Xs) ->
    N = erlang:system_info(schedulers_online),
    pmap(F, lists:enumerate(Xs), N, #{}, 0, #{}).

pmap(_, [], _, Running, _, Done) when map_size(Running) =:= 0 ->
    [V || {_, V} <- lists:sort(maps:to_list(Done))];
pmap(F, [{I, X} | Rest], N, Running, Active, Done) when Active < N ->
    Self = self(),
    {Pid, Ref} = spawn_monitor(fun() -> Self ! {self(), pq_quad_result, F(X)} end),
    pmap(F, Rest, N, Running#{Pid => {I, Ref}}, Active + 1, Done);
pmap(F, Todo, N, Running, Active, Done) ->
    receive
        {Pid, pq_quad_result, V} when is_map_key(Pid, Running) ->
            {I, Ref} = maps:get(Pid, Running),
            erlang:demonitor(Ref, [flush]),
            pmap(F, Todo, N, maps:remove(Pid, Running), Active - 1, Done#{I => V});
        {'DOWN', _, process, Pid, Reason} when is_map_key(Pid, Running) ->
            [exit(P, kill) || P <- maps:keys(Running), P =/= Pid],
            error({pq_quad, {integrand_failed, Reason}})
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
