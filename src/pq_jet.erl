%% The Taylor-series method for nonlinear systems of ordinary differential equations in arbitrary precision (decision
%% record SY2, phase 6): y_i' = F_i(x, y), the F_i expressions of x, the y_i and constants with +, *, integer and real
%% powers, exp, log, sin, cos, sinh, cosh, sqrt.
%%
%% At each point the solution's Taylor coefficients come from automatic differentiation on truncated series (jets):
%% every node of the expressions gets its coefficients order by order by the recurrences of its operation (products by
%% convolution; exp, log, powers, sin/cos, sinh/cosh by their differential equations), and y_(i, k+1) = F_(i, k) / (k+1).
%% The step is chosen from the last coefficients so that the truncation is below the working precision.
%% (R. Barrio, Performance of the Taylor series method for ODEs/DAEs, Appl. Math. Comput. 163 (2005) 525.)
-module(pq_jet).
-export([integrate/6, integrate/7, diff/2]).

%% integrate(Rhs, X, Ys, X0, Y0, X1, Opts) -- Rhs: [expression] (one per y_i); X: the independent symbol; Ys: the state
%% symbols; X0, X1: start and end (X1 may be below X0); Y0: the start values (numbers or mp). Opts: digits => D
%% (default 30), order => K (default from the precision). {ok, [mp]}, the state at X1.
-spec integrate([term()], atom(), [atom()], term(), [term()], term()) -> {ok, [pq_mp:mp()]} | {error, term()}.
integrate(Rhs, X, Ys, X0, Y0, X1) -> integrate(Rhs, X, Ys, X0, Y0, X1, #{}).

-spec integrate([term()], atom(), [atom()], term(), [term()], term(), map()) -> {ok, [pq_mp:mp()]} | {error, term()}.

integrate(Rhs, X, Ys, X0, Y0, X1, Opts) ->
    Digits = maps:get(digits, Opts, 30),
    P = ceil(Digits * 3.3219280948873622) + 24,
    W = P + 32,
    K = maps:get(order, Opts, max(20, ceil(Digits * 0.9))),
    try
        Prog = compile(Rhs, X, Ys),
        Xa = num(X0, W), Xb = num(X1, W),
        Dir = case pq_mp:cmp(Xb, Xa) of -1 -> -1; _ -> 1 end,
        Ys1 = run(Prog, Xa, Xb, [num(V, W) || V <- Y0], Dir, K, W, P, 0),
        {ok, [pq_mp:from(V, P) || V <- Ys1]}
    catch
        throw:{jet, Why} -> {error, {symbolic_unsupported, #{jet => Why}}};
        error:badarith -> {error, {non_finite, jet}}
    end.

num({mp, _, _} = M, W) -> pq_mp:from(M, W);
num(S, W) when is_list(S) -> pq_mp:from_string(S, W);
num(V, W) -> pq_mp:from(V, W).

%%% The expressions as a program: a list of nodes in evaluation order, shared subexpressions once

compile(Rhs, X, Ys) ->
    Index = maps:from_list([{Y, I} || {I, Y} <- lists:enumerate(Ys)]),
    {Outs, St} = lists:mapfoldl(fun(E, S) -> node(E, X, Index, S) end, #{nodes => [], memo => #{}}, Rhs),
    #{nodes => lists:reverse(maps:get(nodes, St)), outs => Outs, n => length(Ys)}.

%% a node: {Id, Op}; ops: {x}, {y, I}, {const, Value}, {add, [Ids]}, {mul, A, B}, {div, A, B}, {pow, A, Real},
%% {exp, A}, {log, A}, {sin, A}, {cos, A}, {sinh, A}, {cosh, A} (sin/cos and sinh/cosh computed in pairs)
node(E, X, Index, St) ->
    case maps:find(E, maps:get(memo, St)) of
        {ok, Id} -> {Id, St};
        error ->
            {Op, St2} = op(E, X, Index, St),
            Id = length(maps:get(nodes, St2)),
            {Id, St2#{nodes => [{Id, Op} | maps:get(nodes, St2)], memo => (maps:get(memo, St2))#{E => Id}}}
    end.

op({sym, X}, X, _, St) -> {{x}, St};
op({sym, Y}, _, Index, St) ->
    case maps:find(Y, Index) of
        {ok, I} -> {{y, I}, St};
        error -> throw({jet, {unbound, Y}})
    end;
op(N, _, _, St) when is_integer(N) -> {{const, N}, St};
op(F, _, _, St) when is_float(F) -> {{const, F}, St};
op({rat, A, B}, _, _, St) -> {{const, {A, B}}, St};
op({const, pi}, _, _, St) -> {{const, pi}, St};
op({add, L}, X, Index, St) ->
    {Ids, St2} = lists:mapfoldl(fun(T, S) -> node(T, X, Index, S) end, St, L),
    {{add, Ids}, St2};
op({mul, [A]}, X, Index, St) -> op(A, X, Index, St);
op({mul, [A | Rest]}, X, Index, St) ->
    {Ia, St1} = node(A, X, Index, St),
    {Ib, St2} = node({mul, Rest}, X, Index, St1),
    {{mul, Ia, Ib}, St2};
op({pow, A, K}, X, Index, St) when is_integer(K), K < 0 ->
    {Ia, St1} = node({pow, A, -K}, X, Index, St),
    {I1, St2} = node(1, X, Index, St1),
    {{'div', I1, Ia}, St2};
op({pow, A, 1}, X, Index, St) -> op(A, X, Index, St);
op({pow, A, 2}, X, Index, St) ->
    {Ia, St1} = node(A, X, Index, St),
    {{mul, Ia, Ia}, St1};
op({pow, A, K}, X, Index, St) when is_integer(K), K > 2 ->
    {Ia, St1} = node(A, X, Index, St),
    {Ib, St2} = node({pow, A, K - 1}, X, Index, St1),
    {{mul, Ia, Ib}, St2};
op({pow, A, {rat, P, Q}}, X, Index, St) ->
    {Ia, St1} = node(A, X, Index, St),
    {{pow, Ia, {P, Q}}, St1};
op({pow, A, R}, X, Index, St) when is_float(R) ->
    {Ia, St1} = node(A, X, Index, St),
    {{pow, Ia, R}, St1};
op({pow, A, E}, X, Index, St) ->
    % a non-constant exponent: exp(E log A)
    op({exp, {mul, [E, {log, A}]}}, X, Index, St);
op({sqrt, A}, X, Index, St) -> op({pow, A, {rat, 1, 2}}, X, Index, St);
op({F, A}, X, Index, St) when F =:= exp; F =:= log; F =:= sin; F =:= cos; F =:= sinh; F =:= cosh ->
    {Ia, St1} = node(A, X, Index, St),
    {{F, Ia}, St1};
op(E, _, _, _) -> throw({jet, {unsupported, case E of T when is_tuple(T) -> element(1, T); _ -> E end}}).

%%% One Taylor step

run(Prog, X, Xb, Y, Dir, K, W, P, Steps) ->
    case Steps > 100000 of true -> throw({jet, too_many_steps}); false -> ok end,
    Rest = pq_mp:sub(Xb, X, W),
    case pq_mp:to_float(Rest) * Dir =< 0 of
        true -> Y;
        false ->
            Coefs = jets(Prog, X, Y, K, W),      % per state variable: [a_0 .. a_K]
            H0 = step(Coefs, K, P),
            H = case abs(pq_mp:to_float(Rest)) =< H0 of
                    true -> Rest;
                    false -> pq_mp:from(Dir * H0, W)
                end,
            Y2 = [horner(Cs, H, W) || Cs <- Coefs],
            run(Prog, pq_mp:add(X, H, W), Xb, Y2, Dir, K, W, P, Steps + 1)
    end.

horner(Cs, H, W) -> lists:foldl(fun(C, Acc) -> pq_mp:add(pq_mp:mul(Acc, H, W), C, W) end, {mp, 0, 0}, lists:reverse(Cs)).

%% the step: the last two coefficients of every variable, |a_k| h^k below 2^-P relative to the solution's size
step(Coefs, K, P) ->
    Hs = [begin
              Scale = max(1.0e-300, lists:max([abs(pq_mp:to_float(C)) || C <- lists:sublist(Cs, 2)])),
              Hk = [estimate(lists:nth(J + 1, Cs), J, Scale, P) || J <- [K - 1, K]],
              lists:min(Hk)
          end || Cs <- Coefs],
    0.7 * lists:min(Hs).

estimate(C, J, Scale, P) ->
    A = abs(to_f(C)),
    case A == 0.0 of
        true -> 1.0e300;
        false -> math:exp((math:log(Scale) - P * 0.6931471805599453 - math:log(A)) / J)
    end.

%% a float of an mp of any size (to_float would overflow): the logarithm route for very large or small values
to_f({mp, 0, _}) -> 0.0;
to_f({mp, M, E}) ->
    L = length(integer_to_list(abs(M), 2)),
    Mf = (M bsr max(0, L - 60)) * 1.0,
    Ex = E + max(0, L - 60),
    case Ex > 1000 orelse Ex < -1000 of
        true -> math:exp(math:log(abs(Mf)) + Ex * 0.6931471805599453) * case Mf < 0 of true -> -1; false -> 1 end;
        false -> Mf * math:pow(2, Ex)
    end.

%% the Taylor coefficients of the solution through x = X0, y(X0) = Y0, to order K
jets(#{nodes := Nodes, outs := Outs, n := N}, X0, Y0, K, W) ->
    % y_i coefficients, grown order by order; node coefficients likewise
    Yc0 = maps:from_list([{I, #{0 => V}} || {I, V} <- lists:enumerate(Y0)]),
    Yc = lists:foldl(fun(Ord, Ycs) ->
                             Vals = eval_nodes(Nodes, X0, Ycs, Ord, W),
                             % y_(i, Ord+1) = F_(i, Ord) / (Ord + 1)
                             maps:map(fun(I, Ci) ->
                                              Fi = maps:get(Ord, maps:get(lists:nth(I, Outs), Vals)),
                                              Ci#{Ord + 1 => pq_mp:divide(Fi, pq_mp:from(Ord + 1, W), W)}
                                      end, Ycs)
                     end, Yc0, lists:seq(0, K - 1)),
    erase({?MODULE, cache}),
    [[maps:get(J, maps:get(I, Yc)) || J <- lists:seq(0, K)] || I <- lists:seq(1, N)].

%% every node's coefficients up to order Ord (those below Ord cached between orders of one step)
eval_nodes(Nodes, X0, Ycs, Ord, W) ->
    Cache0 = case get({?MODULE, cache}) of undefined -> #{}; C -> C end,
    Cache = lists:foldl(fun({Id, Op}, Acc) ->
                                Prev = maps:get(Id, Acc, #{}),
                                V = coef(Op, Ord, Acc, Prev, X0, Ycs, W),
                                Acc#{Id => V}
                        end, Cache0, Nodes),
    put({?MODULE, cache}, Cache),
    Cache.

%% the node's coefficient map with order Ord added
coef({x}, Ord, _, Prev, X0, _, W) ->
    Prev#{Ord => case Ord of 0 -> pq_mp:from(X0, W); 1 -> {mp, 1, 0}; _ -> {mp, 0, 0} end};
coef({y, I}, Ord, _, Prev, _, Ycs, _) -> Prev#{Ord => maps:get(Ord, maps:get(I, Ycs))};
coef({const, V}, Ord, _, Prev, _, _, W) ->
    Prev#{Ord => case Ord of 0 -> case V of pi -> pq_mp:pi(W); _ -> pq_mp:from(V, W) end; _ -> {mp, 0, 0} end};
coef({add, Ids}, Ord, Acc, Prev, _, _, W) ->
    Prev#{Ord => lists:foldl(fun(Id, S) -> pq_mp:add(S, c(Acc, Id, Ord), W) end, {mp, 0, 0}, Ids)};
coef({mul, A, B}, Ord, Acc, Prev, _, _, W) ->
    Prev#{Ord => conv(Acc, A, B, Ord, W)};
coef({'div', A, B}, Ord, Acc, Prev, _, _, W) ->
    % w = a / b: w_k = (a_k - sum_(j=1..k) b_j w_(k-j)) / b_0
    S = lists:foldl(fun(J, Sm) -> pq_mp:add(Sm, pq_mp:mul(c(Acc, B, J), maps:get(Ord - J, Prev), W), W) end,
                    {mp, 0, 0}, lists:seq(1, Ord)),
    Prev#{Ord => pq_mp:divide(pq_mp:sub(c(Acc, A, Ord), S, W), c(Acc, B, 0), W)};
coef({pow, A, R}, Ord, Acc, Prev, _, _, W) ->
    U0 = c(Acc, A, 0),
    Rm = pq_mp:from(R, W),
    case Ord of
        0 -> Prev#{0 => power(U0, Rm, W)};
        _ ->
            % w_k = (1/(k u_0)) sum_(j=0..k-1) (r (k-j) - j) u_(k-j) w_j
            S = lists:foldl(fun(J, Sm) ->
                                    F = pq_mp:sub(pq_mp:mul(Rm, pq_mp:from(Ord - J, W), W), pq_mp:from(J, W), W),
                                    pq_mp:add(Sm, pq_mp:mul(F, pq_mp:mul(c(Acc, A, Ord - J), maps:get(J, Prev), W), W), W)
                            end, {mp, 0, 0}, lists:seq(0, Ord - 1)),
            Prev#{Ord => pq_mp:divide(S, pq_mp:mul(pq_mp:from(Ord, W), U0, W), W)}
    end;
coef({exp, A}, Ord, Acc, Prev, _, _, W) ->
    case Ord of
        0 -> Prev#{0 => pq_mp:exp(c(Acc, A, 0), W)};
        _ -> Prev#{Ord => dsum(Acc, A, Prev, Ord, 1, W)}         % w_k = (1/k) sum j u_j w_(k-j)
    end;
coef({log, A}, Ord, Acc, Prev, _, _, W) ->
    U0 = c(Acc, A, 0),
    case Ord of
        0 -> Prev#{0 => pq_mp:log(U0, W)};
        _ ->
            % w_k = (u_k - (1/k) sum_(j=1..k-1) j w_j u_(k-j)) / u_0
            S = lists:foldl(fun(J, Sm) -> pq_mp:add(Sm, pq_mp:mul(pq_mp:from(J, W), pq_mp:mul(maps:get(J, Prev), c(Acc, A, Ord - J), W), W), W) end,
                            {mp, 0, 0}, lists:seq(1, Ord - 1)),
            Prev#{Ord => pq_mp:divide(pq_mp:sub(c(Acc, A, Ord), pq_mp:divide(S, pq_mp:from(Ord, W), W), W), U0, W)}
    end;
coef({F, A}, Ord, Acc, Prev, _, _, W) when F =:= sin; F =:= cos; F =:= sinh; F =:= cosh ->
    % the pair (s, c) with s' = c u', c' = -+ s u'; the partner's coefficients kept under '$pair'
    Pair0 = maps:get('$pair', Prev, #{}),
    case Ord of
        0 ->
            U0 = c(Acc, A, 0),
            {S0, C0} = case F of
                           sin -> {pq_mp:sin(U0, W), pq_mp:cos(U0, W)};
                           cos -> {pq_mp:cos(U0, W), pq_mp:sin(U0, W)};
                           sinh -> {sinh(U0, W), cosh(U0, W)};
                           cosh -> {cosh(U0, W), sinh(U0, W)}
                       end,
            Prev#{0 => S0, '$pair' => #{0 => C0}};
        _ ->
            % own' = sgn1 partner u', partner' = sgn2 own u'
            {Sg1, Sg2} = case F of sin -> {1, -1}; cos -> {-1, 1}; sinh -> {1, 1}; cosh -> {1, 1} end,
            Own = dsum(Acc, A, Pair0, Ord, Sg1, W),
            Partner = dsum(Acc, A, Prev, Ord, Sg2, W),
            Prev#{Ord => Own, '$pair' => Pair0#{Ord => Partner}}
    end.

%% (Sg/k) sum_(j=1..k) j u_j w_(k-j), w's coefficients in Map
dsum(Acc, A, Map, K, Sg, W) ->
    S = lists:foldl(fun(J, Sm) -> pq_mp:add(Sm, pq_mp:mul(pq_mp:from(J, W), pq_mp:mul(c(Acc, A, J), maps:get(K - J, Map), W), W), W) end,
                    {mp, 0, 0}, lists:seq(1, K)),
    pq_mp:mul(pq_mp:from({Sg, K}, W), S, W).

conv(Acc, A, B, K, W) ->
    lists:foldl(fun(J, S) -> pq_mp:add(S, pq_mp:mul(c(Acc, A, J), c(Acc, B, K - J), W), W) end, {mp, 0, 0}, lists:seq(0, K)).

c(Acc, Id, K) -> maps:get(K, maps:get(Id, Acc)).

power(U, {mp, M, E}, W) when E >= 0 ->
    K = M bsl E,
    case K >= 0 of
        true -> ipow(U, K, W);
        false -> pq_mp:divide({mp, 1, 0}, ipow(U, -K, W), W)
    end;
power(U, R, W) -> pq_mp:exp(pq_mp:mul(R, pq_mp:log(U, W), W), W).

ipow(_, 0, _) -> {mp, 1, 0};
ipow(X, K, W) -> pq_mp:mul(X, ipow(X, K - 1, W), W).

sinh(X, W) -> pq_mp:ldexp(pq_mp:sub(pq_mp:exp(X, W), pq_mp:exp(pq_mp:neg(X), W), W), -1).
cosh(X, W) -> pq_mp:ldexp(pq_mp:add(pq_mp:exp(X, W), pq_mp:exp(pq_mp:neg(X), W), W), -1).

%%% The derivative of an expression in X (for the right-hand sides, which may contain log, which the engine's
%%% canonical form does not)

-spec diff(term(), atom()) -> term().
diff({sym, X}, X) -> 1;
diff({sym, _}, _) -> 0;
diff(N, _) when is_number(N) -> 0;
diff({rat, _, _}, _) -> 0;
diff({const, _}, _) -> 0;
diff({add, L}, X) -> {add, [diff(T, X) || T <- L]};
diff({mul, []}, _) -> 0;
diff({mul, [F]}, X) -> diff(F, X);
diff({mul, [F | Rest]}, X) -> {add, [{mul, [diff(F, X) | Rest]}, {mul, [F, diff({mul, Rest}, X)]}]};
diff({pow, F, K}, X) when is_integer(K); is_float(K) -> {mul, [K, {pow, F, K - 1}, diff(F, X)]};
diff({pow, F, {rat, P, Q}}, X) -> {mul, [{rat, P, Q}, {pow, F, {rat, P - Q, Q}}, diff(F, X)]};
diff({pow, F, G}, X) -> diff({exp, {mul, [G, {log, F}]}}, X);
diff({sqrt, F}, X) -> diff({pow, F, {rat, 1, 2}}, X);
diff({exp, F}, X) -> {mul, [{exp, F}, diff(F, X)]};
diff({log, F}, X) -> {mul, [{pow, F, -1}, diff(F, X)]};
diff({sin, F}, X) -> {mul, [{cos, F}, diff(F, X)]};
diff({cos, F}, X) -> {mul, [-1, {sin, F}, diff(F, X)]};
diff({sinh, F}, X) -> {mul, [{cosh, F}, diff(F, X)]};
diff({cosh, F}, X) -> {mul, [{sinh, F}, diff(F, X)]}.
