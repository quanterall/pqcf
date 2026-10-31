%% Chi-square goodness of fit for the statistical tests of FS 12.4: bins with expected count below 5 pooled (the
%% smallest first, until the pool reaches 5), p-value from the regularized upper incomplete gamma function.
-module(pq_chi2).
-export([sf/2, p_value/2, counts/1]).

%% P(chi2_k > X)
-spec sf(float(), pos_integer()) -> float().
sf(X, _K) when X =< 0 -> 1.0;
sf(X, K) -> q(K / 2, X / 2).

%% Probs: #{Outcome => p} (a missing mass goes into its own bin); Counts: #{Outcome => n}
-spec p_value(#{term() => float()}, #{term() => non_neg_integer()}) -> {float(), map()}.
p_value(Probs, Counts) ->
    N = lists:sum(maps:values(Counts)),
    Rest = max(0.0, 1.0 - lists:sum(maps:values(Probs))),
    Outside = lists:sum([C || K := C <- Counts, not maps:is_key(K, Probs)]),
    Bins0 = [{P * N, maps:get(K, Counts, 0)} || K := P <- Probs] ++ [{Rest * N, Outside} || Rest * N > 0 orelse Outside > 0],
    Bins = pool(lists:sort(Bins0)),
    Stat = lists:sum([(O - E) * (O - E) / E || {E, O} <- Bins]),
    Dof = length(Bins) - 1,
    {sf(Stat, Dof), #{stat => Stat, dof => Dof, shots => N}}.

-spec counts([term()]) -> #{term() => pos_integer()}.
counts(Samples) -> lists:foldl(fun(S, M) -> maps:update_with(S, fun(V) -> V + 1 end, 1, M) end, #{}, Samples).

%% ascending by expected count: merge the smallest bins until each has at least 5 expected
pool([{E1, O1}, {E2, O2} | T]) when E1 < 5 -> pool(lists:sort([{E1 + E2, O1 + O2} | T]));
pool([{E, _O} = B]) when E < 5 -> [B];
pool(Bins) -> Bins.

%%% Q(a, x), Numerical Recipes 6.2: series below a + 1, continued fraction (modified Lentz) above

q(A, X) when X < A + 1 -> 1.0 - p_series(A, X);
q(A, X) -> cf(A, X).

p_series(A, X) ->
    {Sum, _} = lists:foldl(fun(_, {S, {Ap, Del}}) ->
                                   Ap2 = Ap + 1, Del2 = Del * X / Ap2, {S + Del2, {Ap2, Del2}}
                           end, {1 / A, {A, 1 / A}}, lists:seq(1, 10000)),
    Sum * math:exp(-X + A * math:log(X) - lgamma(A)).

cf(A, X) ->
    Tiny = 1.0e-300,
    B0 = X + 1 - A,
    {H, _, _, _} = lists:foldl(fun(I, {H, B, C, D}) ->
                                       An = -I * (I - A),
                                       B2 = B + 2,
                                       D1 = An * D + B2, D2 = case abs(D1) < Tiny of true -> Tiny; false -> D1 end,
                                       C1 = B2 + An / C, C2 = case abs(C1) < Tiny of true -> Tiny; false -> C1 end,
                                       Dd = 1 / D2,
                                       {H * Dd * C2, B2, C2, Dd}
                               end, {1 / B0, B0, 1 / Tiny, 1 / B0}, lists:seq(1, 10000)),
    math:exp(-X + A * math:log(X) - lgamma(A)) * H.

%% Lanczos (g = 7, n = 9)
lgamma(X) when X < 0.5 -> math:log(math:pi() / abs(math:sin(math:pi() * X))) - lgamma(1 - X);
lgamma(X0) ->
    Cs = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059,
          12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7],
    X = X0 - 1,
    {A, _} = lists:foldl(fun(C, {S, I}) -> {S + C / (X + I), I + 1} end, {hd(Cs), 1}, tl(Cs)),
    T = X + 7.5,
    0.5 * math:log(2 * math:pi()) + (X + 0.5) * math:log(T) - T + math:log(A).
