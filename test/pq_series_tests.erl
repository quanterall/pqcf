%% pq_series (decision record SY2, phase 2) against SymPy 1.14 and mpmath (test/oracle/gt_series.py,
%% test/fixtures/gt_series.json): every term's coefficient, the parameters at the fixture's rationals, the atoms
%% (log 2, pi, Euler's gamma, zeta(k), Gamma(q), ...) evaluated in floating point.
-module(pq_series_tests).
-include_lib("eunit/include/eunit.hrl").

fixture() ->
    Dir = filename:dirname(code:which(?MODULE)),
    Paths = [filename:join([Dir, "..", "test", "fixtures", "gt_series.json"]),
             filename:join([Dir, "..", "..", "..", "..", "test", "fixtures", "gt_series.json"]),
             "test/fixtures/gt_series.json"],
    [P | _] = [P || P <- Paths, filelib:is_regular(P)],
    {ok, B} = file:read_file(P),
    json:decode(B).

sympy_test_() ->
    #{<<"cases">> := Cases} = fixture(),
    [{binary_to_list(Name), fun() -> check(C) end} || #{<<"name">> := Name} = C <- Cases].

check(#{<<"erl">> := Erl, <<"n">> := N, <<"at">> := At0, <<"params">> := Ps, <<"terms">> := Want}) ->
    {ok, Toks, _} = erl_scan:string(binary_to_list(Erl) ++ "."),
    {ok, E} = erl_parse:parse_term(Toks),
    At = case At0 of
             0 -> 0;
             <<"infinity">> -> infinity;
             S -> binary_to_integer(S)
         end,
    Bind = maps:from_list([{binary_to_atom(K), P / Q} || K := [P, Q] <- Ps]),
    {ok, {Got, _Order}} = pq_series:terms(E, x, N, #{at => At}),
    GotMap = maps:from_list([{{expo(Ex), M}, num(C, Bind)} || {Ex, M, C} <- Got]),
    % pq_series returns the terms of exponent below N; SymPy also keeps x^N log(x)^m (not O(x^N)): compared below
    Below = fun(P, Q) -> case At of infinity -> -P / Q < N; _ -> P / Q < N end end,
    WantMap = maps:from_list([{{P / Q, M}, binary_to_float_or_int(Re)} || [P, Q, M, Re, _Im] <- Want, Below(P, Q)]),
    Keys = lists:usort(maps:keys(GotMap) ++ maps:keys(WantMap)),
    [begin
         G = maps:get(K, GotMap, 0.0),
         W = maps:get(K, WantMap, 0.0),
         ?assert(abs(G - W) =< 1.0e-10 * max(1.0, abs(W)), {K, got, G, want, W})
     end || K <- Keys].

expo(K) when is_integer(K) -> float(K);
expo({rat, P, Q}) -> P / Q.

binary_to_float_or_int(B) ->
    S = binary_to_list(B),
    try list_to_float(S) catch error:badarg -> float(list_to_integer(S)) end.

%% a coefficient in floating point
num(K, _) when is_integer(K) -> float(K);
num({rat, P, Q}, _) -> P / Q;
num({sym, X}, B) -> maps:get(X, B);
num({add, L}, B) -> lists:sum([num(T, B) || T <- L]);
num({mul, L}, B) -> lists:foldl(fun(T, Acc) -> Acc * num(T, B) end, 1.0, L);
num({pow, X, K}, B) when is_integer(K) -> math:pow(num(X, B), K);
num({pow, X, Y}, B) -> math:pow(num(X, B), num(Y, B));
num({log, X}, B) -> math:log(num(X, B));
num({exp, X}, B) -> math:exp(num(X, B));
num({sin, X}, B) -> math:sin(num(X, B));
num({cos, X}, B) -> math:cos(num(X, B));
num({const, pi}, _) -> math:pi();
num({const, euler_gamma}, _) -> 0.57721566490153286061;
num({zeta, K}, _) -> zeta(K);
num({gamma, X}, B) -> gamma(num(X, B)).

%% zeta(k), k >= 2: the first 50 terms and the Euler-Maclaurin tail
zeta(K) ->
    N = 50,
    lists:sum([math:pow(I, -K) || I <- lists:seq(1, N - 1)])
        + math:pow(N, 1 - K) / (K - 1) + 0.5 * math:pow(N, -K) + K * math:pow(N, -K - 1) / 12
        - K * (K + 1) * (K + 2) * math:pow(N, -K - 3) / 720.

%% Gamma by Lanczos (g = 7, n = 9), with the reflection formula below 1/2
gamma(Z) when Z < 0.5 -> math:pi() / (math:sin(math:pi() * Z) * gamma(1 - Z));
gamma(Z0) ->
    Z = Z0 - 1,
    Cs = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059,
          12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7],
    X = lists:foldl(fun({I, C}, Acc) -> Acc + C / (Z + I) end, hd(Cs), lists:zip(lists:seq(1, 8), tl(Cs))),
    T = Z + 7.5,
    math:sqrt(2 * math:pi()) * math:pow(T, Z + 0.5) * math:exp(-T) * X.

api_test() ->
    X = {sym, x},
    % the expression form, and a refusal
    {ok, S} = pq_series:series({exp, X}, x, 3, #{}),
    ?assertEqual({add, [{mul, [1]}, {mul, [1, {pow, X, 1}]}, {mul, [{rat, 1, 2}, {pow, X, 2}]}]}, S),
    ?assertMatch({error, {symbolic_unsupported, #{series := _}}}, pq_series:series({exp, {pow, X, -1}}, x, 3, #{})).

limit_test() ->
    X = {sym, x},
    % sin(x)/x -> 1; (1 - cos x)/x^2 -> 1/2; x log x -> 0; 1/x diverges; x/(x+1) -> 1 at infinity
    ?assertEqual({ok, 1}, pq_symbolic:limit({mul, [{sin, X}, {pow, X, -1}]}, x, 0)),
    ?assertEqual({ok, {rat, 1, 2}}, pq_symbolic:limit({mul, [{add, [1, {mul, [-1, {cos, X}]}]}, {pow, X, -2}]}, x, 0)),
    ?assertEqual({ok, 0}, pq_symbolic:limit({mul, [X, {log, X}]}, x, 0)),
    ?assertMatch({ok, {diverges, _}}, pq_symbolic:limit({pow, X, -1}, x, 0)),
    ?assertEqual({ok, 1}, pq_symbolic:limit({mul, [X, {pow, {add, [X, 1]}, -1}]}, x, infinity)),
    % (x^2 - 4)/(x - 2) -> 4 at x = 2
    ?assertEqual({ok, 4}, pq_symbolic:limit({mul, [{add, [{pow, X, 2}, -4]}, {pow, {add, [X, -2]}, -1}]}, x, 2)),
    ?assertMatch({ok, {add, _}}, pq_symbolic:series({cos, X}, x, 5)).
