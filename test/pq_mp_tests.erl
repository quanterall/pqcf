%% pq_mp (decision record SY2, phase 2) against mpmath 1.3 at 60 digits (test/oracle/gt_mp.py,
%% test/fixtures/gt_mp.json): every value and zero to 48 significant digits.
-module(pq_mp_tests).
-include_lib("eunit/include/eunit.hrl").

-define(P, 180).   % bits: 54 digits

fixture() ->
    Dir = filename:dirname(code:which(?MODULE)),
    Paths = [filename:join([Dir, "..", "test", "fixtures", "gt_mp.json"]),
             filename:join([Dir, "..", "..", "..", "..", "test", "fixtures", "gt_mp.json"]),
             "test/fixtures/gt_mp.json"],
    [P | _] = [P || P <- Paths, filelib:is_regular(P)],
    {ok, B} = file:read_file(P),
    json:decode(B).

term(S) ->
    {ok, Toks, _} = erl_scan:string(binary_to_list(S) ++ "."),
    {ok, T} = erl_parse:parse_term(Toks),
    T.

%% a decimal string "d.ddd[e±k]" as an exact fraction
decimal(B) ->
    S = binary_to_list(B),
    {Mant, Exp} = case string:split(S, "e") of
                      [M, E] -> {M, list_to_integer(E)};
                      [M] -> {M, 0}
                  end,
    {Sign, Digits} = case Mant of [$- | R] -> {-1, R}; R -> {1, R} end,
    {Int, Frac} = case string:split(Digits, ".") of [I, F] -> {I, F}; [I] -> {I, ""} end,
    N = Sign * list_to_integer(Int ++ Frac),
    E2 = Exp - length(Frac),
    case E2 >= 0 of
        true -> {N * pow(10, E2), 1};
        false -> {N, pow(10, -E2)}
    end.

pow(_, 0) -> 1;
pow(B, K) -> B * pow(B, K - 1).

close(Got, Want) ->
    W = pq_mp:from(decimal(Want), ?P + 40),
    G = pq_mp:from(Got, ?P + 40),
    Rel = abs(pq_mp:to_float(pq_mp:divide(pq_mp:sub(G, W, ?P + 40), W, ?P + 40))),
    Rel < 1.0e-48.

values_test_() ->
    #{<<"values">> := Vs} = fixture(),
    [{binary_to_list(E), fun() ->
                                 {ok, S} = pq_mp:evaluate(term(E), #{}, 52),
                                 ?assert(close(decimal(list_to_binary(S)), V), {E, got, S, want, V})
                         end} || #{<<"erl">> := E, <<"value">> := V} <- Vs].

zeros_test_() ->
    #{<<"zeros">> := Zs} = fixture(),
    [{binary_to_list(Nu) ++ "/" ++ integer_to_list(K),
      fun() ->
              Z = pq_mp:besselj_zero(term(Nu), K, ?P),
              ?assert(close(Z, V), {Nu, K, got, pq_mp:to_string(Z, 52), want, V})
      end} || #{<<"nu">> := Nu, <<"k">> := K, <<"value">> := V} <- Zs].

%% the coefficients pq_series returns are expressions pq_mp evaluates: Y_0's constant term, (2/pi)(gamma - log 2),
%% against SymPy's 30 digits (gt_series.json: -0.0738042951086872252743439279485)
series_coefficient_test() ->
    {ok, {[{0, 0, C} | _], _}} = pq_series:terms({bessely, 0, {sym, x}}, x, 2, #{}),
    {ok, S} = pq_mp:evaluate(C, #{}, 30),
    ?assertEqual("-7.38042951086872252743439279485e-2", S).
