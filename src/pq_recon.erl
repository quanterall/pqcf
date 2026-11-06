%% Reconstruction of a polynomial from its values over prime fields (decision record R1, phase 2): Zippel's sparse
%% interpolation, then the rationals by Chinese remaindering and rational reconstruction.
%%
%% The black box is an evaluator Eval(P, Points) -> [Value | bad] (Points: lists of residues, one per variable, in
%% Vars' order): pq_ff over a compiled program on this node, or any executor that distributes the points (phase 3).
%% Every random choice comes from {Seed, prime, stage, attempt}, so the result is the same however the points are
%% evaluated (D1). Probabilistic as every sparse interpolation is: degrees are found by early termination (Terminate
%% further points that add nothing), and the result is checked at a point no stage used, in a prime no stage used.
-module(pq_recon).
-export([polynomial/3, polynomial/4]).

-define(TERMINATE, 2).
-define(MAX_DEGREE, 1000).

%% The polynomial in Vars (atoms) that the program computes (its first output), as pq_sym_poly's map.
-spec polynomial(tuple(), [atom()], map()) -> {ok, map()} | {error, term()}.
polynomial({slp, 1, Inputs, _, _}, Vars, _Opts) when Inputs =/= Vars -> {error, {bad_shape, {variables, Inputs}}};
polynomial(Program, Vars, Opts) ->
    case pq_ff:prepare(Program) of
        {ok, FF} ->
            Eval = fun(P, Points) ->
                           case pq_ff:eval(FF, P, Points) of
                               {ok, Rows} -> [case R of bad -> bad; [V | _] -> V end || R <- Rows];
                               bad_prime -> [bad || _ <- Points]
                           end
                   end,
            polynomial(Eval, Vars, Opts, ok);
        Error -> Error
    end.

-spec polynomial(fun(), [atom()], map(), ok) -> {ok, map()} | {error, term()}.
polynomial(Eval, Vars, Opts, ok) ->
    Seed = maps:get(seed, Opts, 1),
    try
        P1 = pq_ff:prime(1),
        {Support, C1} = zippel(Eval, P1, length(Vars), Seed),
        Rat = rationals(Eval, Support, [{P1, C1}], 2, none, Seed),
        %% the check: a point no stage used, in a prime no stage used
        ok = check(Eval, Support, Rat, Seed),
        {ok, to_poly(Vars, Support, Rat)}
    catch
        throw:{recon, Why} -> {error, {reconstruction, Why}}
    end.

%%% Zippel over one prime: {Support, Coefficients}, the support a list of exponent lists (one per variable)

zippel(Eval, P, N, Seed) ->
    Anchors = rnd_list(P, N, {Seed, P, anchors}),
    %% stage 0: no variable interpolated yet, the support is the empty monomial
    lists:foldl(fun(K, {S, _C}) -> stage(Eval, P, N, K, S, Anchors, Seed) end, {[[]], [none]}, lists:seq(1, N)).

%% Stage K: from the support S in x_1 .. x_{K-1}, the support in x_1 .. x_K. For x_K = b_0, b_1, ...: the
%% coefficients of S at that b (a transposed Vandermonde system over the points r^0, r^1, ... in x_1 .. x_{K-1}),
%% then each coefficient interpolated in x_K by Newton until Terminate further b add nothing to any of them.
stage(Eval, P, _N, K, S, Anchors, Seed) ->
    Rest = lists:nthtail(K, Anchors),
    {R, Nodes} = nodes(P, K - 1, S, {Seed, P, K}),
    Newton0 = [{[], []} || _ <- S],                    % per monomial: {Xs, Newton coefficients}
    grow(Eval, P, K, S, R, Nodes, Rest, Newton0, 0, 0, Seed).

grow(_Eval, P, K, S, _R, _Nodes, _Rest, Newton, _J, Quiet, _Seed) when Quiet >= ?TERMINATE ->
    support(P, S, Newton, K);
grow(_Eval, _P, _K, _S, _R, _Nodes, _Rest, _Newton, J, _Quiet, _Seed) when J > ?MAX_DEGREE + ?TERMINATE ->
    throw({recon, degree_bound});
grow(Eval, P, K, S, R, Nodes, Rest, Newton, J, Quiet, Seed) ->
    B = point_value(P, {Seed, P, K, b, J}),
    Points = [Ri ++ [B | Rest] || Ri <- powers(R, length(S), P)],
    Values = Eval(P, Points),
    case lists:member(bad, Values) of
        true -> grow(Eval, P, K, S, R, Nodes, Rest, Newton, J + 1, Quiet, Seed);   % another b
        false ->
            Coeffs = solve(P, Nodes, Values),
            Newton2 = [newton_add(P, Nw, B, C) || {Nw, C} <- lists:zip(Newton, Coeffs)],
            AllZero = lists:all(fun({_, Cs}) -> lists:last(Cs) =:= 0 end, Newton2),
            grow(Eval, P, K, S, R, Nodes, Rest, Newton2, J + 1, case AllZero of true -> Quiet + 1; false -> 0 end, Seed)
    end.

%% random r in x_1 .. x_{K-1} whose monomial values on the support are distinct and non-zero (nodes of the system)
nodes(P, K1, S, Key) -> nodes(P, K1, S, Key, 0).
nodes(_P, _K1, _S, _Key, 20) -> throw({recon, no_distinct_nodes});
nodes(P, K1, S, Key, Try) ->
    R = rnd_list(P, K1, {Key, Try}),
    Nodes = [mono_value(P, R, M) || M <- S],
    case length(lists:usort(Nodes)) =:= length(Nodes) andalso not lists:member(0, Nodes) of
        true -> {R, Nodes};
        false -> nodes(P, K1, S, Key, Try + 1)
    end.

mono_value(P, R, M) -> lists:foldl(fun({Ri, E}, A) -> A * powm(Ri, E, P) rem P end, 1, lists:zip(R, M)).

%% the new support: each monomial of S times x_K^e for every non-zero coefficient of its interpolant
support(P, S, Newton, _K) ->
    Pairs = lists:append([[{M ++ [E], C} || {E, C} <- lists:enumerate(0, monomial_form(P, Nw)), C =/= 0]
                          || {M, Nw} <- lists:zip(S, Newton)]),
    case Pairs of
        [] -> {[], []};
        _ -> lists:unzip(Pairs)
    end.

%%% Newton interpolation mod P

%% {Xs, Cs} with the point (X, Y) added: the new coefficient (Y - p(X)) / prod (X - x_i)
newton_add(P, {Xs, Cs}, X, Y) ->
    {Val, Prod} = lists:foldl(fun({Xi, Ci}, {V, Pr}) -> {(V + Ci * Pr) rem P, Pr * mod(X - Xi, P) rem P} end,
                              {0, 1}, lists:zip(Xs, Cs)),
    C = mod(Y - Val, P) * pq_ff:inv(Prod, P) rem P,
    {Xs ++ [X], Cs ++ [C]}.

%% the Newton form as monomial coefficients, from the constant term up, trailing zeros dropped
monomial_form(P, {Xs, Cs}) ->
    %% p = c0 + (x - x0)(c1 + (x - x1)(c2 + ...)): from the innermost out
    Coeffs = lists:foldr(fun({Xi, Ci}, Acc) ->
                                 Shifted = [0 | Acc],                          % x * Acc
                                 Scaled = [mod(-Xi * A, P) || A <- Acc] ++ [0], % -xi * Acc
                                 Sum = [(A + Bv) rem P || {A, Bv} <- lists:zip(Shifted, Scaled)],
                                 [(hd(Sum) + Ci) rem P | tl(Sum)]
                         end, [], lists:zip(Xs, Cs)),
    lists:reverse(lists:dropwhile(fun(C) -> C =:= 0 end, lists:reverse(Coeffs))).

%%% The rationals: the same support in further primes, until the reconstruction is stable

rationals(_Eval, _Support, _Residues, K, _Prev, _Seed) when K > 64 ->
    throw({recon, too_many_primes});
rationals(Eval, Support, Residues, K, Prev, Seed) ->
    Rat = reconstruct(Residues),
    case Rat =/= false andalso Rat =:= Prev of
        true -> Rat;
        false ->
            P = pq_ff:prime(K),
            C = with_support(Eval, P, Support, Seed),
            rationals(Eval, Support, [{P, C} | Residues], K + 1, Rat, Seed)
    end.

%% the coefficients of a known support in the prime P: one transposed Vandermonde system
with_support(_Eval, _P, [], _Seed) -> [];
with_support(Eval, P, Support, Seed) -> with_support(Eval, P, Support, Seed, 0).
with_support(_Eval, _P, _Support, _Seed, 20) -> throw({recon, bad_points});
with_support(Eval, P, Support, Seed, Try) ->
    N = length(hd(Support)),
    {R, Nodes} = nodes(P, N, Support, {Seed, P, support, Try}),
    Values = Eval(P, powers(R, length(Support), P)),
    case lists:member(bad, Values) of
        true -> with_support(Eval, P, Support, Seed, Try + 1);
        false ->
            solve(P, Nodes, Values)
    end.

%% every coefficient from its residues: CRT, then Wang; false while one of them has no reconstruction
reconstruct(Residues) ->
    [{_, C0} | _] = Residues,
    Combined = [lists:foldl(fun({P, Cs}, Acc) ->
                                    Ci = lists:nth(I, Cs),
                                    case Acc of none -> {Ci, P}; _ -> pq_ff:crt(Acc, {Ci, P}) end
                            end, none, Residues) || I <- lists:seq(1, length(C0))],
    Rats = [pq_ff:ratrecon(A, M) || {A, M} <- Combined],
    case lists:member(false, Rats) of true -> false; false -> Rats end.

%% the reconstruction at a fresh point, in a fresh prime
check(_Eval, [], _Rat, _Seed) -> ok;
check(Eval, Support, Rat, Seed) ->
    P = pq_ff:prime(65),
    X = rnd_list(P, length(hd(Support)), {Seed, check}),
    case Eval(P, [X]) of
        [bad] -> ok;
        [V] ->
            Want = lists:foldl(fun({M, {Nu, De}}, A) -> (A + mod(Nu, P) * pq_ff:inv(De, P) rem P * mono_value(P, X, M)) rem P end,
                               0, lists:zip(Support, Rat)),
            case V =:= Want of true -> ok; false -> throw({recon, check_failed}) end
    end.

to_poly(Vars, Support, Rat) ->
    maps:from_list([{{[{{sym, X}, {E, 1}} || {X, E} <- lists:zip(Vars, M), E =/= 0], []}, {q, R, {0, 1}}}
                    || {M, R} <- lists:zip(Support, Rat)]).

%% the points r^1, ..., r^Count (componentwise), each from the one before: not r^0 = (1, ..., 1), where a
%% fraction such as 1 / (x - y) is never defined
powers(R, Count, P) ->
    {Pts, _} = lists:mapfoldl(fun(_, Cur) -> {Cur, [C * Ri rem P || {C, Ri} <- lists:zip(Cur, R)]} end,
                              R, lists:seq(1, Count)),
    Pts.

%% the coefficients c from values_i = sum_j c_j v_j^i (i = 1 .. n): the system from i = 0 gives c_j v_j
solve(P, Nodes, Values) ->
    case pq_nif:ff_tvandermonde(P, Nodes, Values) of
        {ok, Cv} -> [C * pq_ff:inv(V, P) rem P || {C, V} <- lists:zip(Cv, Nodes)];
        singular -> throw({recon, singular})
    end.

%%% Deterministic randomness and residues

rnd_list(P, N, Key) ->
    S0 = rand:seed_s(exsss, erlang:phash2(Key, 1 bsl 32) + 1),
    {L, _} = lists:mapfoldl(fun(_, S) -> {V, S2} = rand:uniform_s(1 bsl 62, S), {V rem (P - 1) + 1, S2} end, S0, lists:seq(1, N)),
    L.

point_value(P, Key) -> hd(rnd_list(P, 1, Key)).

powm(_, 0, P) -> 1 rem P;
powm(B, E, P) when E band 1 =:= 1 -> B * powm(B * B rem P, E bsr 1, P) rem P;
powm(B, E, P) -> powm(B * B rem P, E bsr 1, P).

mod(A, M) -> ((A rem M) + M) rem M.
