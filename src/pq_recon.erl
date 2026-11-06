%% Reconstruction of a polynomial from its values over prime fields (decision record R1, phases 2 and 3): Zippel's
%% sparse interpolation, then the rationals by Chinese remaindering and rational reconstruction; the work as tasks
%% for an executor (pq_rexec: this node, or Erlang nodes), the program and supports in the store (pq_rstore).
%%
%% The black box is a compiled program (pq_symbolic:compile/3, exact_constants), evaluated over prime fields by the
%% tasks. Every random choice comes from {Seed, prime, stage, attempt}, and the columns are taken in order however
%% many run at once, so the result is the same for every executor, node count and batch (D1). Probabilistic as every
%% sparse interpolation is: degrees are found by early termination (Terminate further points that add nothing), and
%% the result is checked at a point no stage used, in a prime no stage used.
-module(pq_recon).
-export([polynomial/3]).

-define(TERMINATE, 2).
-define(MAX_DEGREE, 1000).

%% The polynomial in Vars (the program's inputs, atoms) that the program computes (its first output), as
%% pq_sym_poly's map. Opts: seed (1); executor => local (the default) | {nodes, Nodes} (pq_rexec), and its workers,
%% per_node, task_timeout; batch, the columns computed ahead at a time (default: the executor's slots).
-spec polynomial(tuple(), [atom()], map()) -> {ok, map()} | {error, term()}.
polynomial({slp, 1, Inputs, _, _}, Vars, _Opts) when Inputs =/= Vars -> {error, {bad_shape, {variables, Inputs}}};
polynomial(Program, Vars, Opts) ->
    case pq_ff:prepare(Program) of
        {ok, _} ->
            Seed = maps:get(seed, Opts, 1),
            Exec = maps:get(executor, Opts, local),
            Run = fun(Tasks) -> pq_rexec:map(Tasks, Exec, Opts) end,
            Batch = maps:get(batch, Opts, batch(Exec, Opts)),
            Job = #{key => pq_rstore:put(Program), run => Run, batch => Batch, seed => Seed, n => length(Vars)},
            try
                P1 = pq_ff:prime(1),
                {Support, C1} = zippel(Job, P1),
                Rat = rationals(Job, Support, [{P1, C1}], 2, none),
                ok = check(Job, Support, Rat),
                {ok, to_poly(Vars, Support, Rat)}
            catch
                throw:{recon, Why} -> {error, {reconstruction, Why}}
            end;
        Error -> Error
    end.

batch(local, Opts) -> maps:get(workers, Opts, erlang:system_info(dirty_cpu_schedulers_online));
batch({nodes, Nodes}, Opts) -> max(1, length(Nodes)) * maps:get(per_node, Opts, erlang:system_info(dirty_cpu_schedulers_online)).

%%% Zippel over one prime: {Support, Coefficients}, the support a list of exponent lists (one per variable)

zippel(#{n := N, seed := Seed} = Job, P) ->
    Anchors = rnd_list(P, N, {Seed, P, anchors}),
    %% stage 0: no variable interpolated yet, the support is the empty monomial
    lists:foldl(fun(K, {S, _C}) -> stage(Job, P, K, S, Anchors) end, {[[]], [none]}, lists:seq(1, N)).

%% Stage K: from the support S in x_1 .. x_{K-1}, the support in x_1 .. x_K. For x_K = b_0, b_1, ...: the
%% coefficients of S at that b (one column: a transposed Vandermonde system over the points r^1, r^2, ... in
%% x_1 .. x_{K-1}), then each coefficient interpolated in x_K by Newton until Terminate further b add nothing to
%% any of them. The columns are computed Batch at a time, in parallel, and taken in order.
stage(#{seed := Seed} = Job, P, K, S, Anchors) ->
    Rest = lists:nthtail(K, Anchors),
    {R, _Nodes} = nodes(P, K - 1, S, {Seed, P, K}),
    SKey = pq_rstore:put(lists:append(S)),
    Newton0 = [{[], []} || _ <- S],                    % per monomial: {Xs, Newton coefficients}
    columns(Job, P, K, S, {R, SKey}, Rest, Newton0, 0, 0).

columns(#{key := Key, run := Run, batch := Batch, seed := Seed} = Job, P, K, S, {R, SKey}, Rest, Newton, J, Quiet) ->
    J > ?MAX_DEGREE andalso throw({recon, degree_bound}),
    Bs = [point_value(P, {Seed, P, K, b, I}) || I <- lists:seq(J, J + Batch - 1)],
    Results = Run([{column, Key, P, R, length(S), SKey, [B | Rest]} || B <- Bs]),
    case take(P, lists:zip(Bs, Results), Newton, Quiet) of
        {done, Newton2} -> support(P, S, Newton2);
        {more, Newton2, Quiet2} -> columns(Job, P, K, S, {R, SKey}, Rest, Newton2, J + Batch, Quiet2)
    end.

%% the columns in order: each adds a point to every Newton interpolant, until Terminate in a row add nothing
take(_P, _, Newton, Quiet) when Quiet >= ?TERMINATE -> {done, Newton};
take(_P, [], Newton, Quiet) -> {more, Newton, Quiet};
take(P, [{_B, bad} | Rest], Newton, Quiet) -> take(P, Rest, Newton, Quiet);    % a b where the program divides by 0
take(_P, [{_B, singular} | _], _Newton, _Quiet) -> throw({recon, singular});
take(_P, [{_B, bad_prime} | _], _Newton, _Quiet) -> throw({recon, bad_prime});
take(P, [{B, {ok, Coeffs}} | Rest], Newton, Quiet) ->
    Newton2 = [newton_add(P, Nw, B, C) || {Nw, C} <- lists:zip(Newton, Coeffs)],
    AllZero = lists:all(fun({_, Cs}) -> lists:last(Cs) =:= 0 end, Newton2),
    take(P, Rest, Newton2, case AllZero of true -> Quiet + 1; false -> 0 end).

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
support(P, S, Newton) ->
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

rationals(_Job, _Support, _Residues, K, _Prev) when K > 64 -> throw({recon, too_many_primes});
rationals(Job, Support, Residues, K, Prev) ->
    Rat = reconstruct(Residues),
    case Rat =/= false andalso Rat =:= Prev of
        true -> Rat;
        false ->
            P = pq_ff:prime(K),
            C = with_support(Job, P, Support),
            rationals(Job, Support, [{P, C} | Residues], K + 1, Rat)
    end.

%% the coefficients of a known support in the prime P: one column, over every variable
with_support(_Job, _P, []) -> [];
with_support(Job, P, Support) -> with_support(Job, P, Support, 0).
with_support(_Job, _P, _Support, 20) -> throw({recon, bad_points});
with_support(#{key := Key, run := Run, seed := Seed, n := N} = Job, P, Support, Try) ->
    {R, _} = nodes(P, N, Support, {Seed, P, support, Try}),
    SKey = pq_rstore:put(lists:append(Support)),
    case Run([{column, Key, P, R, length(Support), SKey, []}]) of
        [{ok, Cs}] -> Cs;
        [bad] -> with_support(Job, P, Support, Try + 1);
        [singular] -> with_support(Job, P, Support, Try + 1);
        [Other] -> throw({recon, Other})
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
check(_Job, [], _Rat) -> ok;
check(#{key := Key, run := Run, seed := Seed, n := N}, Support, Rat) ->
    P = pq_ff:prime(65),
    X = rnd_list(P, N, {Seed, check}),
    case Run([{values, Key, P, [X]}]) of
        [{ok, [bad]}] -> ok;
        [{ok, [[V | _]]}] ->
            Want = lists:foldl(fun({M, {Nu, De}}, A) -> (A + mod(Nu, P) * pq_ff:inv(De, P) rem P * mono_value(P, X, M)) rem P end,
                               0, lists:zip(Support, Rat)),
            case V =:= Want of true -> ok; false -> throw({recon, check_failed}) end;
        Other -> throw({recon, {check, Other}})
    end.

to_poly(Vars, Support, Rat) ->
    maps:from_list([{{[{{sym, X}, {E, 1}} || {X, E} <- lists:zip(Vars, M), E =/= 0], []}, {q, R, {0, 1}}}
                    || {M, R} <- lists:zip(Support, Rat)]).

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
