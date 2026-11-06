%% Reconstruction from values over prime fields (decision record R1): polynomials (phases 2, 3) and rational
%% functions (phase 4), by Zippel's sparse interpolation -- one engine for one output or many -- then the rationals
%% by Chinese remaindering and rational reconstruction; the work as tasks for an executor (pq_rexec: this node, or
%% Erlang nodes), the program and supports in the store (pq_rstore).
%%
%% The black box is a compiled program (pq_symbolic:compile/3, exact_constants), evaluated over prime fields by the
%% tasks. Every random choice comes from {Seed, prime, stage, attempt}, and the columns are taken in order however
%% many run at once, so the result is the same for every executor, node count and batch (D1). Probabilistic as every
%% sparse interpolation is: degrees are found by early termination (Terminate further points that add nothing), and
%% the result is checked at a point no stage used, in a prime no stage used.
-module(pq_recon).
-export([polynomial/3, rational/3]).

-define(TERMINATE, 2).
-define(MAX_DEGREE, 1000).
-define(CHECK_PRIME, 65).

%%% Polynomials

%% The polynomial in Vars (the program's inputs, atoms) that the program computes (its first output), as
%% pq_sym_poly's map. Opts: seed (1); executor => local (the default) | {nodes, Nodes} (pq_rexec), and its workers,
%% per_node, slots, task_timeout; batch, the columns computed ahead at a time (default: the executor's slots).
-spec polynomial(tuple(), [atom()], map()) -> {ok, map()} | {error, term()}.
polynomial({slp, 1, Inputs, _, _}, Vars, _Opts) when Inputs =/= Vars -> {error, {bad_shape, {variables, Inputs}}};
polynomial(Program, Vars, Opts) ->
    with_job(Program, length(Vars), Opts,
             fun(Job) ->
                     Key = maps:get(key, Job),
                     Poly = Job#{outs => 1,
                                 task => fun(P, R, N, SKey, Suffix) -> {column, Key, P, R, N, SKey, Suffix} end},
                     [Coeffs] = solve(Poly),
                     {Support, Rat} = Coeffs,
                     ok = check_poly(Poly, Support, Rat),
                     to_poly(Vars, Support, Rat)
             end).

%%% Rational functions

%% The rational function in Vars that the program computes (its first output), as {Numerator, Denominator}:
%% pq_sym_poly maps, coprime, the denominator's coefficient of its last monomial (in Erlang's order of exponent
%% lists) 1. Along the lines x = s + t (1, z_2, .., z_n), for a random shift s, the function is n(t) / d(t) with
%% d(0) = 1; the coefficient of t^r in n and in d is a homogeneous polynomial of degree r in z, found by Zippel in
%% z_2 .. z_n (z_1 = 1), all of them from the same points. The shift is then taken out degree by degree from the top
%% (N_k = n_k - sum_{j > k} (s . grad)^(j - k) N_j / (j - k)!). Opts as for polynomial/3.
-spec rational(tuple(), [atom()], map()) -> {ok, {map(), map()}} | {error, term()}.
rational({slp, 1, Inputs, _, _}, Vars, _Opts) when Inputs =/= Vars -> {error, {bad_shape, {variables, Inputs}}};
rational(_Program, [], _Opts) -> {error, {bad_shape, no_variables}};
rational(Program, Vars, Opts) ->
    N = length(Vars),
    with_job(Program, N - 1, Opts,
             fun(#{key := Key, seed := Seed} = Job) ->
                     Shift = [1 + V rem (1 bsl 20) || V <- rnd_list(1 bsl 62, N, {Seed, shift})],
                     P1 = pq_ff:prime(1),
                     {Dn, Dd} = degrees(Job, P1, Shift),
                     Task = fun(P, R, NM, SKey, Suffix) ->
                                    {line_column, Key, P, [S rem P || S <- Shift], R, NM, SKey, Suffix, ts(Job, P, Dn + Dd + 1), Dn, Dd}
                            end,
                     Rat = Job#{outs => Dn + Dd + 1, task => Task},
                     {Support, Coeffs} = case solve(Rat) of [Sc] -> Sc end,
                     {Num, Den} = unshift(Vars, Shift, Support, Coeffs, Dn, Dd),
                     ok = check_rat(Rat, Vars, Num, Den),
                     {to_map(Vars, Num), to_map(Vars, Den)}
             end).

ts(#{seed := Seed}, P, Count) ->
    distinct(P, Count, {Seed, P, ts}).

distinct(P, Count, Key) ->
    L = lists:usort(rnd_list(P, Count + 8, Key)),
    case length(L) >= Count of
        true -> lists:sublist(L, Count);
        false -> distinct(P, Count, {Key, again})
    end.

%%% The job

with_job(Program, N, Opts, Body) ->
    case pq_ff:prepare(Program) of
        {ok, _} ->
            Exec = maps:get(executor, Opts, local),
            Job = #{key => pq_rstore:put(Program), run => fun(Tasks) -> pq_rexec:map(Tasks, Exec, Opts) end,
                    batch => maps:get(batch, Opts, batch(Exec, Opts)), seed => maps:get(seed, Opts, 1), n => N},
            try {ok, Body(Job)}
            catch throw:{recon, Why} -> {error, {reconstruction, Why}}
            end;
        Error -> Error
    end.

batch(local, Opts) -> maps:get(workers, Opts, erlang:system_info(dirty_cpu_schedulers_online));
batch({nodes, Nodes}, Opts) -> max(1, length(Nodes)) * maps:get(per_node, Opts, erlang:system_info(dirty_cpu_schedulers_online)).

%% {Support, [Coefficients per output]} over the rationals: Zippel over the first prime, the support in further
%% primes until the reconstruction is stable
solve(Job) ->
    P1 = pq_ff:prime(1),
    {Support, C1} = zippel(Job, P1),
    [{Support, rationals(Job, Support, [{P1, C1}], 2, none)}].

%%% Zippel over one prime: {Support, [Coefficients per output]}, the support a list of exponent lists

zippel(#{n := 0} = Job, P) -> {[[]], with_support(Job, P, [[]])};
zippel(#{n := N, seed := Seed, outs := Outs} = Job, P) ->
    Anchors = rnd_list(P, N, {Seed, P, anchors}),
    %% stage 0: no variable interpolated yet, the support is the empty monomial
    lists:foldl(fun(K, {S, _C}) -> stage(Job, P, K, S, Anchors) end, {[[]], lists:duplicate(Outs, none)}, lists:seq(1, N)).

%% Stage K: from the support S in x_1 .. x_{K-1}, the support in x_1 .. x_K. For x_K = b_0, b_1, ...: the
%% coefficients of S at that b in every output (one column: transposed Vandermonde systems over the points r^1,
%% r^2, ... in x_1 .. x_{K-1}), then each interpolated in x_K by Newton until Terminate further b add nothing to any
%% of them. The columns are computed Batch at a time, in parallel, and taken in order.
stage(#{seed := Seed, outs := Outs} = Job, P, K, S, Anchors) ->
    Rest = lists:nthtail(K, Anchors),
    {R, _Nodes} = nodes(P, K - 1, S, {Seed, P, K}),
    SKey = pq_rstore:put(lists:append(S)),
    Newton0 = [[{[], []} || _ <- S] || _ <- lists:seq(1, Outs)],   % per output, per monomial: {Xs, Newton coefficients}
    columns(Job, P, K, S, {R, SKey}, Rest, Newton0, 0, 0).

columns(#{run := Run, batch := Batch, seed := Seed, task := Task} = Job, P, K, S, {R, SKey}, Rest, Newton, J, Quiet) ->
    J > ?MAX_DEGREE andalso throw({recon, degree_bound}),
    Bs = [point_value(P, {Seed, P, K, b, I}) || I <- lists:seq(J, J + Batch - 1)],
    Results = Run([Task(P, R, length(S), SKey, [B | Rest]) || B <- Bs]),
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
take(P, [{B, Result} | Rest], Newton, Quiet) ->
    Outputs = outputs(Result),
    Newton2 = [[newton_add(P, Nw, B, C) || {Nw, C} <- lists:zip(NwO, Cs)] || {NwO, Cs} <- lists:zip(Newton, Outputs)],
    AllZero = lists:all(fun(NwO) -> lists:all(fun({_, Cs}) -> lists:last(Cs) =:= 0 end, NwO) end, Newton2),
    take(P, Rest, Newton2, case AllZero of true -> Quiet + 1; false -> 0 end).

%% a column task's result as coefficients per output
outputs({ok, [L | _] = Lists}) when is_list(L) -> Lists;
outputs({ok, []}) -> [[]];
outputs({ok, Cs}) -> [Cs];
outputs(Other) -> throw({recon, Other}).

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

%% the new support: each monomial of S times x_K^e wherever an output's interpolant has a non-zero coefficient,
%% with every output's coefficient there (0 where it has none)
support(P, S, Newton) ->
    Forms = [[monomial_form(P, Nw) || Nw <- NwO] || NwO <- Newton],      % per output, per monomial
    PerMono = [[lists:nth(Mi, FO) || FO <- Forms] || Mi <- lists:seq(1, length(S))],
    Rows = lists:append([[{M ++ [E], [coef(F, E) || F <- Fs]}
                          || E <- lists:usort([E || F <- Fs, {E, C} <- lists:enumerate(0, F), C =/= 0])]
                         || {M, Fs} <- lists:zip(S, PerMono)]),
    {[M || {M, _} <- Rows], transpose([Cs || {_, Cs} <- Rows], length(Newton))}.

coef(F, E) when E < length(F) -> lists:nth(E + 1, F);
coef(_, _) -> 0.

transpose([], Outs) -> lists:duplicate(Outs, []);
transpose(Rows, _) -> [[lists:nth(O, Row) || Row <- Rows] || O <- lists:seq(1, length(hd(Rows)))].

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
with_support(#{outs := Outs}, _P, []) -> lists:duplicate(Outs, []);
with_support(Job, P, Support) -> with_support(Job, P, Support, 0).
with_support(_Job, _P, _Support, 20) -> throw({recon, bad_points});
with_support(#{run := Run, seed := Seed, n := N, task := Task} = Job, P, Support, Try) ->
    {R, _} = nodes(P, N, Support, {Seed, P, support, Try}),
    SKey = pq_rstore:put(lists:append(Support)),
    case Run([Task(P, R, length(Support), SKey, [])]) of
        [bad] -> with_support(Job, P, Support, Try + 1);
        [singular] -> with_support(Job, P, Support, Try + 1);
        [Result] -> outputs(Result)
    end.

%% every coefficient of every output from its residues: CRT, then Wang; false while one has no reconstruction
reconstruct(Residues) ->
    [{_, C0} | _] = Residues,
    Flat = [{P, lists:append(Cs)} || {P, Cs} <- Residues],
    [{_, F0} | _] = Flat,
    Combined = [lists:foldl(fun({P, Cs}, Acc) ->
                                    Ci = lists:nth(I, Cs),
                                    case Acc of none -> {Ci, P}; _ -> pq_ff:crt(Acc, {Ci, P}) end
                            end, none, Flat) || I <- lists:seq(1, length(F0))],
    Rats = [pq_ff:ratrecon(A, M) || {A, M} <- Combined],
    case lists:member(false, Rats) of
        true -> false;
        false -> split(Rats, [length(C) || C <- C0])
    end.

split([], []) -> [];
split(L, [N | Ns]) -> {A, B} = lists:split(N, L), [A | split(B, Ns)].

%%% Degrees along a line, by Thiele's continued fraction

%% the degrees of numerator and denominator in t along a random line s + t (1, z): the continued fraction grows
%% until Terminate further points are predicted, then its last convergent gives them
degrees(#{key := Key, run := Run, seed := Seed, n := M}, P, Shift) ->
    Z = rnd_list(P, M, {Seed, P, line}),
    Ts = rnd_list(P, 2 * ?MAX_DEGREE + 2 * ?TERMINATE, {Seed, P, line_ts}),
    Point = fun(T) -> [(S + T * Zi) rem P || {S, Zi} <- lists:zip(Shift, [1 | Z])] end,
    thiele(Run, Key, P, Point, Ts, {[], []}, 0).

thiele(_Run, _Key, _P, _Point, [], _State, _Quiet) -> throw({recon, degree_bound});
thiele(_Run, _Key, P, _Point, _Ts, {Xs, As}, Quiet) when Quiet >= ?TERMINATE -> convergent_degrees(P, Xs, As);
thiele(Run, Key, P, Point, Ts, State, Quiet) ->
    {Batch, Rest} = lists:split(min(16, length(Ts)), Ts),
    [{ok, Rows}] = Run([{values, Key, P, [Point(T) || T <- Batch]}]),
    {State2, Quiet2} = lists:foldl(fun({_T, bad}, Acc) -> Acc;
                                      (_, {S, Q}) when Q >= ?TERMINATE -> {S, Q};
                                      ({T, [Y | _]}, {S, Q}) ->
                                          case predict(P, S, T) of
                                              Y -> {S, Q + 1};
                                              _ -> {thiele_add(P, S, T, Y), 0}
                                          end
                                   end, {State, Quiet}, lists:zip(Batch, Rows)),
    thiele(Run, Key, P, Point, Rest, State2, Quiet2).

%% the continued fraction a0 + (t - x0) / (a1 + (t - x1) / (a2 + ...)) at t; none where it isn't defined
predict(_P, {[], []}, _T) -> none;
predict(P, {Xs, As}, T) ->
    [Last | RevA] = lists:reverse(As),
    RevX = tl(lists:reverse(Xs)),
    lists:foldl(fun(_, none) -> none;
                   ({Xi, Ai}, F) -> case F of 0 -> none; _ -> (Ai + mod(T - Xi, P) * pq_ff:inv(F, P)) rem P end
                end, Last, lists:zip(RevX, RevA)).

%% the point added: its reciprocal differences; a point the fraction cannot reach is skipped
thiele_add(P, {Xs, As} = S, T, Y) ->
    R = lists:foldl(fun(_, skip) -> skip;
                       ({Xi, Ai}, Ri) -> case mod(Ri - Ai, P) of 0 -> skip; D -> mod(T - Xi, P) * pq_ff:inv(D, P) rem P end
                    end, Y, lists:zip(Xs, As)),
    case R of skip -> S; _ -> {Xs ++ [T], As ++ [R]} end.

%% the degrees of the last convergent P_k / Q_k: P_k = a_k P_{k-1} + (t - x_{k-1}) P_{k-2}, likewise Q
convergent_degrees(P, Xs, As) ->
    {Pk, Qk, _, _} = lists:foldl(fun({I, A}, {P1, Q1, P2, Q2}) ->
                                         case I of
                                             0 -> {[A], [1], [1], [0]};
                                             _ ->
                                                 X = lists:nth(I, Xs),
                                                 Lin = [mod(-X, P), 1],
                                                 {padd(P, pscale(P, P1, A), pmul(P, Lin, P2)), padd(P, pscale(P, Q1, A), pmul(P, Lin, Q2)), P1, Q1}
                                         end
                                 end, {[0], [1], [1], [0]}, lists:enumerate(0, As)),
    {pdeg(Pk), pdeg(Qk)}.

padd(P, A, B) -> [(X + Y) rem P || {X, Y} <- lists:zip(pad(A, length(B)), pad(B, length(A)))].
pscale(P, A, C) -> [X * C rem P || X <- A].
pmul(P, A, B) -> lists:foldl(fun({I, X}, Acc) -> padd(P, Acc, lists:duplicate(I, 0) ++ [X * Y rem P || Y <- B]) end, [0], lists:enumerate(0, A)).
pad(L, N) -> L ++ lists:duplicate(max(0, N - length(L)), 0).
pdeg(L) -> length(lists:dropwhile(fun(C) -> C =:= 0 end, lists:reverse(L))) - 1.

%%% From the coefficients along the lines to N / D

%% n_r, d_r (polynomials in z_2 .. z_n with z_1 = 1) made homogeneous of degree r again, then the shift taken out
%% from the top degree down; both normalised so that D's coefficient of its last monomial is 1
unshift(Vars, Shift, Support, Coeffs, Dn, _Dd) ->
    {NCs, DCs} = lists:split(Dn + 1, Coeffs),
    Homog = fun(R, Cs) ->
                    maps:from_list([{homog(R, M), C} || {M, C} <- lists:zip(Support, Cs), C =/= {0, 1}])
            end,
    NComps = [Homog(R, Cs) || {R, Cs} <- lists:enumerate(0, NCs)],
    DComps = [#{lists:duplicate(length(Vars), 0) => {1, 1}} | [Homog(R, Cs) || {R, Cs} <- lists:enumerate(1, DCs)]],
    Num = untaylor(Shift, NComps),
    Den = untaylor(Shift, DComps),
    {_, Lead} = lists:max(maps:to_list(Den)),
    Inv = pq_sym_exact:rinv(Lead),
    {pscale_q(Num, Inv), pscale_q(Den, Inv)}.

homog(R, M) ->
    E1 = R - lists:sum(M),
    E1 >= 0 orelse throw({recon, degree_mismatch}),
    [E1 | M].

%% N = sum_k N_k with N_k = n_k - sum_{j > k} (s . grad)^(j - k) N_j / (j - k)!, from the top down; Comps are the
%% homogeneous parts n_0 .. n_R of N(x + s) (exponent list => rational)
untaylor(Shift, Comps) ->
    R = length(Comps) - 1,
    {Parts, _} = lists:foldl(fun(K, {Acc, Stored}) ->
                                     %% D^(j-k) N_j / (j-k)! from D^(j-k-1) N_j / (j-k-1)!
                                     Stored2 = [{J, pscale_q(sgrad(Shift, T), {1, J - K})} || {J, T} <- Stored],
                                     Nk = lists:foldl(fun({_, T}, A) -> psub_q(A, T) end, lists:nth(K + 1, Comps), Stored2),
                                     {[Nk | Acc], [{K, Nk} | Stored2]}
                             end, {[], []}, lists:seq(R, 0, -1)),
    lists:foldl(fun padd_q/2, #{}, Parts).

%% (s . grad) P = sum_i s_i dP/dx_i
sgrad(Shift, P) ->
    maps:fold(fun(E, C, Acc) ->
                      lists:foldl(fun({I, {Ei, Si}}, A) when Ei > 0 ->
                                          E2 = setnth(I, E, Ei - 1),
                                          padd_q(A, #{E2 => pq_sym_exact:rmul(C, {Ei * Si, 1})});
                                     (_, A) -> A
                                  end, Acc, lists:enumerate(lists:zip(E, Shift)))
              end, #{}, P).

setnth(I, L, V) -> {A, [_ | B]} = lists:split(I - 1, L), A ++ [V | B].

padd_q(A, B) ->
    maps:fold(fun(E, C, Acc) ->
                      S = pq_sym_exact:radd(maps:get(E, Acc, {0, 1}), C),
                      case S of {0, _} -> maps:remove(E, Acc); _ -> Acc#{E => S} end
              end, A, B).
psub_q(A, B) -> padd_q(A, maps:map(fun(_, C) -> pq_sym_exact:rneg(C) end, B)).
pscale_q(A, F) -> maps:map(fun(_, C) -> pq_sym_exact:rmul(C, F) end, A).

%%% Checks at a fresh point, in a fresh prime

check_poly(_Job, [], _Rat) -> ok;
check_poly(#{key := Key, run := Run, seed := Seed, n := N}, Support, [Rat]) ->
    P = pq_ff:prime(?CHECK_PRIME),
    X = rnd_list(P, N, {Seed, check}),
    case Run([{values, Key, P, [X]}]) of
        [{ok, [bad]}] -> ok;
        [{ok, [[V | _]]}] ->
            case V =:= qeval(P, X, lists:zip(Support, Rat)) of true -> ok; false -> throw({recon, check_failed}) end;
        Other -> throw({recon, {check, Other}})
    end.

check_rat(#{key := Key, run := Run, seed := Seed}, Vars, Num, Den) ->
    P = pq_ff:prime(?CHECK_PRIME),
    X = rnd_list(P, length(Vars), {Seed, check}),
    case Run([{values, Key, P, [X]}]) of
        [{ok, [bad]}] -> ok;
        [{ok, [[V | _]]}] ->
            Dv = qeval(P, X, maps:to_list(Den)),
            Dv =/= 0 andalso (V =/= qeval(P, X, maps:to_list(Num)) * pq_ff:inv(Dv, P) rem P)
                andalso throw({recon, check_failed}),
            ok;
        Other -> throw({recon, {check, Other}})
    end.

qeval(P, X, Terms) ->
    lists:foldl(fun({M, {Nu, De}}, A) -> (A + mod(Nu, P) * pq_ff:inv(De, P) rem P * mono_value(P, X, M)) rem P end, 0, Terms).

to_poly(Vars, Support, [Rat]) -> to_map(Vars, maps:from_list(lists:zip(Support, Rat))).

to_map(Vars, Terms) ->
    maps:from_list([{{[{{sym, X}, {E, 1}} || {X, E} <- lists:zip(Vars, M), E =/= 0], []}, {q, R, {0, 1}}}
                    || {M, R} <- maps:to_list(Terms)]).

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
