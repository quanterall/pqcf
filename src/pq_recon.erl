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
-export([polynomial/3, rational/3, linear/4]).

-define(TERMINATE, 2).
%% points per line for guessing a denominator (Thiele: degrees of numerator and denominator up to about 45 together)
-define(GUESS_POINTS, 48).
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
    with_job(Program, length(Vars) - 1, Opts, fun(Job) -> rational_job(Job, Vars) end).

rational_job(#{seed := Seed} = Job, Vars) ->
    Shift = [1 + V rem (1 bsl 20) || V <- rnd_list(1 bsl 62, length(Vars), {Seed, shift})],
    P1 = pq_ff:prime(1),
    case degrees(Job, P1, Shift) of
        {-1, _} -> {#{}, #{{[], []} => {q, {1, 1}, {0, 1}}}};       % the zero function: 0 / 1
        {Dn, Dd} -> rational_job(Job, Vars, Shift, Dn, Dd)
    end.

rational_job(#{bb := BB} = Job, Vars, Shift, Dn, Dd) ->
    Task = fun(P, R, NM, SKey, Suffix) ->
                   {line_column, BB, P, [S rem P || S <- Shift], R, NM, SKey, Suffix, ts(Job, P, Dn + Dd + 1), Dn, Dd}
           end,
    Rat = Job#{outs => Dn + Dd + 1, task => Task},
    {Support, Coeffs} = case solve(Rat) of [Sc] -> Sc end,
    {Num, Den} = unshift(Vars, Shift, Support, Coeffs, Dn, Dd),
    ok = check_rat(Rat, Vars, Num, Den),
    {to_map(Vars, Num), to_map(Vars, Den)}.

%%% Linear systems

%% a row as [{Column, Entry}] by increasing column, the zero entries left out; bad if it doesn't fit N columns
sparse_row(Row, N) when is_list(Row), length(Row) =:= N ->
    [{C, E} || {C, E} <- lists:enumerate(0, Row), E =/= 0];
sparse_row(Row, N) when is_map(Row) ->
    case lists:all(fun(C) -> is_integer(C) andalso C >= 0 andalso C < N end, maps:keys(Row)) of
        true -> [{C, E} || {C, E} <- lists:sort(maps:to_list(Row)), E =/= 0];
        false -> bad
    end;
sparse_row(_, _) -> bad.

%% The solution of A u = B, A a list of rows -- each a list of N entries, or a map Column => entry (0-based, the
%% missing entries 0: the form for large sparse systems) -- B a list of right-hand-side entries -- or of rows of M
%% entries, M right-hand sides at once (one per master integral, say), each row a list or a map Rhs => entry with M
%% given as Opts' rhs -- all expressions in Vars: for every component in Opts'
%% components (default all, 0-based), the reduced rational function {N, D} of each right-hand side's solution, as
%% [{Component, {N, D}}] for one right-hand side, [{Component, [{Rhs, {N, D}}]}] for several. Each is reconstructed
%% from sparse eliminations mod p at the sample points (pq_nif:ff_linsys/7), never a symbolic solve; with several
%% right-hand sides, one elimination per point tests them all for 0 (at three points, in the first prime) and only
%% the others are reconstructed (R1 phase 5). Opts as for rational/3, plus the expression limits
%% (max_expr_nodes, ..) for compiling A and B.
%% Only the non-zero entries of A and B are compiled (a dense row's literal zeros are dropped); the elimination
%% carries the right-hand sides as sparse columns, so thousands of them cost little where rows hold a few.
-spec linear([[term()] | #{non_neg_integer() => term()}], [term()] | [[term()] | #{non_neg_integer() => term()}], [atom()], map()) ->
          {ok, list()} | {error, term()}.
linear(A, B, Vars, Opts) ->
    N = length(B),
    Multi = B =/= [] andalso (is_list(hd(B)) orelse is_map(hd(B))),
    Bs = case Multi of true -> B; false -> [[Bi] || Bi <- B] end,
    M = case Bs of [] -> 0; [R | _] when is_map(R) -> maps:get(rhs, Opts, 0); [R | _] -> length(R) end,
    Sparse = [sparse_row(Row, N) || Row <- A],
    SparseB = [sparse_row(Row, M) || Row <- Bs],
    case Vars =/= [] andalso length(A) =:= N andalso not lists:member(bad, Sparse ++ SparseB) andalso M > 0 of
        false -> {error, {bad_shape, linear_system}};
        true ->
            Pattern = {[[C || {C, _} <- Row] || Row <- Sparse], [[C || {C, _} <- Row] || Row <- SparseB]},
            %% a large system needs the expression limits raised: Opts' limit keys reach compile
            Limits = maps:with([max_poly_terms, max_expand_power, max_expr_nodes, sym_timeout_ms, sym_max_heap_words], Opts),
            Entries = [E || Row <- Sparse ++ SparseB, {_, E} <- Row],
            case pq_symbolic:compile(Entries, Vars, Limits#{exact_constants => true}) of
                {ok, Prog} ->
                    Wanted = maps:get(components, Opts, lists:seq(0, N - 1)),
                    Results = [{W, component(Prog, Vars, {N, Pattern}, M, W, Opts)} || W <- Wanted],
                    case [E || {_, {error, _} = E} <- Results] of
                        [] when Multi -> {ok, [{W, R} || {W, {ok, R}} <- Results]};
                        [] -> {ok, [{W, NumDen} || {W, {ok, [{0, NumDen}]}} <- Results]};
                        [E | _] -> E
                    end;
                Error -> Error
            end
    end.

%% one component for every right-hand side: the zero test at three points, then a reconstruction per non-zero one
component(Prog, Vars, NP, M, W, Opts) ->
    Zero = {#{}, #{{[], []} => {q, {1, 1}, {0, 1}}}},
    case zero_columns(Prog, Vars, NP, M, W, Opts) of
        {ok, Zeros} ->
            NonZero = [C || C <- lists:seq(0, M - 1), not lists:member(C, Zeros)],
            Guesses = case maps:get(guess, Opts, true) of
                          true -> guess_denominators(Prog, Vars, NP, M, W, NonZero, Opts);
                          false -> #{}
                      end,
            Cols = [case lists:member(C, Zeros) of
                        true -> {C, {ok, Zero}};
                        false -> {C, column(Prog, Vars, NP, M, W, C, maps:get(C, Guesses, none), Opts)}
                    end || C <- lists:seq(0, M - 1)],
            case [E || {_, {error, _} = E} <- Cols] of
                [] -> {ok, [{C, R} || {C, {ok, R}} <- Cols]};
                [E | _] -> E
            end;
        Error -> Error
    end.

%% One right-hand side's solution: with a guessed denominator D, F * D reconstructed as a polynomial (far fewer
%% points than a rational function) and the factors it still contains divided out; without one, or if that fails
%% its final check, as a rational function
column(Prog, Vars, NP, M, W, C, none, Opts) ->
    with_job(Prog, length(Vars) - 1, Opts, {linsys, NP, M, W, C}, fun(Job) -> rational_job(Job, Vars) end);
column(Prog, Vars, {N, {PA, PB}} = NP, M, W, C, {ok, Factors}, Opts) ->
    Scale = [{Cs, C0, E} || {{Cs, C0}, E} <- Factors],
    case with_job(Prog, length(Vars), Opts, {linsys, {N, {PA, PB, Scale}}, M, W, C}, fun(Job) -> scaled_job(Job, Vars, Factors) end) of
        {ok, _} = Ok -> Ok;
        {error, _} -> column(Prog, Vars, NP, M, W, C, none, Opts)
    end.

%% F * D as a polynomial (Zippel, as polynomial/3, on the scaled black box), then N / D reduced: each factor of D
%% divided out of the numerator as often as it goes; D normalised as unshift/6 does (its last monomial's coefficient 1)
scaled_job(#{bb := BB} = Job, Vars, Factors) ->
    Poly = Job#{key => BB, outs => 1, task => fun(P, R, NM, SKey, Suffix) -> {column, BB, P, R, NM, SKey, Suffix} end},
    [{Support, Rat}] = solve(Poly),
    ok = check_poly(Poly, Support, Rat),
    Num0 = maps:from_list([{S, C} || {S, C} <- lists:zip(Support, hd(Rat)), C =/= {0, 1}]),
    {Num, Den} = lists:foldl(fun({{Cs, C0}, E}, {Nu, De}) ->
                                     F = linear_poly(Cs, C0),
                                     {Nu2, Left} = divide_times(Nu, F, E),
                                     {Nu2, lists:foldl(fun(_, A) -> pmul_q(A, F) end, De, lists:seq(1, Left))}
                             end, {Num0, #{lists:duplicate(length(Vars), 0) => {1, 1}}}, Factors),
    {_, Lead} = lists:max(maps:to_list(Den)),
    Inv = pq_sym_exact:rinv(Lead),
    {to_map(Vars, pscale_q(Num, Inv)), to_map(Vars, pscale_q(Den, Inv))}.

linear_poly(Cs, C0) ->
    N = length(Cs),
    Unit = fun(I) -> [case J of I -> 1; _ -> 0 end || J <- lists:seq(1, N)] end,
    maps:from_list([{Unit(I), {C, 1}} || {I, C} <- lists:enumerate(Cs), C =/= 0] ++ [{lists:duplicate(N, 0), {C0, 1}} || C0 =/= 0]).

%% {P / F^k, E - k}: F divided out of P as often as it goes, at most E times
divide_times(P, _F, 0) -> {P, 0};
divide_times(P, F, E) ->
    case map_size(P) > 0 andalso pq_simplify:divide(P, F) of
        {ok, Q} -> divide_times(Q, F, E - 1);
        _ -> {P, E}
    end.

pmul_q(A, B) ->
    lists:foldl(fun(T, Acc) -> padd_q(Acc, T) end, #{},
                [#{[X + Y || {X, Y} <- lists:zip(Ka, Kb)] => pq_sym_exact:rmul(Va, Vb)} || Ka := Va <- A, Kb := Vb <- B]).

%% Guessing denominators (as in the reduction programs' denominator guessing): every non-zero right-hand side's
%% univariate denominator on one random line per variable (Thiele, all right-hand sides from the same points,
%% through the linear system's cache), matched against an alphabet of candidate factors (pq_simplify:alphabet/1):
%% each factor's multiplicity is how often its root on the line divides the denominator there. A guess stands only if
%% the factors found account for the whole univariate denominator on every line; otherwise the column has none.
%% #{Column => {ok, [{{Coefficients, Constant}, Power}]} | none}
guess_denominators(_Prog, _Vars, _NP, _M, _W, [], _Opts) -> #{};
guess_denominators(Prog, Vars, NP, M, W, Cols, Opts) ->
    Alpha = [linear_form(A, Vars) || A <- pq_simplify:alphabet(Vars)],
    case with_job(Prog, length(Vars) - 1, Opts, {linsys, NP, M, W, 0},
                  fun(#{bb := BB, run := Run, seed := Seed}) ->
                          P = pq_ff:prime(1),
                          Lines = [begin
                                       Base = rnd_list(P, length(Vars), {Seed, guess, I}),
                                       Ts = distinct(P, ?GUESS_POINTS, {Seed, guess_t, I}),
                                       [{ok, Rows}] = Run([{linsys_all, BB, P, [setnth(I, Base, T) || T <- Ts]}]),
                                       {I, Base, lists:zip(Ts, Rows)}
                                   end || I <- lists:seq(1, length(Vars))],
                          maps:from_list([{C, guess_column(P, Alpha, Lines, C)} || C <- Cols])
                  end) of
        {ok, Guesses} -> Guesses;
        {error, _} -> #{}
    end.

guess_column(P, Alpha, Lines, C) ->
    PerLine = [line_factors(P, Alpha, I, Base, [{T, lists:nth(C + 1, Row)} || {T, Row} <- Pts, Row =/= bad]) || {I, Base, Pts} <- Lines],
    case lists:member(none, PerLine) of
        true -> none;
        false ->
            Max = lists:foldl(fun(Ms, Acc) -> maps:fold(fun(F, E, A) -> A#{F => max(E, maps:get(F, A, 0))} end, Acc, Ms) end, #{}, PerLine),
            {ok, [{F, E} || F := E <- Max, E > 0]}
    end.

%% the alphabet's factors in the univariate denominator along variable I through Base: #{Factor => Power} or none
line_factors(P, Alpha, I, Base, Pts) ->
    case thiele_fit(P, Pts, {[], []}, 0) of
        none -> none;
        {Xs, As} ->
            {_, Q} = convergent(P, Xs, As),
            Dq = pdeg(Q),
            Ms = maps:from_list([{F, root_power(P, Q, mod(-B * pq_ff:inv(A, P), P))}
                                 || {Cs, C0} = F <- Alpha,
                                    A <- [mod(lists:nth(I, Cs), P)], A =/= 0,
                                    B <- [mod(C0 + lists:sum([Cj * Xj || {J, Cj, Xj} <- lists:zip3(lists:seq(1, length(Cs)), Cs, Base), J =/= I]), P)]]),
            case lists:sum(maps:values(Ms)) =:= Dq of
                true -> Ms;
                false -> none
            end
    end.

%% Thiele's continued fraction through the points until it predicts ?TERMINATE + 1 in a row
thiele_fit(_P, _Pts, _State, Quiet) when Quiet > ?TERMINATE -> none;
thiele_fit(_P, [], _State, _Quiet) -> none;
thiele_fit(P, [{T, Y} | Rest], State, Quiet) ->
    case predict(P, State, T) of
        Y when Quiet >= ?TERMINATE -> State;
        Y -> thiele_fit(P, Rest, State, Quiet + 1);
        _ -> thiele_fit(P, Rest, thiele_add(P, State, T, Y), 0)
    end.

%% how often (x - R) divides the polynomial Q (coefficients from x^0 up, mod P)
root_power(P, Q, R) ->
    case pdeg(Q) < 1 of
        true -> 0;
        false ->
            %% synthetic division: Q = (x - R) Q2 + Rem
            [Lead | Rev] = lists:reverse(lists:sublist(Q, pdeg(Q) + 1)),
            {Q2Rev, Rem} = lists:foldl(fun(Ck, {Acc, Carry}) -> V = mod(Ck + Carry * R, P), {[Carry | Acc], V} end, {[], Lead}, Rev),
            case Rem of
                0 -> 1 + root_power(P, Q2Rev, R);
                _ -> 0
            end
    end.

linear_form(Ast, Vars) ->
    Zero = maps:from_list([{V, 0} || V <- Vars]),
    {ok, C0} = pq_symbolic:evaluate(Ast, Zero, #{exact => true}),
    {[begin {ok, X} = pq_symbolic:evaluate(Ast, Zero#{V => 1}, #{exact => true}), X - C0 end || V <- Vars], C0}.

%% the right-hand sides whose solution component W is 0 at three random points of the first prime (a rational
%% function that isn't 0 vanishes there with probability about its degree / p, 10^-17)
zero_columns(_Prog, _Vars, _NP, 1, _W, _Opts) -> {ok, []};
zero_columns(Prog, Vars, NP, M, W, Opts) ->
    with_job(Prog, length(Vars) - 1, Opts, {linsys, NP, M, W, 0},
             fun(#{bb := BB, run := Run, seed := Seed}) ->
                     P = pq_ff:prime(1),
                     Points = [rnd_list(P, length(Vars), {Seed, zero, I}) || I <- [1, 2, 3]],
                     [{ok, Rows}] = Run([{linsys_all, BB, P, Points}]),
                     Good = [Row || Row <- Rows, Row =/= bad],
                     Good =/= [] orelse throw({recon, singular_system}),
                     [C || C <- lists:seq(0, M - 1), lists:all(fun(Row) -> lists:nth(C + 1, Row) =:= 0 end, Good)]
             end).

ts(#{seed := Seed}, P, Count) ->
    distinct(P, Count, {Seed, P, ts}).

%% the first Count distinct values of the key's stream, in its order: a prefix of any longer request, so that
%% reconstructions of different degrees share their sample points (and a linear system's cache serves them all)
distinct(P, Count, Key) ->
    L = dedup(rnd_list(P, 2 * Count + 8, Key), #{}),
    case length(L) >= Count of
        true -> lists:sublist(L, Count);
        false -> distinct(P, Count, {Key, again})
    end.

dedup([], _) -> [];
dedup([X | Xs], Seen) when is_map_key(X, Seen) -> dedup(Xs, Seen);
dedup([X | Xs], Seen) -> [X | dedup(Xs, Seen#{X => true})].

%%% The job

with_job(Program, N, Opts, Body) -> with_job(Program, N, Opts, program, Body).

%% Kind: program (the program's first output) or {linsys, {Size, Pattern}, Rhs, Want, Col} (R1 phase 5)
with_job(Program, N, Opts, Kind, Body) ->
    case pq_ff:prepare(Program) of
        {ok, _} ->
            Exec = maps:get(executor, Opts, local),
            Key = pq_rstore:put(Program),
            BB = case Kind of
                     program -> Key;
                     {linsys, {Size, Pattern}, Rhs, Want, Col} -> {linsys, Key, pq_rstore:put(Pattern), Size, Rhs, Want, Col}
                 end,
            Job = #{key => Key, bb => BB, run => fun(Tasks) -> pq_rexec:map(Tasks, Exec, Opts) end,
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
    columns(Job, P, K, S, {R, SKey}, Rest, Newton0, 0, 0, min(maps:get(batch, Job), ?TERMINATE + 2)).

%% a stage needs about degree + 1 + Terminate columns, unknown in advance: the first round computes Terminate + 2,
%% each further round twice as many, up to the executor's slots -- so a small answer doesn't pay for columns
%% computed ahead and never used (20 cores were slower than one on the IBP benchmark before), and a large one still
%% fills every slot. The columns are taken in order, so the result is the same for any round sizes.
columns(#{run := Run, batch := Batch, seed := Seed, task := Task} = Job, P, K, S, {R, SKey}, Rest, Newton, J, Quiet, Round) ->
    J > ?MAX_DEGREE andalso throw({recon, degree_bound}),
    Bs = [point_value(P, {Seed, P, K, b, I}) || I <- lists:seq(J, J + Round - 1)],
    Results = Run([Task(P, R, length(S), SKey, [B | Rest]) || B <- Bs]),
    case take(P, lists:zip(Bs, Results), Newton, Quiet) of
        {done, Newton2} -> support(P, S, Newton2);
        {more, Newton2, Quiet2} -> columns(Job, P, K, S, {R, SKey}, Rest, Newton2, J + Round, Quiet2, min(Batch, 2 * Round))
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
degrees(#{bb := Key, run := Run, seed := Seed, n := M}, P, Shift) ->
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
    {Pk, Qk} = convergent(P, Xs, As),
    {pdeg(Pk), pdeg(Qk)}.

%% the continued fraction's numerator and denominator polynomials (coefficients from x^0 up, mod P)
convergent(P, Xs, As) ->
    {Pk, Qk, _, _} = lists:foldl(fun({I, A}, {P1, Q1, P2, Q2}) ->
                                         case I of
                                             0 -> {[A], [1], [1], [0]};
                                             _ ->
                                                 X = lists:nth(I, Xs),
                                                 Lin = [mod(-X, P), 1],
                                                 {padd(P, pscale(P, P1, A), pmul(P, Lin, P2)), padd(P, pscale(P, Q1, A), pmul(P, Lin, Q2)), P1, Q1}
                                         end
                                 end, {[0], [1], [1], [0]}, lists:enumerate(0, As)),
    {Pk, Qk}.

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

check_rat(#{bb := Key, run := Run, seed := Seed}, Vars, Num, Den) ->
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
