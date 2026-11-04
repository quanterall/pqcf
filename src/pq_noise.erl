%% Noise (FS 5.8): loss channels inserted into a circuit. Loss is exact on Gaussian states and goes through
%% density matrices on Fock states (FR-N03, FR-F05).
%%
%% Phase 4 (FS 15 item 4, docs/decisions/N1-noise-models.md): partial distinguishability of the photons of boson
%% sampling (visibility => M, or a Gram matrix), and dark counts of the detectors (dark_counts => Nu: Poisson counts
%% of mean Nu per detector and shot, added to photon counts, a click for threshold detectors).
-module(pq_noise).
-include("pq.hrl").
-export([insert_loss/3, detector_efficiency/2]).
-export([shared_probabilities/2, visibility_gram/2, gram/2, dark_means/2, dark_rows/4, dark_distribution/4]).

%% FR-N02: after every instruction that is not a measurement, a loss on each mode it touches; Etas maps modes to
%% transmissivities (numeric or symbolic), modes without an entry get none. Opts is reserved.
-spec insert_loss(#circuit{}, #{non_neg_integer() => pq_symbolic:sym_expr()}, map()) -> {ok, #circuit{}}.
insert_loss(#circuit{instrs = Is} = C, Etas, _Opts) ->
    New = lists:append([[I | [{loss, M, maps:get(M, Etas)} || M <- touched(I), maps:is_key(M, Etas)]]
                        || I <- Is]),
    {ok, C#circuit{instrs = New}}.

%% FR-N04: a loss of transmissivity Eta (or Etas per mode) immediately before each measurement, on its modes
-spec detector_efficiency(#circuit{}, pq_symbolic:sym_expr() | #{non_neg_integer() => pq_symbolic:sym_expr()}) ->
          {ok, #circuit{}}.
detector_efficiency(#circuit{instrs = Is} = C, Eta) ->
    EtaOf = fun(M) when is_map(Eta) -> maps:get(M, Eta, none); (_) -> Eta end,
    New = lists:append([case measured(I) of
                            [] -> [I];
                            Ms -> [{loss, M, EtaOf(M)} || M <- Ms, EtaOf(M) =/= none] ++ [I]
                        end || I <- Is]),
    {ok, C#circuit{instrs = New}}.

touched({beamsplitter, A, B, _, _}) -> [A, B];
touched({unitary, Ms, _}) -> Ms;
touched({custom, _, Args}) -> maps:get(modes, Args, []);
touched({Tag, _}) when Tag =:= measure_pnr; Tag =:= measure_threshold -> [];
touched({measure_homodyne, _, _}) -> [];
touched({loss, _, _}) -> [];
touched(I) -> [element(2, I)].

measured({measure_pnr, Ms}) -> Ms;
measured({measure_threshold, Ms}) -> Ms;
measured({measure_homodyne, M, _}) -> [M];
measured(_) -> [].

%%% Partial distinguishability

%% visibility => M | [M_i] (per input mode): photon a of mode i is in the shared internal state with probability
%% p_i = sqrt(M_i), so two photons with the same M show the Hong-Ou-Mandel visibility M
-spec shared_probabilities(pos_integer(), number() | [number()]) -> {ok, [float()]} | {error, pq_core:reason()}.
shared_probabilities(N, M) when is_number(M) -> shared_probabilities(N, lists:duplicate(N, M));
shared_probabilities(N, Ms) when is_list(Ms), length(Ms) =:= N ->
    case lists:all(fun(X) -> is_number(X) andalso X >= 0 andalso X =< 1 end, Ms) of
        true -> {ok, [math:sqrt(X) || X <- Ms]};
        false -> {error, {bad_shape, #{visibility => Ms, expected => <<"numbers in [0, 1]">>}}}
    end;
shared_probabilities(N, Ms) -> {error, {bad_shape, #{visibility => Ms, modes => N}}}.

%% the Gram matrix of that model over the photons in input-mode order: 1 on the diagonal, sqrt(p_a p_b) elsewhere
-spec visibility_gram([non_neg_integer()], number() | [number()]) -> {ok, pq_core:packed_cmatrix()} | {error, pq_core:reason()}.
visibility_gram(In, M) ->
    case shared_probabilities(length(In), M) of
        {ok, Ps} ->
            Ph = lists:append([lists:duplicate(K, P) || {K, P} <- lists:zip(In, Ps)]),
            Rows = [[case A =:= B of true -> 1.0; false -> math:sqrt(Pa * Pb) end
                     || {B, Pb} <- lists:enumerate(Ph)] || {A, Pa} <- lists:enumerate(Ph)],
            pq_core:pack({cmat, Rows});
        E -> E
    end.

%% gram => Rows (n x n numbers or {complex, Re, Im}, the photons in input-mode order): Hermitian with a unit
%% diagonal; positivity shows as the distribution summing to 1 (checked by the caller)
-spec gram([non_neg_integer()], term()) -> {ok, pq_core:packed_cmatrix()} | {error, pq_core:reason()}.
gram(In, Rows) ->
    N = lists:sum(In),
    C = fun({complex, R, I}) -> {float(R), float(I)}; (X) when is_number(X) -> {float(X), 0.0}; (_) -> bad end,
    Ok = is_list(Rows) andalso length(Rows) =:= N andalso
        lists:all(fun(R) -> is_list(R) andalso length(R) =:= N andalso not lists:member(bad, [C(X) || X <- R]) end, Rows),
    case Ok of
        false -> {error, {bad_shape, #{gram => Rows, expected => {N, N}}}};
        true ->
            G = [[C(X) || X <- R] || R <- Rows],
            At = fun(A, B) -> lists:nth(B, lists:nth(A, G)) end,
            Herm = lists:all(fun({A, B}) ->
                                     {R1, I1} = At(A, B), {R2, I2} = At(B, A),
                                     abs(R1 - R2) < 1.0e-12 andalso abs(I1 + I2) < 1.0e-12
                             end, [{A, B} || A <- lists:seq(1, N), B <- lists:seq(1, N)]),
            Diag = lists:all(fun(A) -> {R, I} = At(A, A), abs(R - 1) < 1.0e-12 andalso abs(I) < 1.0e-12 end, lists:seq(1, N)),
            case {Herm andalso Diag, Herm andalso Diag andalso psd(G)} of
                {true, true} -> pq_core:pack({cmat, [[{complex, R, I} || {R, I} <- Row] || Row <- G]});
                {true, false} -> {error, {bad_shape, #{gram => not_positive_semidefinite}}};
                {false, _} -> {error, {bad_shape, #{gram => not_hermitian_unit_diagonal}}}
            end
    end.

%% positive semidefinite: a Cholesky factorisation G = L L^dagger whose pivots are all >= -1e-12 (a pivot below
%% 1e-12 is a zero column, the vectors being linearly dependent)
psd(G) ->
    N = length(G),
    At = fun(A, B) -> lists:nth(B + 1, lists:nth(A + 1, G)) end,
    psd(0, N, At, #{}).

psd(J, N, _At, _L) when J =:= N -> true;
psd(J, N, At, L) ->
    Lj = fun(I, K) -> maps:get({I, K}, L, {0.0, 0.0}) end,
    {Gr, _} = At(J, J),
    D = Gr - lists:sum([R * R + I * I || K <- lists:seq(0, J - 1), {R, I} <- [Lj(J, K)]]),
    if
        D < -1.0e-12 -> false;
        D < 1.0e-12 ->
            %% a zero pivot: the rest of the column must vanish too
            Col = [cdot_rest(At(I, J), I, J, Lj) || I <- lists:seq(J + 1, N - 1)],
            lists:all(fun({R, I}) -> abs(R) < 1.0e-9 andalso abs(I) < 1.0e-9 end, Col) andalso psd(J + 1, N, At, L);
        true ->
            S = math:sqrt(D),
            L2 = lists:foldl(fun(I, Acc) -> {R, Im} = cdot_rest(At(I, J), I, J, Lj), Acc#{{I, J} => {R / S, Im / S}} end,
                             L#{{J, J} => {S, 0.0}}, lists:seq(J + 1, N - 1)),
            psd(J + 1, N, At, L2)
    end.

%% G[i][j] - sum_k L[i][k] conj(L[j][k])
cdot_rest({Gr, Gi}, I, J, Lj) ->
    lists:foldl(fun(K, {R, Im}) ->
                        {Ar, Ai} = Lj(I, K), {Br, Bi} = Lj(J, K),
                        {R - (Ar * Br + Ai * Bi), Im - (Ai * Br - Ar * Bi)}
                end, {Gr, Gi}, lists:seq(0, J - 1)).

%%% Dark counts

%% dark_counts => Nu | #{Mode => Nu}: the mean per measured column (modes without an entry: none)
-spec dark_means(number() | map(), [non_neg_integer()]) -> {ok, [float()]} | {error, pq_core:reason()}.
dark_means(Nu, Modes) ->
    Means = [case Nu of
                 #{} -> maps:get(M, Nu, 0.0);
                 _ -> Nu
             end || M <- Modes],
    case lists:all(fun(X) -> is_number(X) andalso X >= 0 end, Means) of
        true -> {ok, [float(X) || X <- Means]};
        false -> {error, {bad_shape, #{dark_counts => Nu, expected => <<"non-negative means">>}}}
    end.

%% the rows of shots Start, Start + 1, ... with dark counts: shot k draws one uniform per counting column (in column
%% order) from stream (dark_seed(Seed), k), so the noise does not depend on chunking or on the sampler's own
%% streams (D1); a PNR column gets n + Poisson(Nu) by inversion, a threshold column clicks when it did or with
%% probability 1 - exp(-Nu); homodyne columns (floats) are untouched
-spec dark_rows([list()], non_neg_integer(), map(), [pnr | threshold | homodyne]) -> [list()].
dark_rows(Rows, Start, #{seed := Seed, dark_means := Means}, Kinds) ->
    DSeed = dark_seed(Seed),
    Cols = length([K || K <- Kinds, K =/= homodyne]),
    [begin
         {ok, Us} = pq_nif:rng_u64(DSeed, Start + I, Cols),
         noisy(Row, Kinds, Means, [(U bsr 11) * math:pow(2, -53) || U <- Us])
     end || {I, Row} <- lists:enumerate(0, Rows)].

dark_seed(Seed) -> (Seed bxor 16#6461726B636F756E) band 16#FFFFFFFFFFFFFFFF.

noisy([], [], [], _) -> [];
noisy([X | Row], [homodyne | Ks], [_ | Ms], Us) -> [X | noisy(Row, Ks, Ms, Us)];
noisy([X | Row], [threshold | Ks], [Nu | Ms], [U | Us]) ->
    [case X > 0 orelse U < (1 - math:exp(-Nu)) of true -> 1; false -> 0 end | noisy(Row, Ks, Ms, Us)];
noisy([X | Row], [pnr | Ks], [Nu | Ms], [U | Us]) -> [X + poisson(Nu, U) | noisy(Row, Ks, Ms, Us)].

poisson(Nu, _U) when Nu == 0 -> 0;
poisson(Nu, U) -> poisson(Nu, U, 0, math:exp(-Nu), math:exp(-Nu)).

poisson(_Nu, U, K, _P, F) when U < F; K >= 10000 -> K;
poisson(Nu, U, K, P, F) ->
    P2 = P * Nu / (K + 1),
    poisson(Nu, U, K + 1, P2, F + P2).

%% the exact distribution with dark counts: threshold patterns exactly; PNR patterns with at most K extra counts
%% per column, K the smallest with a Poisson tail below 1e-13 (the rest reported as truncated mass)
-spec dark_distribution(#{[non_neg_integer()] => float()}, [float()], pnr | threshold, map()) ->
          {#{[non_neg_integer()] => float()}, float()}.
dark_distribution(Map, Means, threshold, _Opts) ->
    Q = [(1 - math:exp(-Nu)) || Nu <- Means],
    Out = maps:fold(fun(P, Pr, Acc) ->
                            lists:foldl(fun({Pat, W}, A) -> maps:update_with(Pat, fun(V) -> V + Pr * W end, Pr * W, A) end,
                                        Acc, click_branches(P, Q))
                    end, #{}, Map),
    {Out, 0.0};
dark_distribution(Map, Means, pnr, _Opts) ->
    Extra = [[{K, pois(Nu, K)} || K <- lists:seq(0, kmax(Nu))] || Nu <- Means],
    Out = maps:fold(fun(P, Pr, Acc) ->
                            lists:foldl(fun({Add, W}, A) ->
                                                Pat = [X + Y || {X, Y} <- lists:zip(P, Add)],
                                                maps:update_with(Pat, fun(V) -> V + Pr * W end, Pr * W, A)
                                        end, Acc, combos(Extra))
                    end, #{}, Map),
    {Out, max(0.0, lists:sum(maps:values(Map)) - lists:sum(maps:values(Out)))}.

click_branches([], []) -> [{[], 1.0}];
click_branches([X | P], [Q | Qs]) ->
    Rest = click_branches(P, Qs),
    case X > 0 of
        true -> [{[1 | R], W} || {R, W} <- Rest];
        false -> [{[1 | R], W * Q} || {R, W} <- Rest, Q > 0] ++ [{[0 | R], W * (1 - Q)} || {R, W} <- Rest]
    end.

combos([]) -> [{[], 1.0}];
combos([L | Ls]) -> [{[K | R], W * Wk} || {K, Wk} <- L, {R, W} <- combos(Ls)].

pois(Nu, K) when Nu == 0 -> case K of 0 -> 1.0; _ -> 0.0 end;
pois(Nu, K) -> math:exp(K * math:log(Nu) - Nu - lgamma1(K)).

lgamma1(K) -> lists:sum([math:log(J) || J <- lists:seq(1, K)]).

kmax(Nu) when Nu == 0 -> 0;
kmax(Nu) -> kmax(Nu, 0, pois(Nu, 0)).
kmax(_Nu, K, F) when 1 - F < 1.0e-13; K >= 200 -> K;
kmax(Nu, K, F) -> kmax(Nu, K + 1, F + pois(Nu, K + 1)).
