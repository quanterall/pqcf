%% States (FS 3.4, 4.5, FR-G07, FR-G09, A.12): the Fock basis index of FS 3.4 (normative for every backend and
%% test), Gaussian state preparation, xxpp/xpxp conversion, reduced states and the physicality check.
-module(pq_state).
-include("pq.hrl").
-export([basis_size/2, index/3, pattern/3,
         vacuum/2, coherent/3, squeezed_vacuum/3, thermal/3,
         to_xpxp/1, from_xpxp/1, from_xpxp/2, reduce/2, check_physical/1,
         mean_cov/1, from_mean_cov/3,
         fock_state/2, fock_superposition/2, fock_dm/1]).

%%% Fock basis (FS 3.4): all occupation tuples with total at most Cutoff, ordered by total photon number, then
%%% lexicographically ascending from mode 0. D = C(Cutoff + m, m).

-spec basis_size(pos_integer(), non_neg_integer()) -> pos_integer().
basis_size(NModes, Cutoff) -> binomial(Cutoff + NModes, NModes).

-spec index([non_neg_integer()], pos_integer(), non_neg_integer()) ->
          {ok, non_neg_integer()} | {error, pq_core:reason()}.
index(Pattern, NModes, Cutoff) when length(Pattern) =:= NModes ->
    N = lists:sum(Pattern),
    case N =< Cutoff andalso lists:all(fun(X) -> is_integer(X) andalso X >= 0 end, Pattern) of
        true -> {ok, below(N, NModes) + rank(Pattern, N)};
        false -> {error, {index_out_of_range, #{pattern => Pattern, cutoff => Cutoff}}}
    end;
index(Pattern, NModes, _Cutoff) -> {error, {bad_shape, #{pattern => Pattern, n_modes => NModes}}}.

-spec pattern(non_neg_integer(), pos_integer(), non_neg_integer()) ->
          {ok, [non_neg_integer()]} | {error, pq_core:reason()}.
pattern(I, NModes, Cutoff) ->
    case is_integer(I) andalso I >= 0 andalso I < basis_size(NModes, Cutoff) of
        true ->
            N = total_at(I, NModes, 0),
            {ok, unrank(I - below(N, NModes), N, NModes)};
        false -> {error, {index_out_of_range, I}}
    end.

%% the number of patterns with total below N in M modes: C(N - 1 + M, M)
below(0, _M) -> 0;
below(N, M) -> binomial(N - 1 + M, M).

%% patterns in K modes with total S
count(0, 0) -> 1;
count(0, _) -> 0;
count(K, S) -> binomial(S + K - 1, K - 1).

rank([_], _N) -> 0;
rank([P | Rest], N) ->
    lists:sum([count(length(Rest), N - A) || A <- lists:seq(0, P - 1)]) + rank(Rest, N - P).

total_at(I, M, N) ->
    case below(N + 1, M) > I of
        true -> N;
        false -> total_at(I, M, N + 1)
    end.

unrank(_R, N, 1) -> [N];
unrank(R, N, M) -> first(R, N, M, 0).

first(R, N, M, A) ->
    C = count(M - 1, N - A),
    case R < C of
        true -> [A | unrank(R, N - A, M - 1)];
        false -> first(R - C, N, M, A + 1)
    end.

binomial(_N, 0) -> 1;
binomial(N, K) when K > N -> 0;
binomial(N, K) -> binomial(N, K, 1, 1).
binomial(_N, 0, Num, Den) -> Num div Den;
binomial(N, K, Num, Den) -> binomial(N - 1, K - 1, Num * N, Den * K).

%%% Gaussian states (FR-G09)

-spec vacuum(pos_integer(), number()) -> #gaussian_state{}.
vacuum(N, Hbar) ->
    from_mean_cov(lists:duplicate(2 * N, 0.0), diag(lists:duplicate(2 * N, Hbar / 2)), Hbar).

%% the vacuum displaced by Alphas (FS 3.3): x = sqrt(2 hbar) Re alpha, p = sqrt(2 hbar) Im alpha
-spec coherent(pos_integer(), [number() | {complex, number(), number()}], number()) -> #gaussian_state{}.
coherent(N, Alphas, Hbar) when length(Alphas) =:= N ->
    K = math:sqrt(2 * Hbar),
    Cs = [cplx(A) || A <- Alphas],
    from_mean_cov([K * Re || {Re, _} <- Cs] ++ [K * Im || {_, Im} <- Cs], diag(lists:duplicate(2 * N, Hbar / 2)), Hbar).

%% the vacuum squeezed by Zs, z = R e^(i Phi) (A.4): per mode (hbar/2) S S^T
-spec squeezed_vacuum(pos_integer(), [number() | {complex, number(), number()}], number()) -> #gaussian_state{}.
squeezed_vacuum(N, Zs, Hbar) when length(Zs) =:= N ->
    Blocks = [begin
                  {Re, Im} = cplx(Z),
                  R = math:sqrt(Re * Re + Im * Im), Phi = math:atan2(Im, Re),
                  C = math:cosh(R), S = math:sinh(R),
                  {A, B, D} = {C - math:cos(Phi) * S, -math:sin(Phi) * S, C + math:cos(Phi) * S},
                  %% S S^T for S = [[A, B], [B, D]]
                  {A * A + B * B, A * B + B * D, B * B + D * D}
              end || Z <- Zs],
    Cov = [[cov_entry(I, J, N, Blocks, Hbar) || J <- lists:seq(0, 2 * N - 1)] || I <- lists:seq(0, 2 * N - 1)],
    from_mean_cov(lists:duplicate(2 * N, 0.0), Cov, Hbar).

cov_entry(I, J, N, Blocks, Hbar) ->
    {Mi, Qi} = {I rem N, I div N},
    {Mj, Qj} = {J rem N, J div N},
    case Mi =:= Mj of
        false -> 0.0;
        true ->
            {Xx, Xp, Pp} = lists:nth(Mi + 1, Blocks),
            Hbar / 2 * case {Qi, Qj} of {0, 0} -> Xx; {1, 1} -> Pp; _ -> Xp end
    end.

%% covariance (2 nbar + 1)(hbar/2) I per mode, zero mean
-spec thermal(pos_integer(), [number()], number()) -> #gaussian_state{}.
thermal(N, MeanPhotons, Hbar) when length(MeanPhotons) =:= N ->
    V = [(2 * Nb + 1) * Hbar / 2 || Nb <- MeanPhotons],
    from_mean_cov(lists:duplicate(2 * N, 0.0), diag(V ++ V), Hbar).

cplx({complex, Re, Im}) -> {float(num(Re)), float(num(Im))};
cplx(X) -> {float(num(X)), 0.0}.

num({rat, N, D}) -> N / D;
num(X) -> X.

diag(Vs) ->
    N = length(Vs),
    [[case I =:= J of true -> float(lists:nth(I, Vs)); false -> 0.0 end || J <- lists:seq(1, N)] || I <- lists:seq(1, N)].

%% a state from lists (mean of 2m, covariance rows), packed
-spec from_mean_cov([number()], [[number()]], number()) -> #gaussian_state{}.
from_mean_cov(Mean, Cov, Hbar) ->
    {ok, M} = pq_core:pack({rvec, Mean}),
    {ok, C} = pq_core:pack({rmat, Cov}),
    #gaussian_state{n_modes = length(Mean) div 2, hbar = Hbar, mean = M, cov = C}.

%% {Mean, CovRows} as lists of floats
-spec mean_cov(#gaussian_state{}) -> {[float()], [[float()]]}.
mean_cov(#gaussian_state{mean = M, cov = C}) ->
    {ok, {rvec, Mean}} = pq_core:unpack(M),
    {ok, {rmat, Cov}} = pq_core:unpack(C),
    {Mean, Cov}.

%%% FR-G07

%% {Mean, Cov} in xpxp order (x_0, p_0, x_1, p_1, ...), as lists
-spec to_xpxp(#gaussian_state{}) -> {[float()], [[float()]]}.
to_xpxp(#gaussian_state{n_modes = N} = S) ->
    {Mean, Cov} = mean_cov(S),
    P = xpxp_order(N),
    {[lists:nth(I + 1, Mean) || I <- P], [[lists:nth(J + 1, lists:nth(I + 1, Cov)) || J <- P] || I <- P]}.

%% FS 7: with the configured hbar
-spec from_xpxp({[number()], [[number()]]}) -> #gaussian_state{}.
from_xpxp(MC) -> from_xpxp(MC, pq_core:get(hbar)).

-spec from_xpxp({[number()], [[number()]]}, number()) -> #gaussian_state{}.
from_xpxp({Mean, Cov}, Hbar) ->
    N = length(Mean) div 2,
    %% the xxpp index of position k in xpxp order, inverted
    Inv = [K || {_, K} <- lists:sort(lists:zip(xpxp_order(N), lists:seq(0, 2 * N - 1)))],
    from_mean_cov([lists:nth(K + 1, Mean) || K <- Inv], [[lists:nth(L + 1, lists:nth(K + 1, Cov)) || L <- Inv] || K <- Inv], Hbar).

xpxp_order(N) -> lists:append([[J, J + N] || J <- lists:seq(0, N - 1)]).

%% the reduced state on Modes (partial trace), in the order of Modes
-spec reduce(#gaussian_state{}, [non_neg_integer()]) -> {ok, #gaussian_state{}} | {error, pq_core:reason()}.
reduce(#gaussian_state{hbar = H, mean = M, cov = C}, Modes) ->
    case pq_nif:gaussian_marginal(M, C, Modes) of
        {ok, {M2, C2}} -> {ok, #gaussian_state{n_modes = length(Modes), hbar = H, mean = M2, cov = C2}};
        {error, _} = E -> E
    end.

%% A.12: Cov + i (hbar/2) Omega positive semidefinite within physicality_tol
-spec check_physical(#gaussian_state{}) -> ok | {error, pq_core:reason()}.
check_physical(#gaussian_state{hbar = H, cov = C}) ->
    case pq_nif:gaussian_physicality(C, float(H)) of
        {ok, L} ->
            case L >= -pq_core:get(physicality_tol) of
                true -> ok;
                false -> {error, {non_physical_state, #{min_eigenvalue => L}}}
            end;
        {error, _} = E -> E
    end.

%%% FR-F01: Fock states from occupation patterns, superpositions within the cutoff, density matrices from pure states

%% the basis state |Pattern> with the given cutoff (n_modes = length(Pattern))
-spec fock_state([non_neg_integer()], non_neg_integer()) -> {ok, #fock_state{}} | {error, pq_core:reason()}.
fock_state(Pattern, Cutoff) -> fock_superposition([{1, Pattern}], Cutoff).

%% sum of Coeff |Pattern>, normalised; Coeff a number or {complex, Re, Im}; repeated patterns add up
-spec fock_superposition([{number() | {complex, number(), number()}, [non_neg_integer()]}], non_neg_integer()) ->
          {ok, #fock_state{}} | {error, pq_core:reason()}.
fock_superposition([{_, P0} | _] = Terms, Cutoff) when is_integer(Cutoff), Cutoff >= 0 ->
    N = length(P0),
    try
        N >= 1 orelse throw({bad_shape, #{pattern => P0}}),
        Acc = lists:foldl(fun({C, P}, M) ->
                                  length(P) =:= N orelse throw({bad_shape, #{pattern => P, n_modes => N}}),
                                  lists:all(fun(X) -> is_integer(X) andalso X >= 0 end, P) orelse throw({bad_shape, #{pattern => P}}),
                                  lists:sum(P) =< Cutoff orelse throw({unsupported_input, #{pattern => P, cutoff => Cutoff}}),
                                  {ok, I} = index(P, N, Cutoff),
                                  {Re, Im} = coeff(C),
                                  maps:update_with(I, fun({A, B}) -> {A + Re, B + Im} end, {Re, Im}, M)
                          end, #{}, Terms),
        Norm = math:sqrt(lists:sum([A * A + B * B || {A, B} <- maps:values(Acc)])),
        Norm > 0 orelse throw({bad_shape, zero_norm}),
        D = basis_size(N, Cutoff),
        Amps = << <<(element(1, Z) / Norm):64/float-little, (element(2, Z) / Norm):64/float-little>>
                  || K <- lists:seq(0, D - 1), Z <- [maps:get(K, Acc, {0.0, 0.0})] >>,
        {ok, #fock_state{n_modes = N, cutoff = Cutoff, amps = {cvec, D, Amps}}}
    catch
        throw:{Cat, _} = R when is_atom(Cat) -> {error, R}
    end;
fock_superposition(Terms, Cutoff) -> {error, {bad_shape, #{terms => Terms, cutoff => Cutoff}}}.

coeff({complex, Re, Im}) when is_number(Re), is_number(Im) -> {float(Re), float(Im)};
coeff(X) when is_number(X) -> {float(X), 0.0};
coeff(X) -> throw({bad_shape, #{coefficient => X}}).

%% |psi><psi| of a pure state, or the mixture sum w_i |psi_i><psi_i| (weights normalised, states of one shape)
-spec fock_dm(#fock_state{} | [{number(), #fock_state{}}]) -> {ok, #fock_dm{}} | {error, pq_core:reason()}.
fock_dm(#fock_state{} = S) -> {ok, pq_fock_dm:from_pure(S)};
fock_dm([{_, #fock_state{n_modes = N, cutoff = C}} | _] = Mix) ->
    W = lists:sum([X || {X, _} <- Mix]),
    case lists:all(fun({X, #fock_state{n_modes = N1, cutoff = C1}}) -> is_number(X) andalso X >= 0 andalso N1 =:= N andalso C1 =:= C;
                      (_) -> false
                   end, Mix) andalso W > 0 of
        false -> {error, {bad_shape, mixture}};
        true ->
            Rhos = [{X / W, element(2, pq_core:unpack((pq_fock_dm:from_pure(S))#fock_dm.rho))} || {X, S} <- Mix],
            D = basis_size(N, C),
            Sum = [[lists:foldl(fun({Wi, {cmat, R}}, {complex, A, B}) ->
                                        {complex, Ai, Bi} = lists:nth(J, lists:nth(I, R)),
                                        {complex, A + Wi * Ai, B + Wi * Bi}
                                end, {complex, 0.0, 0.0}, Rhos) || J <- lists:seq(1, D)] || I <- lists:seq(1, D)],
            {ok, Rho} = pq_core:pack({cmat, Sum}),
            {ok, #fock_dm{n_modes = N, cutoff = C, rho = Rho}}
    end;
fock_dm(Other) -> {error, {bad_shape, #{fock_dm => Other}}}.

