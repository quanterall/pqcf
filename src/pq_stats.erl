%% Statistics (FS 5.12): mean_photon_number (FR-X01), g2 (FR-X02), tv_distance (FR-X03), histogram (FR-X04),
%% fidelity (FR-X05).
-module(pq_stats).
-include("pq.hrl").
-export([mean_photon_number/2, g2/2, tv_distance/2, histogram/2, fidelity/2]).

%% Gaussian state: sum over Modes of (mu_x^2 + mu_p^2 + Sigma_xx + Sigma_pp) / (2 hbar) - 1/2; Fock state: sum n p(n)
%% over Modes; a result: the sample mean of each listed column.
-spec mean_photon_number(#gaussian_state{} | #fock_state{} | #result{}, [non_neg_integer()]) ->
          {ok, float() | [float()]} | {error, pq_core:reason()}.
mean_photon_number(#gaussian_state{n_modes = N, hbar = H} = S, Modes) ->
    case valid(Modes, N) of
        true ->
            {Mu, Cov} = pq_state:mean_cov(S),
            At = fun(I) -> lists:nth(I + 1, Mu) end,
            V = fun(I) -> lists:nth(I + 1, lists:nth(I + 1, Cov)) end,
            {ok, lists:sum([(At(J) * At(J) + At(J + N) * At(J + N) + V(J) + V(J + N)) / (2 * H) - 0.5 || J <- Modes])};
        false -> {error, {mode_out_of_range, Modes}}
    end;
mean_photon_number(#fock_state{n_modes = N, cutoff = C, amps = A}, Modes) ->
    case valid(Modes, N) andalso pq_nif:fock_pnr_probabilities(A, Modes, N, C) of
        {ok, PP} -> {ok, lists:sum([lists:sum(P) * Pr || {P, Pr} <- PP])};
        false -> {error, {mode_out_of_range, Modes}};
        {error, _} = E -> E
    end;
mean_photon_number(#result{measured = Measured, samples = Samples}, Modes) ->
    Cols = [string:str(Measured, [M]) || M <- Modes],
    case lists:member(0, Cols) orelse Samples =:= [] of
        true -> {error, {mode_out_of_range, Modes}};
        false ->
            N = length(Samples),
            {ok, [lists:sum([lists:nth(Ci, Row) || Row <- Samples]) / N || Ci <- Cols]}
    end.

valid(Modes, N) -> Modes =/= [] andalso lists:all(fun(M) -> is_integer(M) andalso M >= 0 andalso M < N end, Modes)
                   andalso length(lists:usort(Modes)) =:= length(Modes).

%% FR-X02: <n(n-1)> / <n>^2 on one mode, from a Fock state or density matrix, a Gaussian state (Wick's theorem on
%% its moments, exactly) or a result's samples
-spec g2(#fock_state{} | #fock_dm{} | #gaussian_state{} | #result{}, non_neg_integer()) ->
          {ok, float()} | {error, pq_core:reason()}.
g2(#result{} = R, Mode) ->
    case column(R, Mode) of
        {ok, Ns} -> ratio(lists:sum([N * (N - 1) || N <- Ns]) / length(Ns), lists:sum(Ns) / length(Ns));
        E -> E
    end;
g2(#gaussian_state{n_modes = N, hbar = H} = S, Mode) when Mode < N ->
    %% a = alpha + d: <n> = Nc + |alpha|^2, <a+a+aa> = 2 Nc^2 + |Mc|^2 + 4 |alpha|^2 Nc + 2 Re(conj(alpha)^2 Mc) + |alpha|^4
    {Mu, Cov} = pq_state:mean_cov(S),
    At = fun(I, J) -> lists:nth(J + 1, lists:nth(I + 1, Cov)) end,
    {Vx, Vp, C} = {At(Mode, Mode), At(Mode + N, Mode + N), At(Mode, Mode + N)},
    K = math:sqrt(2 * H),
    {Ar, Ai} = {lists:nth(Mode + 1, Mu) / K, lists:nth(Mode + N + 1, Mu) / K},
    A2 = Ar * Ar + Ai * Ai,
    Nc = (Vx + Vp) / (2 * H) - 0.5,
    {Mr, Mi} = {(Vx - Vp) / (2 * H), 2 * C / (2 * H)},
    %% conj(alpha)^2 = (Ar^2 - Ai^2) - 2 i Ar Ai
    ReTerm = (Ar * Ar - Ai * Ai) * Mr + 2 * Ar * Ai * Mi,
    ratio(2 * Nc * Nc + (Mr * Mr + Mi * Mi) + 4 * A2 * Nc + 2 * ReTerm + A2 * A2, Nc + A2);
g2(Source, Mode) ->
    case distribution(Source, Mode) of
        {ok, PP} ->
            ratio(lists:sum([P * N * (N - 1) || {[N], P} <- PP]), lists:sum([P * N || {[N], P} <- PP]));
        E -> E
    end.

ratio(_, M) when M == 0 -> {error, {undefined_statistic, zero_mean}};
ratio(A, M) -> {ok, A / (M * M)}.

distribution(#fock_state{n_modes = N, cutoff = C, amps = A}, Mode) when Mode < N ->
    pq_nif:fock_pnr_probabilities(A, [Mode], N, C);
distribution(#fock_dm{n_modes = N, cutoff = C, rho = R}, Mode) when Mode < N ->
    pq_nif:fock_pnr_probabilities(R, [Mode], N, C);
distribution(_, Mode) -> {error, {mode_out_of_range, Mode}}.

column(#result{measured = Ms, samples = S}, Mode) ->
    case string:str(Ms, [Mode]) of
        0 -> {error, {mode_out_of_range, Mode}};
        _ when S =:= [] -> {error, {undefined_statistic, no_samples}};
        I -> {ok, [lists:nth(I, Row) || Row <- S]}
    end.

%% FR-X04: #{Pattern => Count} over the given columns
-spec histogram(#result{}, [non_neg_integer()]) -> {ok, #{list() => pos_integer()}} | {error, pq_core:reason()}.
histogram(#result{measured = Ms, samples = S}, Modes) ->
    Is = [string:str(Ms, [M]) || M <- Modes],
    case lists:member(0, Is) of
        true -> {error, {mode_out_of_range, Modes}};
        false ->
            {ok, lists:foldl(fun(Row, Acc) ->
                                     K = [lists:nth(I, Row) || I <- Is],
                                     maps:update_with(K, fun(X) -> X + 1 end, 1, Acc)
                             end, #{}, S)}
    end.

%% FR-X03: 1/2 sum |p - q| over the union of patterns; a result becomes its normalised histogram over its columns
-spec tv_distance(map() | #result{}, map() | #result{}) -> {ok, float()} | {error, pq_core:reason()}.
tv_distance(P, Q) ->
    case {probs(P), probs(Q)} of
        {{ok, A}, {ok, B}} ->
            Keys = lists:usort(maps:keys(A) ++ maps:keys(B)),
            {ok, 0.5 * lists:sum([abs(maps:get(K, A, 0.0) - maps:get(K, B, 0.0)) || K <- Keys])};
        {{error, _} = E, _} -> E;
        {_, E} -> E
    end.

probs(#result{measured = Ms, samples = S} = R) when S =/= [] ->
    {ok, H} = histogram(R, Ms),
    N = length(S),
    {ok, maps:map(fun(_, C) -> C / N end, H)};
probs(#result{}) -> {error, {undefined_statistic, no_samples}};
probs(M) when is_map(M) -> {ok, M}.

%% FR-X05: |<a|b>|^2 for two pure Fock states, <psi|rho|psi> for a pure state and a density matrix
-spec fidelity(#fock_state{} | #fock_dm{}, #fock_state{} | #fock_dm{}) -> {ok, float()} | {error, pq_core:reason()}.
fidelity(#fock_dm{} = R, #fock_state{} = P) -> fidelity(P, R);
fidelity(A, B) ->
    case {dims(A), dims(B)} of
        {D, D} -> fid(A, B);
        {Da, Db} -> {error, {bad_shape, #{a => Da, b => Db}}}
    end.

dims(#fock_state{n_modes = N, cutoff = C}) -> {N, C};
dims(#fock_dm{n_modes = N, cutoff = C}) -> {N, C}.

fid(#fock_state{amps = A}, #fock_state{amps = B}) ->
    {ok, {cvec, X}} = pq_core:unpack(A),
    {ok, {cvec, Y}} = pq_core:unpack(B),
    {Re, Im} = lists:foldl(fun({{complex, Ar, Ai}, {complex, Br, Bi}}, {S, T}) -> {S + Ar * Br + Ai * Bi, T + Ar * Bi - Ai * Br} end,
                           {0.0, 0.0}, lists:zip(X, Y)),
    {ok, Re * Re + Im * Im};
fid(#fock_state{amps = A}, #fock_dm{rho = R}) ->
    {ok, {cvec, X}} = pq_core:unpack(A),
    {ok, {cmat, Rows}} = pq_core:unpack(R),
    %% sum_ij conj(x_i) rho_ij x_j
    Re = lists:sum([begin
                        {complex, Xr, Xi} = Xi0,
                        lists:sum([begin
                                       {complex, Rr, Ri} = Rij, {complex, Yr, Yi} = Xj,
                                       %% conj(x_i) * rho_ij * x_j, real part
                                       Pr = Rr * Yr - Ri * Yi, Pi = Rr * Yi + Ri * Yr,
                                       Xr * Pr + Xi * Pi
                                   end || {Rij, Xj} <- lists:zip(Row, X)])
                    end || {Xi0, Row} <- lists:zip(X, Rows)]),
    {ok, Re};
fid(#fock_dm{}, #fock_dm{}) -> {error, {unsupported_input, mixed_mixed_fidelity}}.
