%% Statistics (FS 5.12). Phase 1: mean_photon_number (FR-X01).
-module(pq_stats).
-include("pq.hrl").
-export([mean_photon_number/2]).

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
