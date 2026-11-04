%% The density-matrix Fock backend (FS 5.5, FR-F05): rho is D x D in the basis of FS 3.4 (Decision 3), within
%% max_dm_dim. Gates act as G rho G^dagger; loss by the Kraus sum of A.11 (exact); displacement and squeezing
%% truncated, their cumulative trace deficit a warning or, under strict, an error (FS 3.5).
-module(pq_fock_dm).
-behaviour(pq_backend).
-compile({no_auto_import, [apply/3]}).
-include("pq.hrl").
-export([name/0, capabilities/0, init/3, apply/3, sample/3, final/1, from_pure/1, measure/3]).

name() -> fock_dm.

capabilities() ->
    #{beamsplitter => exact, phaseshift => exact, unitary => exact, displacement => truncated,
      squeezing => truncated, kerr => exact, loss => exact, thermal_loss => truncated, measure_pnr => midcircuit,
      measure_threshold => terminal, measure_homodyne => unsupported}.

init(Input, #circuit{n_modes = N} = Circ, Opts) ->
    case pq_fock:cutoff(Input, Opts) of
        {ok, Cut} ->
            D = pq_state:basis_size(N, Cut),
            Max = pq_core:get(max_dm_dim, Opts),
            case D =< Max of
                false -> {error, {dimension_limit, D, Max}};
                true ->
                    case pq_fock:init(Input, Circ, Opts#{max_fock_dim => Max}) of
                        {ok, {#fock_state{} = S, Def}} -> {ok, {from_pure(S), Def}};
                        E -> E
                    end
            end;
        E -> E
    end.

final({S, _Deficit}) -> S.

%% |psi><psi|
-spec from_pure(#fock_state{}) -> #fock_dm{}.
from_pure(#fock_state{n_modes = N, cutoff = C, amps = A}) ->
    {ok, {cvec, V}} = pq_core:unpack(A),
    {ok, Rho} = pq_core:pack({cmat, [[mul(X, conj(Y)) || Y <- V] || X <- V]}),
    #fock_dm{n_modes = N, cutoff = C, rho = Rho}.

mul({complex, A, B}, {complex, C, D}) -> {complex, A * C - B * D, A * D + B * C}.
conj({complex, A, B}) -> {complex, A, -B}.

apply(I, {#fock_dm{n_modes = N, cutoff = C, rho = R} = S, Def}, Opts) ->
    Res = case I of
              {beamsplitter, M1, M2, T, P} -> pq_nif:fock_dm_apply_beamsplitter(R, T, P, M1, M2, N, C);
              {phaseshift, M, P} -> pq_nif:fock_dm_apply_phaseshift(R, P, M, N, C);
              {kerr, M, K} -> pq_nif:fock_dm_apply_kerr(R, K, M, N, C);
              {unitary, Ms, U} -> pq_nif:fock_dm_apply_passive(R, U, Ms, N, C);
              {displacement, M, Rr, P} -> pq_nif:fock_dm_apply_displacement(R, Rr, P, M, N, C);
              {squeezing, M, Rr, P} -> pq_nif:fock_dm_apply_squeezing(R, Rr, P, M, N, C);
              {loss, M, Eta} -> pq_nif:fock_dm_apply_loss(R, float(Eta), M, N, C);
              {thermal_loss, M, Eta, Nth} -> pq_nif:fock_dm_apply_thermal_loss(R, float(Eta), float(Nth), M, N, C);
              _ -> {error, {unsupported_instruction, element(1, I), fock_dm}}
          end,
    case Res of
        {ok, {Rho, D}} -> pq_fock:truncation(S#fock_dm{rho = Rho}, Def + D, Opts);
        {ok, Rho} -> {ok, {S#fock_dm{rho = Rho}, Def}, []};
        {error, _} = E -> E
    end.

sample(Ms, St, Opts) -> pq_fock:sample(Ms, St, Opts).

measure(I, St, Opts) -> pq_fock:measure(I, St, Opts).
