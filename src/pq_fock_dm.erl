%% The density-matrix Fock backend (FS 5.5, FR-F05): rho is D x D in the basis of FS 3.4 (Decision 3), within
%% max_dm_dim. Gates act as G rho G^dagger; loss by the Kraus sum of A.11 (exact); displacement and squeezing
%% truncated, their cumulative trace deficit a warning or, under strict, an error (FS 3.5).
%%
%% dm_storage => blocks (Phase 4, FS 15 item 6; docs/decisions/B1-dm-blocks.md): inside the backend rho is kept as
%% photon-number blocks ({dm_blocks, Blocks} in the rho field), zero blocks not stored, within max_dm_entries
%% instead of max_dm_dim; the entries equal the dense ones bit for bit. Decision 3: the public state (final/1) is the
%% full matrix, within max_dm_dim.
-module(pq_fock_dm).
-behaviour(pq_backend).
-compile({no_auto_import, [apply/3]}).
-include("pq.hrl").
-export([name/0, capabilities/0, init/3, apply/3, sample/3, final/1, from_pure/1, measure/3]).
-export([promote/2, pnr_probabilities/2]).

name() -> fock_dm.

capabilities() ->
    #{beamsplitter => exact, phaseshift => exact, unitary => exact, displacement => truncated,
      squeezing => truncated, kerr => exact, loss => exact, thermal_loss => truncated, measure_pnr => midcircuit,
      measure_threshold => terminal, measure_homodyne => unsupported}.

init(Input, #circuit{n_modes = N} = Circ, Opts) ->
    case pq_fock:cutoff(Input, Opts) of
        {ok, Cut} ->
            D = pq_state:basis_size(N, Cut),
            Blocks = pq_core:get(dm_storage, Opts) =:= blocks,
            Max = case Blocks of
                      true -> pq_core:get(max_fock_dim, Opts);
                      false -> pq_core:get(max_dm_dim, Opts)
                  end,
            case D =< Max of
                false -> {error, {dimension_limit, D, Max}};
                true ->
                    case pq_fock:init(Input, Circ, Opts#{max_fock_dim => Max}) of
                        {ok, {#fock_state{} = S, Def}} -> promote_checked(S, Def, Opts);
                        E -> E
                    end
            end;
        E -> E
    end.

final({#fock_dm{rho = {dm_blocks, B}, n_modes = N, cutoff = C} = S, _Deficit}) ->
    D = pq_state:basis_size(N, C),
    Max = pq_core:get(max_dm_dim),
    case D =< Max of
        true ->
            {ok, Rho} = pq_nif:fock_dmb_to_dense(B, N, C),
            S#fock_dm{rho = Rho};
        false -> {error, {dimension_limit, D, Max}}
    end;
final({S, _Deficit}) -> S.

%% a pure state as the backend's density matrix: blocks under dm_storage => blocks, else the full matrix
-spec promote(#fock_state{}, map()) -> #fock_dm{}.
promote(#fock_state{n_modes = N, cutoff = C, amps = A} = S, Opts) ->
    case pq_core:get(dm_storage, Opts) of
        blocks ->
            {ok, B} = pq_nif:fock_dmb_from(A, N, C),
            #fock_dm{n_modes = N, cutoff = C, rho = {dm_blocks, B}};
        dense -> from_pure(S)
    end.

promote_checked(S, Def, Opts) ->
    Dm = promote(S, Opts),
    case entries_ok(Dm, Opts) of
        ok -> {ok, {Dm, Def}};
        E -> E
    end.

entries_ok(#fock_dm{rho = {dm_blocks, B}}, Opts) ->
    Entries = lists:sum([byte_size(Bin) div 16 || {_, _, Bin} <- B]),
    Max = pq_core:get(max_dm_entries, Opts),
    case Entries =< Max of
        true -> ok;
        false -> {error, {dimension_limit, #{dm_entries => Entries, max_dm_entries => Max}}}
    end;
entries_ok(_, _) -> ok.

%% FR-F03 on the backend's state, either storage
-spec pnr_probabilities(#fock_dm{}, [non_neg_integer()]) -> {ok, [{[non_neg_integer()], float()}]} | {error, term()}.
pnr_probabilities(#fock_dm{rho = {dm_blocks, B}, n_modes = N, cutoff = C}, Modes) -> pq_nif:fock_dmb_pnr_probabilities(B, Modes, N, C);
pnr_probabilities(#fock_dm{rho = R, n_modes = N, cutoff = C}, Modes) -> pq_nif:fock_pnr_probabilities(R, Modes, N, C).

%% |psi><psi|
-spec from_pure(#fock_state{}) -> #fock_dm{}.
from_pure(#fock_state{n_modes = N, cutoff = C, amps = A}) ->
    {ok, {cvec, V}} = pq_core:unpack(A),
    {ok, Rho} = pq_core:pack({cmat, [[mul(X, conj(Y)) || Y <- V] || X <- V]}),
    #fock_dm{n_modes = N, cutoff = C, rho = Rho}.

mul({complex, A, B}, {complex, C, D}) -> {complex, A * C - B * D, A * D + B * C}.
conj({complex, A, B}) -> {complex, A, -B}.

apply(I, {#fock_dm{n_modes = N, cutoff = C, rho = {dm_blocks, B}} = S, Def}, Opts) ->
    Op = case I of
             {beamsplitter, M1, M2, T, P} -> {beamsplitter, T, P, M1, M2};
             {phaseshift, M, P} -> {phaseshift, P, M};
             {kerr, M, K} -> {kerr, K, M};
             {unitary, Ms, U} -> {passive, U, Ms};
             {displacement, M, Rr, P} -> {displacement, Rr, P, M};
             {squeezing, M, Rr, P} -> {squeezing, Rr, P, M};
             {loss, M, Eta} -> {loss, float(Eta), M};
             {thermal_loss, M, Eta, Nth} -> {thermal_loss, float(Eta), float(Nth), M};
             _ -> none
         end,
    case Op of
        none -> {error, {unsupported_instruction, element(1, I), fock_dm}};
        _ ->
            case pq_nif:fock_dmb_apply(B, Op, N, C) of
                {ok, {B2, D, _Entries}} ->
                    S2 = S#fock_dm{rho = {dm_blocks, B2}},
                    %% the trace deficit is judged after the truncated operations only, as on the dense storage
                    Truncated = lists:member(element(1, I), [displacement, squeezing, thermal_loss]),
                    case {entries_ok(S2, Opts), Truncated} of
                        {ok, true} -> pq_fock:truncation(S2, Def + D, Opts);
                        {ok, false} -> {ok, {S2, Def}, []};
                        {E, _} -> E
                    end;
                {error, _} = E -> E
            end
    end;
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
