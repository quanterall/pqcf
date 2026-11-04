%% The pure-state Fock backend (FS 5.5): amplitudes in the basis of FS 3.4 under a cutoff. Exact: beamsplitter,
%% phase shift, passive unitary (through the Clements mesh), Kerr. Truncated: displacement and squeezing, their
%% cumulative norm deficit a warning above truncation_tol, an error under strict (FS 3.5). Loss needs a density
%% matrix (FR-F05, Phase 2). Terminal PNR and threshold sampling; mid-circuit PNR in Phase 2.
-module(pq_fock).
-behaviour(pq_backend).
-compile({no_auto_import, [apply/3]}).
-include("pq.hrl").
-export([name/0, capabilities/0, capabilities/1, init/3, apply/3, sample/3, final/1, measure/3]).
%% shared with pq_fock_dm
-export([cutoff/2, truncation/3, basis/4]).

name() -> fock.

capabilities() ->
    #{beamsplitter => exact, phaseshift => exact, unitary => exact, displacement => truncated,
      squeezing => truncated, kerr => exact, loss => unsupported, thermal_loss => unsupported, measure_pnr => midcircuit,
      measure_threshold => terminal, measure_homodyne => unsupported}.

%% FR-F05: with promote_dm, loss promotes the state to a density matrix
capabilities(Opts) ->
    case maps:get(promote_dm, Opts, false) orelse maps:get(noise, Opts, none) =:= trajectory of
        true -> (capabilities())#{loss => exact, thermal_loss => truncated};
        false -> capabilities()
    end.

%% the backend state: the public state and the cumulative norm deficit
init(Input, #circuit{n_modes = N}, Opts) ->
    case cutoff(Input, Opts) of
        {ok, Cut} ->
            D = pq_state:basis_size(N, Cut),
            Max = pq_core:get(max_fock_dim, Opts),
            case D =< Max of
                false -> {error, {dimension_limit, D, Max}};
                true -> start(Input, N, Cut, D)
            end;
        E -> E
    end.

cutoff({fock_state, #fock_state{cutoff = C}}, _Opts) -> {ok, C};
cutoff(Input, Opts) ->
    Photons = case Input of {fock, Occ} -> lists:sum(Occ); _ -> 0 end,
    case maps:find(cutoff, Opts) of
        {ok, C} when is_integer(C), C >= Photons -> {ok, C};
        {ok, C} -> {error, {unsupported_input, #{cutoff => C, photons => Photons}}};
        error when Photons > 0 -> {ok, Photons};
        error -> {error, {unsupported_input, cutoff_required}}
    end.

start(vacuum, N, Cut, D) -> {ok, {basis(N, Cut, D, lists:duplicate(N, 0)), 0.0}};
start({fock, Occ}, N, Cut, D) -> {ok, {basis(N, Cut, D, Occ), 0.0}};
start({fock_state, #fock_state{n_modes = N} = S}, N, _Cut, _D) -> {ok, {S, 0.0}};
start(Other, _N, _Cut, _D) -> {error, {unsupported_input, Other}}.

basis(N, Cut, D, Occ) ->
    {ok, I} = pq_state:index(Occ, N, Cut),
    Amps = << <<(case K =:= I of true -> 1.0; false -> 0.0 end):64/float-little, 0.0:64/float-little>>
              || K <- lists:seq(0, D - 1) >>,
    #fock_state{n_modes = N, cutoff = Cut, amps = {cvec, D, Amps}}.

%% a state promoted to a density matrix is presented by pq_fock_dm (Decision 3: the full matrix)
final({#fock_dm{}, _} = St) -> pq_fock_dm:final(St);
final({S, _Deficit}) -> S.

%% FR-F05: loss promotes the pure state to a density matrix when promote_dm is set; from then on pq_fock_dm acts
%% FR-F05 trajectory mode (sampling only): a lost-photon count drawn from the shot's stream, the state kept pure
apply({loss, M, Eta}, {#fock_state{n_modes = N, cutoff = C, amps = A} = S, Def}, #{noise := trajectory, stream_ref := _} = Opts) ->
    case pq_nif:fock_apply_loss_trajectory(A, float(Eta), M, N, C, maps:with([stream_ref], Opts)) of
        {ok, {Amps, _Lost}} -> {ok, {S#fock_state{amps = Amps}, Def}, []};
        {error, _} = E -> E
    end;
apply({thermal_loss, _, _, _} = I, {#fock_state{} = S, Def}, Opts) ->
    case maps:get(promote_dm, Opts, false) of
        true -> pq_fock_dm:apply(I, {pq_fock_dm:promote(S, Opts), Def}, Opts);
        false -> {error, {unsupported_instruction, thermal_loss, fock}}
    end;
apply({loss, _, _} = I, {#fock_state{amps = A} = S, Def}, Opts) ->
    case maps:get(promote_dm, Opts, false) of
        true -> pq_fock_dm:apply(I, {pq_fock_dm:promote(S#fock_state{amps = A}, Opts), Def}, Opts);
        false -> {error, {unsupported_instruction, loss, fock}}
    end;
apply(I, {#fock_dm{}, _} = St, Opts) -> pq_fock_dm:apply(I, St, Opts);
apply(I, {#fock_state{n_modes = N, cutoff = C, amps = A} = S, Def}, Opts) ->
    R = case I of
            {beamsplitter, M1, M2, T, P} -> pq_nif:fock_apply_beamsplitter(A, T, P, M1, M2, N, C);
            {phaseshift, M, P} -> pq_nif:fock_apply_phaseshift(A, P, M, N, C);
            {kerr, M, K} -> pq_nif:fock_apply_kerr(A, K, M, N, C);
            {unitary, Ms, U} -> pq_nif:fock_apply_passive(A, U, Ms, N, C);
            {displacement, M, Rr, P} -> pq_nif:fock_apply_displacement(A, Rr, P, M, N, C);
            {squeezing, M, Rr, P} -> pq_nif:fock_apply_squeezing(A, Rr, P, M, N, C);
            _ -> {error, {unsupported_instruction, element(1, I), fock}}
        end,
    case R of
        {ok, {Amps, D}} -> truncation(S#fock_state{amps = Amps}, Def + D, Opts);
        {ok, Amps} -> {ok, {S#fock_state{amps = Amps}, Def}, []};
        {error, _} = E -> E
    end.

truncation(S, Def, Opts) ->
    Tol = pq_core:get(truncation_tol, Opts),
    Info = #{deficit => Def, instr_index => maps:get(instr_index, Opts, undefined)},
    case Def > Tol of
        false -> {ok, {S, Def}, []};
        true ->
            case pq_core:get(strict, Opts) of
                true -> {error, {truncation, Info}};
                false -> {ok, {S, Def}, [{truncation, Info}]}
            end
    end.

%% terminal measurements: the joint PNR distribution of every measured mode, sampled (one uniform per shot);
%% threshold columns report n > 0
sample(Ms, {#fock_dm{} = S, _}, Opts) -> sample_probs(Ms, fun(Modes) -> pq_fock_dm:pnr_probabilities(S, Modes) end, S#fock_dm.cutoff, Opts);
sample(Ms, {#fock_state{n_modes = N, cutoff = C, amps = A}, _}, Opts) -> sample_state(Ms, A, N, C, Opts).

sample_state(Ms, A, N, C, Opts) ->
    sample_probs(Ms, fun(Modes) -> pq_nif:fock_pnr_probabilities(A, Modes, N, C) end, C, Opts).

sample_probs(Ms, ProbsOf, C, Opts) ->
    Cols = lists:append([[{M, kind(K)} || M <- Modes] || {K, Modes} <- Ms]),
    Modes = [M || {M, _} <- Cols],
    case ProbsOf(Modes) of
        {ok, PatProbs} ->
            {Pats, Probs} = lists:unzip(PatProbs),
            {ok, P} = pq_core:pack({rvec, Probs}),
            {ok, Idx} = pq_nif:sample_categorical(P, maps:get(shots, Opts, 1), maps:with([seed, stream_start, stream_ref], Opts)),
            PatT = list_to_tuple(Pats),
            Kinds = [Kd || {_, Kd} <- Cols],
            Samples = [[case Kd of threshold -> min(X, 1); pnr -> X end
                        || {X, Kd} <- lists:zip(element(I + 1, PatT), Kinds)] || I <- Idx],
            Kind = case lists:usort(Kinds) of [One] -> One; _ -> mixed end,
            Extra = case Kind of mixed -> #{column_kinds => Kinds}; _ -> #{} end,
            {ok, Kind, Modes, Samples, Extra#{cutoff => C}};
        {error, _} = E -> E
    end.

%% FR-M07: PNR on Modes of one shot: the outcome drawn from the shot's stream (one uniform) or fixed, the state
%% projected onto it and renormalised (the measured modes keep their outcome; nothing is renumbered)
measure({measure_pnr, Modes}, {#fock_dm{rho = {dm_blocks, B}, n_modes = N, cutoff = C} = S, Def}, Opts) ->
    %% photon-number blocks (FS 15 item 6): the same draw, the projection on the blocks
    Out = case maps:find(fix_outcome, Opts) of
              {ok, O} -> O;
              error ->
                  {ok, PP} = pq_fock_dm:pnr_probabilities(S, Modes),
                  {Pats, Probs} = lists:unzip(PP),
                  {ok, P} = pq_core:pack({rvec, Probs}),
                  {ok, [I]} = pq_nif:sample_categorical(P, 1, maps:with([stream_ref], Opts)),
                  lists:nth(I + 1, Pats)
          end,
    case pq_nif:fock_dmb_measure_pnr(B, Modes, Out, N, C) of
        {ok, {Prob, _}} when Prob =< 0.0 -> {error, {postselection_exhausted, #{fixed_outcome_probability => 0.0}}};
        {ok, {_, B2}} -> {ok, Out, {S#fock_dm{rho = {dm_blocks, B2}}, Def}};
        {error, _} = E -> E
    end;
measure({measure_pnr, Modes}, {S, Def}, Opts) ->
    {Kind, N, C, State} = case S of
                              #fock_state{n_modes = N0, cutoff = C0, amps = A} -> {pure, N0, C0, A};
                              #fock_dm{n_modes = N0, cutoff = C0, rho = R} -> {dm, N0, C0, R}
                          end,
    Outcome = case maps:find(fix_outcome, Opts) of
                  {ok, O} -> {ok, O};
                  error ->
                      {ok, PP} = pq_nif:fock_pnr_probabilities(State, Modes, N, C),
                      {Pats, Probs} = lists:unzip(PP),
                      {ok, P} = pq_core:pack({rvec, Probs}),
                      {ok, [I]} = pq_nif:sample_categorical(P, 1, maps:with([stream_ref], Opts)),
                      {ok, lists:nth(I + 1, Pats)}
              end,
    {ok, Out} = Outcome,
    case pq_nif:fock_measure_pnr(State, Modes, Out, N, C) of
        {ok, {Prob, _}} when Prob =< 0.0 -> {error, {postselection_exhausted, #{fixed_outcome_probability => 0.0}}};
        {ok, {_, New}} when Kind =:= pure -> {ok, Out, {S#fock_state{amps = New}, Def}};
        {ok, {_, New}} -> {ok, Out, {S#fock_dm{rho = New}, Def}};
        {error, _} = E -> E
    end;
measure(I, _S, _Opts) -> {error, {unsupported_instruction, element(1, I), #{backend => fock, placement => midcircuit}}}.

kind(measure_pnr) -> pnr;
kind(measure_threshold) -> threshold.
