%% Adjoint (reverse-mode) gradients on the pure-state Fock backend (Phase 4, FS 15 item 3; decision record
%% docs/decisions/AD1-adjoint-gradients.md). One forward and one backward pass in the kernel
%% (pq_nif:fock_adjoint_gradient/5) give the derivative of the observable with respect to every gate parameter;
%% the chain rule (pq_gradient_chain) turns those into derivatives with respect to the circuit's symbols.
%%
%% Gates: beamsplitter, phase shift, Kerr, numeric passive unitaries, displacement and squeezing (exact for the
%% truncated gates the backend applies). The state lives under the cutoff of the Fock backend (FS 3.4), so the
%% gradient is the one of the truncated model: for a vacuum or Gaussian-free circuit it equals the Gaussian
%% backend's within the truncated mass. With loss or thermal loss the same reverse pass runs on density matrices
%% (pq_nif:fock_dm_adjoint_gradient/5, within max_dm_dim): the loss transmissivity is differentiated, thermal loss
%% is passed through (its parameters refused). Mid-circuit conditioning and feed-forward are refused
%% (differentiate_unsupported); use mode numeric or parameter_shift there.
-module(pq_gradient_adjoint).
-behaviour(pq_gradient).
-include("pq.hrl").
-export([supports/1, gradient/4]).

%% the measurements are ignored, as by pq_gradient:value/3
supports(K) -> gate_kind(K) orelse lists:member(K, [measure_pnr, measure_threshold, measure_homodyne]).

gate_kind(K) -> lists:member(K, [beamsplitter, phaseshift, kerr, unitary, displacement, squeezing, loss, thermal_loss]).

-spec gradient(#circuit{}, pq_gradient:observable(), map(), map()) -> {ok, #{atom() => float()}} | {error, pq_core:reason()}.
gradient(C, Obs, Binding, Opts) ->
    Syms = maps:get(symbols, Opts, lists:sort(maps:keys(Binding))),
    case maps:get(backend, Opts, fock) of
        fock ->
            case pq_gradient_chain:prepare(C, Obs, Binding, Syms) of
                {ok, P} -> run(P, Opts);
                E -> E
            end;
        Other -> {error, {differentiate_unsupported, #{mode => adjoint, backend => Other}}}
    end.

run(#{lowered := #circuit{instrs = Is, n_modes = N} = Low, terms := Terms, symbols := Syms, coefs := Coefs} = P, Opts) ->
    Lossy = lists:any(fun(I) -> lists:member(element(1, I), [loss, thermal_loss]) end, Is),
    Thermal = [K || {K, {{thermal_loss, _, _, _}, Ds}} <- lists:enumerate(0, lists:zip(Is, Coefs)), lists:any(fun(D) -> map_size(D) > 0 end, Ds)],
    case {[element(1, I) || I <- Is, not gate_kind(element(1, I))], Thermal} of
        {[], [K | _]} -> {error, {differentiate_unsupported, #{mode => adjoint, instr_index => K, instruction => thermal_loss}}};
        {[], []} ->
            case pq_fock:init(Low#circuit.input, Low, Opts) of
                {ok, {#fock_state{amps = A, cutoff = Cut}, _}} ->
                    Gates = [gate(I) || I <- Is],
                    KTerms = [kterm(W, B) || {W, _, B} <- Terms],
                    D = pq_state:basis_size(N, Cut),
                    MaxDm = pq_core:get(max_dm_dim, Opts),
                    R = case Lossy of
                            false -> pq_nif:fock_adjoint_gradient(A, Gates, KTerms, N, Cut);
                            true when D > MaxDm -> {error, {dimension_limit, D, MaxDm}};
                            true -> pq_nif:fock_dm_adjoint_gradient(A, Gates, KTerms, N, Cut)
                        end,
                    case R of
                        {ok, {Values, Grads}} -> {ok, pq_gradient_chain:combine(P, Grads, Values, Syms)};
                        {error, _} = E -> E
                    end;
                {error, _} = E -> E
            end;
        {Kinds, _} -> {error, {differentiate_unsupported, #{mode => adjoint, instructions => lists:usort(Kinds)}}}
    end.

gate({beamsplitter, M1, M2, T, P}) -> {beamsplitter, T, P, M1, M2};
gate({phaseshift, M, P}) -> {phaseshift, P, M};
gate({kerr, M, K}) -> {kerr, K, M};
gate({unitary, Ms, U}) -> {passive, U, Ms};
gate({displacement, M, R, P}) -> {displacement, R, P, M};
gate({squeezing, M, R, P}) -> {squeezing, R, P, M};
gate({loss, M, Eta}) -> {loss, float(Eta), M};
gate({thermal_loss, M, Eta, Nth}) -> {thermal_loss, float(Eta), float(Nth), M}.

kterm(W, {probability, Modes, Pattern}) -> {probability, W, Modes, Pattern};
kterm(W, {mean_photon_number, Modes}) -> {photons, W, Modes}.
