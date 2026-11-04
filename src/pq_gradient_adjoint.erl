%% Adjoint (reverse-mode) gradients on the pure-state Fock backend (Phase 4, FS 15 item 3; decision record
%% docs/decisions/AD1-adjoint-gradients.md). One forward and one backward pass in the kernel
%% (pq_nif:fock_adjoint_gradient/5) give the derivative of the observable with respect to every gate parameter;
%% the chain rule (pq_gradient_chain) turns those into derivatives with respect to the circuit's symbols.
%%
%% Gates: beamsplitter, phase shift, Kerr, numeric passive unitaries, displacement and squeezing (exact for the
%% truncated gates the backend applies). The state lives under the cutoff of the Fock backend (FS 3.4), so the
%% gradient is the one of the truncated model: for a vacuum or Gaussian-free circuit it equals the Gaussian
%% backend's within the truncated mass. Loss, mid-circuit conditioning and feed-forward are refused
%% (differentiate_unsupported); use mode numeric or parameter_shift there.
-module(pq_gradient_adjoint).
-behaviour(pq_gradient).
-include("pq.hrl").
-export([supports/1, gradient/4]).

%% the measurements are ignored, as by pq_gradient:value/3
supports(K) -> gate_kind(K) orelse lists:member(K, [measure_pnr, measure_threshold, measure_homodyne]).

gate_kind(K) -> lists:member(K, [beamsplitter, phaseshift, kerr, unitary, displacement, squeezing]).

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

run(#{lowered := #circuit{instrs = Is, n_modes = N} = Low, terms := Terms, symbols := Syms} = P, Opts) ->
    case [element(1, I) || I <- Is, not gate_kind(element(1, I))] of
        [] ->
            case pq_fock:init(Low#circuit.input, Low, Opts) of
                {ok, {#fock_state{amps = A, cutoff = Cut}, _}} ->
                    Gates = [gate(I) || I <- Is],
                    KTerms = [kterm(W, B) || {W, _, B} <- Terms],
                    case pq_nif:fock_adjoint_gradient(A, Gates, KTerms, N, Cut) of
                        {ok, {Values, Grads}} -> {ok, pq_gradient_chain:combine(P, Grads, Values, Syms)};
                        {error, _} = E -> E
                    end;
                {error, _} = E -> E
            end;
        Kinds -> {error, {differentiate_unsupported, #{mode => adjoint, instructions => lists:usort(Kinds)}}}
    end.

gate({beamsplitter, M1, M2, T, P}) -> {beamsplitter, T, P, M1, M2};
gate({phaseshift, M, P}) -> {phaseshift, P, M};
gate({kerr, M, K}) -> {kerr, K, M};
gate({unitary, Ms, U}) -> {passive, U, Ms};
gate({displacement, M, R, P}) -> {displacement, R, P, M};
gate({squeezing, M, R, P}) -> {squeezing, R, P, M}.

kterm(W, {probability, Modes, Pattern}) -> {probability, W, Modes, Pattern};
kterm(W, {mean_photon_number, Modes}) -> {photons, W, Modes}.
