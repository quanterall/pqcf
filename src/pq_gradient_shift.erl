%% Parameter-shift rules (Phase 4, FS 15 item 3; decision record docs/decisions/AD1-adjoint-gradients.md): exact
%% derivatives from evaluations of the circuit at shifted parameter values, gate by gate, with the chain rule of
%% pq_gradient_chain for the circuit's symbols.
%%
%% A gate parameter p enters the observable as a trigonometric polynomial f(p) = sum_(k <= R) a_k cos(k w p) +
%% b_k sin(k w p) when the photon number is bounded by P (the Fock backends' cutoff, or the photons of a Fock input).
%% The equidistant rule (Wierichs, Izaac, Wang and Lin, Quantum 6, 677 (2022), eq. 34) is then exact:
%%   f'(p) = w sum_(mu = 1..2R) f(p + x_mu / w) (-1)^(mu - 1) / (4 R sin^2(x_mu / 2)),  x_mu = (2 mu - 1) pi / (2 R).
%% Frequencies, from the spectra of the generators (fock.rs):
%%   phase shift (w 1, R P)        e^(i phi n), n <= P
%%   beamsplitter theta (w 1, R 2P) eigenvalues N, N-2, ..., -N of -i G in the sector N <= P; later gates couple
%%                                  sectors N and N' of the two modes, so lambda - lambda' takes every integer
%%                                  up to N + N' (odd ones too: measured on a three-mode circuit)
%%   beamsplitter phi (w 1, R 2P)   U[j][k] e^(i phi (k - j)), |k - j| <= N, in the ket and the bra
%%   Kerr (w 1, R P^2)             e^(i kappa n^2)
%%   displacement phi (w 1, R 2P), squeezing phi (w 1, R P)   <m|G|n> e^(i phi (m - n)) and e^(i phi (m - n)/2)
%% The magnitudes r of displacement and squeezing and the transmissivity of loss are not periodic: refused
%% (differentiate_unsupported); so are Gaussian-backend circuits, whose photon number is unbounded.
-module(pq_gradient_shift).
-behaviour(pq_gradient).
-include("pq.hrl").
-export([supports/1, gradient/4, rule/2]).

%% the kinds that may appear in the circuit; which of their parameters may be shifted is checked per parameter
supports(K) -> lists:member(K, [beamsplitter, phaseshift, kerr, unitary, displacement, squeezing, loss,
                                measure_pnr, measure_threshold, measure_homodyne]).

-spec gradient(#circuit{}, pq_gradient:observable(), map(), map()) -> {ok, #{atom() => float()}} | {error, pq_core:reason()}.
gradient(C, Obs, Binding, Opts0) ->
    Syms = maps:get(symbols, Opts0, lists:sort(maps:keys(Binding))),
    Opts = maps:without([bindings, mode, strategy, symbols], Opts0),
    try
        P = ok(pq_gradient_chain:prepare(C, Obs, Binding, Syms)),
        #{lowered := Low, coefs := Coefs, terms := Terms} = P,
        Bound = photon_bound(Low, Opts),
        NumObs = {weighted_sum, [{W, B} || {W, _, B} <- Terms]},
        F = fun(Circ) -> ok(pq_gradient:value(Circ, NumObs, Opts)) end,
        Grads = [[case maps:size(D) of
                      0 -> 0.0;
                      _ -> shift_derivative(F, Low, K, J, I, Bound)
                  end || {J, D} <- lists:enumerate(0, Ds)]
                 || {K, {I, Ds}} <- lists:enumerate(0, lists:zip(Low#circuit.instrs, Coefs))],
        Values = case lists:any(fun({_, DW, _}) -> map_size(DW) > 0 end, Terms) of
                     true -> [ok(pq_gradient:value(Low, B, Opts)) || {_, _, B} <- Terms];
                     false -> [0.0 || _ <- Terms]
                 end,
        {ok, pq_gradient_chain:combine(P, Grads, Values, Syms)}
    catch
        throw:{error, _} = E -> E
    end.

%% the largest photon number the state can hold
photon_bound(#circuit{input = In} = Low, Opts) ->
    case pq_api:backend_for(Low, Opts) of
        gaussian -> fail({differentiate_unsupported, #{mode => parameter_shift, backend => gaussian,
                                                       reason => <<"the photon number is unbounded">>}});
        _ -> ok(pq_fock:cutoff(In, Opts))
    end.

%% {W, R} for parameter J of instruction I under the photon bound P
-spec rule({atom(), non_neg_integer()}, non_neg_integer()) -> {pos_integer(), non_neg_integer()} | none.
rule({phaseshift, 0}, P) -> {1, P};
rule({beamsplitter, 0}, P) -> {1, 2 * P};
rule({beamsplitter, 1}, P) -> {1, 2 * P};
rule({kerr, 0}, P) -> {1, P * P};
rule({displacement, 1}, P) -> {1, 2 * P};
rule({squeezing, 1}, P) -> {1, P};
rule(_, _) -> none.

shift_derivative(F, Low, K, J, I, Bound) ->
    case rule({element(1, I), J}, Bound) of
        none -> fail({differentiate_unsupported, #{mode => parameter_shift, instr_index => K, instruction => element(1, I),
                                                   parameter => J, reason => <<"not periodic">>}});
        {_, 0} -> 0.0;
        {W, R} ->
            P0 = lists:nth(J + 1, pq_gradient_chain:gate_params(I)),
            lists:sum([begin
                           X = (2 * Mu - 1) * math:pi() / (2 * R),
                           Sign = case Mu rem 2 of 1 -> 1; 0 -> -1 end,
                           Circ = Low#circuit{instrs = set_nth(K, set_param(I, J, P0 + X / W), Low#circuit.instrs)},
                           W * Sign * F(Circ) / (4 * R * math:pow(math:sin(X / 2), 2))
                       end || Mu <- lists:seq(1, 2 * R)])
    end.

set_param({beamsplitter, A, B, _, P}, 0, V) -> {beamsplitter, A, B, V, P};
set_param({beamsplitter, A, B, T, _}, 1, V) -> {beamsplitter, A, B, T, V};
set_param({phaseshift, M, _}, 0, V) -> {phaseshift, M, V};
set_param({kerr, M, _}, 0, V) -> {kerr, M, V};
set_param({displacement, M, R, _}, 1, V) -> {displacement, M, R, V};
set_param({squeezing, M, R, _}, 1, V) -> {squeezing, M, R, V}.

set_nth(K, X, L) ->
    {A, [_ | B]} = lists:split(K, L),
    A ++ [X | B].

ok({ok, V}) -> V;
ok({error, _} = E) -> throw(E).

-spec fail(term()) -> no_return().
fail(R) -> throw({error, R}).
