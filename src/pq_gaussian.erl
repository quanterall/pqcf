%% The Gaussian backend (FS 5.4): mean and covariance, xxpp, maps of A.2-A.6 applied one instruction at a time.
%% Phase 1: terminal homodyne sampling (FR-G05). Photon counting (5.7) and mid-circuit homodyne (FR-G08) come in
%% Phase 2.
-module(pq_gaussian).
-behaviour(pq_backend).
-compile({no_auto_import, [apply/3]}).
-include("pq.hrl").
-export([name/0, capabilities/0, init/3, apply/3, sample/3]).

name() -> gaussian.

capabilities() ->
    #{beamsplitter => exact, phaseshift => exact, unitary => exact, displacement => exact, squeezing => exact,
      kerr => unsupported, loss => exact, measure_pnr => terminal, measure_threshold => terminal,
      measure_homodyne => midcircuit}.

hbar(#circuit{hbar = undefined}) -> pq_core:get(hbar);
hbar(#circuit{hbar = H}) -> H.

init(vacuum, #circuit{n_modes = N} = C, _Opts) -> {ok, pq_state:vacuum(N, hbar(C))};
init({gaussian, #gaussian_state{} = S}, _C, Opts) ->
    case maps:get(check_physical, Opts, true) of
        true ->
            case pq_state:check_physical(S) of
                ok -> {ok, S};
                E -> E
            end;
        false -> {ok, S}
    end;
init({fock, _}, _C, _Opts) -> {error, {unsupported_input, fock}};
init(Other, _C, _Opts) -> {error, {unsupported_input, Other}}.

apply({beamsplitter, A, B, T, P}, S, _Opts) ->
    passive(S, [[{complex, math:cos(T), 0.0}, {complex, -math:cos(-P) * math:sin(T), -math:sin(-P) * math:sin(T)}],
                [{complex, math:cos(P) * math:sin(T), math:sin(P) * math:sin(T)}, {complex, math:cos(T), 0.0}]], [A, B]);
apply({phaseshift, M, P}, S, _Opts) -> passive(S, [[{complex, math:cos(P), math:sin(P)}]], [M]);
apply({unitary, Ms, U}, #gaussian_state{mean = Mu, cov = C} = S, _Opts) ->
    result(S, pq_nif:gaussian_apply_passive(Mu, C, U, Ms));
apply({displacement, M, R, P}, #gaussian_state{hbar = H} = S, _Opts) ->
    K = math:sqrt(2 * H) * R,
    symplectic(S, [[1, 0], [0, 1]], [K * math:cos(P), K * math:sin(P)], [M]);
apply({squeezing, M, R, P}, S, _Opts) ->
    C = math:cosh(R), Sh = math:sinh(R),
    symplectic(S, [[C - math:cos(P) * Sh, -math:sin(P) * Sh], [-math:sin(P) * Sh, C + math:cos(P) * Sh]], [0, 0], [M]);
apply({loss, M, Eta}, #gaussian_state{hbar = H, mean = Mu, cov = C} = S, _Opts) ->
    result(S, pq_nif:gaussian_apply_loss(Mu, C, M, float(Eta), float(H)));
apply(I, _S, _Opts) -> {error, {unsupported_instruction, element(1, I), gaussian}}.

passive(#gaussian_state{mean = Mu, cov = C} = S, Rows, Modes) ->
    {ok, U} = pq_core:pack({cmat, Rows}),
    result(S, pq_nif:gaussian_apply_passive(Mu, C, U, Modes)).

symplectic(#gaussian_state{mean = Mu, cov = C} = S, Rows, D, Modes) ->
    {ok, Sm} = pq_core:pack({rmat, Rows}),
    {ok, Dv} = pq_core:pack({rvec, D}),
    result(S, pq_nif:gaussian_apply_symplectic(Mu, C, Sm, Dv, Modes)).

result(S, {ok, {Mu, C}}) -> {ok, S#gaussian_state{mean = Mu, cov = C}, []};
result(_S, {error, _} = E) -> E.

%% terminal homodyne on every listed mode: each mode rotated by -Phi into x, the measured modes' marginal sampled
sample(Ms, S, Opts) ->
    case lists:all(fun(I) -> element(1, I) =:= measure_homodyne end, Ms) of
        false -> {error, {unsupported_instruction, element(1, hd([I || I <- Ms, element(1, I) =/= measure_homodyne])),
                          #{backend => gaussian, available => phase_2}}};
        true ->
            Modes = [M || {measure_homodyne, M, _} <- Ms],
            Rot = lists:foldl(fun({measure_homodyne, M, P}, {ok, St, _}) -> apply({phaseshift, M, -P}, St, Opts);
                                 (_, E) -> E
                              end, {ok, S, []}, Ms),
            case Rot of
                {ok, St, _} ->
                    {ok, Red} = pq_state:reduce(St, Modes),
                    Shots = maps:get(shots, Opts, 1),
                    NifOpts = maps:with([seed, stream_start, cancel], Opts),
                    case pq_nif:sample_gaussian(Red#gaussian_state.mean, Red#gaussian_state.cov, Shots, {homodyne, 0.0},
                                                float(Red#gaussian_state.hbar), NifOpts) of
                        {ok, {Packed, _}} ->
                            {ok, {rmat, Rows}} = pq_core:unpack(Packed),
                            {ok, homodyne, Modes, Rows, #{}};
                        {error, _} = E -> E
                    end;
                E -> E
            end
    end.
