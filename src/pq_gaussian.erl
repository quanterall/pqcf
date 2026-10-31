%% The Gaussian backend (FS 5.4): mean and covariance, xxpp, maps of A.2-A.6 applied one instruction at a time.
%% Phase 1: terminal homodyne sampling (FR-G05). Photon counting (5.7) and mid-circuit homodyne (FR-G08) come in
%% Phase 2.
-module(pq_gaussian).
-behaviour(pq_backend).
-compile({no_auto_import, [apply/3]}).
-include("pq.hrl").
-export([name/0, capabilities/0, init/3, apply/3, sample/3, measure/3]).
%% FR-G04: symbolic evolution and the conversions between symbolic and numeric states
-export([symbolic_evolve/2, to_symbolic/1, lower/2]).

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

%%% FR-G04, SY 6.2: symbolic states hold mean and covariance as canonical matrices (2m x 1 and 2m x 2m, xxpp)

%% The state before any measurement of a Gaussian circuit with symbolic parameters, from vacuum or a Gaussian input
%% (numeric or symbolic); n_modes at most max_symbolic_modes. Terminal measurements are ignored; sampling needs
%% lower/2 first.
-spec symbolic_evolve(#circuit{}, map()) -> {ok, #gaussian_state{}} | {error, pq_core:reason()}.
symbolic_evolve(#circuit{n_modes = N, input = In, instrs = Is} = C, Opts) ->
    Measurement = fun(I) -> lists:member(element(1, I), [measure_pnr, measure_threshold, measure_homodyne]) end,
    {Body, Tail} = lists:splitwith(fun(I) -> not Measurement(I) end, Is),
    H = hbar(C),
    Start = case In of
                vacuum -> {ok, lists:duplicate(2 * N, 0),
                           [[case I =:= J of true -> {mul, [H, {rat, 1, 2}]}; false -> 0 end || J <- lists:seq(1, 2 * N)]
                            || I <- lists:seq(1, 2 * N)]};
                {gaussian, #gaussian_state{mean = {matrix, M}, cov = {matrix, Cv}}} -> {ok, [X || [X] <- M], Cv};
                {gaussian, #gaussian_state{} = S0} -> {Mv, Cv} = pq_state:mean_cov(S0), {ok, Mv, Cv};
                Other -> {error, {unsupported_input, Other}}
            end,
    case {Start, lists:all(Measurement, Tail)} of
        {_, false} -> {error, {unsupported_instruction, {midcircuit_measurement, symbolic_gaussian}}};
        {{error, _} = E, _} -> E;
        {{ok, Mean0, Cov0}, true} ->
            case pq_symbolic:gaussian_evolve(C#circuit{instrs = Body}, Opts#{mean => Mean0, cov => Cov0, hbar => H}) of
                {ok, {Mean, Cov}} -> {ok, #gaussian_state{n_modes = N, hbar = H, mean = Mean, cov = Cov}};
                E -> E
            end
    end.

%% a numeric state as a symbolic one (its entries become float constants)
-spec to_symbolic(#gaussian_state{}) -> #gaussian_state{}.
to_symbolic(#gaussian_state{mean = {matrix, _}} = S) -> S;
to_symbolic(#gaussian_state{} = S) ->
    {Mv, Cv} = pq_state:mean_cov(S),
    S#gaussian_state{mean = {matrix, [[X] || X <- Mv]}, cov = {matrix, Cv}}.

%% a symbolic state evaluated under Binding: numeric, ready for sampling; complex entries are an error
-spec lower(#gaussian_state{}, map()) -> {ok, #gaussian_state{}} | {error, pq_core:reason()}.
lower(#gaussian_state{mean = {matrix, _} = M, cov = {matrix, _} = Cv, hbar = H}, Binding) ->
    case {pq_symbolic:evaluate(M, Binding), pq_symbolic:evaluate(Cv, Binding)} of
        {{ok, Mp}, {ok, Cp}} ->
            {ok, {cmat, MRows}} = pq_core:unpack(Mp),
            {ok, {cmat, CRows}} = pq_core:unpack(Cp),
            try
                {ok, pq_state:from_mean_cov([real_entry(X) || [X] <- MRows], [[real_entry(X) || X <- R] || R <- CRows], H)}
            catch
                throw:{non_physical_state, _} = R -> {error, R}
            end;
        {{error, _} = E, _} -> E;
        {_, {error, _} = E} -> E
    end;
lower(#gaussian_state{} = S, _Binding) -> {ok, S}.

real_entry({complex, Re, Im}) when abs(Im) =< 1.0e-12 * (1 + abs(Re)) -> Re;
real_entry({complex, _, _} = Z) -> throw({non_physical_state, #{complex_entry => Z}}).

%% FR-G08: homodyne on one mode of one shot, the rest conditioned (A.7); the measured mode is then reset to vacuum
%% in place, so no mode is renumbered (no later instruction may use it, FR-M06)
measure({measure_homodyne, M, Phi}, #gaussian_state{n_modes = N, hbar = H, mean = Mu, cov = C} = S, Opts) when N >= 2 ->
    NifOpts = maps:merge(maps:with([stream_ref], Opts),
                         case maps:find(fix_outcome, Opts) of {ok, Y0} -> #{fix_outcome => float(Y0)}; error -> #{} end),
    case pq_nif:gaussian_measure_homodyne(Mu, C, M, float(Phi), float(H), NifOpts) of
        {ok, {Y, Mu2, C2}} ->
            {ok, {rvec, Mv}} = pq_core:unpack(Mu2),
            {ok, {rmat, Cv}} = pq_core:unpack(C2),
            {ok, [Y], reinsert(S, M, Mv, Cv)};
        {error, _} = E -> E
    end;
measure({measure_homodyne, _, _}, _S, _Opts) -> {error, {unsupported_input, homodyne_on_the_last_mode}};
measure(I, _S, _Opts) -> {error, {unsupported_instruction, element(1, I), #{backend => gaussian, placement => midcircuit}}}.

%% put back mode M as vacuum, uncorrelated, into a state of N - 1 modes (xxpp)
reinsert(#gaussian_state{n_modes = N, hbar = H} = S, M, Mv, Cv) ->
    Old = fun(I) -> case I rem N of M -> none; J -> (I div N) * (N - 1) + (case J > M of true -> J - 1; false -> J end) end end,
    Mean = [case Old(I) of none -> 0.0; K -> lists:nth(K + 1, Mv) end || I <- lists:seq(0, 2 * N - 1)],
    Cov = [[case {Old(I), Old(J)} of
                {none, _} when I =:= J -> H / 2;
                {none, _} -> 0.0;
                {_, none} -> 0.0;
                {A, B} -> lists:nth(B + 1, lists:nth(A + 1, Cv))
            end || J <- lists:seq(0, 2 * N - 1)] || I <- lists:seq(0, 2 * N - 1)],
    New = pq_state:from_mean_cov(Mean, Cov, H),
    S#gaussian_state{mean = New#gaussian_state.mean, cov = New#gaussian_state.cov}.

%% terminal homodyne on every listed mode: each mode rotated by -Phi into x, the measured modes' marginal sampled
sample(Ms, S, Opts) ->
    Kinds = lists:usort([element(1, I) || I <- Ms]),
    case Kinds of
        [measure_homodyne] -> sample_homodyne(Ms, S, Opts);
        [K] when K =:= measure_pnr; K =:= measure_threshold -> sample_counting(K, Ms, S, Opts);
        _ -> {error, {unsupported_instruction, mixed_detectors, #{backend => gaussian}}}
    end.

%% photon counting (FS 5.7): the measured modes' marginal, sampled by the chain rule; truncated_mass recorded
sample_counting(K, Ms, S, Opts) ->
    Modes = lists:append([Mo || {_, Mo} <- Ms]),
    {ok, Red} = pq_state:reduce(S, Modes),
    Det = case K of measure_pnr -> pnr; measure_threshold -> threshold end,
    NifOpts = maps:merge(maps:with([seed, stream_start, cancel, stream_ref], Opts),
                         #{max_photons => pq_core:get(max_photons, Opts),
                           max_per_mode => maps:get(max_per_mode, Opts, pq_core:get(max_photons, Opts)),
                           max_clicks => pq_core:get(max_clicks, Opts)}),
    case pq_nif:sample_gaussian(Red#gaussian_state.mean, Red#gaussian_state.cov, maps:get(shots, Opts, 1), Det,
                                float(Red#gaussian_state.hbar), NifOpts) of
        {ok, {Packed, #{truncated_mass := T}}} ->
            {ok, {imat, Rows}} = pq_core:unpack(Packed),
            Tol = pq_core:get(truncation_tol, Opts),
            W = [{gbs_truncation, #{truncated_mass => T}} || T > Tol],
            case W =/= [] andalso pq_core:get(strict, Opts) of
                true -> {error, {gbs_truncation, #{truncated_mass => T}}};
                false -> {ok, Det, Modes, Rows, #{truncated_mass => T, extra_warnings => W}}
            end;
        {error, _} = E -> E
    end.

sample_homodyne(Ms, S, Opts) ->
    case true of
        true ->
            Modes = [M || {measure_homodyne, M, _} <- Ms],
            Rot = lists:foldl(fun({measure_homodyne, M, P}, {ok, St, _}) -> apply({phaseshift, M, -P}, St, Opts);
                                 (_, E) -> E
                              end, {ok, S, []}, Ms),
            case Rot of
                {ok, St, _} ->
                    {ok, Red} = pq_state:reduce(St, Modes),
                    Shots = maps:get(shots, Opts, 1),
                    NifOpts = maps:with([seed, stream_start, cancel, stream_ref], Opts),
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
