%% Convention converters (FS 3.3, Decision 1, Phase 3): circuits to and from Perceval and Piquasso, as plain data
%% (component lists) that a script turns into the library's objects. The conventions are measured and asserted in
%% test/oracle/converter_conventions.py (perceval 1.3.1, piquasso 8.0.2):
%%
%%   Perceval  BS.Rx(T, tl, bl, tr, br) = diag(e^(i tr), e^(i br)) B(T/2, pi/2) diag(e^(i tl), e^(i bl)), on adjacent
%%             modes; our B(theta, phi) = BS.Rx(2 theta, tl = tr = 0, bl = pi/2 - phi, br = phi - pi/2); PS(phi) =
%%             e^(i phi). Passive circuits only; a beamsplitter between non-adjacent modes becomes a Unitary over
%%             the span.
%%   Piquasso  Beamsplitter, Phaseshifter, Squeezing, Displacement, Kerr equal ours; its covariance is twice ours.
%%
%% Circuits are lowered first (numeric parameters).
-module(pq_convert).
-include("pq.hrl").
-export([to_perceval/2, from_perceval/1, to_piquasso/2, from_piquasso/1, perceval_unitary/2]).

-type component() :: {atom(), [non_neg_integer()], map()}.

%%% Perceval

-spec to_perceval(#circuit{}, map()) -> {ok, #{n_modes := pos_integer(), input := term(), components := [component()]}}
                                        | {error, pq_core:reason()}.
to_perceval(C, Binding) ->
    case pq_circuit:lower(C, Binding) of
        {ok, #circuit{n_modes = N, input = In, instrs = Is}} ->
            Body = [I || I <- Is, not measurement(I)],
            case [element(1, I) || I <- Body, not lists:member(element(1, I), [beamsplitter, phaseshift, unitary])] of
                [] -> {ok, #{n_modes => N, input => perceval_input(In), components => lists:append([pcvl(I, N) || I <- Body])}};
                [K | _] -> {error, {unsupported_instruction, K, perceval}}
            end;
        E -> E
    end.

perceval_input({fock, Occ}) -> Occ;
perceval_input(vacuum) -> vacuum;
perceval_input(_) -> unsupported.

pcvl({beamsplitter, A, B, T, P}, _N) when B =:= A + 1 ->
    [{'BS.Rx', [A, B], #{theta => 2 * T, phi_tl => 0.0, phi_tr => 0.0, phi_bl => math:pi() / 2 - P, phi_br => P - math:pi() / 2}}];
pcvl({beamsplitter, A, B, T, P}, _N) ->      % not on (m, m+1): its 2 x 2 matrix as a Unitary over the span
    span_unitary([A, B], [[cos(T), cneg(cexp(-P, math:sin(T)))], [cexp(P, math:sin(T)), cos(T)]]);
pcvl({phaseshift, M, P}, _N) -> [{'PS', [M], #{phi => P}}];
pcvl({unitary, Ms, Packed}, _N) ->
    {ok, {cmat, Rows}} = pq_core:unpack(Packed),
    span_unitary(Ms, [[{Re, Im} || {complex, Re, Im} <- R] || R <- Rows]).

cos(T) -> {math:cos(T), 0.0}.
cexp(P, S) -> {math:cos(P) * S, math:sin(P) * S}.
cneg({A, B}) -> {-A, -B}.

%% a k x k unitary on arbitrary modes as a Unitary on the consecutive span min..max (identity elsewhere)
span_unitary(Ms, G) ->
    Lo = lists:min(Ms), Hi = lists:max(Ms),
    Span = lists:seq(Lo, Hi),
    Pos = maps:from_list(lists:zip(Ms, lists:seq(1, length(Ms)))),
    Entry = fun(R, Cc) ->
                    case {maps:find(R, Pos), maps:find(Cc, Pos)} of
                        {{ok, I}, {ok, J}} -> lists:nth(J, lists:nth(I, G));
                        _ when R =:= Cc -> {1.0, 0.0};
                        _ -> {0.0, 0.0}
                    end
            end,
    [{'Unitary', Span, #{matrix => [[Entry(R, Cc) || Cc <- Span] || R <- Span]}}].

%% our circuit from Perceval components (BS.Rx, PS, Unitary)
-spec from_perceval(#{n_modes := pos_integer(), components := [component()], input => term()}) ->
          {ok, #circuit{}} | {error, pq_core:reason()}.
from_perceval(#{n_modes := N, components := Cs} = D) ->
    Input = case maps:get(input, D, vacuum) of Occ when is_list(Occ) -> {fock, Occ}; _ -> vacuum end,
    try
        Instrs = lists:append([from_pcvl(X) || X <- Cs]),
        {ok, C} = pq_circuit:new(N, Input),
        {ok, C#circuit{instrs = Instrs}}
    catch
        throw:{error, _} = E -> E
    end.

from_pcvl({'BS.Rx', [A, B], P}) when B =:= A + 1 ->
    G = fun(K) -> maps:get(K, P, 0.0) end,
    [{phaseshift, A, G(phi_tl)}, {phaseshift, B, G(phi_bl)}, {beamsplitter, A, B, G(theta) / 2, math:pi() / 2},
     {phaseshift, A, G(phi_tr)}, {phaseshift, B, G(phi_br)}];
from_pcvl({'PS', [M], #{phi := P}}) -> [{phaseshift, M, P}];
from_pcvl({'Unitary', Ms, #{matrix := Rows}}) ->
    {ok, U} = pq_core:pack({cmat, [[{complex, Re, Im} || {Re, Im} <- R] || R <- Rows]}),
    [{unitary, Ms, U}];
from_pcvl(X) -> throw({error, {unsupported_instruction, X, perceval}}).

%% the unitary of Perceval components with the measured formulas (for tests and checks)
-spec perceval_unitary(pos_integer(), [component()]) -> pq_core:packed_cmatrix().
perceval_unitary(N, Cs) ->
    {ok, C} = from_perceval(#{n_modes => N, components => Cs}),
    pq_sampling:passive_unitary(C).

measurement(I) -> lists:member(element(1, I), [measure_pnr, measure_threshold, measure_homodyne]).

%%% Piquasso

-spec to_piquasso(#circuit{}, map()) -> {ok, map()} | {error, pq_core:reason()}.
to_piquasso(C, Binding) ->
    case pq_circuit:lower(C, Binding) of
        {ok, #circuit{n_modes = N, input = In, instrs = Is, hbar = H}} ->
            try
                {ok, #{d => N, hbar => case H of undefined -> pq_core:get(hbar); _ -> H end,
                       input => piq_input(In), instructions => [piq(I) || I <- Is]}}
            catch
                throw:{error, _} = E -> E
            end;
        E -> E
    end.

piq_input(vacuum) -> vacuum;
piq_input({fock, Occ}) -> {number_state, Occ};
piq_input({gaussian, S}) ->
    {Mean, Cov} = pq_state:mean_cov(S),
    %% Piquasso's covariance is twice ours; both xxpp here
    {gaussian, #{mean => Mean, cov => [[2 * X || X <- R] || R <- Cov], order => xxpp}}.

piq({beamsplitter, A, B, T, P}) -> {'Beamsplitter', [A, B], #{theta => T, phi => P}};
piq({phaseshift, M, P}) -> {'Phaseshifter', [M], #{phi => P}};
piq({squeezing, M, R, P}) -> {'Squeezing', [M], #{r => R, phi => P}};
piq({displacement, M, R, P}) -> {'Displacement', [M], #{r => R, phi => P}};
piq({kerr, M, K}) -> {'Kerr', [M], #{xi => K}};
piq({unitary, Ms, Packed}) ->
    {ok, {cmat, Rows}} = pq_core:unpack(Packed),
    {'Interferometer', Ms, #{matrix => [[{Re, Im} || {complex, Re, Im} <- R] || R <- Rows]}};
piq({measure_pnr, Ms}) -> {'ParticleNumberMeasurement', Ms, #{}};
piq({measure_threshold, Ms}) -> {'ThresholdMeasurement', Ms, #{}};
piq({measure_homodyne, M, P}) -> {'HomodyneMeasurement', [M], #{phi => P}};
piq(I) -> throw({error, {unsupported_instruction, element(1, I), piquasso}}).

-spec from_piquasso(map()) -> {ok, #circuit{}} | {error, pq_core:reason()}.
from_piquasso(#{d := N, instructions := Is} = D) ->
    H = maps:get(hbar, D, pq_core:get(hbar)),
    Input = case maps:get(input, D, vacuum) of
                vacuum -> vacuum;
                {number_state, Occ} -> {fock, Occ};
                {gaussian, #{mean := Mean, cov := Cov, order := xxpp}} ->
                    {gaussian, pq_state:from_mean_cov(Mean, [[X / 2 || X <- R] || R <- Cov], H)}
            end,
    try
        Instrs = [from_piq(X) || X <- Is],
        {ok, C} = pq_circuit:new(N, Input),
        {ok, C#circuit{instrs = Instrs, hbar = H}}
    catch
        throw:{error, _} = E -> E
    end.

from_piq({'Beamsplitter', [A, B], #{theta := T, phi := P}}) -> {beamsplitter, A, B, T, P};
from_piq({'Phaseshifter', [M], #{phi := P}}) -> {phaseshift, M, P};
from_piq({'Squeezing', [M], #{r := R, phi := P}}) -> {squeezing, M, R, P};
from_piq({'Displacement', [M], #{r := R, phi := P}}) -> {displacement, M, R, P};
from_piq({'Kerr', [M], #{xi := K}}) -> {kerr, M, K};
from_piq({'Interferometer', Ms, #{matrix := Rows}}) ->
    {ok, U} = pq_core:pack({cmat, [[{complex, Re, Im} || {Re, Im} <- R] || R <- Rows]}),
    {unitary, Ms, U};
from_piq({'ParticleNumberMeasurement', Ms, _}) -> {measure_pnr, Ms};
from_piq({'ThresholdMeasurement', Ms, _}) -> {measure_threshold, Ms};
from_piq({'HomodyneMeasurement', [M], #{phi := P}}) -> {measure_homodyne, M, P};
from_piq(X) -> throw({error, {unsupported_instruction, X, piquasso}}).
