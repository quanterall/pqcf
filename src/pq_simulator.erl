%% Running a circuit on a backend (FS 5.3, 9.1): validate, check capabilities, lower once, evolve instruction by
%% instruction, then the terminal measurements. Phase 1: single-process execution; the job server (pq_job) takes
%% over for parallel shots.
-module(pq_simulator).
-include("pq.hrl").
-export([evolve/3, run/3]).

-spec evolve(module(), #circuit{}, map()) -> {ok, term(), [{atom(), term()}]} | {error, pq_core:reason()}.
evolve(Backend, C, Opts) ->
    case prepare(Backend, C, Opts) of
        {ok, Low} ->
            Unitary = [I || I <- Low#circuit.instrs, not measurement(I)],
            case evolve_instrs(Backend, Low, Unitary, Opts) of
                {ok, S, W} -> {ok, public(Backend, S), W};
                E -> E
            end;
        E -> E
    end.

public(Backend, S) ->
    case erlang:function_exported(Backend, final, 1) of
        true -> Backend:final(S);
        false -> S
    end.

-spec run(module(), #circuit{}, map()) -> {ok, #result{}} | {error, pq_core:reason()}.
run(Backend, C, Opts) ->
    Seed = maps:get(seed, Opts, random_seed()),
    Shots = maps:get(shots, Opts, 1),
    case prepare(Backend, C, Opts) of
        {ok, Low} ->
            {Body, Ms} = lists:splitwith(fun(I) -> not measurement(I) end, Low#circuit.instrs),
            Mid = [K || {K, I} <- lists:enumerate(length(Body), Ms), not measurement(I)],
            case evolve_instrs(Backend, Low, Body, Opts) of
                _ when Mid =/= [] ->
                    %% mid-circuit measurement with later instructions: Phase 2 (FR-G08, FR-M07)
                    {error, {unsupported_instruction, midcircuit_measurement, #{instr_index => hd(Mid), available => phase_2}}};
                {ok, State, Warnings} when Ms =/= [] ->
                    case Backend:sample(Ms, State, Opts#{seed => Seed, shots => Shots, stream_start => 0}) of
                        {ok, Kind, Modes, Samples, Extra0} ->
                            {More, Extra} = {maps:get(extra_warnings, Extra0, []), maps:remove(extra_warnings, Extra0)},
                            Meta = maps:merge(Extra, #{seed => Seed, shots_requested => Shots,
                                                       shots_accepted => length(Samples), acceptance_rate => 1.0,
                                                       backend => Backend:name(),
                                                       warnings => Warnings ++ maps:get(warnings, Low#circuit.meta, []) ++ More}),
                            {ok, #result{kind = Kind, measured = Modes, samples = Samples, metadata = Meta}};
                        {error, _} = E -> E
                    end;
                {ok, _, _} -> {error, {unsupported_input, no_measurement}};
                E -> E
            end;
        E -> E
    end.

prepare(Backend, C, Opts) ->
    case pq_circuit:validate(C) of
        ok ->
            case pq_circuit:check_capabilities(C, Backend) of
                ok -> pq_circuit:lower(C, maps:get(bindings, Opts, #{}));
                E -> E
            end;
        E -> E
    end.

evolve_instrs(Backend, Low, Instrs, Opts) ->
    case Backend:init(Low#circuit.input, Low, Opts) of
        {ok, S0} ->
            lists:foldl(fun({K, I}, {ok, S, W}) ->
                                case Backend:apply(I, S, Opts#{instr_index => K}) of
                                    {ok, S2, W2} -> {ok, S2, W ++ W2};
                                    E -> E
                                end;
                           (_, E) -> E
                        end, {ok, S0, []}, lists:enumerate(0, Instrs));
        E -> E
    end.

measurement(I) -> lists:member(element(1, I), [measure_pnr, measure_threshold, measure_homodyne]).

%% a missing seed is drawn (not with rand: FS 8.3) and recorded in the result
random_seed() ->
    <<S:64>> = crypto:strong_rand_bytes(8),
    S.
