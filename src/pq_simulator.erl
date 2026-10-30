%% Running a circuit on a backend (FS 5.3, 9.1, FR-M05-M07): validate, check capabilities, lower once, evolve
%% instruction by instruction, measure.
%%
%% Without mid-circuit measurements the backend samples the terminal measurements of all shots in one call. With
%% them, everything before the first measurement is evolved once and the rest runs shot by shot (through the job
%% server for parallel shots): shot k consumes its stream (Seed, k) in instruction order -- each mid-circuit
%% measurement's draw, then the terminal sampler's (FS 9.3). Measured modes stay in the state (no renumbering).
%%
%% Post-selection (FR-M05): postselect => predicate() (FS 4.7) on the samples; shots_mode attempts (default: the
%% accepted ones of `shots` attempts) or accepted (attempts until `shots` are accepted, at most max_attempts); the
%% attempt index is the stream index, so results do not depend on chunking.
-module(pq_simulator).
-include("pq.hrl").
-export([evolve/3, run/3, shots/3, eval/3]).

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
run(Backend, C, Opts0) ->
    Seed = maps:get(seed, Opts0, random_seed()),
    Shots = maps:get(shots, Opts0, 1),
    Opts = Opts0#{seed => Seed},
    try
        Low = ok(prepare(Backend, C, Opts)),
        {Body, Rest} = lists:splitwith(fun(I) -> not measurement(I) end, Low#circuit.instrs),
        Rest =/= [] orelse throw({error, {unsupported_input, no_measurement}}),
        Mid = lists:any(fun(I) -> not measurement(I) end, Rest),
        Mid andalso (erlang:function_exported(Backend, measure, 3) orelse
                     throw({error, {unsupported_instruction, midcircuit_measurement, Backend:name()}})),
        {S0, Warnings} = ok2(evolve_instrs(Backend, Low, Body, Opts)),
        Spec = #{kind => circuit, backend => Backend, state => S0, rest => Rest, opts => Opts, seed => Seed},
        Rows = fun(Start, N) -> rows(Spec, Mid, Start, N, Opts) end,
        Measured = lists:append([modes(I) || I <- Rest, measurement(I)]),
        {Samples, Attempts, Extra} = postselected(Rows, Shots, Opts#{measured_columns => Measured}),
        Kinds = maps:get(column_kinds, Extra, column_kinds(Rest, Mid)),
        Kind = case lists:usort(Kinds) of [K] -> K; _ -> mixed end,
        Meta0 = maps:merge(maps:without([column_kinds, extra_warnings], Extra),
                           #{seed => Seed, shots_requested => Shots, shots_accepted => length(Samples),
                             attempts => Attempts, acceptance_rate => length(Samples) / max(1, Attempts),
                             backend => Backend:name(),
                             warnings => Warnings ++ maps:get(warnings, Low#circuit.meta, []) ++
                                 maps:get(extra_warnings, Extra, [])}),
        Meta = case Kind of mixed -> Meta0#{column_kinds => Kinds}; _ -> Meta0 end,
        {ok, #result{kind = Kind, measured = Measured, samples = Samples, metadata = Meta}}
    catch
        throw:{error, _} = E -> E
    end.

ok({ok, V}) -> V;
ok({error, _} = E) -> throw(E).

ok2({ok, A, B}) -> {A, B};
ok2({error, _} = E) -> throw(E).

%% the samples of shots Start .. Start + N - 1, and what the backend adds to the metadata
rows(#{backend := Backend, state := S0, rest := Rest} = Spec, false, Start, N, Opts) ->
    _ = Spec,
    case Backend:sample(Rest, S0, Opts#{shots => N, stream_start => Start}) of
        {ok, _Kind, _Modes, Samples, Extra} -> {Samples, Extra};
        {error, _} = E -> throw(E)
    end;
rows(Spec, true, Start, N, Opts) ->
    Parallel = maps:get(parallel, Opts, true) andalso N >= pq_core:get(parallel_threshold, Opts),
    R = case Parallel of
            true -> pq_parallel:run_job(Spec, (maps:with([max_inflight, max_retries, max_chunk, chunk, timeout], Opts))#{shots => N, offset => Start});
            false -> shots(Spec, Start, N)
        end,
    case R of
        {ok, Samples} -> {Samples, #{}};
        {error, _} = E -> throw(E)
    end.

%% shot by shot (a worker's chunk): every shot from the state before its first measurement
-spec shots(map(), non_neg_integer(), non_neg_integer()) -> {ok, [list()]} | {error, term()}.
shots(#{backend := Backend, state := S0, rest := Rest, opts := Opts, seed := Seed}, First, N) ->
    try
        {ok, [shot(Backend, S0, Rest, Opts, pq_nif:stream_new(Seed, K), []) || K <- lists:seq(First, First + N - 1)]}
    catch
        throw:{error, _} = E -> E
    end.

shot(_Backend, _S, [], _Opts, _R, Acc) -> Acc;
shot(Backend, S, [I | T] = Is, Opts, R, Acc) ->
    case measurement(I) of
        true ->
            case lists:all(fun measurement/1, T) of
                true ->
                    %% the terminal block: one sample from the same stream
                    case Backend:sample(Is, S, Opts#{shots => 1, stream_ref => R}) of
                        {ok, _, _, [Row], _} -> Acc ++ Row;
                        {error, _} = E -> throw(E)
                    end;
                false ->
                    case Backend:measure(I, S, Opts#{stream_ref => R}) of
                        {ok, Out, S2} -> shot(Backend, S2, T, Opts, R, Acc ++ Out);
                        {error, _} = E -> throw(E)
                    end
            end;
        false ->
            case Backend:apply(I, S, Opts) of
                {ok, S2, _} -> shot(Backend, S2, T, Opts, R, Acc);
                {error, _} = E -> throw(E)
            end
    end.

column_kinds(Rest, _Mid) ->
    lists:append([[kind(I) || _ <- modes(I)] || I <- Rest, measurement(I)]).

kind({measure_pnr, _}) -> pnr;
kind({measure_threshold, _}) -> threshold;
kind({measure_homodyne, _, _}) -> homodyne.

modes({measure_homodyne, M, _}) -> [M];
modes({_, Ms}) -> Ms.

%%% Post-selection (FR-M05, FS 4.7)

postselected(Rows, Shots, Opts) ->
    case maps:find(postselect, Opts) of
        error ->
            {S, Extra} = Rows(0, Shots),
            {S, Shots, Extra};
        {ok, Pred} ->
            Accept = fun(Samples) -> [X || X <- Samples, eval(Pred, X, Opts)] end,
            _ = Pred,
            case pq_core:get(shots_mode, Opts) of
                attempts ->
                    {S, Extra} = Rows(0, Shots),
                    {Accept(S), Shots, Extra};
                accepted ->
                    Max = case pq_core:get(max_attempts, Opts) of derived -> 1000 * Shots; M -> M end,
                    accepted(Rows, Accept, Shots, Max, 0, [], #{})
            end
    end.

accepted(_Rows, _Accept, Shots, _Max, A, Acc, Extra) when length(Acc) >= Shots ->
    {lists:sublist(Acc, Shots), A, Extra};
accepted(_Rows, _Accept, Shots, Max, A, Acc, _Extra) when A >= Max ->
    _ = Shots,
    throw({error, {postselection_exhausted, #{attempts => A, accepted => length(Acc)}}});
accepted(Rows, Accept, Shots, Max, A, Acc, _Extra) ->
    B = min(max(Shots, 64), Max - A),
    {S, Extra} = Rows(A, B),
    accepted(Rows, Accept, Shots, Max, A + B, Acc ++ Accept(S), Extra).

%% a predicate on one sample; modes refer to the circuit's modes (the result's measured columns, passed in Opts as
%% the instruction order is known only at run time)
-spec eval(term(), list(), map()) -> boolean().
eval(Pred, Sample, Opts) -> ev(Pred, col_map(Opts), Sample).

col_map(#{measured_columns := Cols}) -> Cols;
col_map(_) -> undefined.

ev({'and', Ps}, M, S) -> lists:all(fun(P) -> ev(P, M, S) end, Ps);
ev({'or', Ps}, M, S) -> lists:any(fun(P) -> ev(P, M, S) end, Ps);
ev({'not', P}, M, S) -> not ev(P, M, S);
ev({local_fun, F}, _M, S) -> F(S);
ev({eq, Mo, N}, M, S) -> v(Mo, M, S) == N;
ev({ge, Mo, N}, M, S) -> v(Mo, M, S) >= N;
ev({le, Mo, N}, M, S) -> v(Mo, M, S) =< N;
ev({in, Mo, Ns}, M, S) -> lists:member(v(Mo, M, S), Ns);
ev({total_eq, Ms, N}, M, S) -> lists:sum([v(Mo, M, S) || Mo <- Ms]) == N;
ev({total_ge, Ms, N}, M, S) -> lists:sum([v(Mo, M, S) || Mo <- Ms]) >= N;
ev({total_le, Ms, N}, M, S) -> lists:sum([v(Mo, M, S) || Mo <- Ms]) =< N.

v(Mode, undefined, S) -> lists:nth(Mode + 1, S);
v(Mode, Cols, S) -> lists:nth(string:str(Cols, [Mode]), S).

%%% Preparation and evolution

prepare(Backend, C, Opts) ->
    case pq_circuit:validate(C) of
        ok ->
            case pq_circuit:check_capabilities(C, Backend, Opts) of
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
