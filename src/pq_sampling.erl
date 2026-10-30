%% Sampling (FS 5.6, 7.1). Phase 1: boson sampling of a passive circuit with a Fock input and PNR detection on
%% every mode, through the job server for parallel shots (FR-B04, FS 9.1).
-module(pq_sampling).
-include("pq.hrl").
-export([sample/3, passive_unitary/1]).

-spec sample(boson, #circuit{}, map()) -> {ok, #result{}} | {error, pq_core:reason()}.
sample(boson, C, Opts) ->
    Seed = maps:get(seed, Opts, random_seed()),
    Shots = maps:get(shots, Opts, 1),
    case prepare(C, Opts) of
        {ok, U, In, Measured} ->
            Alg = maps:get(algorithm, Opts, clifford_clifford),
            Spec = #{kind => boson, unitary => U, input => In, seed => Seed, algorithm => Alg,
                     max_enumeration => pq_core:get(max_enumeration, Opts), crash => maps:get(crash, Opts, [])},
            Parallel = maps:get(parallel, Opts, true) andalso Shots >= pq_core:get(parallel_threshold, Opts),
            R = case Parallel of
                    %% brute_force builds the whole distribution per call: one chunk (A-5)
                    true when Alg =:= brute_force -> pq_parallel:run_job(Spec, (job_opts(Opts))#{shots => Shots, chunk => Shots});
                    true -> pq_parallel:run_job(Spec, (job_opts(Opts))#{shots => Shots});
                    false -> pq_worker:run_chunk({Spec, 0, Shots})
                end,
            case R of
                {ok, Rows} ->
                    Cols = [[lists:nth(M + 1, Row) || M <- Measured] || Row <- Rows],
                    {ok, #result{kind = pnr, measured = Measured, samples = Cols,
                                 metadata = #{seed => Seed, shots_requested => Shots, shots_accepted => length(Rows),
                                              acceptance_rate => 1.0, backend => boson, algorithm => Alg,
                                              warnings => []}}};
                {error, _} = E -> E
            end;
        E -> E
    end.

job_opts(Opts) -> maps:with([max_inflight, max_retries, max_chunk, chunk, timeout], Opts).

%% the lowered circuit's single-photon unitary, its Fock input and the measured modes
prepare(C, Opts) ->
    case pq_circuit:lower(C, maps:get(bindings, Opts, #{})) of
        {ok, #circuit{input = {fock, In}, instrs = Is}} ->
            {Body, Ms} = lists:splitwith(fun(I) -> element(1, I) =/= measure_pnr end, Is),
            Measured = lists:append([Modes || {measure_pnr, Modes} <- Ms]),
            Bad = [element(1, I) || I <- Body ++ Ms,
                                    not lists:member(element(1, I), [beamsplitter, phaseshift, unitary, measure_pnr])],
            case {Bad, Measured} of
                {[], [_ | _]} -> {ok, passive_unitary(C#circuit{instrs = Body, input = {fock, In}}), In, Measured};
                {[], []} -> {error, {unsupported_input, no_measurement}};
                {[K | _], _} -> {error, {unsupported_instruction, K, boson}}
            end;
        {ok, _} -> {error, {unsupported_input, boson_sampling_needs_a_fock_input}};
        E -> E
    end.

%% U = G_k .. G_1 of the lowered passive instructions (FS 3.3 single-photon matrices), packed
-spec passive_unitary(#circuit{}) -> pq_core:packed_cmatrix().
passive_unitary(#circuit{n_modes = N, instrs = Is}) ->
    Id = [[case I =:= J of true -> {1.0, 0.0}; false -> {0.0, 0.0} end || J <- lists:seq(1, N)] || I <- lists:seq(1, N)],
    U = lists:foldl(fun apply_rows/2, list_to_tuple([list_to_tuple(R) || R <- Id]), Is),
    {ok, P} = pq_core:pack({cmat, [[{complex, Re, Im} || {Re, Im} <- tuple_to_list(R)] || R <- tuple_to_list(U)]}),
    P.

%% U <- G U: G acts on the rows of the instruction's modes
apply_rows({phaseshift, M, P}, U) -> rows(U, [M], [[{math:cos(P), math:sin(P)}]]);
apply_rows({beamsplitter, A, B, T, P}, U) ->
    C = math:cos(T), S = math:sin(T),
    rows(U, [A, B], [[{C, 0.0}, {-math:cos(P) * S, math:sin(P) * S}], [{math:cos(P) * S, math:sin(P) * S}, {C, 0.0}]]);
apply_rows({unitary, Ms, Packed}, U) ->
    {ok, {cmat, G}} = pq_core:unpack(Packed),
    rows(U, Ms, [[{Re, Im} || {complex, Re, Im} <- Row] || Row <- G]).

rows(U, Ms, G) ->
    Old = [element(M + 1, U) || M <- Ms],
    New = [list_to_tuple([lists:foldl(fun({Gij, Row}, {Sr, Si}) ->
                                              {Ar, Ai} = element(Col, Row),
                                              {Gr, Gi} = Gij,
                                              {Sr + Gr * Ar - Gi * Ai, Si + Gr * Ai + Gi * Ar}
                                      end, {0.0, 0.0}, lists:zip(GRow, Old))
                          || Col <- lists:seq(1, tuple_size(element(1, U)))]) || GRow <- G],
    lists:foldl(fun({M, R}, Acc) -> setelement(M + 1, Acc, R) end, U, lists:zip(Ms, New)).

random_seed() ->
    <<S:64>> = crypto:strong_rand_bytes(8),
    S.
