%% A worker (FS 2.2): runs one chunk of shots and reports to its job, then exits. Temporary under pq_worker_sup;
%% the job decides on retries.
-module(pq_worker).
-export([start_link/3, run_chunk/1]).

start_link(Job, Chunk, Attempt) ->
    {ok, proc_lib:spawn_link(fun() -> Job ! {chunk_done, self(), element(1, Chunk), run(Chunk, Attempt)} end)}.

run({Id, First, N, Spec}, Attempt) ->
    %% failure injection for tests (FS 12.5): crash on the listed {ChunkId, Attempt}
    lists:member({Id, Attempt}, maps:get(crash, Spec, [])) andalso exit(injected_crash),
    run_chunk({Spec, First, N}).

%% {Spec, FirstShot, NShots} -> {ok, Samples} | {error, Reason}; shot k uses stream (Seed, k)
-spec run_chunk({map(), non_neg_integer(), non_neg_integer()}) -> {ok, [list()]} | {error, term()}.
run_chunk({#{kind := circuit} = Spec, First, N}) ->
    pq_simulator:shots(Spec, First, N);
run_chunk({#{kind := boson, unitary := U, input := In, seed := Seed} = Spec, First, N}) ->
    Opts = maps:merge(maps:with([algorithm, max_enumeration, cancel, loss, burn_in, thin, chains], Spec),
                      #{seed => Seed, stream_start => First}),
    case pq_nif:sample_boson(U, In, N, Opts) of
        {ok, {Packed, _}} ->
            {ok, {imat, Rows}} = pq_core:unpack(Packed),
            {ok, Rows};
        {error, _} = E -> E
    end.
