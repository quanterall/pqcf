%% A worker (FS 2.2): runs one chunk of shots and reports to its job, then exits. Temporary under pq_worker_sup;
%% the job decides on retries.
-module(pq_worker).
-export([start_link/3, run_chunk/1, run_cached/3]).

%% The worker links to its job before anything else, so the job always sees the real exit reason (linking from
%% the job after the start raced with an early crash: noproc).
start_link(Job, Chunk, Attempt) ->
    {ok, proc_lib:spawn_link(fun() ->
                                     link(Job),
                                     Job ! {chunk_done, self(), element(1, Chunk), run(Chunk, Attempt)}
                             end)}.

run({Id, First, N, Spec}, Attempt) ->
    %% failure injection for tests (FS 12.5): crash on the listed {ChunkId, Attempt}
    lists:member({Id, Attempt}, maps:get(crash, Spec, [])) andalso exit(injected_crash),
    case maps:get(target, Spec, node()) of
        Node when Node =:= node() -> run_chunk({maps:remove(target, Spec), First, N});
        Node -> remote(Node, maps:without([target, cancel], Spec), First, N)
    end.

%% a chunk on another node (FS 9.5): erpc with only the spec's hash; the spec itself only when the node lacks it.
%% A lost node or a crash there exits this proxy, and the job retries the chunk elsewhere.
remote(Node, Spec, First, N) ->
    Hash = pq_spec_cache:hash(Spec),
    case erpc:call(Node, ?MODULE, run_cached, [Hash, First, N]) of
        spec_missing ->
            ok = erpc:call(Node, pq_spec_cache, put, [Hash, Spec]),
            erpc:call(Node, ?MODULE, run_cached, [Hash, First, N]);
        R -> R
    end.

%% on the worker node: the cached spec's chunk, or spec_missing
-spec run_cached(binary(), non_neg_integer(), non_neg_integer()) -> {ok, [list()]} | {error, term()} | spec_missing.
run_cached(Hash, First, N) ->
    case pq_spec_cache:get(Hash) of
        {ok, #{halt_node := Node}} when Node =:= node() -> erlang:halt(1);       % the FR-P07 test hook
        {ok, Spec} -> run_chunk({Spec, First, N});
        error -> spec_missing
    end.

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
