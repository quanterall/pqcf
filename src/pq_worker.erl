%% A worker (FS 2.2): runs one chunk of shots and reports to its job, then exits. Temporary under pq_worker_sup;
%% the job decides on retries.
-module(pq_worker).
-export([start_link/3, run_chunk/1, run_cached/4, cancel_remote/1]).

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
%% A lost node or a crash there exits this proxy, and the job retries the chunk elsewhere. The proxy traps exits:
%% when its job stops (cancelled, timed out, failed), it sets the chunk's cancel token on that node, so the remote
%% kernel stops within its polling interval instead of running to the end of the chunk (FS 9.4).
remote(Node, Spec, First, N) ->
    process_flag(trap_exit, true),
    Hash = pq_spec_cache:hash(Spec),
    Key = make_ref(),
    case await(Node, Key, erpc:send_request(Node, ?MODULE, run_cached, [Hash, First, N, Key])) of
        spec_missing ->
            ok = erpc:call(Node, pq_spec_cache, put, [Hash, Spec]),
            await(Node, Key, erpc:send_request(Node, ?MODULE, run_cached, [Hash, First, N, Key]));
        R -> R
    end.

await(Node, Key, Req) ->
    receive
        {'EXIT', _From, Reason} ->
            erpc:cast(Node, ?MODULE, cancel_remote, [Key]),
            exit(Reason);
        Msg ->
            case erpc:check_response(Msg, Req) of
                {response, R} -> R;
                no_response -> await(Node, Key, Req)
            end
    end.

%% on the worker node: the cached spec's chunk under a cancel token of this node, or spec_missing
-spec run_cached(binary(), non_neg_integer(), non_neg_integer(), reference()) ->
          {ok, [list()] | pq_core:packed_imatrix()} | {ok, [list()], map()} | {error, term()} | spec_missing.
run_cached(Hash, First, N, Key) ->
    case pq_spec_cache:get(Hash) of
        {ok, #{halt_node := Node}} when Node =:= node() -> erlang:halt(1);       % the FR-P07 test hook
        {ok, Spec} ->
            Token = pq_nif:cancel_token_new(),
            true = ets:insert(pq_remote_cancel, {Key, Token}),
            try run_chunk({Spec#{cancel => Token}, First, N})
            after ets:delete(pq_remote_cancel, Key)
            end;
        error -> spec_missing
    end.

-spec cancel_remote(reference()) -> ok.
cancel_remote(Key) ->
    case ets:lookup(pq_remote_cancel, Key) of
        [{_, Token}] -> _ = pq_nif:cancel_token_set(Token), ok;
        [] -> ok
    end.

%% {Spec, FirstShot, NShots} -> {ok, Samples} | {error, Reason}; shot k uses stream (Seed, k)
-spec run_chunk({map(), non_neg_integer(), non_neg_integer()}) -> {ok, [list()] | pq_core:packed_imatrix()} | {ok, [list()], map()} | {error, term()}.
run_chunk({#{kind := circuit, terminal := true} = Spec, First, N}) ->
    pq_simulator:terminal(Spec, First, N);
run_chunk({#{kind := circuit} = Spec, First, N}) ->
    pq_simulator:shots(Spec, First, N);
run_chunk({#{kind := boson, unitary := U, input := In, seed := Seed} = Spec, First, N}) ->
    Opts = maps:merge(maps:with([algorithm, max_enumeration, cancel, loss, shared, burn_in, thin, chains], Spec),
                      #{seed => Seed, stream_start => First}),
    %% packed: the rows travel to the caller as a reference-counted binary, made lists there (pq_parallel:to_rows/1)
    case pq_nif:sample_boson(U, In, N, Opts) of
        {ok, {{imat, _, _, _} = Packed, _}} -> {ok, Packed};
        {error, _} = E -> E
    end.
