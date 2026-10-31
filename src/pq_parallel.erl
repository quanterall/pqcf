%% Shot-parallel execution (FS 5.10, 9.1-9.5): on this node, or with distributed => true on every registered node
%% that passes the handshake (pq_node_registry), chunks spread over them and retried elsewhere when a node fails.
%% Before dispatch the job's memory is estimated against max_job_memory (FR-P08).
-module(pq_parallel).
-export([run_job/2, chunks/3, job_opts/1, cancel/1, status/1, register_node/0]).

%% the run options that reach the job server, for every caller (the simulator passed only five, so distributed and
%% isolate_kernels were dropped for circuits other than boson sampling)
-spec job_opts(map()) -> map().
job_opts(Opts) ->
    maps:with([max_inflight, max_retries, max_chunk, chunk, timeout, distributed, nodes, max_job_memory, auto_throttle,
               isolate_kernels, isolated_nodes, halt_first_isolated, notify, partial_results], Opts).

%% FS 7, 9.4: a job's handle is the pid sent by notify => Pid as {pq_job, Job}
-spec cancel(pid()) -> ok.
cancel(Job) -> pq_job:cancel(Job).

-spec status(pid()) -> {ok, map()}.
status(Job) -> pq_job:status(Job).

-spec register_node() -> term().
register_node() -> pq_node_registry:register_node().

%% Spec: what a worker runs (pq_worker:run_chunk/1), including the seed; Opts: shots, max_inflight, max_retries,
%% max_chunk, chunk (a fixed chunk size), timeout, partial_results, distributed, notify, max_job_memory,
%% auto_throttle. Returns the samples in shot order.
-spec run_job(map(), map()) -> {ok, [list()]} | {ok, [list()], [{pos_integer(), map()}]} | {error, pq_core:reason()}.
run_job(Spec, Opts) ->
    case maps:get(isolate_kernels, Opts, false) of
        true -> isolated(Spec, Opts);
        false -> run_job(Spec, Opts, maps:get(distributed, Opts, false))
    end.

%% FR-P07: the job's chunks run on peer nodes it starts (isolate_kernels nodes, default 2) and stops afterwards; a
%% hard fault in a kernel kills only its node, and the job retries the chunk on another (FS 9.4)
isolated(Spec, Opts) ->
    case is_alive() of
        false -> {error, {unsupported_input, isolate_kernels_needs_distribution}};
        true ->
            Paths = lists:append([["-pa", P] || P <- code:get_path()]),
            Peers = [begin
                         {ok, Pid, Node} = peer:start(#{name => peer:random_name(pq_isolated), args => Paths}),
                         {ok, _} = erpc:call(Node, application, ensure_all_started, [pq]),
                         {Pid, Node}
                     end || _ <- lists:seq(1, maps:get(isolated_nodes, Opts, 2))],
            Nodes = [N || {_, N} <- Peers],
            %% test hook (FS 12.5): the first isolated node halts inside its first chunk, like a crashed kernel
            Spec2 = case maps:get(halt_first_isolated, Opts, false) of
                        true -> Spec#{halt_node => hd(Nodes)};
                        false -> Spec
                    end,
            try run_job(Spec2, Opts#{nodes => Nodes}, nodes)
            after [catch peer:stop(P) || {P, _} <- Peers]
            end
    end.

run_job(Spec, Opts, Mode) ->
    Shots = maps:get(shots, Opts),
    {Nodes, _Hetero} = case Mode of
                           nodes -> {maps:get(nodes, Opts), false};
                           true ->
                               case pq_node_registry:eligible_nodes(Opts) of
                                   {[], H} -> {[node()], H};
                                   NH -> NH
                               end;
                           false -> {[node()], false}
                       end,
    Workers = maps:get(max_inflight, Opts, erlang:system_info(dirty_cpu_schedulers_online) * length(Nodes)),
    Size = case maps:find(chunk, Opts) of
               {ok, C} -> C;
               error -> max(1, min(pq_core:get(max_chunk, Opts), ceil_div(Shots, 4 * Workers)))
           end,
    Offset = maps:get(offset, Opts, 0),
    %% FR-P08: the spec held per worker plus the spec cache, against max_job_memory (default half the memory);
    %% auto_throttle lowers max_inflight to fit instead of refusing
    PerWorker = erlang:external_size(Spec),
    Cache = pq_core:get(spec_cache_bytes),
    Limit = maps:get(max_job_memory, Opts, memory_half()),
    Inflight = case maps:get(auto_throttle, Opts, false) of
                   true when PerWorker * Workers + Cache > Limit -> max(1, min(Workers, (Limit - Cache) div PerWorker));
                   _ -> Workers
               end,
    Estimate = PerWorker * Inflight + Cache,
    JobOpts = maps:merge(maps:with([notify, partial_results], Opts),
                         #{chunks => [{Id, F + Offset, N} || {Id, F, N} <- chunks(Shots, Size, 0)], max_inflight => Inflight,
                           max_retries => pq_core:get(max_retries, Opts), timeout => maps:get(timeout, Opts, infinity),
                           nodes => Nodes}),
    case Shots of
        _ when Estimate > Limit -> {error, {dimension_limit, Estimate, Limit}};
        0 -> {ok, []};
        %% the pq application not started and nothing distributed: the whole job as one chunk here (the same
        %% samples, D1), so the library also works without its supervision tree
        _ when Mode =:= false ->
            case whereis(pq_job_sup) of
                undefined -> inline(Spec, Offset, Shots);
                _ -> start(Spec, JobOpts)
            end;
        _ -> start(Spec, JobOpts)
    end.

inline(Spec, First, N) ->
    case pq_worker:run_chunk({Spec, First, N}) of
        {ok, Rows, Extra} -> {ok, Rows, [{N, Extra}]};
        R -> R
    end.

start(Spec, JobOpts) ->
    Ref = make_ref(),
    {ok, _Job} = pq_job_sup:start_job(Spec, JobOpts, {self(), Ref}),
    receive
        {Ref, Reply} -> Reply
    end.

%% [{ChunkId, FirstShot, NShots}] covering 0..Shots-1
-spec chunks(non_neg_integer(), pos_integer(), non_neg_integer()) -> [{non_neg_integer(), non_neg_integer(), pos_integer()}].
chunks(Shots, Size, Id) when Id * Size >= Shots -> [];
chunks(Shots, Size, Id) -> [{Id, Id * Size, min(Size, Shots - Id * Size)} | chunks(Shots, Size, Id + 1)].

ceil_div(A, B) -> (A + B - 1) div B.

memory_half() -> total_memory_from_proc() div 2.

total_memory_from_proc() ->
    case file:read_file("/proc/meminfo") of
        {ok, <<"MemTotal:", Rest/binary>>} ->
            [Kb | _] = string:lexemes(binary_to_list(Rest), " \n"),
            list_to_integer(Kb) * 1024;
        _ -> 1 bsl 40
    end.
