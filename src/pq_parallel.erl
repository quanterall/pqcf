%% Shot-parallel execution on this node (FS 5.10, 9.1-9.4). Distribution across nodes (pg, erpc, handshake, spec
%% cache) comes in Phase 2.
-module(pq_parallel).
-export([run_job/2, chunks/3]).

%% Spec: what a worker runs (pq_worker:run_chunk/1), including the seed; Opts: shots, max_inflight, max_retries,
%% max_chunk, chunk (a fixed chunk size), timeout. Returns the samples in shot order.
-spec run_job(map(), map()) -> {ok, [list()]} | {error, pq_core:reason()}.
run_job(Spec, Opts) ->
    Shots = maps:get(shots, Opts),
    Workers = maps:get(max_inflight, Opts, erlang:system_info(dirty_cpu_schedulers_online)),
    Size = case maps:find(chunk, Opts) of
               {ok, C} -> C;
               error -> max(1, min(pq_core:get(max_chunk, Opts), ceil_div(Shots, 4 * Workers)))
           end,
    JobOpts = #{chunks => chunks(Shots, Size, 0), max_inflight => Workers,
                max_retries => pq_core:get(max_retries, Opts), timeout => maps:get(timeout, Opts, infinity)},
    case Shots of
        0 -> {ok, []};
        _ ->
            Ref = make_ref(),
            Self = self(),
            {ok, _Job} = pq_job_sup:start_job(Spec, JobOpts, {Self, Ref}),
            receive
                {Ref, Reply} -> Reply
            end
    end.

%% [{ChunkId, FirstShot, NShots}] covering 0..Shots-1
-spec chunks(non_neg_integer(), pos_integer(), non_neg_integer()) -> [{non_neg_integer(), non_neg_integer(), pos_integer()}].
chunks(Shots, Size, Id) when Id * Size >= Shots -> [];
chunks(Shots, Size, Id) -> [{Id, Id * Size, min(Size, Shots - Id * Size)} | chunks(Shots, Size, Id + 1)].

ceil_div(A, B) -> (A + B - 1) div B.
