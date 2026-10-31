%% One run's job (FS 2.2, 9.1-9.4): splits the shots into chunks, dispatches at most max_inflight at a time to
%% temporary workers, retries a failed chunk up to max_retries, merges in shot order and replies to the caller.
%% Cancellation sets the run's cancel token (kernels poll it); a timeout cancels and fails with timeout, or with
%% partial_results => true with the samples finished so far and their shot (stream) indices (FS 9.4). status/1
%% reports progress; with notify => Pid the job sends {pq_job, self()} on start, the handle for cancel/status.
%% Telemetry events of FS 9.6 through pq_telemetry. Chunks may report metadata, {ok, Samples, Extra}: the reply is then
%% {ok, Samples, [{NShots, Extra}]} in shot order.
-module(pq_job).
-behaviour(gen_statem).
-export([start_link/3, cancel/1, status/1]).
-export([init/1, callback_mode/0, running/3, terminate/3]).

-record(d, {spec, from, queue = [], inflight = #{}, attempts = #{}, done = #{}, max_inflight, max_retries,
            token, chunks, nodes = [node()], partial = false, started, started_chunks = #{}, extras = #{}}).

start_link(Spec, Opts, From) -> gen_statem:start_link(?MODULE, {Spec, Opts, From}, []).

cancel(Job) -> gen_statem:cast(Job, cancel).

%% {ok, #{state := running, ...}} with the progress, or {ok, #{state := done}} once the job has replied
-spec status(pid()) -> {ok, map()}.
status(Job) ->
    try gen_statem:call(Job, status) of
        S -> {ok, S}
    catch
        exit:{R, _} when R =:= noproc; R =:= normal -> {ok, #{state => done}}
    end.

callback_mode() -> state_functions.

init({Spec, Opts, From}) ->
    process_flag(trap_exit, true),
    Token = pq_nif:cancel_token_new(),
    Chunks = maps:get(chunks, Opts),
    D = #d{spec = Spec#{cancel => Token}, from = From, queue = Chunks, chunks = length(Chunks),
           max_inflight = maps:get(max_inflight, Opts), max_retries = maps:get(max_retries, Opts), token = Token,
           nodes = maps:get(nodes, Opts, [node()]), partial = maps:get(partial_results, Opts, false),
           started = erlang:monotonic_time()},
    _ = case maps:get(notify, Opts, undefined) of
            Pid when is_pid(Pid) -> Pid ! {pq_job, self()};
            _ -> ok
        end,
    pq_telemetry:execute([pq, job, start], #{system_time => erlang:system_time()},
                         #{job => self(), chunks => length(Chunks), nodes => D#d.nodes, max_inflight => D#d.max_inflight}),
    Actions = case maps:get(timeout, Opts, infinity) of
                  infinity -> [];
                  T -> [{{timeout, run}, T, expired}]
              end,
    {ok, running, dispatch(D), Actions}.

running({call, From}, status, D) ->
    {keep_state_and_data, [{reply, From, #{state => running, chunks => D#d.chunks, done => map_size(D#d.done),
                                             inflight => map_size(D#d.inflight), queued => length(D#d.queue),
                                             retries => lists:sum(maps:values(D#d.attempts)), nodes => D#d.nodes}}]};
running(info, {chunk_done, Pid, Id, {ok, Samples, Extra}}, D) ->
    %% a chunk that also reports metadata (terminal sampling): kept with its shot count for the caller to merge
    {Id, First, N} = maps:get(Pid, D#d.inflight),
    running(info, {chunk_done, Pid, Id, {ok, Samples}}, D#d{extras = (D#d.extras)#{First => {N, Extra}}});
running(info, {chunk_done, Pid, Id, {ok, Samples}}, D) ->
    {Id, First, N} = maps:get(Pid, D#d.inflight),
    chunk_event(stop, Id, D, #{shots => N}),
    finish_or_continue(dispatch(D#d{inflight = maps:remove(Pid, D#d.inflight), done = (D#d.done)#{First => Samples}}), Id);
running(info, {chunk_done, Pid, Id, {error, cancelled}}, D) ->
    stop_with({error, cancelled}, D#d{inflight = maps:remove(Pid, D#d.inflight)}, Id);
running(info, {chunk_done, Pid, Id, {error, R}}, D) ->
    stop_with({error, {chunk_failed, Id, R}}, D#d{inflight = maps:remove(Pid, D#d.inflight)}, Id);
running(info, {'EXIT', Pid, normal}, _D) when is_pid(Pid) -> keep_state_and_data;
running(info, {'EXIT', Pid, Reason}, D) ->
    case maps:take(Pid, D#d.inflight) of
        {{Id, First, N}, In} ->
            A = maps:get(Id, D#d.attempts, 0) + 1,
            case A > D#d.max_retries of
                true -> stop_with({error, {chunk_failed, Id, Reason}}, D#d{inflight = In}, Id);
                false ->
                    chunk_event(retry, Id, D, #{attempt => A, reason => Reason}),
                    {keep_state, dispatch(D#d{inflight = In, queue = [{Id, First, N} | D#d.queue],
                                              attempts = (D#d.attempts)#{Id => A}})}
            end;
        error -> keep_state_and_data
    end;
running(cast, cancel, D) ->
    pq_nif:cancel_token_set(D#d.token),
    stop_with({error, cancelled}, D, cancel);
running({timeout, run}, expired, #d{partial = true, done = Done} = D) ->
    pq_nif:cancel_token_set(D#d.token),
    Sorted = lists:sort(maps:to_list(Done)),
    Partial = #{samples => lists:append([S || {_, S} <- Sorted]),
                shot_indices => lists:append([lists:seq(F, F + length(S) - 1) || {F, S} <- Sorted])},
    stop_with({error, {timeout, #{partial => Partial}}}, D, timeout);
running({timeout, run}, expired, D) ->
    pq_nif:cancel_token_set(D#d.token),
    stop_with({error, timeout}, D, timeout).

dispatch(#d{queue = [{Id, First, N} | Rest], inflight = In, max_inflight = Max} = D) when map_size(In) < Max ->
    Attempt = maps:get(Id, D#d.attempts, 0),
    %% a node per chunk, another one on each retry (FS 9.4)
    Node = lists:nth((Id + Attempt) rem length(D#d.nodes) + 1, D#d.nodes),
    {ok, Pid} = pq_worker_sup:start_worker(self(), {Id, First, N, (D#d.spec)#{target => Node}}, Attempt),
    pq_telemetry:execute([pq, chunk, start], #{system_time => erlang:system_time()},
                         #{job => self(), chunk => Id, first_shot => First, shots => N, node => Node, attempt => Attempt}),
    dispatch(D#d{queue = Rest, inflight = In#{Pid => {Id, First, N}},
                 started_chunks = (D#d.started_chunks)#{Id => erlang:monotonic_time()}});
dispatch(D) -> D.

finish_or_continue(#d{queue = [], inflight = In, done = Done, from = From} = D, _Id) when map_size(In) =:= 0 ->
    %% merged in shot order (FS 9.3, FR-P03)
    Samples = lists:append([S || {_, S} <- lists:sort(maps:to_list(Done))]),
    Reply = case map_size(D#d.extras) of
                0 -> {ok, Samples};
                _ -> {ok, Samples, [E || {_, E} <- lists:sort(maps:to_list(D#d.extras))]}
            end,
    gen_statem:reply(From, Reply),
    job_stop(D, ok),
    {stop, normal};
finish_or_continue(D, _Id) -> {keep_state, D}.

stop_with(Reply, D, _Id) ->
    pq_nif:cancel_token_set(D#d.token),
    %% not kill: a proxy of a remote chunk traps this and cancels the kernel on its node
    [exit(P, {shutdown, cancelled}) || P <- maps:keys(D#d.inflight)],
    gen_statem:reply(D#d.from, Reply),
    {error, R} = Reply,
    job_stop(D, element(1, if is_tuple(R) -> R; true -> {R} end)),
    {stop, normal}.

job_stop(D, Result) ->
    pq_telemetry:execute([pq, job, stop], #{duration => erlang:monotonic_time() - D#d.started},
                         #{job => self(), result => Result, chunks => D#d.chunks,
                           retries => lists:sum(maps:values(D#d.attempts))}).

chunk_event(stop, Id, D, Meta) ->
    T0 = maps:get(Id, D#d.started_chunks, erlang:monotonic_time()),
    pq_telemetry:execute([pq, chunk, stop], #{duration => erlang:monotonic_time() - T0}, Meta#{job => self(), chunk => Id});
chunk_event(retry, Id, _D, Meta) ->
    pq_telemetry:execute([pq, chunk, retry], #{count => 1}, Meta#{job => self(), chunk => Id}).

terminate(normal, _State, _D) -> ok;
terminate(Reason, _State, #d{} = D) ->
    pq_telemetry:execute([pq, job, exception], #{duration => erlang:monotonic_time() - D#d.started},
                         #{job => self(), reason => Reason}).
