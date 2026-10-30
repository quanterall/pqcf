%% One run's job (FS 2.2, 9.1-9.4): splits the shots into chunks, dispatches at most max_inflight at a time to
%% temporary workers, retries a failed chunk up to max_retries, merges in shot order and replies to the caller.
%% Cancellation sets the run's cancel token (kernels poll it); a timeout cancels and fails with timeout.
-module(pq_job).
-behaviour(gen_statem).
-export([start_link/3, cancel/1]).
-export([init/1, callback_mode/0, running/3, terminate/3]).

-record(d, {spec, from, queue = [], inflight = #{}, attempts = #{}, done = #{}, max_inflight, max_retries,
            token, chunks}).

start_link(Spec, Opts, From) -> gen_statem:start_link(?MODULE, {Spec, Opts, From}, []).

cancel(Job) -> gen_statem:cast(Job, cancel).

callback_mode() -> state_functions.

init({Spec, Opts, From}) ->
    process_flag(trap_exit, true),
    Token = pq_nif:cancel_token_new(),
    Chunks = maps:get(chunks, Opts),
    D = #d{spec = Spec#{cancel => Token}, from = From, queue = Chunks, chunks = length(Chunks),
           max_inflight = maps:get(max_inflight, Opts), max_retries = maps:get(max_retries, Opts), token = Token},
    Actions = case maps:get(timeout, Opts, infinity) of
                  infinity -> [];
                  T -> [{{timeout, run}, T, expired}]
              end,
    {ok, running, dispatch(D), Actions}.

running(info, {chunk_done, Pid, Id, {ok, Samples}}, D) ->
    {Id, First, _} = maps:get(Pid, D#d.inflight),
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
                    telemetry_retry(Id, Reason),
                    {keep_state, dispatch(D#d{inflight = In, queue = [{Id, First, N} | D#d.queue],
                                              attempts = (D#d.attempts)#{Id => A}})}
            end;
        error -> keep_state_and_data
    end;
running(cast, cancel, D) ->
    pq_nif:cancel_token_set(D#d.token),
    stop_with({error, cancelled}, D, cancel);
running({timeout, run}, expired, D) ->
    pq_nif:cancel_token_set(D#d.token),
    stop_with({error, timeout}, D, timeout).

dispatch(#d{queue = [{Id, First, N} | Rest], inflight = In, max_inflight = Max} = D) when map_size(In) < Max ->
    Attempt = maps:get(Id, D#d.attempts, 0),
    {ok, Pid} = pq_worker_sup:start_worker(self(), {Id, First, N, D#d.spec}, Attempt),
    dispatch(D#d{queue = Rest, inflight = In#{Pid => {Id, First, N}}});
dispatch(D) -> D.

finish_or_continue(#d{queue = [], inflight = In, done = Done, from = From}, _Id) when map_size(In) =:= 0 ->
    %% merged in shot order (FS 9.3, FR-P03)
    Samples = lists:append([S || {_, S} <- lists:sort(maps:to_list(Done))]),
    gen_statem:reply(From, {ok, Samples}),
    {stop, normal};
finish_or_continue(D, _Id) -> {keep_state, D}.

stop_with(Reply, D, _Id) ->
    pq_nif:cancel_token_set(D#d.token),
    [exit(P, kill) || P <- maps:keys(D#d.inflight)],
    gen_statem:reply(D#d.from, Reply),
    {stop, normal}.

telemetry_retry(_Id, _Reason) -> ok.

terminate(_Reason, _State, _D) -> ok.
