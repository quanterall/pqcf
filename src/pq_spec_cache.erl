%% The node-local cache of run specifications (FS 9.5): a job's spec is sent to a node once, chunk requests carry
%% only its hash. Bounded by spec_cache_bytes; the oldest entries go first.
-module(pq_spec_cache).
-behaviour(gen_server).
-export([start_link/0, put/2, get/1, hash/1]).
-export([init/1, handle_call/3, handle_cast/2]).

-define(T, pq_spec_cache).

start_link() -> gen_server:start_link({local, ?MODULE}, ?MODULE, [], []).

-spec hash(map()) -> binary().
hash(Spec) -> crypto:hash(sha256, term_to_binary(maps:remove(cancel, Spec), [deterministic])).

-spec put(binary(), map()) -> ok.
put(Hash, Spec) -> gen_server:call(?MODULE, {put, Hash, Spec}).

-spec get(binary()) -> {ok, map()} | error.
get(Hash) ->
    case ets:lookup(?T, Hash) of
        [{_, Spec, _, _}] -> {ok, Spec};
        [] -> error
    end.

init([]) ->
    ?T = ets:new(?T, [named_table, set, protected, {read_concurrency, true}]),
    %% the cancel tokens of chunks running here for jobs on other nodes (pq_worker:run_cached/4)
    pq_remote_cancel = ets:new(pq_remote_cancel, [named_table, set, public]),
    {ok, #{bytes => 0, seq => 0}}.

handle_call({put, Hash, Spec}, _From, #{bytes := B, seq := Q} = S) ->
    case ets:member(?T, Hash) of
        true -> {reply, ok, S};
        false ->
            Size = erlang:external_size(Spec),
            ets:insert(?T, {Hash, Spec, Size, Q}),
            {reply, ok, evict(S#{bytes := B + Size, seq := Q + 1})}
    end.

handle_cast(_, S) -> {noreply, S}.

evict(#{bytes := B} = S) ->
    Max = pq_core:get(spec_cache_bytes),
    case B > Max andalso ets:info(?T, size) > 1 of
        false -> S;
        true ->
            [{H, _, Sz, _} | _] = lists:keysort(4, ets:tab2list(?T)),
            ets:delete(?T, H),
            evict(S#{bytes := B - Sz})
    end.
