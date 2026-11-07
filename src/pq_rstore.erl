%% The minimal store of reconstruction jobs (decision record R1, phase 3): blobs (programs, supports, results) keyed
%% by the SHA-256 of their deterministic external form, so the same data has the same key on every node. In memory
%% (ETS), and also on disk when the application environment sets store_dir. A worker node fetches a blob it lacks
%% from the coordinator once, by key, and keeps it; programs are also kept prepared for the field evaluator
%% (pq_ff), so a node prepares each program once.
%%
%% The tables belong to this server, under pq_sup; on a node where the pq application doesn't run, the first call
%% starts it, so a single node needs nothing set up.
-module(pq_rstore).
-compile({no_auto_import, [get/1, put/1]}).
-behaviour(gen_server).
-export([start_link/0, put/1, get/1, fetch/2, prepared/2, key/1]).
-export([init/1, handle_call/3, handle_cast/2]).

-define(BLOBS, pq_rstore_blobs).
-define(PREPARED, pq_rstore_prepared).

%% under pq_sup; a store already started on first use (before the application) is adopted
start_link() ->
    case whereis(?MODULE) of
        undefined -> gen_server:start_link({local, ?MODULE}, ?MODULE, [], []);
        Pid -> link(Pid), {ok, Pid}
    end.

init([]) ->
    _ = ets:new(?BLOBS, [named_table, public, set, {read_concurrency, true}]),
    _ = ets:new(?PREPARED, [named_table, public, set, {read_concurrency, true}]),
    {ok, #{}}.

handle_call(ping, _From, S) -> {reply, ok, S}.
handle_cast(_, S) -> {noreply, S}.

%% the tables exist: the one init creates last is checked, and after a start -- or a start that another caller
%% got to first, whose init may still be running -- a call returns only once init is done (two tasks arriving
%% together on a fresh node raced here)
ensure() ->
    case ets:whereis(?PREPARED) of
        undefined ->
            case gen_server:start({local, ?MODULE}, ?MODULE, [], []) of
                {ok, _} -> ok;
                {error, {already_started, _}} -> ok
            end,
            ok = gen_server:call(?MODULE, ping, 60000);
        _ -> ok
    end.

%% the key of a term: SHA-256 of its deterministic external form
-spec key(term()) -> binary().
key(Term) -> crypto:hash(sha256, term_to_binary(Term, [deterministic])).

%% stores the term; its key
-spec put(term()) -> binary().
put(Term) ->
    ensure(),
    Bin = term_to_binary(Term, [deterministic]),
    Key = crypto:hash(sha256, Bin),
    ets:insert(?BLOBS, {Key, Term}),
    case application:get_env(pq, store_dir) of
        {ok, Dir} ->
            File = filename:join(Dir, binary:encode_hex(Key, lowercase)),
            case filelib:is_regular(File) of
                true -> ok;
                false -> ok = filelib:ensure_dir(File), ok = file:write_file(File, Bin)
            end;
        undefined -> ok
    end,
    Key.

-spec get(binary()) -> {ok, term()} | not_found.
get(Key) ->
    ensure(),
    case ets:lookup(?BLOBS, Key) of
        [{_, Term}] -> {ok, Term};
        [] ->
            case application:get_env(pq, store_dir) of
                {ok, Dir} ->
                    case file:read_file(filename:join(Dir, binary:encode_hex(Key, lowercase))) of
                        {ok, Bin} ->
                            %% a file that doesn't hash to its name is ignored
                            case crypto:hash(sha256, Bin) of
                                Key -> Term = binary_to_term(Bin, [safe]), ets:insert(?BLOBS, {Key, Term}), {ok, Term};
                                _ -> not_found
                            end;
                        _ -> not_found
                    end;
                undefined -> not_found
            end
    end.

%% the blob, from this node or else from Home (kept here afterwards)
-spec fetch(binary(), node()) -> {ok, term()} | not_found.
fetch(Key, Home) ->
    case get(Key) of
        {ok, _} = Found -> Found;
        not_found when Home =:= node() -> not_found;
        not_found ->
            case erpc:call(Home, ?MODULE, get, [Key], 60000) of
                {ok, Term} = Found ->
                    %% stored only if it really is the blob asked for
                    case key(Term) of
                        Key -> ets:insert(?BLOBS, {Key, Term}), Found;
                        _ -> not_found
                    end;
                not_found -> not_found
            end
    end.

%% the black box of a task prepared for pq_ff, once per node: a program by its key, or
%% {linsys, Key, PatternKey, N, M, Want, Col}, the program's outputs solved as A u = B at each point (A's sparsity
%% pattern a blob of its own), u[Want] for right-hand side Col its value (R1 phase 5)
-spec prepared(binary() | {linsys, binary(), binary(), pos_integer(), pos_integer(), non_neg_integer(), non_neg_integer()},
               node()) -> {ok, tuple()} | {error, term()}.
prepared({linsys, Key, PKey, N, M, Want, Col} = BB, Home) ->
    ensure(),
    case ets:lookup(?PREPARED, BB) of
        [{_, FF}] -> {ok, FF};
        [] ->
            case {prepared(Key, Home), fetch(PKey, Home)} of
                {{ok, {ff, Ref, Inputs, _}}, {ok, Pattern}} ->
                    case pq_nif:ff_linsys(Ref, N, M, Want, Col, Pattern) of
                        {ok, R} -> FF = {ff, R, Inputs, 1}, ets:insert(?PREPARED, {BB, FF}), {ok, FF};
                        Error -> Error
                    end;
                {{ok, _}, not_found} -> {error, {not_found, PKey}};
                {Error, _} -> Error
            end
    end;
prepared(Key, Home) ->
    ensure(),
    case ets:lookup(?PREPARED, Key) of
        [{_, FF}] -> {ok, FF};
        [] ->
            case fetch(Key, Home) of
                {ok, Program} ->
                    case pq_ff:prepare(Program) of
                        {ok, FF} -> ets:insert(?PREPARED, {Key, FF}), {ok, FF};
                        Error -> Error
                    end;
                not_found -> {error, {not_found, Key}}
            end
    end.
