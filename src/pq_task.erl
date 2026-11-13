%% The task behaviour: how an application gives PQCF work to distribute (decision record P1). A task is plain data,
%% {pq_task, Module, Args}, and its result a pure function of it, so the result is the same on whichever node runs it
%% and however the work is spread (D1). Module implements this behaviour:
%%
%%   run(Args, Ctx) -> Result                  the work itself; Ctx: #{home => node()} (the coordinating node, for
%%                                             fetching data from its store, pq_rstore)
%%   resources(Args) -> #{cpu => N, gpu => N, memory_mb => N, data => [Key]}     (optional) what the task needs,
%%                                             for scheduling; default #{}
%%   describe(Args) -> #{atom() => term()}     (optional) what the task is, for its provenance record; default #{}
%%
%% A node runs only the task modules it allows: PQCF's own (pq_ff_task), those in the application environment
%% (pq, task_modules) and those registered with allow/1. A task for another module fails with
%% {task_not_allowed, Module} on that node: an executor never runs code it was not told to.
%%
%% PQCF's own tasks, the finite-field reconstruction's {column, ...} etc. tuples, are pq_ff_task's arguments;
%% the bare tuples are still accepted and mean {pq_task, pq_ff_task, Tuple}.
-module(pq_task).
-export([normalise/1, run/2, run_traced/2, allow/1, allowed/0, resources/1]).
-export_type([task/0, provenance/0]).

-callback run(Args :: term(), Ctx :: #{home := node()}) -> term().
-callback resources(Args :: term()) -> map().
-callback describe(Args :: term()) -> map().
-optional_callbacks([resources/1, describe/1]).

-type task() :: {pq_task, module(), term()} | tuple().
-type provenance() :: #{task := binary(), module := module(), code := binary(), describe := map(), node := node(),
                        pq_version := binary(), started := integer(), seconds := float()}.

-define(BUILTIN, [pq_ff_task]).
-define(KEY, {?MODULE, allowed}).

%% {Module, Args} of a task
-spec normalise(task()) -> {module(), term()}.
normalise({pq_task, Mod, Args}) when is_atom(Mod) -> {Mod, Args};
normalise(Tuple) when is_tuple(Tuple) -> {pq_ff_task, Tuple}.

-spec allow(module()) -> ok.
allow(Mod) when is_atom(Mod) ->
    persistent_term:put(?KEY, lists:usort([Mod | persistent_term:get(?KEY, [])])).

-spec allowed() -> [module()].
allowed() ->
    lists:usort(?BUILTIN ++ application:get_env(pq, task_modules, []) ++ persistent_term:get(?KEY, [])).

%% The task's result, on this node.
-spec run(task(), node()) -> term().
run(Task, Home) ->
    {Mod, Args} = normalise(Task),
    lists:member(Mod, allowed()) orelse error({task_not_allowed, Mod}),
    Mod:run(Args, #{home => Home}).

%% The task's result and its provenance record: what ran (the task's hash, the module and its code's checksum, the
%% module's own description), where (node, PQCF's version) and when (start in milliseconds since the epoch, wall
%% time in seconds). The result is the same as run/2's; only the record carries times.
-spec run_traced(task(), node()) -> {term(), provenance()}.
run_traced(Task, Home) ->
    {Mod, Args} = normalise(Task),
    Started = erlang:system_time(millisecond),
    T0 = erlang:monotonic_time(),
    Result = run(Task, Home),
    Secs = erlang:convert_time_unit(erlang:monotonic_time() - T0, native, microsecond) / 1.0e6,
    {Result, #{task => hex(crypto:hash(sha256, term_to_binary(Task, [deterministic]))),
               module => Mod, code => hex(Mod:module_info(md5)), describe => optional(Mod, describe, Args),
               node => node(), pq_version => version(), started => Started, seconds => Secs}}.

-spec resources(task()) -> map().
resources(Task) ->
    {Mod, Args} = normalise(Task),
    optional(Mod, resources, Args).

optional(Mod, F, Args) ->
    _ = code:ensure_loaded(Mod),
    case erlang:function_exported(Mod, F, 1) of
        true -> Mod:F(Args);
        false -> #{}
    end.

version() ->
    case application:get_key(pq, vsn) of
        {ok, V} -> list_to_binary(V);
        undefined -> _ = application:load(pq), case application:get_key(pq, vsn) of {ok, V} -> list_to_binary(V); _ -> <<"unknown">> end
    end.

hex(Bin) -> binary:encode_hex(Bin, lowercase).
