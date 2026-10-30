%% The backend behaviour (FS 5.3).
-module(pq_backend).

-type capability() :: exact | truncated | terminal | midcircuit | unsupported.
-export_type([capability/0]).

-callback name() -> atom().
-callback capabilities() -> #{atom() => capability()}.
-callback init(Input :: term(), Circuit :: tuple(), Opts :: map()) -> {ok, State :: term()} | {error, term()}.
-callback apply(Instruction :: tuple(), State :: term(), Opts :: map()) ->
              {ok, State :: term(), Warnings :: [{atom(), term()}]} | {error, term()}.
%% the terminal measurements of a run: samples, one row per shot, one column per measured mode, and metadata
%% to add to the result (for example column_kinds)
-callback sample(Measurements :: [tuple()], State :: term(), Opts :: map()) ->
              {ok, Kind :: atom(), Modes :: [non_neg_integer()], Samples :: [list()], Meta :: map()} | {error, term()}.
%% optional: the public state of evolve/3 (when the backend keeps more, such as a cumulative norm deficit)
-callback final(State :: term()) -> term().
-optional_callbacks([final/1]).
