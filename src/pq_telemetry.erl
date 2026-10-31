%% Observability (FS 9.6): events through the `telemetry` library when it is loaded, nothing otherwise. A soft
%% dependency, so a numerical library does not pull it in: applications that want the events add `telemetry`.
%% Events: [pq, job, start | stop | exception], [pq, chunk, start | stop | retry].
-module(pq_telemetry).
-export([execute/3]).

-spec execute([atom()], map(), map()) -> ok.
execute(Event, Measurements, Meta) ->
    case available() of
        true -> apply(telemetry, execute, [Event, Measurements, Meta]);
        false -> ok
    end.

%% checked on every call (cheap once loaded), so starting telemetry after pq still works
available() ->
    erlang:function_exported(telemetry, execute, 3) orelse
        (code:ensure_loaded(telemetry) =:= {module, telemetry} andalso erlang:function_exported(telemetry, execute, 3)).
