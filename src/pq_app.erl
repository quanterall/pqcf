-module(pq_app).
-behaviour(application).
-export([start/2, stop/1]).

start(_Type, _Args) ->
    case pq_nif:check_version() of
        ok -> pq_sup:start_link();
        {error, _} = Error -> Error
    end.

stop(_State) ->
    ok.
