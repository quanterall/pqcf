%% The top supervisor (FS 2.2). Its children (pq_config, pq_node_registry, pq_job_sup, pq_worker_sup) arrive with
%% the modules that need them.
-module(pq_sup).
-behaviour(supervisor).
-export([start_link/0, init/1]).

start_link() ->
    supervisor:start_link({local, ?MODULE}, ?MODULE, []).

init([]) ->
    {ok, {#{strategy => one_for_one, intensity => 5, period => 10}, []}}.
