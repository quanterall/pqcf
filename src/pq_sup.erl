%% The top supervisor (FS 2.2): pq_job_sup and pq_worker_sup; pq_config and pq_node_registry arrive with
%% distribution (Phase 2).
-module(pq_sup).
-behaviour(supervisor).
-export([start_link/0, init/1]).

start_link() ->
    supervisor:start_link({local, ?MODULE}, ?MODULE, []).

init([]) ->
    Sup = fun(M) -> #{id => M, start => {M, start_link, []}, type => supervisor} end,
    {ok, {#{strategy => one_for_one, intensity => 5, period => 10}, [Sup(pq_worker_sup), Sup(pq_job_sup)]}}.
