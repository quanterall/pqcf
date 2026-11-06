%% The top supervisor (FS 2.2): the pg scope pq, the node registry, the spec cache, the reconstruction store (R1),
%% pq_worker_sup and pq_job_sup.
-module(pq_sup).
-behaviour(supervisor).
-export([start_link/0, init/1]).

start_link() ->
    supervisor:start_link({local, ?MODULE}, ?MODULE, []).

init([]) ->
    Sup = fun(M) -> #{id => M, start => {M, start_link, []}, type => supervisor} end,
    Worker = fun(M) -> #{id => M, start => {M, start_link, []}, type => worker} end,
    Pg = #{id => pg_pq, start => {pg, start_link, [pq]}, type => worker},
    {ok, {#{strategy => one_for_one, intensity => 5, period => 10},
          [Pg, Worker(pq_node_registry), Worker(pq_spec_cache), Worker(pq_rstore), Sup(pq_worker_sup), Sup(pq_job_sup)]}}.
