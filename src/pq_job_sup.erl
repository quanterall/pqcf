-module(pq_job_sup).
-behaviour(supervisor).
-export([start_link/0, init/1, start_job/3]).

start_link() -> supervisor:start_link({local, ?MODULE}, ?MODULE, []).

init([]) ->
    {ok, {#{strategy => simple_one_for_one, intensity => 100, period => 1},
          [#{id => pq_job, start => {pq_job, start_link, []}, restart => temporary, type => worker}]}}.

start_job(Spec, Opts, From) -> supervisor:start_child(?MODULE, [Spec, Opts, From]).
