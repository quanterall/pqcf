-module(pq_worker_sup).
-behaviour(supervisor).
-export([start_link/0, init/1, start_worker/3]).

start_link() -> supervisor:start_link({local, ?MODULE}, ?MODULE, []).

init([]) ->
    {ok, {#{strategy => simple_one_for_one, intensity => 1000, period => 1},
          [#{id => pq_worker, start => {pq_worker, start_link, []}, restart => temporary, type => worker}]}}.

start_worker(Job, Chunk, Attempt) -> supervisor:start_child(?MODULE, [Job, Chunk, Attempt]).
