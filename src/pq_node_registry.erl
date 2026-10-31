%% Node registration and the handshake (FS 9.5, FR-P02): nodes join pg scope `pq`, group `workers`; a job takes a
%% node only if its pq version, NIF {Major, Minor} and OTP release equal the job's node's, and its target triple
%% too unless allow_heterogeneous_arch is set (then only D2 holds and results say heterogeneous => true).
-module(pq_node_registry).
-behaviour(gen_server).
-export([start_link/0, register_node/0, unregister_node/0, info/0, eligible_nodes/1]).
-export([init/1, handle_call/3, handle_cast/2]).

start_link() -> gen_server:start_link({local, ?MODULE}, ?MODULE, [], []).

%% join the pool of worker nodes
-spec register_node() -> ok.
register_node() -> gen_server:call(?MODULE, register).

-spec unregister_node() -> ok.
unregister_node() -> gen_server:call(?MODULE, unregister).

%% this node's handshake data
-spec info() -> map().
info() ->
    {Ma, Mi, Pa, Triple} = pq_nif:nif_version(),
    {ok, Vsn} = application:get_key(pq, vsn),
    #{pq => Vsn, nif => {Ma, Mi, Pa}, triple => Triple, otp => erlang:system_info(otp_release)}.

%% the registered nodes that pass the handshake against this node, this node first if registered; and whether any
%% of them has another target triple
-spec eligible_nodes(map()) -> {[node()], boolean()}.
eligible_nodes(Opts) ->
    Me = info(),
    Hetero = pq_core:get(allow_heterogeneous_arch, Opts),
    Nodes = lists:usort([node(P) || P <- pg:get_members(pq, workers)]),
    Checked = [{N, I} || N <- Nodes, {ok, I} <- [remote_info(N)], compatible(Me, I, Hetero)],
    Ordered = [N || {N, _} <- Checked, N =:= node()] ++ [N || {N, _} <- Checked, N =/= node()],
    {Ordered, lists:any(fun({_, I}) -> maps:get(triple, I) =/= maps:get(triple, Me) end, Checked)}.

remote_info(N) when N =:= node() -> {ok, info()};
remote_info(N) ->
    try {ok, erpc:call(N, ?MODULE, info, [], 5000)}
    catch _:_ -> error
    end.

compatible(#{pq := V, nif := {Ma, Mi, _}, otp := O, triple := T}, #{pq := V, nif := {Ma, Mi, _}, otp := O, triple := T2}, Hetero) ->
    T =:= T2 orelse Hetero;
compatible(_, _, _) -> false.

init([]) -> {ok, #{}}.

handle_call(register, _From, S) ->
    ok = pg:join(pq, workers, self()),
    {reply, ok, S};
handle_call(unregister, _From, S) ->
    _ = pg:leave(pq, workers, self()),
    {reply, ok, S}.

handle_cast(_, S) -> {noreply, S}.
