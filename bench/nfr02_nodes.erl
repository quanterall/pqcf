%% NFR-02 across nodes (FS 11: at least 0.7 across 4 nodes): the same 100,000 boson-sampling shots on 1 peer node
%% and on 4, each node with S dirty schedulers (default 2), the calling node only dispatching; efficiency
%% T1 / (4 * T4), best of 5. Here the 4 nodes are peers on one host (the provisional reference): it measures the
%% distribution layer (erpc, the spec cache, merging), not a network.
%%   erl -sname main -pa _build/default/lib/*/ebin -pa <dir of this beam> -noshell -s nfr02_nodes main
-module(nfr02_nodes).
-include("../include/pq.hrl").
-export([main/0]).

main() ->
    {ok, _} = application:ensure_all_started(pq),
    S = list_to_integer(os:getenv("NODE_SCHEDULERS", "2")),
    Paths = lists:append([["-pa", P] || P <- code:get_path()]),
    Peers = [begin
                 {ok, Pid, Node} = peer:start(#{name => list_to_atom("w" ++ integer_to_list(K)),
                                                args => Paths ++ ["+S", io_lib:format("~b:~b", [S, S]),
                                                                  "+SDcpu", io_lib:format("~b:~b", [S, S])]}),
                 {ok, _} = erpc:call(Node, application, ensure_all_started, [pq]),
                 ok = erpc:call(Node, pq_node_registry, register_node, []),
                 {Pid, Node}
             end || K <- lists:seq(1, 4)],
    Nodes = [N || {_, N} <- Peers],
    [run(Name, C, Nodes, S) || {Name, C} <- [{"6 photons in 12 modes (Phase 1)", circuit(6, 12)},
                                             {"10 photons in 20 modes", circuit(10, 20)}]],
    [peer:stop(P) || {P, _} <- Peers],
    halt().

circuit(Ph, M) ->
    {ok, C0} = pq_circuit:new(M, {fock, lists:duplicate(Ph, 1) ++ lists:duplicate(M - Ph, 0)}),
    Layer = fun(Off) -> [{beamsplitter, J, J + 1, 0.3 + 0.1 * J, 0.2 * J} || J <- lists:seq(Off, M - 2, 2)] end,
    C0#circuit{instrs = lists:append([Layer(L rem 2) || L <- lists:seq(0, M - 1)]) ++ [{measure_pnr, lists:seq(0, M - 1)}]}.

run(Name, C, Nodes, S) ->
    Time = fun(Ns) ->
                   Ts = [begin
                             {T, {ok, #result{samples = Smp}}} =
                                 timer:tc(pq_sampling, sample, [boson, C, #{shots => 100000, seed => 1, distributed => nodes,
                                                                            nodes => Ns, max_inflight => S * length(Ns)}]),
                             {T / 1.0e6, erlang:phash2(Smp)}
                         end || _ <- lists:seq(1, 5)],
                   {lists:min([T || {T, _} <- Ts]), element(2, hd(Ts))}
           end,
    _ = Time(Nodes),                                     % warm-up: the spec reaches every node's cache
    {T1, H1} = Time([hd(Nodes)]),
    {T4, H4} = Time(Nodes),
    io:format("~s: 1 node ~.2f s, 4 nodes ~.2f s, efficiency ~.2f, samples identical: ~p~n",
              [Name, T1, T4, T1 / (4 * T4), H1 =:= H4]).
