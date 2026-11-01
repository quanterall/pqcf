%% NFR-02 benchmark (FS 11): parallel efficiency T1 / (p * Tp) of 100,000 boson-sampling shots on one node.
%% Configuration (fixed in Phase 1): 6 photons in 12 modes, a fixed interferometer of nearest-neighbour
%% beamsplitters, clifford_clifford, seed 1. Run: erl -pa _build/default/lib/*/ebin -noshell -s nfr02 main
-module(nfr02).
-include("../include/pq.hrl").
-export([main/0]).

main() ->
    {ok, _} = application:ensure_all_started(pq),
    {ok, C0} = pq_circuit:new(12, {fock, lists:duplicate(6, 1) ++ lists:duplicate(6, 0)}),
    Layer = fun(Off) -> [{beamsplitter, J, J + 1, 0.3 + 0.1 * J, 0.2 * J} || J <- lists:seq(Off, 10, 2)] end,
    C = C0#circuit{instrs = lists:append([Layer(L rem 2) || L <- lists:seq(0, 11)]) ++ [{measure_pnr, lists:seq(0, 11)}]},
    Max = erlang:system_info(dirty_cpu_schedulers_online),
    Ps = lists:usort([P || P <- [1, 2, 4, 8, 16], P =< Max] ++ [Max]),
    %% the best of 5 runs per p: the minimum is the least disturbed by other load on a shared machine
    Run = fun(P) ->
                  Ts = [begin
                            {T, {ok, #result{samples = S}}} =
                                timer:tc(pq_sampling, sample, [boson, C, #{shots => 100000, seed => 1, max_inflight => P}]),
                            {T / 1.0e6, erlang:phash2(S)}
                        end || _ <- lists:seq(1, 5)],
                  {lists:min([T || {T, _} <- Ts]), element(2, hd(Ts))}
          end,
    _ = Run(Max),                                        % warm-up
    Results = [{P, Run(P)} || P <- Ps],
    {_, {T1, H1}} = hd(Results),
    io:format("dirty CPU schedulers online: ~p~n", [Max]),
    [io:format("p = ~2b: ~.2f s, efficiency ~.2f, samples identical to p = 1: ~p~n", [P, T, T1 / (P * T), H =:= H1])
     || {P, {T, H}} <- Results],
    halt().
