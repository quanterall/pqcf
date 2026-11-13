%% G1, PQCF's side: each graph of specs.json encoded as a Gaussian boson sampler (squeezing r_i on every mode, then
%% the interferometer U, then photon counting on every mode), sampled with PQCF's Gaussian backend, and the samples
%% post-selected on exactly k photons, at most one per mode: the k-vertex subgraphs. Writes gbs.json.
-module(g1_pqcf).
-export([main/1]).
-include("../../include/pq.hrl").

main(Dir) ->
    {ok, Bin} = file:read_file(filename:join(Dir, "specs.json")),
    Specs = json:decode(Bin),
    Out = [run(S) || S <- Specs],
    ok = file:write_file(filename:join(Dir, "gbs.json"), json:encode(Out)).

run(#{<<"name">> := Name, <<"modes">> := N, <<"r">> := Rs, <<"U">> := Rows, <<"k">> := K, <<"target">> := Target,
      <<"max_shots">> := MaxShots, <<"seed">> := Seed, <<"max_photons">> := Max} = Spec) ->
    {ok, U} = pq_core:pack({cmat, [[{complex, float(Re), float(Im)} || [Re, Im] <- Row] || Row <- Rows]}),
    Modes = lists:seq(0, N - 1),
    Instrs = [{squeezing, M, float(R), 0.0} || {M, R} <- lists:zip(Modes, Rs)]
        ++ [{unitary, Modes, U}, {measure_pnr, Modes}],
    {ok, C0} = pq_circuit:new(N),
    C = C0#circuit{instrs = Instrs},
    T0 = erlang:monotonic_time(),
    %% "collisions": keep every pattern with K photons, as photon counts (G3); otherwise collision-free ones, as the
    %% modes that clicked (G1)
    Keep = case maps:get(<<"collisions">>, Spec, false) of
               true -> fun(P) -> lists:sum(P) =:= K end;
               false -> fun(P) -> lists:sum(P) =:= K andalso lists:max(P) =< 1 end
           end,
    Form = case maps:get(<<"collisions">>, Spec, false) of
               true -> fun(P) -> P end;
               false -> fun(P) -> [I - 1 || {I, 1} <- lists:zip(lists:seq(1, N), P)] end
           end,
    {Kept, Shots, Totals, Trunc} = batches(C, {Keep, Form}, K, Target, MaxShots, Seed, Max, 0, [], 0, #{}, []),
    Secs = erlang:convert_time_unit(erlang:monotonic_time() - T0, native, microsecond) / 1.0e6,
    #{<<"name">> => Name, <<"k">> => K, <<"shots">> => Shots, <<"seconds">> => Secs,
      <<"truncated_mass">> => lists:max(Trunc),
      <<"photon_totals">> => maps:from_list([{integer_to_binary(T), X} || {T, X} <- maps:to_list(Totals)]),
      <<"samples">> => lists:sublist(lists:reverse(Kept), Target)}.

%% batches of 100,000 shots, each one call to the sampler (one prefix cache) with its own seed, until Target
%% k-photon collision-free samples are kept or MaxShots are spent
-define(BATCH, 100000).
batches(_C, _F, _K, Target, MaxShots, _Seed, _Max, _B, Kept, Shots, Totals, Trunc)
  when length(Kept) >= Target; Shots >= MaxShots ->
    {Kept, Shots, Totals, [0.0 | Trunc]};
batches(C, {Keep, Form} = F, K, Target, MaxShots, Seed, Max, B, Kept, Shots, Totals, Trunc) ->
    {ok, #result{samples = Samples, metadata = Meta}} =
        pq_simulator:run(pq_gaussian, C, #{shots => ?BATCH, seed => Seed * 100000 + B, max_photons => Max,
                                           parallel => false}),
    Totals1 = lists:foldl(fun(P, Acc) -> maps:update_with(lists:sum(P), fun(X) -> X + 1 end, 1, Acc) end,
                          Totals, Samples),
    New = [Form(P) || P <- Samples, Keep(P)],
    batches(C, F, K, Target, MaxShots, Seed, Max, B + 1, lists:reverse(New) ++ Kept, Shots + ?BATCH, Totals1,
            [maps:get(truncated_mass, Meta, 0.0) | Trunc]).
