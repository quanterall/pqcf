%% PQCF's side of the photonic head-to-head (README.md): every task of inputs.json through PQCF's public kernels and
%% its Fock backend, timed; writes pqcf.json in the format of peers.py.
-module(photonic_pqcf).
-export([main/1]).
-include("../../include/pq.hrl").

%% a size slower than this is the last one run for that task
-define(BUDGET, 120.0).

main(Dir) ->
    {ok, Bin} = file:read_file(filename:join(Dir, "inputs.json")),
    In = json:decode(Bin),
    Out = #{<<"tool">> => <<"PQCF">>, <<"version">> => version(),
            <<"permanent">> => sized(maps:get(<<"permanent">>, In), <<"n">>, fun permanent/1),
            <<"torontonian">> => sized(maps:get(<<"torontonian">>, In), <<"k">>, fun torontonian/1),
            <<"boson_sampling">> => sized(maps:get(<<"boson_sampling">>, In), <<"photons">>, fun boson/1),
            <<"boson_check">> => boson_check(maps:get(<<"boson_check">>, In)),
            <<"fock">> => fock(maps:get(<<"fock">>, In))},
    ok = file:write_file(filename:join(Dir, "pqcf.json"), json:encode(Out)).

version() ->
    {ok, V} = application:get_key(pq, vsn),
    list_to_binary(V).

%% each size in order until one exceeds the budget
%% one untimed call on the smallest case first, as for the peers (peers.py)
sized([C0 | _] = Cases, Key, F) -> _ = F(C0), sized(Cases, Key, F, []).
sized([], _, _, Acc) -> lists:reverse(Acc);
sized([C | Rest], Key, F, Acc) ->
    R = (F(C))#{<<"size">> => maps:get(Key, C)},
    case maps:get(<<"seconds">>, R) > ?BUDGET of
        true -> lists:reverse([R | Acc]);
        false -> sized(Rest, Key, F, [R | Acc])
    end.

%% the median of repeated runs: once if it takes a second or more, else up to 7 times within about two seconds
timed(F) ->
    {T1, V} = once(F),
    Ts = repeat(F, T1, [T1]),
    {median(Ts), V}.
repeat(_, _, Ts) when length(Ts) >= 7 -> Ts;
repeat(F, T1, Ts) when T1 * (length(Ts) + 1) < 2.0, T1 < 1.0 -> {T, _} = once(F), repeat(F, T1, [T | Ts]);
repeat(_, _, Ts) -> Ts.
once(F) ->
    T0 = erlang:monotonic_time(),
    V = F(),
    {erlang:convert_time_unit(erlang:monotonic_time() - T0, native, microsecond) / 1.0e6, V}.
median(Ts) -> lists:nth((length(Ts) + 1) div 2, lists:sort(Ts)).

cmat(Rows) ->
    {ok, M} = pq_core:pack({cmat, [[{complex, float(Re), float(Im)} || [Re, Im] <- Row] || Row <- Rows]}),
    M.

permanent(#{<<"matrix">> := Rows}) ->
    M = cmat(Rows),
    {T, {ok, {{Re, Im}, _}}} = timed(fun() -> pq_nif:permanent(M, #{threads => 1}) end),
    #{<<"seconds">> => T, <<"value">> => [Re, Im]}.

torontonian(#{<<"matrix">> := Rows}) ->
    M = cmat(Rows),
    {T, {ok, {V, _}}} = timed(fun() -> pq_nif:torontonian(M, #{threads => 1}) end),
    #{<<"seconds">> => T, <<"value">> => V}.

input(Photons, Modes) -> lists:duplicate(Photons, 1) ++ lists:duplicate(Modes - Photons, 0).

boson(#{<<"photons">> := N, <<"modes">> := M, <<"shots">> := Shots, <<"unitary">> := Rows}) ->
    U = cmat(Rows),
    {T, {ok, _}} = once(fun() -> pq_nif:sample_boson(U, input(N, M), Shots, #{seed => 1, threads => 1}) end),
    #{<<"seconds">> => T, <<"shots">> => Shots, <<"samples_per_second">> => Shots / T}.

%% the samples of the statistics check, as counts per output pattern
boson_check(#{<<"photons">> := N, <<"modes">> := M, <<"shots">> := Shots, <<"unitary">> := Rows}) ->
    {ok, {S, _}} = pq_nif:sample_boson(cmat(Rows), input(N, M), Shots, #{seed => 1}),
    {ok, {imat, Patterns}} = pq_core:unpack(S),
    counts(Patterns).

counts(Patterns) ->
    Count = lists:foldl(fun(P, Acc) -> maps:update_with(P, fun(C) -> C + 1 end, 1, Acc) end, #{}, Patterns),
    [#{<<"pattern">> => P, <<"count">> => C} || {P, C} <- lists:sort(maps:to_list(Count))].

gate([<<"phaseshift">>, M, P]) -> {phaseshift, M, float(P)};
gate([<<"beamsplitter">>, M1, M2, T, P]) -> {beamsplitter, M1, M2, float(T), float(P)};
gate([<<"squeezing">>, M, R, P]) -> {squeezing, M, float(R), float(P)};
gate([<<"displacement">>, M, R, P]) -> {displacement, M, float(R), float(P)};
gate([<<"kerr">>, M, K]) -> {kerr, M, float(K)}.

%% the CVQNN layer on vacuum at each total-photon cutoff; the mean photon number of every mode as the answer
fock(#{<<"modes">> := N, <<"gates">> := Gates, <<"total_cutoffs">> := Cuts}) ->
    {ok, C0} = pq_circuit:new(N, {fock, lists:duplicate(N, 0)}),
    C = C0#circuit{instrs = [gate(G) || G <- Gates]},
    [begin
         Opts = #{cutoff => Cut, max_fock_dim => 100000000},
         {T, {ok, S, _}} = timed(fun() -> pq_api:evolve(fock, C, Opts) end),
         #fock_state{amps = {cvec, D, _} = A} = S,
         Mean = [begin
                     {ok, PP} = pq_nif:fock_pnr_probabilities(A, [I], N, Cut),
                     lists:sum([K * P || {[K], P} <- PP])
                 end || I <- lists:seq(0, N - 1)],
         {ok, All} = pq_nif:fock_pnr_probabilities(A, [0], N, Cut),
         #{<<"cutoff">> => Cut, <<"truncation">> => <<"total">>, <<"dimension">> => D, <<"seconds">> => T,
           <<"mean_photons">> => Mean, <<"norm">> => lists:sum([P || {_, P} <- All])}
     end || Cut <- Cuts].
