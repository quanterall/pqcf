%% The public facade (FS 7, 7.1). The backend is chosen from the circuit unless `backend` is given: gaussian when the
%% input is vacuum or Gaussian and every instruction is Gaussian-supported (FS 5.3), otherwise fock, or fock_dm when
%% the circuit has loss; `boson` (pq_sampling) for passive circuits with a Fock input.
-module(pq_api).
-include("pq.hrl").
-export([circuit/1, circuit/2, run/2, run/3, sample/2, sample/3, evolve/2, evolve/3, probabilities/2, run_batch/3,
         backend_for/2]).

-type result(T) :: {ok, T} | {error, pq_core:reason()}.

-spec circuit(pos_integer()) -> result(#circuit{}).
circuit(N) -> pq_circuit:new(N).

-spec circuit(pos_integer(), pq_circuit:input_spec()) -> result(#circuit{}).
circuit(N, Input) -> pq_circuit:new(N, Input).

-spec run(#circuit{}, map()) -> result(#result{}).
run(C, Opts) ->
    case backend_for(C, Opts) of
        boson -> pq_sampling:sample(boson, C, Opts);
        B -> pq_simulator:run(backend_module(B), C, Opts)
    end.

%% FS 7: the three-argument forms take the backend first, like pq_simulator and pq_sampling
-spec run(gaussian | fock | fock_dm | boson, #circuit{}, map()) -> result(#result{}).
run(Backend, C, Opts) -> run(C, Opts#{backend => Backend}).

-spec sample(gaussian | fock | fock_dm | boson, #circuit{}, map()) -> result(#result{}).
sample(Backend, C, Opts) -> sample(C, Opts#{backend => Backend}).

-spec evolve(gaussian | fock | fock_dm, #circuit{}, map()) -> {ok, term(), [{atom(), term()}]} | {error, pq_core:reason()}.
evolve(Backend, C, Opts) -> evolve(C, Opts#{backend => Backend}).

%% the same as run/2, with shots required
-spec sample(#circuit{}, map()) -> result(#result{}).
sample(C, #{shots := _} = Opts) -> run(C, Opts);
sample(_C, _Opts) -> {error, {unsupported_input, shots_required}}.

-spec evolve(#circuit{}, map()) -> {ok, term(), [{atom(), term()}]} | {error, pq_core:reason()}.
evolve(C, Opts) ->
    case backend_for(C, Opts) of
        boson -> {error, {unsupported_input, #{backend => boson, evolve => no_state}}};
        B -> pq_simulator:evolve(backend_module(B), C, Opts)
    end.

%% FS 7.1: exact probabilities of the PNR (or threshold) outcomes of the measured modes, within max_enumeration
-spec probabilities(#circuit{}, map()) -> {ok, #{[non_neg_integer()] => float()}, map()} | {error, pq_core:reason()}.
probabilities(#circuit{instrs = Is} = C, Opts) ->
    Ms = [I || I <- Is, lists:member(element(1, I), [measure_pnr, measure_threshold])],
    Modes = lists:append([Mo || {_, Mo} <- Ms]),
    Kind = case lists:usort([element(1, I) || I <- Ms]) of
               [measure_threshold] -> threshold;
               _ -> pnr
           end,
    Unitary = C#circuit{instrs = [I || I <- Is, not lists:member(element(1, I), [measure_pnr, measure_threshold, measure_homodyne])]},
    case Modes of
        [] -> {error, {unsupported_input, no_measurement}};
        _ -> probs(backend_for(C, Opts), Unitary, Modes, Kind, Opts)
    end.

probs(B, C, Modes, Kind, Opts) when B =:= fock; B =:= fock_dm ->
    case pq_simulator:evolve(backend_module(B), C, Opts) of
        {ok, S, W} ->
            {State, N, Cut} = case S of
                                  #fock_state{amps = A, n_modes = N0, cutoff = C0} -> {A, N0, C0};
                                  #fock_dm{rho = R, n_modes = N0, cutoff = C0} -> {R, N0, C0}
                              end,
            {ok, PP} = pq_nif:fock_pnr_probabilities(State, Modes, N, Cut),
            Map = lists:foldl(fun({P, Pr}, Acc) ->
                                      K = case Kind of threshold -> [min(X, 1) || X <- P]; pnr -> P end,
                                      maps:update_with(K, fun(V) -> V + Pr end, Pr, Acc)
                              end, #{}, PP),
            {ok, Map, #{backend => B, cutoff => Cut, warnings => W}};
        E -> E
    end;
probs(gaussian, C, Modes, Kind, Opts) ->
    case pq_simulator:evolve(pq_gaussian, C, Opts) of
        {ok, S, W} ->
            {ok, Red} = pq_state:reduce(S, Modes),
            K = length(Modes),
            MaxN = pq_core:get(max_photons, Opts),
            Pats = case Kind of
                       threshold -> [P || P <- patterns_upto(K, K), lists:max([0 | P]) =< 1];
                       pnr -> patterns_upto(K, MaxN)
                   end,
            Max = pq_core:get(max_enumeration, Opts),
            case length(Pats) =< Max of
                false -> {error, {dimension_limit, length(Pats), Max}};
                true ->
                    Map = maps:from_list([{P, element(2, pq_nif:gaussian_probability(Red#gaussian_state.mean, Red#gaussian_state.cov,
                                                                                    float(Red#gaussian_state.hbar), P, Kind))}
                                          || P <- Pats]),
                    Missing = 1 - lists:sum(maps:values(Map)),
                    {ok, Map, #{backend => gaussian, truncated_mass => max(Missing, 0.0), warnings => W}}
            end;
        E -> E
    end;
probs(boson, C, Modes, _Kind, Opts) ->
    case pq_circuit:lower(C, maps:get(bindings, Opts, #{})) of
        {ok, #circuit{input = {fock, In}, n_modes = N} = L} ->
            U = pq_sampling:passive_unitary(L),
            Pats = patterns_exact(N, lists:sum(In)),
            Max = pq_core:get(max_enumeration, Opts),
            case length(Pats) =< Max of
                false -> {error, {dimension_limit, length(Pats), Max}};
                true ->
                    Full = [{P, element(1, element(2, pq_nif:bs_probability(U, In, P, #{})))} || P <- Pats],
                    Map = lists:foldl(fun({P, Pr}, Acc) ->
                                              Key = [lists:nth(M + 1, P) || M <- Modes],
                                              maps:update_with(Key, fun(V) -> V + Pr end, Pr, Acc)
                                      end, #{}, Full),
                    {ok, Map, #{backend => boson}}
            end;
        {ok, _} -> {error, {unsupported_input, boson_sampling_needs_a_fock_input}};
        E -> E
    end.

%% FR-V02: one run per binding, lowered per binding, dispatched in parallel; results in binding order
-spec run_batch(#circuit{}, [map()], map()) -> result([#result{}]).
run_batch(C, Bindings, Opts) ->
    Self = self(),
    Refs = [begin
                Ref = make_ref(),
                spawn_link(fun() -> Self ! {Ref, run(C, Opts#{bindings => B})} end),
                Ref
            end || B <- Bindings],
    Results = [receive {Ref, R} -> R end || Ref <- Refs],
    case [E || {error, _} = E <- Results] of
        [] -> {ok, [R || {ok, R} <- Results]};
        [E | _] -> E
    end.

%% FS 7.1's choice of backend
-spec backend_for(#circuit{}, map()) -> gaussian | fock | fock_dm | boson.
backend_for(#circuit{input = In, instrs = Is}, Opts) ->
    case maps:find(backend, Opts) of
        {ok, B} -> B;
        error ->
            Kinds = [element(1, I) || I <- Is],
            Gaussian = (In =:= vacuum orelse element(1, In) =:= gaussian) andalso
                lists:all(fun(K) -> maps:get(K, pq_gaussian:capabilities(), unsupported) =/= unsupported end, Kinds),
            if
                Gaussian -> gaussian;
                true ->
                    case lists:member(loss, Kinds) of
                        true -> fock_dm;
                        false -> fock
                    end
            end
    end.

backend_module(gaussian) -> pq_gaussian;
backend_module(fock) -> pq_fock;
backend_module(fock_dm) -> pq_fock_dm.

patterns_upto(K, Max) -> lists:append([patterns_exact(K, N) || N <- lists:seq(0, Max)]).

patterns_exact(1, N) -> [[N]];
patterns_exact(K, N) -> [[A | R] || A <- lists:seq(0, N), R <- patterns_exact(K - 1, N - A)].
