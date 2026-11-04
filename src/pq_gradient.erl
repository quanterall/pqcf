%% Gradients (FR-V03, Decision 6): gradient(Circuit, Observable, Binding, Opts) with mode numeric (central finite
%% differences, always available; step, default 1e-6) or symbolic (passive circuits with a Fock input: the
%% observable as a canonical expression in the entries of the symbolic unitary, differentiated and evaluated;
%% anything else is differentiate_unsupported). Phase 4 adds adjoint (one forward and one backward pass on the
%% pure-state Fock backend, pq_gradient_adjoint) and parameter_shift (exact shift rules on the Fock backends,
%% pq_gradient_shift). A strategy module implementing this behaviour may be passed as strategy => Module; it must
%% say for which instruction kinds its rule holds.
-module(pq_gradient).
-include("pq.hrl").
-export([gradient/4, symbolic_gradient/3, value/3, supports/1]).
-export_type([observable/0]).

-callback supports(InstructionKind :: atom()) -> boolean().
-callback gradient(#circuit{}, Observable :: term(), Binding :: map(), Opts :: map()) ->
              {ok, #{atom() => float()}} | {error, term()}.

-type observable() :: {probability, [non_neg_integer()], [non_neg_integer()]} | {mean_photon_number, [non_neg_integer()]}
                    | {weighted_sum, [{pq_symbolic:sym_expr(), observable()}]}.

%% central differences hold for every instruction
supports(_) -> true.

%% the derivatives with respect to the symbols of Binding (or Opts symbols => [atom()]), at Binding
-spec gradient(#circuit{}, observable(), map(), map()) -> {ok, #{atom() => float()}} | {error, pq_core:reason()}.
gradient(C, Obs, Binding, Opts) ->
    Syms = maps:get(symbols, Opts, lists:sort(maps:keys(Binding))),
    case maps:find(strategy, Opts) of
        {ok, Mod} when Mod =/= ?MODULE ->
            case lists:all(fun(I) -> Mod:supports(element(1, I)) end, C#circuit.instrs) of
                true -> Mod:gradient(C, Obs, Binding, Opts);
                false -> {error, {differentiate_unsupported, {strategy, Mod}}}
            end;
        _ ->
            case maps:get(mode, Opts, numeric) of
                numeric -> numeric(C, Obs, Binding, Syms, Opts);
                adjoint -> pq_gradient_adjoint:gradient(C, Obs, Binding, Opts#{symbols => Syms});
                parameter_shift -> pq_gradient_shift:gradient(C, Obs, Binding, Opts#{symbols => Syms});
                symbolic ->
                    case symbolic_gradient(C, Obs, Syms) of
                        {ok, Exprs} ->
                            Vals = [{S, pq_symbolic:evaluate(E, maps:merge(numbers(C), Binding), #{domains => C#circuit.symbols})}
                                    || S := E <- Exprs],
                            case [R || {_, {error, _} = R} <- Vals] of
                                [] -> {ok, maps:from_list([{S, re(V)} || {S, {ok, V}} <- Vals])};
                                [E | _] -> E
                            end;
                        E -> E
                    end
            end
    end.

numbers(#circuit{params = P}) -> maps:filter(fun(_, V) -> is_number(V) end, P).

re(X) when is_float(X) -> X;
re({complex, R, _}) -> R.

numeric(C, Obs, Binding, Syms, Opts) ->
    H = maps:get(step, Opts, 1.0e-6),
    try
        {ok, maps:from_list([{S, (num_value(C, Obs, shift(Binding, S, H), Opts) -
                                  num_value(C, Obs, shift(Binding, S, -H), Opts)) / (2 * H)} || S <- Syms])}
    catch
        throw:{error, _} = E -> E
    end.

shift(B, S, H) -> B#{S => float(value_of(maps:get(S, B))) + H}.

value_of({rat, N, D}) -> N / D;
value_of(X) -> X.

num_value(C, Obs, B, Opts) ->
    case value(C, Obs, Opts#{bindings => B}) of
        {ok, V} -> V;
        {error, _} = E -> throw(E)
    end.

%% the observable's value on the circuit's output state (FS 4.8) under Opts bindings
-spec value(#circuit{}, observable(), map()) -> {ok, float()} | {error, pq_core:reason()}.
value(C, {weighted_sum, Terms}, Opts) ->
    B = maps:get(bindings, Opts, #{}),
    Vals = [{pq_symbolic:evaluate(pq_symbolic:substitute(W, maps:merge(numbers(C), B)), #{}), value(C, O, Opts)} || {W, O} <- Terms],
    case [E || {{error, _} = E, _} <- Vals] ++ [E || {_, {error, _} = E} <- Vals] of
        [] -> {ok, lists:sum([re(W) * V || {{ok, W}, {ok, V}} <- Vals])};
        [E | _] -> E
    end;
value(#circuit{instrs = Is} = C, {probability, Modes, Pattern}, Opts) ->
    Body = [I || I <- Is, not measurement(I)],
    case pq_api:probabilities(C#circuit{instrs = Body ++ [{measure_pnr, Modes}]}, Opts) of
        {ok, Map, _} -> {ok, maps:get(Pattern, Map, 0.0)};
        E -> E
    end;
value(#circuit{instrs = Is} = C, {mean_photon_number, Modes}, Opts) ->
    Body = C#circuit{instrs = [I || I <- Is, not measurement(I)]},
    case pq_api:evolve(Body, maps:remove(backend, Opts)) of
        {ok, S, _} -> pq_stats:mean_photon_number(S, Modes);
        E -> E
    end.

measurement(I) -> lists:member(element(1, I), [measure_pnr, measure_threshold, measure_homodyne]).

%% Decision 6: the derivative expressions, for passive circuits with a Fock input, or Gaussian circuits with a
%% Gaussian input (FR-G04)
-spec symbolic_gradient(#circuit{}, observable(), [atom()]) -> {ok, #{atom() => pq_symbolic:sym_expr()}} | {error, pq_core:reason()}.
symbolic_gradient(#circuit{input = {fock, In}, instrs = Is} = C, Obs, Syms) ->
    Body = [I || I <- Is, not measurement(I)],
    case lists:all(fun(I) -> lists:member(element(1, I), [beamsplitter, phaseshift, unitary]) end, Body) of
        false -> {error, {differentiate_unsupported, {instructions, [element(1, I) || I <- Body]}}};
        true ->
            Sub = maps:filter(fun(K, _) -> not lists:member(K, Syms) end, C#circuit.params),
            Bc = C#circuit{instrs = [map_instr(I, Sub) || I <- Body]},
            case pq_symbolic:circuit_unitary(Bc) of
                {ok, U} ->
                    case observable_expr(U, In, Obs) of
                        {ok, E} ->
                            D = [{S, pq_symbolic:differentiate(E, S, #{domains => C#circuit.symbols})} || S <- Syms],
                            case [R || {_, {error, _} = R} <- D] of
                                [] -> {ok, maps:from_list([{S, X} || {S, {ok, X}} <- D])};
                                [R | _] -> R
                            end;
                        Err -> Err
                    end;
                Err -> Err
            end
    end;
symbolic_gradient(#circuit{input = In, instrs = Is} = C, Obs, Syms) when In =:= vacuum; element(1, In) =:= gaussian ->
    %% FR-G04, SY 6.4: a Gaussian input, the observable from the symbolic mean and covariance
    Body = [I || I <- Is, not measurement(I)],
    Gaussian = [beamsplitter, phaseshift, unitary, displacement, squeezing, loss],
    case lists:all(fun(I) -> lists:member(element(1, I), Gaussian) end, Body) of
        false -> {error, {differentiate_unsupported, {instructions, [element(1, I) || I <- Body]}}};
        true ->
            Sub = maps:filter(fun(K, _) -> not lists:member(K, Syms) end, C#circuit.params),
            case pq_gaussian:symbolic_evolve(C#circuit{instrs = [map_instr(I, Sub) || I <- Body]}, #{}) of
                {ok, S} ->
                    case gaussian_expr(S, Obs) of
                        {ok, E} -> differentiate_all(E, Syms, C);
                        Err -> Err
                    end;
                Err -> Err
            end
    end;
symbolic_gradient(_C, Obs, _Syms) -> {error, {differentiate_unsupported, Obs}}.

differentiate_all(E, Syms, C) ->
    D = [{S, pq_symbolic:differentiate(E, S, #{domains => C#circuit.symbols})} || S <- Syms],
    case [R || {_, {error, _} = R} <- D] of
        [] -> {ok, maps:from_list([{S, X} || {S, {ok, X}} <- D])};
        [R | _] -> R
    end.

%% FR-X01 on the symbolic state: per mode (mu_x^2 + mu_p^2 + Sigma_xx + Sigma_pp) / (2 hbar) - 1/2
gaussian_expr(#gaussian_state{n_modes = N, hbar = H, mean = {matrix, M}, cov = {matrix, Cv}}, {mean_photon_number, Modes}) ->
    Mu = fun(K) -> hd(lists:nth(K + 1, M)) end,
    Sg = fun(K) -> lists:nth(K + 1, lists:nth(K + 1, Cv)) end,
    pq_symbolic:canonicalise({add, [{add, [{mul, [{pow, {mul, [2, H]}, -1},
                                                  {add, [{pow, Mu(J), 2}, {pow, Mu(N + J), 2}, Sg(J), Sg(N + J)]}]},
                                           {rat, -1, 2}]} || J <- Modes]});
gaussian_expr(S, {weighted_sum, Terms}) ->
    Es = [{W, gaussian_expr(S, O)} || {W, O} <- Terms],
    case [E || {_, {error, _} = E} <- Es] of
        [] -> pq_symbolic:canonicalise({add, [{mul, [W, E]} || {W, {ok, E}} <- Es]});
        [E | _] -> E
    end;
gaussian_expr(_S, Obs) -> {error, {differentiate_unsupported, Obs}}.

map_instr(I, Sub) when map_size(Sub) =:= 0 -> I;
map_instr({beamsplitter, A, B, T, P}, Sub) -> {beamsplitter, A, B, pq_symbolic:substitute(T, Sub), pq_symbolic:substitute(P, Sub)};
map_instr({phaseshift, M, P}, Sub) -> {phaseshift, M, pq_symbolic:substitute(P, Sub)};
map_instr({displacement, M, R, P}, Sub) -> {displacement, M, pq_symbolic:substitute(R, Sub), pq_symbolic:substitute(P, Sub)};
map_instr({squeezing, M, R, P}, Sub) -> {squeezing, M, pq_symbolic:substitute(R, Sub), pq_symbolic:substitute(P, Sub)};
map_instr({loss, M, Eta}, Sub) -> {loss, M, pq_symbolic:substitute(Eta, Sub)};
map_instr(I, _Sub) -> I.

observable_expr(U, In, {probability, Modes, Pattern}) ->
    {matrix, Rows} = U,
    N = length(Rows),
    Total = lists:sum(In),
    Rest = Total - lists:sum(Pattern),
    Others = [J || J <- lists:seq(0, N - 1), not lists:member(J, Modes)],
    case Rest >= 0 of
        false -> {ok, 0};
        true ->
            Completions = [full(N, Modes, Pattern, Others, C) || C <- compositions(length(Others), Rest)],
            Ps = [pq_symbolic:bs_probability(U, In, T) || T <- Completions],
            case [E || {error, _} = E <- Ps] of
                [] -> pq_symbolic:canonicalise({add, [P || {ok, P} <- Ps]});
                [E | _] -> E
            end
    end;
observable_expr({matrix, Rows}, In, {mean_photon_number, Modes}) ->
    %% sum over j in Modes, i of S_i |U_ji|^2 (exact for linear optics)
    pq_symbolic:canonicalise({add, [{mul, [Si, lists:nth(I, lists:nth(J + 1, Rows)), {conj, lists:nth(I, lists:nth(J + 1, Rows))}]}
                                    || J <- Modes, {I, Si} <- lists:enumerate(1, In), Si > 0]});
observable_expr(U, In, {weighted_sum, Terms}) ->
    Es = [{W, observable_expr(U, In, O)} || {W, O} <- Terms],
    case [E || {_, {error, _} = E} <- Es] of
        [] -> pq_symbolic:canonicalise({add, [{mul, [W, E]} || {W, {ok, E}} <- Es]});
        [E | _] -> E
    end.

full(N, Modes, Pattern, Others, C) ->
    M = maps:from_list(lists:zip(Modes, Pattern) ++ lists:zip(Others, C)),
    [maps:get(J, M) || J <- lists:seq(0, N - 1)].

compositions(0, 0) -> [[]];
compositions(0, _) -> [];
compositions(K, N) -> [[A | R] || A <- lists:seq(0, N), R <- compositions(K - 1, N - A)].
