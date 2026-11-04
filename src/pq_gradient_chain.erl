%% The chain rule shared by the gate-level gradient strategies (pq_gradient_adjoint, pq_gradient_shift): a gradient
%% with respect to the circuit's symbols from derivatives with respect to the gates' numeric parameters.
%%
%% Every parameter expression is expanded the way lowering resolves it (FR-C04): a symbol outside the differentiated
%% ones takes its value from the binding, else its default (which may itself use the differentiated symbols), until
%% only those remain; then d(expression)/d(symbol) is evaluated at the effective binding. Weighted sums are
%% flattened into terms with their weights expanded the same way, so a weight that depends on a symbol contributes
%% dW/ds times its observable's value. Internal: not part of the frozen API.
-module(pq_gradient_chain).
-include("pq.hrl").
-export([prepare/4, gate_params/1, combine/4]).

-type basic() :: {probability, [non_neg_integer()], [non_neg_integer()]} | {mean_photon_number, [non_neg_integer()]}.
-type prepared() :: #{body := [tuple()], coefs := [[#{atom() => float()}]], terms := [{float(), #{atom() => float()}, basic()}],
                      lowered := #circuit{}, symbols := [atom()]}.
-export_type([prepared/0, basic/0]).

%% the circuit without its measurements, lowered at Binding; per instruction and parameter the derivatives of the
%% parameter with respect to Syms; the observable as terms {Weight, dWeight/dSym, Basic}
-spec prepare(#circuit{}, term(), map(), [atom()]) -> {ok, prepared()} | {error, pq_core:reason()}.
prepare(#circuit{instrs = Is, input = In} = C, Obs, Binding, Syms) ->
    try
        Body = [I || I <- Is, not measurement(I)],
        [fail({differentiate_unsupported, {instructions, [K]}}) || I <- Body, K <- [element(1, I)],
                                                                    K =:= custom orelse K =:= feed_forward orelse K =:= conditional],
        InputSyms = pq_circuit:all_symbols(C#circuit{instrs = [], params = #{}}),
        [fail({differentiate_unsupported, {input, S}}) || S <- Syms, lists:member(S, InputSyms), symbolic_input(In)],
        Bc = C#circuit{instrs = Body},
        Env = ok(pq_circuit:effective_binding(Bc, Binding)),
        Point = maps:with(Syms, Env),
        Subst = maps:without(Syms, maps:merge(C#circuit.params, Binding)),
        [fail({differentiate_unsupported, {symbolic_unitary, K}})
         || {K, {unitary, _, U}} <- lists:enumerate(0, Body), element(1, U) =/= cmat,
            lists:any(fun(S) -> lists:member(S, Syms) end, free(expand(U, Subst, Syms)))],
        Low = ok(pq_circuit:lower(Bc, Binding)),
        Coefs = [[derivatives(expand(E, Subst, Syms), Syms, Point, C) || E <- param_exprs(I)] || I <- Body],
        Terms = [{value(W, Point, C), derivatives(W, Syms, Point, C), B}
                 || {W0, B} <- flatten(Obs, 1), W <- [expand(W0, Subst, Syms)]],
        {ok, #{body => Body, coefs => Coefs, terms => Terms, lowered => Low, symbols => Syms}}
    catch
        throw:{error, _} = E -> E
    end.

symbolic_input({gaussian, #gaussian_state{mean = {rvec, _, _}, cov = {rmat, _, _, _}}}) -> false;
symbolic_input({fock_state, #fock_state{amps = {cvec, _, _}}}) -> false;
symbolic_input(_) -> true.

%% the parameter expressions of an instruction, in the order of its fields; a symbolic unitary has none that the
%% gate-level strategies can differentiate
param_exprs({beamsplitter, _, _, T, P}) -> [T, P];
param_exprs({phaseshift, _, P}) -> [P];
param_exprs({kerr, _, K}) -> [K];
param_exprs({displacement, _, R, P}) -> [R, P];
param_exprs({squeezing, _, R, P}) -> [R, P];
param_exprs({loss, _, Eta}) -> [Eta];
param_exprs(_) -> [].

%% the numeric parameters of a lowered instruction, in the same order
-spec gate_params(tuple()) -> [float()].
gate_params({beamsplitter, _, _, T, P}) -> [T, P];
gate_params({phaseshift, _, P}) -> [P];
gate_params({kerr, _, K}) -> [K];
gate_params({displacement, _, R, P}) -> [R, P];
gate_params({squeezing, _, R, P}) -> [R, P];
gate_params({loss, _, Eta}) -> [Eta];
gate_params(_) -> [].

%% the gradient: sum over gates and parameters of dF/dp * dp/ds, plus sum over terms of dW/ds * value
-spec combine(prepared(), [[float()]], [float()], [atom()]) -> #{atom() => float()}.
combine(#{coefs := Coefs, terms := Terms}, GateGrads, TermValues, Syms) ->
    Gates = lists:zip(Coefs, GateGrads),
    maps:from_list([{S, lists:sum([maps:get(S, D, 0.0) * G || {Ds, Gs} <- Gates, {D, G} <- lists:zip(Ds, pad(Gs, Ds))])
                        + lists:sum([maps:get(S, DW, 0.0) * V || {{_, DW, _}, V} <- lists:zip(Terms, TermValues)])}
                    || S <- Syms]).

%% passive gates report no parameter gradients
pad(Gs, Ds) when length(Gs) =:= length(Ds) -> Gs;
pad([], Ds) -> [0.0 || _ <- Ds].

flatten({weighted_sum, Terms}, W) -> lists:append([flatten(O, {mul, [W, Wi]}) || {Wi, O} <- Terms]);
flatten({probability, _, _} = O, W) -> [{W, O}];
flatten({mean_photon_number, _} = O, W) -> [{W, O}];
flatten(O, _) -> fail({differentiate_unsupported, O}).

%% substitute until only the differentiated symbols (and unbound ones) are left; a cycle of defaults is refused
expand(E, Subst, Syms) -> expand(E, Subst, Syms, map_size(Subst) + 1).

expand(E, Subst, Syms, Left) ->
    case [X || X <- free(E), maps:is_key(X, Subst), not lists:member(X, Syms)] of
        [] -> E;
        _ when Left =:= 0 -> fail({differentiate_unsupported, {cyclic_defaults, free(E)}});
        _ -> expand(pq_symbolic:substitute(E, Subst), Subst, Syms, Left - 1)
    end.

free(E) when is_number(E) -> [];
free(E) -> pq_symbolic:free_symbols(E).

derivatives(E, Syms, Point, C) ->
    Fs = free(E),
    maps:from_list([{S, value(ok(pq_symbolic:differentiate(E, S, #{domains => C#circuit.symbols})), Point, C)}
                    || S <- Syms, lists:member(S, Fs)]).

value(E, _Point, _C) when is_number(E) -> float(E);
value(E, Point, C) ->
    case pq_symbolic:evaluate(E, Point, #{domains => C#circuit.symbols}) of
        {ok, X} when is_float(X) -> X;
        {ok, {complex, R, _}} -> R;
        {ok, X} when is_number(X) -> float(X);
        {ok, Other} -> fail({differentiate_unsupported, {not_scalar, Other}});
        {error, _} = Err -> throw(Err)
    end.

measurement(I) -> lists:member(element(1, I), [measure_pnr, measure_threshold, measure_homodyne]).

ok({ok, V}) -> V;
ok({error, _} = E) -> throw(E).

-spec fail(term()) -> no_return().
fail(R) -> throw({error, R}).
