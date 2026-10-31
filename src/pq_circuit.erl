%% Circuits (FS 4.3, 4.4, 5.1): construction, editing, validation, parameters, composition, symbols,
%% serialisation and the gate library. Pure and functional: every operation returns a new circuit.
%% Lowering (FR-C07), the capability check of FS 5.3 and the mesh (FR-C10) arrive with the backends (Phase 1);
%% JSON (FR-C09) is Phase 3.
-module(pq_circuit).
-include("pq.hrl").
-export([new/1, new/2, add/2, insert/3, replace/3, delete/2,
         validate/1, effective_binding/2, compose/2, compose/3,
         free_symbols/1, all_symbols/1, to_binary/1, from_binary/1, lower/2, check_capabilities/2, check_capabilities/3, mesh/1, to_json/1, from_json/1, from_json/2,
         two_mode_squeezing/3, two_mode_squeezing/4, mzi/4]).

-export_type([circuit/0, instruction/0, input_spec/0]).

-type circuit() :: #circuit{}.
-type mode() :: non_neg_integer().
-type param() :: pq_symbolic:sym_expr().
-type instruction() ::
        {beamsplitter, mode(), mode(), param(), param()}
      | {phaseshift, mode(), param()}
      | {displacement, mode(), param(), param()}
      | {squeezing, mode(), param(), param()}
      | {kerr, mode(), param()}
      | {loss, mode(), param()}
      | {unitary, [mode()], pq_symbolic:sym_expr() | pq_core:numeric_matrix()}
      | {measure_pnr, [mode()]}
      | {measure_threshold, [mode()]}
      | {measure_homodyne, mode(), param()}
      | {custom, atom(), map()}.
-type input_spec() :: vacuum | {fock, [non_neg_integer()]} | {gaussian, #gaussian_state{}}
                    | {fock_state, #fock_state{}}.
-type result(T) :: {ok, T} | {error, pq_core:reason()}.

-define(BINARY_VERSION, 1).

%%% FR-C01, FR-C02

-spec new(pos_integer()) -> result(circuit()).
new(N) -> new(N, vacuum).

-spec new(pos_integer(), input_spec()) -> result(circuit()).
new(N, Input) when is_integer(N), N >= 1 ->
    C = #circuit{n_modes = N, input = Input},
    try check_input(C) of
        ok -> {ok, C}
    catch
        throw:{error, _} = E -> E
    end;
new(N, _) -> {error, {bad_shape, #{n_modes => N}}}.

-spec add(circuit(), instruction()) -> circuit().
add(#circuit{instrs = Is} = C, I) -> C#circuit{instrs = Is ++ [I]}.

-spec insert(circuit(), non_neg_integer(), instruction()) -> result(circuit()).
insert(#circuit{instrs = Is} = C, Pos, I) when is_integer(Pos), Pos >= 0, Pos =< length(Is) ->
    {A, B} = lists:split(Pos, Is),
    {ok, C#circuit{instrs = A ++ [I | B]}};
insert(_C, Pos, _I) -> {error, {index_out_of_range, Pos}}.

-spec replace(circuit(), non_neg_integer(), instruction()) -> result(circuit()).
replace(#circuit{instrs = Is} = C, Pos, I) when is_integer(Pos), Pos >= 0, Pos < length(Is) ->
    {A, [_ | B]} = lists:split(Pos, Is),
    {ok, C#circuit{instrs = A ++ [I | B]}};
replace(_C, Pos, _I) -> {error, {index_out_of_range, Pos}}.

-spec delete(circuit(), non_neg_integer()) -> result(circuit()).
delete(#circuit{instrs = Is} = C, Pos) when is_integer(Pos), Pos >= 0, Pos < length(Is) ->
    {A, [_ | B]} = lists:split(Pos, Is),
    {ok, C#circuit{instrs = A ++ B}};
delete(_C, Pos) -> {error, {index_out_of_range, Pos}}.

%%% FR-C03: validation (without the backend capability check)

-spec validate(circuit()) -> ok | {error, pq_core:reason()}.
validate(#circuit{n_modes = N} = C) when is_integer(N), N >= 1 ->
    try
        ok = check_hbar(C),
        ok = check_maps(C),
        ok = check_input(C),
        _ = lists:foldl(fun(I, {K, Measured}) -> {K + 1, check_instr(I, K, C, Measured)} end,
                        {0, #{}}, C#circuit.instrs),
        ok
    catch
        throw:{error, _} = E -> E
    end;
validate(#circuit{n_modes = N}) -> {error, {bad_shape, #{n_modes => N}}};
validate(Other) -> {error, {bad_term, {not_a_circuit, Other}}}.

-spec fail(pq_core:reason()) -> no_return().
fail(R) -> throw({error, R}).

check_hbar(#circuit{hbar = undefined}) -> ok;
check_hbar(#circuit{hbar = H}) when is_number(H), H > 0 -> ok;
check_hbar(#circuit{hbar = H}) -> fail({bad_shape, #{hbar => H}}).

check_maps(#circuit{params = P, symbols = S}) when is_map(P), is_map(S) ->
    maps:foreach(fun(K, V) ->
                         is_atom(K) orelse fail({bad_term, {param_name, K}}),
                         scalar(V, {param, K})
                 end, P),
    maps:foreach(fun(K, D) ->
                         (is_atom(K) andalso lists:member(D, [real, positive, nonnegative, complex]))
                             orelse fail({bad_term, {symbol_domain, K, D}})
                 end, S),
    ok;
check_maps(_) -> fail({bad_term, params_or_symbols}).

check_input(#circuit{input = vacuum}) -> ok;
check_input(#circuit{n_modes = N, input = {fock, Occ}}) ->
    (is_list(Occ) andalso length(Occ) =:= N andalso lists:all(fun(X) -> is_integer(X) andalso X >= 0 end, Occ))
        orelse fail({unsupported_input, {fock, Occ}}),
    ok;
check_input(#circuit{n_modes = N, hbar = H, input = {gaussian, #gaussian_state{} = G}}) ->
    G#gaussian_state.n_modes =:= N orelse fail({bad_shape, #{input_modes => G#gaussian_state.n_modes, n_modes => N}}),
    (H =:= undefined orelse H == G#gaussian_state.hbar) orelse fail({hbar_mismatch, H, G#gaussian_state.hbar}),
    vector_or_sym(G#gaussian_state.mean, rvec, 2 * N, mean),
    matrix_or_sym(G#gaussian_state.cov, rmat, 2 * N, cov),
    ok;
check_input(#circuit{n_modes = N, input = {fock_state, #fock_state{n_modes = M, cutoff = Cut, amps = A}}}) ->
    M =:= N orelse fail({bad_shape, #{input_modes => M, n_modes => N}}),
    (is_integer(Cut) andalso Cut >= 0) orelse fail({bad_shape, #{cutoff => Cut}}),
    vector_or_sym(A, cvec, binomial(Cut + N, N), amps),
    ok;
check_input(#circuit{input = I}) -> fail({unsupported_input, I}).

vector_or_sym({Tag, L, _} = V, Tag, L, _What) ->
    pq_core:validate(V) =:= ok orelse fail({bad_shape, V}), ok;
vector_or_sym({Tag, L, _}, Tag, Want, What) -> fail({bad_shape, #{What => L, expected => Want}});
vector_or_sym({matrix, _} = E, _Tag, Len, What) ->
    pq_symbolic:shape(E) =:= {ok, {Len, 1}} orelse fail({bad_shape, #{What => E, expected => {Len, 1}}}), ok;
vector_or_sym(E, _, Len, What) -> fail({bad_shape, #{What => E, expected => Len}}).

matrix_or_sym({Tag, R, R, _} = M, Tag, R, _What) ->
    pq_core:validate(M) =:= ok orelse fail({bad_shape, M}), ok;
matrix_or_sym(E, _Tag, R, What) ->
    pq_symbolic:shape(E) =:= {ok, {R, R}} orelse fail({bad_shape, #{What => E, expected => {R, R}}}), ok.

binomial(N, K) -> binomial(N, K, 1, 1).
binomial(_N, 0, Num, Den) -> Num div Den;
binomial(N, K, Num, Den) -> binomial(N - 1, K - 1, Num * N, Den * K).

%% one instruction at position K; Measured maps each measured mode to the index of its measurement
check_instr(I, K, C, Measured) ->
    Modes = modes(I, K, C),
    lists:foreach(fun(M) -> (is_integer(M) andalso M >= 0 andalso M < C#circuit.n_modes)
                                orelse fail({mode_out_of_range, #{mode => M, instr_index => K}})
                  end, Modes),
    length(lists:usort(Modes)) =:= length(Modes) orelse fail({duplicate_modes, #{modes => Modes, instr_index => K}}),
    [fail({measured_mode_reuse, M, K}) || M <- Modes, maps:is_key(M, Measured)],
    check_params(I, K, C),
    case is_measurement(I) of
        true -> maps:merge(Measured, maps:from_keys(Modes, K));
        false -> Measured
    end.

modes({beamsplitter, A, B, _, _}, _, _) -> [A, B];
modes({phaseshift, M, _}, _, _) -> [M];
modes({displacement, M, _, _}, _, _) -> [M];
modes({squeezing, M, _, _}, _, _) -> [M];
modes({kerr, M, _}, _, _) -> [M];
modes({loss, M, _}, _, _) -> [M];
modes({unitary, Ms, _}, _K, _) when is_list(Ms), Ms =/= [] -> Ms;
modes({measure_pnr, Ms}, _, _) when is_list(Ms), Ms =/= [] -> Ms;
modes({measure_threshold, Ms}, _, _) when is_list(Ms), Ms =/= [] -> Ms;
modes({measure_homodyne, M, _}, _, _) -> [M];
modes({custom, Name, Args}, K, C) when is_atom(Name), is_map(Args) ->
    Mod = custom_module(Name),
    case Mod:validate(Args, C#circuit.n_modes) of
        ok -> maps:get(modes, Args, []);
        {error, R} -> fail({unsupported_instruction, Name, #{instr_index => K, reason => R}})
    end;
modes(I, K, _) -> fail({unsupported_instruction, I, #{instr_index => K}}).

is_measurement(I) -> lists:member(element(1, I), [measure_pnr, measure_threshold, measure_homodyne]).

check_params({beamsplitter, _, _, T, P}, K, _) -> scalar(T, K), scalar(P, K);
check_params({phaseshift, _, P}, K, _) -> scalar(P, K);
check_params({displacement, _, R, P}, K, _) -> scalar(R, K), scalar(P, K);
check_params({squeezing, _, R, P}, K, _) -> scalar(R, K), scalar(P, K);
check_params({kerr, _, Kappa}, K, _) -> scalar(Kappa, K);
check_params({measure_homodyne, _, P}, K, _) -> scalar(P, K);
check_params({loss, _, Eta}, K, C) ->
    scalar(Eta, K),
    case numeric_value(Eta, C) of
        {ok, V} when is_float(V), V >= 0.0, V =< 1.0 -> ok;
        {ok, V} -> fail({domain_violation, #{param => eta, value => V, instr_index => K}});
        symbolic -> ok                                   % checked after lowering
    end;
check_params({unitary, Ms, U}, K, C) -> check_unitary(U, length(Ms), K, C);
check_params(_, _, _) -> ok.

scalar(E, Where) ->
    case pq_symbolic:shape(E) of
        {ok, scalar} -> ok;
        _ -> fail({bad_shape, #{param => E, at => Where}})
    end.

%% the value of a parameter without free symbols (defaults applied), else symbolic
numeric_value(E, #circuit{params = P, symbols = S}) ->
    Sub = pq_symbolic:substitute(E, P),
    case pq_symbolic:free_symbols(Sub) of
        [] ->
            case pq_symbolic:evaluate(Sub, #{}, #{domains => S}) of
                {ok, V} -> {ok, V};
                {error, R} -> fail(R)
            end;
        _ -> symbolic
    end.

%% a numeric k x k matrix must be unitary within unitarity_tol; a symbolic one only needs the shape
check_unitary({cmat, KK, KK, _} = U, KK, Idx, _C) ->
    pq_core:validate(U) =:= ok orelse fail({bad_shape, #{unitary => U, instr_index => Idx}}),
    {ok, {cmat, Rows}} = case pq_core:unpack(U) of
                             {ok, _} = Ok -> Ok;
                             {error, R} -> fail(R)
                         end,
    Dev = unitarity_deviation(Rows),
    Tol = pq_core:get(unitarity_tol),
    Dev =< Tol orelse fail({not_unitary, #{deviation => Dev, instr_index => Idx}}),
    ok;
check_unitary({cmat, R, Cc, _}, KK, Idx, _C) ->
    fail({bad_shape, #{unitary => {R, Cc}, expected => {KK, KK}, instr_index => Idx}});
check_unitary(E, KK, Idx, _C) ->
    pq_symbolic:shape(E) =:= {ok, {KK, KK}} orelse fail({bad_shape, #{unitary => E, expected => {KK, KK},
                                                                         instr_index => Idx}}),
    ok.

%% the largest entry of U^dagger U - I
unitarity_deviation(Rows) ->
    Cols = transpose(Rows),
    K = length(Rows),
    lists:max([abs_c(sub_c(dot(Ci, Cj), case I =:= J of true -> {1.0, 0.0}; false -> {0.0, 0.0} end))
               || {I, Ci} <- lists:zip(lists:seq(1, K), Cols), {J, Cj} <- lists:zip(lists:seq(1, K), Cols)]).

transpose([[] | _]) -> [];
transpose(Rows) -> [[hd(R) || R <- Rows] | transpose([tl(R) || R <- Rows])].

%% <a, b> = sum conj(a_k) b_k
dot(A, B) ->
    lists:foldl(fun({{complex, Ar, Ai}, {complex, Br, Bi}}, {Sr, Si}) ->
                        {Sr + Ar * Br + Ai * Bi, Si + Ar * Bi - Ai * Br}
                end, {0.0, 0.0}, lists:zip(A, B)).

sub_c({A, B}, {C, D}) -> {A - C, B - D}.
abs_c({A, B}) -> math:sqrt(A * A + B * B).

custom_module(Name) ->
    Mods = [M || M <- pq_core:get(custom_instructions), code:ensure_loaded(M) =:= {module, M},
                 erlang:function_exported(M, name, 0), M:name() =:= Name],
    case Mods of
        [M] -> M;
        [] -> fail({unsupported_instruction, Name, custom});
        _ -> fail({unsupported_instruction, Name, duplicate_custom})
    end.

%%% FR-C04: the effective binding

%% The defaults overridden by Binding, resolved in dependency order; every free symbol of the instructions and
%% of the input must have a value.
-spec effective_binding(circuit(), pq_symbolic:binding()) -> result(#{atom() => term()}).
effective_binding(#circuit{params = P, symbols = S} = C, Binding) ->
    case pq_symbolic:resolve_params(P, Binding, #{domains => S}) of
        {ok, Values} ->
            case [X || X <- instr_symbols(C), not maps:is_key(X, Values)] of
                [] -> {ok, Values};
                Missing -> {error, {unbound_symbols, lists:usort(Missing)}}
            end;
        {error, _} = E -> E
    end.

%%% FR-C07: lowering

%% The circuit under the effective binding with every parameter a float, every unitary a packed matrix, custom
%% instructions expanded, meta.lowered = true and no symbols. Validated before and after.
-spec lower(circuit(), pq_symbolic:binding()) -> result(circuit()).
lower(#circuit{} = C, Binding) ->
    try
        ok = ok_or_fail(validate(C)),
        {ok, Env} = ok_or_fail(effective_binding(C, Binding)),
        Sub = maps:map(fun(_, V) -> number_ast(V) end, Env),
        Instrs = lists:append([lower_instr(I, K, Sub, Env, C) || {K, I} <- lists:enumerate(0, C#circuit.instrs)]),
        Unused = [K || K <- maps:keys(Binding), not lists:member(K, all_symbols(C))],
        Meta = (C#circuit.meta)#{lowered => true, binding => Env,
                                 warnings => [{unused_bindings, lists:sort(Unused)} || Unused =/= []]},
        Low = C#circuit{instrs = Instrs, params = #{}, meta = Meta},
        ok = ok_or_fail(validate(Low)),
        {ok, Low}
    catch
        throw:{error, _} = E -> E
    end.

ok_or_fail(ok) -> ok;
ok_or_fail({ok, _} = R) -> R;
ok_or_fail({error, R}) -> fail(R).

number_ast(X) when is_float(X) -> X;
number_ast({complex, R, I}) -> {complex, R, I}.

lower_instr({custom, Name, Args}, K, _Sub, Env, C) ->
    Mod = custom_module(Name),
    case Mod:expand(Args, Env) of
        {ok, Instrs} ->
            [fail({unsupported_instruction, custom, #{instr_index => K, nested => true}}) ||
                I <- Instrs, element(1, I) =:= custom],
            lists:append([lower_instr(I, K, maps:map(fun(_, V) -> number_ast(V) end, Env), Env, C) || I <- Instrs]);
        {error, R} -> fail({unsupported_instruction, Name, #{instr_index => K, reason => R}})
    end;
lower_instr({unitary, Ms, {cmat, _, _, _} = U}, _K, _Sub, _Env, _C) -> [{unitary, Ms, U}];
lower_instr({unitary, Ms, U}, K, Sub, _Env, C) ->
    case pq_symbolic:evaluate(pq_symbolic:substitute(U, Sub), #{}, #{domains => C#circuit.symbols}) of
        {ok, {cmat, _, _, _} = P} -> [{unitary, Ms, P}];
        {ok, _} -> fail({bad_shape, #{unitary => U, instr_index => K}});
        {error, R} -> fail(R)
    end;
lower_instr(I, K, Sub, _Env, C) ->
    [map_params(I, fun(E) -> real_value(E, K, Sub, C) end)].

real_value(E, K, Sub, C) ->
    case pq_symbolic:evaluate(pq_symbolic:substitute(E, Sub), #{}, #{domains => C#circuit.symbols}) of
        {ok, X} when is_float(X) -> X;
        {ok, {complex, Re, Im}} when abs(Im) =< 1.0e-12 * (1 + abs(Re)) -> Re;
        {ok, V} -> fail({unsupported_input, #{complex_parameter => V, instr_index => K}});
        {error, R} -> fail(R)
    end.

%% FS 5.3 and FR-M06: every instruction kind supported by the backend; a measurement the backend allows only at
%% the end ('terminal') must be followed by measurements only.
-spec check_capabilities(circuit(), module()) -> ok | {error, pq_core:reason()}.
check_capabilities(C, Backend) -> check_capabilities(C, Backend, #{}).

%% with the run's options: a backend may export capabilities/1 (pq_fock: loss is exact under promote_dm)
-spec check_capabilities(circuit(), module(), map()) -> ok | {error, pq_core:reason()}.
check_capabilities(#circuit{instrs = Is}, Backend, Opts) ->
    _ = code:ensure_loaded(Backend),
    Caps = case erlang:function_exported(Backend, capabilities, 1) of
               true -> Backend:capabilities(Opts);
               false -> Backend:capabilities()
           end,
    Name = Backend:name(),
    try
        lists:foreach(
          fun({K, I}) ->
                  Kind = element(1, I),
                  case maps:get(Kind, Caps, unsupported) of
                      unsupported when Kind =/= custom -> fail({unsupported_instruction, Kind, Name});
                      terminal ->
                          Later = lists:nthtail(K + 1, Is),
                          [fail({unsupported_instruction, Kind, #{backend => Name, placement => midcircuit,
                                                                  instr_index => K}})
                           || L <- Later, not is_measurement(L)];
                      _ -> ok
                  end
          end, lists:enumerate(0, Is)),
        ok
    catch
        throw:{error, _} = E -> E
    end.

%%% FR-C05: composition

-spec compose(circuit(), circuit()) -> result(circuit()).
compose(C1, C2) -> compose(C1, C2, #{}).

%% Opts: modes => #{ModeOfC2 => ModeOfC1} to embed a smaller circuit; rename => #{Old => New} for C2's symbols.
-spec compose(circuit(), circuit(), map()) -> result(circuit()).
compose(#circuit{n_modes = N1} = C1, #circuit{n_modes = N2} = C2, Opts) ->
    try
        Rename = maps:get(rename, Opts, #{}),
        C2r = rename(C2, Rename),
        Instrs2 = case maps:find(modes, Opts) of
                      {ok, Map} -> embed(C2r#circuit.instrs, Map, N1);
                      error when N1 =:= N2 -> C2r#circuit.instrs;
                      error -> fail({bad_shape, #{n_modes => {N1, N2}}})
                  end,
        Hbar = case {C1#circuit.hbar, C2#circuit.hbar} of
                   {undefined, H} -> H;
                   {H, undefined} -> H;
                   {H1, H2} when H1 == H2 -> H1;
                   {H1, H2} -> fail({hbar_mismatch, H1, H2})
               end,
        Params = merge_checked(C1#circuit.params, C2r#circuit.params, fun same_expr/2),
        Symbols = merge_checked(C1#circuit.symbols, C2r#circuit.symbols, fun(A, B) -> A =:= B end),
        {ok, C1#circuit{instrs = C1#circuit.instrs ++ Instrs2, hbar = Hbar, params = Params, symbols = Symbols,
                        meta = maps:merge(C2#circuit.meta, C1#circuit.meta)}}
    catch
        throw:{error, _} = E -> E
    end.

merge_checked(A, B, Same) ->
    maps:fold(fun(K, V, Acc) ->
                      case Acc of
                          #{K := V0} -> Same(V0, V) orelse fail({param_conflict, K}), Acc;
                          _ -> Acc#{K => V}
                      end
              end, A, B).

same_expr(A, B) -> A =:= B orelse pq_symbolic:canonicalise(A) =:= pq_symbolic:canonicalise(B).

rename(C, R) when map_size(R) =:= 0 -> C;
rename(#circuit{instrs = Is, params = P, symbols = S, input = In} = C, R) ->
    Sub = maps:map(fun(_, New) -> {sym, New} end, R),
    Key = fun(K) -> maps:get(K, R, K) end,
    C#circuit{instrs = [map_params(I, fun(E) -> pq_symbolic:substitute(E, Sub) end) || I <- Is],
              params = maps:from_list([{Key(K), pq_symbolic:substitute(V, Sub)} || K := V <- P]),
              symbols = maps:from_list([{Key(K), D} || K := D <- S]),
              input = In}.

embed(Instrs, Map, N1) ->
    M = fun(X) ->
                case Map of
                    #{X := Y} when is_integer(Y), Y >= 0, Y < N1 -> Y;
                    _ -> fail({mode_out_of_range, #{mode => X, mode_map => Map}})
                end
        end,
    [map_modes(I, M) || I <- Instrs].

map_modes({beamsplitter, A, B, T, P}, M) -> {beamsplitter, M(A), M(B), T, P};
map_modes({unitary, Ms, U}, M) -> {unitary, [M(X) || X <- Ms], U};
map_modes({measure_pnr, Ms}, M) -> {measure_pnr, [M(X) || X <- Ms]};
map_modes({measure_threshold, Ms}, M) -> {measure_threshold, [M(X) || X <- Ms]};
map_modes({custom, Name, Args}, M) -> {custom, Name, Args#{modes => [M(X) || X <- maps:get(modes, Args, [])]}};
map_modes(I, M) when is_integer(element(2, I)) -> setelement(2, I, M(element(2, I))).

map_params({beamsplitter, A, B, T, P}, F) -> {beamsplitter, A, B, F(T), F(P)};
map_params({unitary, Ms, U}, F) -> {unitary, Ms, case element(1, U) of cmat -> U; _ -> F(U) end};
map_params({Tag, _} = I, _F) when Tag =:= measure_pnr; Tag =:= measure_threshold -> I;
map_params({custom, _, _} = I, _F) -> I;
map_params({Tag, M, X}, F) -> {Tag, M, F(X)};
map_params({Tag, M, X, Y}, F) -> {Tag, M, F(X), F(Y)}.

%%% FR-C06: symbols

%% symbols of the instructions and the input that have no default
-spec free_symbols(circuit()) -> [atom()].
free_symbols(#circuit{params = P} = C) -> [X || X <- instr_symbols(C), not maps:is_key(X, P)].

-spec all_symbols(circuit()) -> [atom()].
all_symbols(#circuit{params = P} = C) ->
    lists:usort(instr_symbols(C) ++ maps:keys(P) ++ lists:append([pq_symbolic:free_symbols(V) || V <- maps:values(P)])).

instr_symbols(#circuit{instrs = Is, input = In}) ->
    Exprs = lists:append([params_of(I) || I <- Is]) ++ input_exprs(In),
    lists:usort(lists:append([pq_symbolic:free_symbols(E) || E <- Exprs])).

params_of({beamsplitter, _, _, T, P}) -> [T, P];
params_of({unitary, _, {cmat, _, _, _}}) -> [];
params_of({unitary, _, U}) -> [U];
params_of({Tag, _}) when Tag =:= measure_pnr; Tag =:= measure_threshold -> [];
params_of({custom, _, _}) -> [];
params_of({_, _, X}) -> [X];
params_of({_, _, X, Y}) -> [X, Y].

input_exprs({gaussian, #gaussian_state{mean = M, cov = Cv}}) -> [E || E <- [M, Cv], not packed(E)];
input_exprs({fock_state, #fock_state{amps = A}}) -> [A || not packed(A)];
input_exprs(_) -> [].

packed(T) -> is_tuple(T) andalso lists:member(element(1, T), [cmat, rmat, imat, cvec, rvec]).

%%% FR-C08: serialisation

-spec to_binary(circuit()) -> binary().
to_binary(#circuit{} = C) -> term_to_binary({pq_circuit, ?BINARY_VERSION, C}, [deterministic]).

%% Untrusted input: decoded with [safe] (no new atoms), checked for funs, pids, ports and references, the record
%% shape checked, then validated.
-spec from_binary(binary()) -> result(circuit()).
from_binary(Bin) when is_binary(Bin) ->
    try binary_to_term(Bin, [safe]) of
        {pq_circuit, ?BINARY_VERSION, C} ->
            case plain(C) andalso is_record(C, circuit) andalso tuple_size(C) =:= record_info(size, circuit) of
                true ->
                    case validate(C) of
                        ok -> {ok, C};
                        {error, _} = E -> E
                    end;
                false -> {error, {bad_term, not_a_circuit}}
            end;
        {pq_circuit, V, _} -> {error, {bad_term, {version, V}}};
        _ -> {error, {bad_term, not_a_circuit}}
    catch
        error:badarg -> {error, {bad_term, decode}}
    end;
from_binary(_) -> {error, {bad_term, not_a_binary}}.

%% no funs, pids, ports or references anywhere in the term
plain(T) when is_function(T); is_pid(T); is_port(T); is_reference(T) -> false;
plain(T) when is_tuple(T) -> lists:all(fun plain/1, tuple_to_list(T));
plain(T) when is_list(T) -> plain_list(T);
plain(T) when is_map(T) -> lists:all(fun plain/1, maps:keys(T) ++ maps:values(T));
plain(_) -> true.

plain_list([]) -> true;
plain_list([H | T]) -> plain(H) andalso plain_list(T);
plain_list(Tail) -> plain(Tail).                        % an improper list's tail

%%% FR-C10: the Clements mesh of a numeric unitary, as instructions on modes 0..k-1: beamsplitters in order,
%%% then one phase shift per mode. Replaying them reproduces U to 1e-12 k (AN 9).

-spec mesh(pq_core:packed_cmatrix()) -> result([instruction()]).
mesh({cmat, K, K, _} = U) ->
    case pq_nif:unitary_decompose(U) of
        {ok, {Bs, Phases}} ->
            {ok, {rvec, Ps}} = pq_core:unpack(Phases),
            {ok, [{beamsplitter, A, B, T, P} || {A, B, T, P} <- Bs] ++
                 [{phaseshift, J, Ph} || {J, Ph} <- lists:enumerate(0, Ps)]};
        {error, _} = E -> E
    end;
mesh(U) -> {error, {bad_shape, #{unitary => U}}}.

%%% FR-C09: JSON interchange (Appendix B), in pq_json

-spec to_json(circuit()) -> result(binary()).
to_json(C) -> pq_json:to_json(C).

-spec from_json(binary()) -> result(circuit()).
from_json(Bin) -> pq_json:from_json(Bin).

-spec from_json(binary(), map()) -> result(circuit()).
from_json(Bin, Opts) -> pq_json:from_json(Bin, Opts).

%%% FR-C11: the gate library, as instruction lists

%% B(pi/4, 0)^dagger (S(z) x S(-z)) B(pi/4, 0), z = R e^(i Phi): the Strawberry Fields decomposition
%% FS 7: phase 0
-spec two_mode_squeezing(mode(), mode(), param()) -> [instruction()].
two_mode_squeezing(M1, M2, R) -> two_mode_squeezing(M1, M2, R, 0).

-spec two_mode_squeezing(mode(), mode(), param(), param()) -> [instruction()].
two_mode_squeezing(M1, M2, R, Phi) ->
    Q = {mul, [{rat, 1, 4}, {const, pi}]},
    [{beamsplitter, M1, M2, Q, 0},
     {squeezing, M1, R, Phi},
     {squeezing, M2, {mul, [-1, R]}, Phi},
     {beamsplitter, M1, M2, {mul, [-1, Q]}, 0}].

%% beamsplitter(pi/4, 0), phaseshift(Theta), beamsplitter(pi/4, 0), phaseshift(Phi)
-spec mzi(mode(), mode(), param(), param()) -> [instruction()].
mzi(M1, M2, Theta, Phi) ->
    Q = {mul, [{rat, 1, 4}, {const, pi}]},
    [{beamsplitter, M1, M2, Q, 0}, {phaseshift, M1, Theta}, {beamsplitter, M1, M2, Q, 0}, {phaseshift, M1, Phi}].
