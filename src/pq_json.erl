%% JSON interchange of circuits (FR-C09, FS Appendix B), with the OTP json module. to_json/1 writes canonical
%% expressions in the parameter subset of B.3 (with amendment B-9's rational exponents); from_json/1,2 reads them
%% back, atom-safely: every symbol name of the document is checked (and counted against max_json_symbols, or
%% matched against allowed_symbols) before any atom is created, then parsed and canonicalised.
-module(pq_json).
-include("pq.hrl").
-export([to_json/1, from_json/1, from_json/2]).

-define(V, 1).

-spec to_json(#circuit{}) -> {ok, binary()} | {error, pq_core:reason()}.
to_json(#circuit{} = C) ->
    try
        Doc = #{<<"pq">> => ?V,
                <<"n_modes">> => C#circuit.n_modes,
                <<"input">> => input_out(C#circuit.input),
                <<"symbols">> => maps:from_list([{atom_to_binary(K), atom_to_binary(D)} || K := D <- C#circuit.symbols]),
                <<"params">> => maps:from_list([{atom_to_binary(K), param_out(V)} || K := V <- C#circuit.params]),
                <<"instructions">> => [instr_out(I) || I <- C#circuit.instrs],
                <<"meta">> => meta_out(C#circuit.meta)},
        Doc2 = case C#circuit.hbar of undefined -> Doc; H -> Doc#{<<"hbar">> => H} end,
        {ok, iolist_to_binary(json:encode(Doc2))}
    catch
        throw:{error, _} = E -> E
    end.

-spec fail(pq_core:reason()) -> no_return().
fail(R) -> throw({error, R}).

input_out(vacuum) -> #{<<"kind">> => <<"vacuum">>};
input_out({fock, Occ}) -> #{<<"kind">> => <<"fock">>, <<"occupation">> => Occ};
input_out({gaussian, #gaussian_state{mean = {rvec, _, _} = M, cov = {rmat, _, _, _} = Cv, hbar = H}}) ->
    {ok, {rvec, Mean}} = pq_core:unpack(M),
    {ok, {rmat, Rows}} = pq_core:unpack(Cv),
    #{<<"kind">> => <<"gaussian">>, <<"mean">> => Mean, <<"cov">> => lists:append(Rows), <<"hbar">> => H};
input_out({fock_state, _}) -> fail({json_input, fock_state});
input_out(_) -> fail({json_input, symbolic_gaussian}).

%% a parameter: a JSON number when it is a plain integer or float, else its canonical text (exact numbers stay exact)
param_out(V) ->
    case pq_symbolic:canonicalise(V) of
        {ok, N} when is_integer(N); is_float(N) -> N;
        {ok, E} ->
            pq_sym_print:in_subset(E) orelse fail({json_input, {expression, E}}),
            pq_sym_print:ascii(E);
        {error, R} -> fail(R)
    end.

instr_out({beamsplitter, A, B, T, P}) -> #{<<"op">> => <<"beamsplitter">>, <<"m1">> => A, <<"m2">> => B, <<"theta">> => param_out(T), <<"phi">> => param_out(P)};
instr_out({phaseshift, M, P}) -> #{<<"op">> => <<"phaseshift">>, <<"mode">> => M, <<"phi">> => param_out(P)};
instr_out({displacement, M, R, P}) -> #{<<"op">> => <<"displacement">>, <<"mode">> => M, <<"r">> => param_out(R), <<"phi">> => param_out(P)};
instr_out({squeezing, M, R, P}) -> #{<<"op">> => <<"squeezing">>, <<"mode">> => M, <<"r">> => param_out(R), <<"phi">> => param_out(P)};
instr_out({kerr, M, K}) -> #{<<"op">> => <<"kerr">>, <<"mode">> => M, <<"kappa">> => param_out(K)};
instr_out({loss, M, E}) -> #{<<"op">> => <<"loss">>, <<"mode">> => M, <<"eta">> => param_out(E)};
instr_out({unitary, Ms, {cmat, _, _, _} = U}) ->
    {ok, {cmat, Rows}} = pq_core:unpack(U),
    #{<<"op">> => <<"unitary">>, <<"modes">> => Ms,
      <<"matrix">> => [#{<<"re">> => Re, <<"im">> => Im} || {complex, Re, Im} <- lists:append(Rows)]};
instr_out({unitary, Ms, {matrix, Rows}}) ->
    #{<<"op">> => <<"unitary">>, <<"modes">> => Ms, <<"matrix">> => [param_out(E) || E <- lists:append(Rows)]};
instr_out({measure_pnr, Ms}) -> #{<<"op">> => <<"measure_pnr">>, <<"modes">> => Ms};
instr_out({measure_threshold, Ms}) -> #{<<"op">> => <<"measure_threshold">>, <<"modes">> => Ms};
instr_out({measure_homodyne, M, P}) -> #{<<"op">> => <<"measure_homodyne">>, <<"mode">> => M, <<"phi">> => param_out(P)};
%% feed-forward (Phase 4): the gate, and each outcome symbol with the mode it is read from
instr_out({feed_forward, I, Map}) ->
    #{<<"op">> => <<"feed_forward">>, <<"gate">> => instr_out(I),
      <<"outcomes">> => maps:from_list([{atom_to_binary(K), M} || K := M <- Map])};
instr_out({conditional, Pred, Is}) ->
    #{<<"op">> => <<"conditional">>, <<"predicate">> => pred_out(Pred), <<"gates">> => [instr_out(I) || I <- Is]};
instr_out({custom, _, _}) -> fail({json_instruction, custom});
instr_out(I) -> fail({json_instruction, element(1, I)}).

%% a predicate (FS 4.7) as data: {"op": "eq" | "ge" | "le", "mode", "value"}, {"op": "in", "mode", "values"},
%% {"op": "total_eq" | "total_ge" | "total_le", "modes", "value"}, {"op": "and" | "or", "args"}, {"op": "not", "arg"}
pred_out({Op, Ps}) when Op =:= 'and'; Op =:= 'or' -> #{<<"op">> => atom_to_binary(Op), <<"args">> => [pred_out(P) || P <- Ps]};
pred_out({'not', P}) -> #{<<"op">> => <<"not">>, <<"arg">> => pred_out(P)};
pred_out({in, M, Vs}) -> #{<<"op">> => <<"in">>, <<"mode">> => M, <<"values">> => Vs};
pred_out({Op, Ms, V}) when is_list(Ms) -> #{<<"op">> => atom_to_binary(Op), <<"modes">> => Ms, <<"value">> => V};
pred_out({Op, M, V}) when Op =:= eq; Op =:= ge; Op =:= le -> #{<<"op">> => atom_to_binary(Op), <<"mode">> => M, <<"value">> => V};
pred_out(P) -> fail({json_instruction, {predicate, P}}).

pred_in(#{<<"op">> := Op} = P) when Op =:= <<"and">>; Op =:= <<"or">> ->
    {binary_to_existing_atom(Op), [pred_in(X) || X <- list(maps:get(<<"args">>, P, undefined))]};
pred_in(#{<<"op">> := <<"not">>, <<"arg">> := A}) -> {'not', pred_in(A)};
pred_in(#{<<"op">> := <<"in">>, <<"mode">> := M, <<"values">> := Vs}) -> {in, int(M), [num(V) || V <- list(Vs)]};
pred_in(#{<<"op">> := Op, <<"modes">> := Ms, <<"value">> := V})
  when Op =:= <<"total_eq">>; Op =:= <<"total_ge">>; Op =:= <<"total_le">> ->
    {binary_to_existing_atom(Op), [int(M) || M <- list(Ms)], num(V)};
pred_in(#{<<"op">> := Op, <<"mode">> := M, <<"value">> := V}) when Op =:= <<"eq">>; Op =:= <<"ge">>; Op =:= <<"le">> ->
    {binary_to_existing_atom(Op), int(M), num(V)};
pred_in(P) -> fail({json_instruction, {predicate, P}}).

meta_out(M) ->
    maps:from_list([begin
                        is_atom(K) orelse fail({json_input, meta}),
                        (is_binary(V) orelse is_number(V) orelse is_boolean(V)) orelse fail({json_input, meta}),
                        {atom_to_binary(K), V}
                    end || K := V <- M]).

%%% Reading

-spec from_json(binary()) -> {ok, #circuit{}} | {error, pq_core:reason()}.
from_json(Bin) -> from_json(Bin, #{}).

%% Opts: allowed_symbols => [atom()] (existing atoms only), max_json_symbols
-spec from_json(binary(), map()) -> {ok, #circuit{}} | {error, pq_core:reason()}.
from_json(Bin, Opts) when is_binary(Bin) ->
    try
        Doc = decode(Bin),
        is_map(Doc) orelse fail({json_parse, 1}),
        case maps:get(<<"pq">>, Doc, undefined) of
            ?V -> ok;
            V -> fail({json_version, V})
        end,
        Atoms = atoms(Doc, Opts),
        A = fun(Name) -> maps:get(Name, Atoms) end,
        N = int(maps:get(<<"n_modes">>, Doc, undefined)),
        Hbar = case maps:find(<<"hbar">>, Doc) of {ok, H} when is_number(H) -> H; error -> undefined; {ok, H} -> fail({json_number, H}) end,
        Syms = maps:from_list([{A(K), domain(D)} || K := D <- maps:get(<<"symbols">>, Doc, #{})]),
        P = fun(V) -> param_in(V, A) end,
        C = #circuit{n_modes = N, hbar = Hbar,
                     input = input_in(maps:get(<<"input">>, Doc, #{<<"kind">> => <<"vacuum">>}), N, Hbar),
                     symbols = Syms,
                     params = maps:from_list([{A(K), P(V)} || K := V <- maps:get(<<"params">>, Doc, #{})]),
                     instrs = [instr_in(I, P, A) || I <- list(maps:get(<<"instructions">>, Doc, []))],
                     meta = maps:from_list([{A(K), meta_value(V)} || K := V <- maps:get(<<"meta">>, Doc, #{})])},
        case pq_circuit:validate(C) of
            ok -> {ok, C};
            {error, _} = Invalid -> Invalid
        end
    catch
        throw:{error, _} = E -> E
    end;
from_json(_, _) -> {error, {json_parse, 1}}.

decode(Bin) ->
    try json:decode(Bin)
    catch error:R -> fail({json_parse, R})
    end.

int(N) when is_integer(N), N >= 0 -> N;
int(V) -> fail({json_number, V}).

list(L) when is_list(L) -> L;
list(V) -> fail({json_input, V}).

domain(<<"real">>) -> real;
domain(<<"positive">>) -> positive;
domain(<<"nonnegative">>) -> nonnegative;
domain(<<"complex">>) -> complex;
domain(D) -> fail({json_symbol, D}).

meta_value(V) when is_binary(V); is_number(V); is_boolean(V) -> V;
meta_value(_) -> fail({json_input, meta}).

%% every name in the document: parameter strings, symbols, params and meta keys; checked, then made atoms
atoms(Doc, Opts) ->
    Strings = [V || V <- maps:values(maps:get(<<"params">>, Doc, #{})), is_binary(V)] ++
              lists:append([param_strings(I) || I <- list(maps:get(<<"instructions">>, Doc, []))]),
    FromStrings = lists:append([case pq_sym_parse:names(S) of
                                    {ok, Ns} -> Ns;
                                    {error, {sym_parse, Pos}} -> fail({json_parse, Pos})
                                end || S <- Strings]),
    Keys = lists:append([maps:keys(maps:get(K, Doc, #{})) || K <- [<<"symbols">>, <<"params">>, <<"meta">>]]) ++
           lists:append([outcome_names(I) || I <- list(maps:get(<<"instructions">>, Doc, []))]),
    Names = lists:usort(FromStrings ++ Keys),
    [valid_name(Nm) orelse fail({json_symbol, Nm}) || Nm <- Names],
    Max = maps:get(max_json_symbols, Opts, pq_core:get(max_json_symbols)),
    length(Names) =< Max orelse fail({json_symbol, {too_many, length(Names), Max}}),
    case maps:find(allowed_symbols, Opts) of
        {ok, Allowed} ->
            maps:from_list([{Nm, existing(Nm, Allowed)} || Nm <- Names]);
        error ->
            maps:from_list([{Nm, binary_to_atom(Nm)} || Nm <- Names])
    end.

existing(Nm, Allowed) ->
    try binary_to_existing_atom(Nm) of
        A -> lists:member(A, Allowed) orelse fail({json_symbol, Nm}), A
    catch error:badarg -> fail({json_symbol, Nm})
    end.

valid_name(<<C, Rest/binary>>) when (C >= $a andalso C =< $z) orelse (C >= $A andalso C =< $Z) ->
    byte_size(Rest) < 64 andalso
        lists:all(fun(X) -> (X >= $a andalso X =< $z) orelse (X >= $A andalso X =< $Z) orelse (X >= $0 andalso X =< $9) orelse X =:= $_ end,
                  binary_to_list(Rest)) andalso
        not lists:member(<<C, Rest/binary>>, [<<"pi">>, <<"e">>, <<"i">>, <<"sin">>, <<"cos">>, <<"exp">>, <<"conj">>, <<"sqrt">>]);
valid_name(_) -> false.

%% a feed-forward gate's outcome symbols (in a conditional's gates too)
outcome_names(#{<<"op">> := <<"feed_forward">>, <<"outcomes">> := O}) when is_map(O) -> maps:keys(O);
outcome_names(#{<<"op">> := <<"conditional">>, <<"gates">> := Gs}) when is_list(Gs) -> lists:append([outcome_names(G) || G <- Gs]);
outcome_names(_) -> [].

%% the parameter strings of one instruction: its parameter fields and string matrix entries (not "op")
param_strings(#{<<"op">> := <<"feed_forward">>} = I) ->
    param_strings(maps:get(<<"gate">>, I, #{}));
param_strings(#{<<"op">> := <<"conditional">>} = I) ->
    lists:append([param_strings(G) || G <- case maps:get(<<"gates">>, I, []) of L when is_list(L) -> L; _ -> [] end]);
param_strings(I) when is_map(I) ->
    [V || K <- [<<"theta">>, <<"phi">>, <<"r">>, <<"kappa">>, <<"eta">>], {ok, V} <- [maps:find(K, I)], is_binary(V)] ++
        [E || E <- case maps:get(<<"matrix">>, I, []) of L when is_list(L) -> L; _ -> [] end, is_binary(E)];
param_strings(_) -> [].

%% a number stays a number; a string is parsed (atoms already made) and canonicalised
param_in(V, _A) when is_number(V) -> V;
param_in(S, A) when is_binary(S) ->
    {ok, Names} = pq_sym_parse:names(S),
    Allowed = [A(N) || N <- Names],
    case pq_sym_parse:parse(S, #{allowed_symbols => Allowed}) of
        {ok, E} ->
            case pq_symbolic:canonicalise(E) of
                {ok, Cn} -> Cn;
                {error, R} -> fail(R)
            end;
        {error, {sym_parse, Pos}} -> fail({json_parse, Pos})
    end;
param_in(V, _A) -> fail({json_number, V}).

input_in(#{<<"kind">> := <<"vacuum">>}, _N, _H) -> vacuum;
input_in(#{<<"kind">> := <<"fock">>, <<"occupation">> := Occ}, _N, _H) -> {fock, [int(X) || X <- list(Occ)]};
input_in(#{<<"kind">> := <<"gaussian">>, <<"mean">> := Mean, <<"cov">> := Cov} = G, N, H) ->
    Hg = maps:get(<<"hbar">>, G, case H of undefined -> pq_core:get(hbar); _ -> H end),
    (H =:= undefined orelse H == Hg) orelse fail({hbar_mismatch, H, Hg}),
    length(Cov) =:= 4 * N * N andalso length(Mean) =:= 2 * N orelse fail({json_input, gaussian}),
    Rows = chunk(Cov, 2 * N),
    S = pq_state:from_mean_cov([num(X) || X <- Mean], [[num(X) || X <- R] || R <- Rows], Hg),
    {gaussian, S};
input_in(#{<<"kind">> := <<"fock_state">>}, _N, _H) -> fail({json_input, fock_state});
input_in(I, _N, _H) -> fail({json_input, I}).

num(X) when is_number(X) -> X;
num(X) -> fail({json_number, X}).

chunk([], _) -> [];
chunk(L, K) -> {A, B} = lists:split(K, L), [A | chunk(B, K)].

instr_in(#{<<"op">> := <<"feed_forward">>, <<"gate">> := G, <<"outcomes">> := O}, P, A) when is_map(O) ->
    {feed_forward, instr_in(G, P, A), maps:from_list([{A(K), int(M)} || K := M <- O])};
instr_in(#{<<"op">> := <<"conditional">>, <<"predicate">> := Pr, <<"gates">> := Gs}, P, A) ->
    {conditional, pred_in(Pr), [instr_in(G, P, A) || G <- list(Gs)]};
instr_in(I, P, _A) -> instr_in(I, P).

instr_in(#{<<"op">> := Op} = I, P) ->
    G = fun(K) -> case maps:find(K, I) of {ok, V} -> V; error -> fail({json_instruction, {Op, missing, K}}) end end,
    case Op of
        <<"beamsplitter">> -> {beamsplitter, int(G(<<"m1">>)), int(G(<<"m2">>)), P(G(<<"theta">>)), P(G(<<"phi">>))};
        <<"phaseshift">> -> {phaseshift, int(G(<<"mode">>)), P(G(<<"phi">>))};
        <<"displacement">> -> {displacement, int(G(<<"mode">>)), P(G(<<"r">>)), P(G(<<"phi">>))};
        <<"squeezing">> -> {squeezing, int(G(<<"mode">>)), P(G(<<"r">>)), P(G(<<"phi">>))};
        <<"kerr">> -> {kerr, int(G(<<"mode">>)), P(G(<<"kappa">>))};
        <<"loss">> -> {loss, int(G(<<"mode">>)), P(G(<<"eta">>))};
        <<"measure_pnr">> -> {measure_pnr, [int(M) || M <- list(G(<<"modes">>))]};
        <<"measure_threshold">> -> {measure_threshold, [int(M) || M <- list(G(<<"modes">>))]};
        <<"measure_homodyne">> -> {measure_homodyne, int(G(<<"mode">>)), P(G(<<"phi">>))};
        <<"unitary">> ->
            Ms = [int(M) || M <- list(G(<<"modes">>))],
            K = length(Ms),
            Entries = list(G(<<"matrix">>)),
            length(Entries) =:= K * K orelse fail({json_instruction, {unitary, size}}),
            {unitary, Ms, matrix_in(Entries, K, P)};
        _ -> fail({json_instruction, Op})
    end;
instr_in(I, _P) -> fail({json_instruction, I}).

%% all {re, im} objects: a packed matrix; any parameter string: a symbolic matrix
matrix_in(Entries, K, P) ->
    case lists:all(fun(#{<<"re">> := R, <<"im">> := I}) -> is_number(R) andalso is_number(I); (_) -> false end, Entries) of
        true ->
            {ok, M} = pq_core:pack({cmat, chunk([{complex, R, I} || #{<<"re">> := R, <<"im">> := I} <- Entries], K)}),
            M;
        false ->
            {matrix, chunk([case E of
                                #{<<"re">> := R, <<"im">> := I} -> {complex, num(R), num(I)};
                                _ -> P(E)
                            end || E <- Entries], K)}
    end.
