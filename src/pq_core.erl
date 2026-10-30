%% Fundamental types, configuration, the error taxonomy and the packed-container codec (FS 2.1, 4.1, 8.2, 10,
%% Appendix E). Pure: callable without the application started.
-module(pq_core).
-compile({no_auto_import, [get/1]}).
-export([config/0, get/1, get/2, defaults/0,
         pack/1, unpack/1, validate/1,
         categories/0, is_reason/1]).

-export_type([mode/0, seed/0, domain/0, exponent/0, packed_cmatrix/0, packed_rmatrix/0, packed_imatrix/0,
              packed_cvector/0, packed_rvector/0, numeric_matrix/0, numeric_vector/0, nif_complex/0,
              reason/0]).

-type mode()     :: non_neg_integer().
-type seed()     :: 0..16#FFFFFFFFFFFFFFFF.
-type domain()   :: real | positive | nonnegative | complex.
-type exponent() :: integer() | {rat, integer(), pos_integer()} | float().

-type packed_cmatrix() :: {cmat, Rows :: pos_integer(), Cols :: pos_integer(), binary()}.
-type packed_rmatrix() :: {rmat, Rows :: pos_integer(), Cols :: pos_integer(), binary()}.
-type packed_imatrix() :: {imat, Rows :: non_neg_integer(), Cols :: pos_integer(), binary()}.
-type packed_cvector() :: {cvec, Len :: pos_integer(), binary()}.
-type packed_rvector() :: {rvec, Len :: pos_integer(), binary()}.
-type packed() :: packed_cmatrix() | packed_rmatrix() | packed_imatrix() | packed_cvector() | packed_rvector().

-type numeric_matrix() :: packed_cmatrix() | packed_rmatrix().
-type numeric_vector() :: packed_cvector() | packed_rvector().
-type nif_complex()    :: {Re :: float(), Im :: float()}.

-type reason() :: atom() | {atom(), term()}.

%% A complex entry for pack/1: a real number, {rat, N, D}, or {complex, Re, Im} with such parts.
-type real_in()    :: number() | {rat, integer(), pos_integer()}.
-type complex_in() :: real_in() | {complex, real_in(), real_in()}.
-type unpacked() :: {cmat, [[{complex, float(), float()}]]} | {rmat, [[float()]]} | {imat, [[non_neg_integer()]]}
                  | {cvec, [{complex, float(), float()}]} | {rvec, [float()]}.

%%% Configuration (Appendix E). Keys that depend on others (max_per_mode, max_attempts, max_inflight, chains) are
%%% resolved where they are used; here they are `derived`.

-spec defaults() -> #{atom() => term()}.
defaults() ->
    #{hbar => 2, unitarity_tol => 1.0e-10, physicality_tol => 1.0e-10, truncation_tol => 1.0e-6,
      strict => false, check_physical => false,
      max_symbolic_modes => 8, max_expr_nodes => 1000000, simplify_max_passes => 100,
      max_poly_terms => 100000, max_expand_power => 8, max_symbolic_photons => 6,
      sym_timeout_ms => 10000, sym_max_heap_words => 20000000,
      max_fock_dim => 5000000, max_dm_dim => 4096, max_enumeration => 1000000,
      max_photons => 16, max_per_mode => derived, max_clicks => 20,
      shots_mode => attempts, max_attempts => derived, parallel_threshold => 64, max_chunk => 10000,
      max_inflight => derived, max_retries => 2, kernel_hard_timeout => 30000,
      spec_cache_bytes => 256 * 1024 * 1024, allow_heterogeneous_arch => false, pq_nif_threads => 1,
      max_json_symbols => 1024, custom_instructions => [],
      burn_in => 100, thin => 10, chains => derived}.

%% The defaults overridden by the application environment of `pq`.
-spec config() -> #{atom() => term()}.
config() ->
    maps:merge(defaults(), maps:with(maps:keys(defaults()), maps:from_list(application:get_all_env(pq)))).

-spec get(atom()) -> term().
get(Key) ->
    case application:get_env(pq, Key) of
        {ok, V} -> V;
        undefined -> maps:get(Key, defaults())
    end.

%% A per-call value from Opts, else the configured one.
-spec get(atom(), map()) -> term().
get(Key, Opts) ->
    case Opts of
        #{Key := V} -> V;
        _ -> get(Key)
    end.

%%% Packed containers (FS 8.2): little-endian binary64, row-major, complex interleaved (re, im); imat u16.
%%% For tests and small inputs only; production paths build binaries directly.

-spec pack({cmat, [[complex_in()]]} | {rmat | rvec, [[real_in()]] | [real_in()]} | {imat, [[non_neg_integer()]]}
           | {cvec, [complex_in()]}) ->
          {ok, packed()} | {error, reason()}.
pack({Kind, Rows}) when (Kind =:= cmat orelse Kind =:= rmat orelse Kind =:= imat), is_list(Rows) ->
    case shape(Rows) of
        {ok, R, C} when Kind =:= imat; R > 0 ->
            entries(Kind, lists:append(Rows), fun(Bin) -> {Kind, R, C, Bin} end);
        {ok, R, C} -> {error, {bad_shape, #{kind => Kind, rows => R, cols => C}}};
        {error, _} = E -> E
    end;
pack({Kind, Entries}) when (Kind =:= cvec orelse Kind =:= rvec), is_list(Entries), Entries =/= [] ->
    entries(Kind, Entries, fun(Bin) -> {Kind, length(Entries), Bin} end);
pack(Other) ->
    {error, {bad_shape, #{input => kind_of(Other)}}}.

shape([First | _] = Rows) when is_list(First), First =/= [] ->
    C = length(First),
    case lists:all(fun(R) -> is_list(R) andalso length(R) =:= C end, Rows) of
        true -> {ok, length(Rows), C};
        false -> {error, {bad_shape, #{rows => length(Rows), ragged => true}}}
    end;
shape([]) -> {ok, 0, 0};
shape(_) -> {error, {bad_shape, #{rows => not_a_list_of_rows}}}.

entries(Kind, Es, Wrap) ->
    try << <<(entry(Kind, E))/binary>> || E <- Es >> of
        Bin -> {ok, Wrap(Bin)}
    catch
        throw:Reason -> {error, Reason}
    end.

entry(cmat, E) -> complex_bin(E);
entry(cvec, E) -> complex_bin(E);
entry(rmat, E) -> <<(real(E)):64/float-little>>;
entry(rvec, E) -> <<(real(E)):64/float-little>>;
entry(imat, N) when is_integer(N), N >= 0, N =< 65535 -> <<N:16/little>>;
entry(imat, N) when is_integer(N), N > 65535 -> throw({photon_number_overflow, N});
entry(imat, N) -> throw({bad_shape, #{entry => N}}).

complex_bin({complex, Re, Im}) -> <<(real(Re)):64/float-little, (real(Im)):64/float-little>>;
complex_bin(X) -> <<(real(X)):64/float-little, 0.0:64/float-little>>.

real(X) when is_integer(X) -> float(X);
real(X) when is_float(X) -> X;
real({rat, N, D}) when is_integer(N), is_integer(D), D > 0 -> N / D;
real(X) -> throw({bad_shape, #{entry => X}}).

-spec unpack(packed()) -> {ok, unpacked()} | {error, reason()}.
unpack(P) ->
    case validate(P) of
        ok -> decode(P);
        {error, _} = E -> E
    end.

decode({cmat, _R, C, Bin}) -> rows(cmat, C, cdecode(Bin, []));
decode({rmat, _R, C, Bin}) -> rows(rmat, C, rdecode(Bin, []));
decode({imat, _R, C, Bin}) -> {ok, {imat, chunk([N || <<N:16/little>> <= Bin], C)}};
decode({cvec, _L, Bin}) -> finite(cvec, cdecode(Bin, []));
decode({rvec, _L, Bin}) -> finite(rvec, rdecode(Bin, [])).

rows(Kind, C, Es) ->
    case finite(Kind, Es) of
        {ok, {Kind, L}} -> {ok, {Kind, chunk(L, C)}};
        E -> E
    end.

finite(_Kind, non_finite) -> {error, {non_finite, unpack}};
finite(Kind, L) -> {ok, {Kind, L}}.

%% Erlang has no NaN or infinity: the float match fails on them, which is the non_finite check.
cdecode(<<>>, Acc) -> lists:reverse(Acc);
cdecode(<<Re:64/float-little, Im:64/float-little, Rest/binary>>, Acc) -> cdecode(Rest, [{complex, Re, Im} | Acc]);
cdecode(_, _) -> non_finite.

rdecode(<<>>, Acc) -> lists:reverse(Acc);
rdecode(<<X:64/float-little, Rest/binary>>, Acc) -> rdecode(Rest, [X | Acc]);
rdecode(_, _) -> non_finite.

chunk([], _C) -> [];
chunk(L, C) -> {Row, Rest} = lists:split(C, L), [Row | chunk(Rest, C)].

%% Shape check of a packed container: the binary size equals entry size times its dimensions.
-spec validate(term()) -> ok | {error, reason()}.
validate({cmat, R, C, B}) when is_integer(R), R > 0, is_integer(C), C > 0, is_binary(B) -> size_is(B, 16 * R * C, {R, C});
validate({rmat, R, C, B}) when is_integer(R), R > 0, is_integer(C), C > 0, is_binary(B) -> size_is(B, 8 * R * C, {R, C});
validate({imat, R, C, B}) when is_integer(R), R >= 0, is_integer(C), C > 0, is_binary(B) -> size_is(B, 2 * R * C, {R, C});
validate({cvec, L, B}) when is_integer(L), L > 0, is_binary(B) -> size_is(B, 16 * L, L);
validate({rvec, L, B}) when is_integer(L), L > 0, is_binary(B) -> size_is(B, 8 * L, L);
validate(Other) -> {error, {bad_shape, #{input => kind_of(Other)}}}.

size_is(B, Expected, Dims) when byte_size(B) =:= Expected -> _ = Dims, ok;
size_is(B, Expected, Dims) -> {error, {bad_shape, #{dims => Dims, bytes => byte_size(B), expected => Expected}}}.

kind_of(T) when is_tuple(T), tuple_size(T) > 0, is_atom(element(1, T)) -> element(1, T);
kind_of(_) -> unknown.

%%% Error taxonomy (FS 10.2, with amendments A-3 and B-7). Every Reason is {Category, Detail} or a bare category.

-spec categories() -> #{atom() => error | warning}.
categories() ->
    Errors = [mode_out_of_range, duplicate_modes, index_out_of_range, bad_shape, not_unitary, not_symplectic,
              non_physical_state, unbound_symbols, param_conflict, unsupported_instruction, unsupported_input,
              measured_mode_reuse, dimension_limit, symbolic_limit, negative_probability, singular_matrix,
              simplify_diverged, differentiate_unsupported, hbar_mismatch, postselection_exhausted,
              predicate_not_distributable, cancelled, timeout, node_down, chunk_failed, nif_version_mismatch,
              nif_panic, photon_number_overflow, bad_term, non_finite, undefined_statistic,
              json_version, json_input, json_instruction, json_parse, json_symbol, json_number,
              sym_parse, domain_violation, param_cycle, symbolic_unsupported, bad_expr],
    Warnings = [unused_bindings, truncation, gbs_truncation],
    maps:merge(maps:from_keys(Errors, error), maps:from_keys(Warnings, warning)).

-spec is_reason(term()) -> boolean().
is_reason({Cat, _Detail}) when is_atom(Cat) -> maps:is_key(Cat, categories());
is_reason(Cat) when is_atom(Cat) -> maps:is_key(Cat, categories());
is_reason(_) -> false.
