%% The only module that calls Rust (FS 2.1). The functions below are replaced by the NIF library at load time.
-module(pq_nif).
-export([nif_version/0, check_version/0, cancel_token_new/0, cancel_token_set/1, rng_u64/3,
         permanent/2, permanent_repeated/4, bs_probability/4]).
-on_load(load/0).

%% The NIF {Major, Minor} this code was written against (FS 8.1).
-define(NIF_MAJOR, 0).
-define(NIF_MINOR, 1).

load() ->
    Dir = case code:priv_dir(pq) of
              {error, bad_name} -> filename:join(filename:dirname(filename:dirname(code:which(?MODULE))), "priv");
              D -> D
          end,
    erlang:load_nif(filename:join(Dir, "pq_nif"), 0).

-spec check_version() -> ok | {error, {nif_version_mismatch, map()}}.
check_version() ->
    case nif_version() of
        {?NIF_MAJOR, ?NIF_MINOR, _, _} -> ok;
        Found -> {error, {nif_version_mismatch, #{expected => {?NIF_MAJOR, ?NIF_MINOR}, found => Found}}}
    end.

-spec nif_version() -> {non_neg_integer(), non_neg_integer(), non_neg_integer(), binary()}.
nif_version() -> erlang:nif_error(nif_not_loaded).

-spec cancel_token_new() -> reference().
cancel_token_new() -> erlang:nif_error(nif_not_loaded).

-spec cancel_token_set(reference()) -> ok.
cancel_token_set(_Token) -> erlang:nif_error(nif_not_loaded).

%% The raw ChaCha8 stream (Seed, Stream) for the golden test GT-R01 (FS 8.3).
-spec rng_u64(non_neg_integer(), non_neg_integer(), non_neg_integer()) ->
          {ok, [non_neg_integer()]} | {error, {dimension_limit, #{count := non_neg_integer(), max := 65536}}}.
rng_u64(_Seed, _Stream, _Count) -> erlang:nif_error(nif_not_loaded).

%%% The permanent family (FS 8.4, AN 3; dirty CPU). Opts: cancel (a token), threads (default 1).
%%% Info carries the scale S of AN 2.4 (amendment A-1): the error bound is about n * eps * S.

-type opts() :: #{cancel => reference(), threads => pos_integer()}.
-type perm_info() :: #{scale := float(), positions := non_neg_integer(), blocks := non_neg_integer()}.
-type kernel_error() :: {error, pq_core:reason()}.

-spec permanent(pq_core:packed_cmatrix(), opts()) -> {ok, {pq_core:nif_complex(), perm_info()}} | kernel_error().
permanent(_M, _Opts) -> erlang:nif_error(nif_not_loaded).

-spec permanent_repeated(pq_core:packed_cmatrix(), [non_neg_integer()], [non_neg_integer()], opts()) ->
          {ok, {pq_core:nif_complex(), perm_info()}} | kernel_error().
permanent_repeated(_M, _RowMult, _ColMult, _Opts) -> erlang:nif_error(nif_not_loaded).

%% FR-B02: probability of OutputOcc from InputOcc through the unitary U (U is not checked here).
-spec bs_probability(pq_core:packed_cmatrix(), [non_neg_integer()], [non_neg_integer()], opts()) ->
          {ok, {float(), #{permanent_scale := float(), positions := non_neg_integer(),
                           blocks := non_neg_integer()}}} | kernel_error().
bs_probability(_U, _InputOcc, _OutputOcc, _Opts) -> erlang:nif_error(nif_not_loaded).
