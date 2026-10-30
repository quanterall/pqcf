%% The only module that calls Rust (FS 2.1). The functions below are replaced by the NIF library at load time.
-module(pq_nif).
-export([nif_version/0, check_version/0, cancel_token_new/0, cancel_token_set/1, rng_u64/3,
         permanent/2, permanent_repeated/4, bs_probability/4, sample_boson/4,
         hafnian/2, hafnian_repeated/3, loop_hafnian/4, torontonian/2]).
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

%%% Boson sampling (FS 8.4, AN 5). Shot i uses stream (seed, stream_start + i), so samples do not depend on how
%%% shots are split between calls. Opts: algorithm (clifford_clifford, default | brute_force), seed (default 0),
%%% stream_start (default 0), max_enumeration (brute_force, default 1,000,000), cancel.

-spec sample_boson(pq_core:packed_cmatrix(), [non_neg_integer()], non_neg_integer(), map()) ->
          {ok, {pq_core:packed_imatrix(), #{algorithm := atom(), shots := non_neg_integer()}}} | kernel_error().
sample_boson(U, InputOcc, NShots, Opts) ->
    case sample_boson_raw(U, InputOcc, NShots, Opts) of
        {error, {probability_sum, Total}} -> error({pq_bug, {probability_sum, Total}});
        Other -> Other
    end.

sample_boson_raw(_U, _InputOcc, _NShots, _Opts) -> erlang:nif_error(nif_not_loaded).

%%% Hafnians and the torontonian (FS 8.4, AN 4; dirty CPU). Info carries the scale S (amendment A-1).

-type haf_info() :: #{scale := float(), subsets := non_neg_integer()}.

-spec hafnian(pq_core:packed_cmatrix(), opts()) -> {ok, {pq_core:nif_complex(), haf_info()}} | kernel_error().
hafnian(_M, _Opts) -> erlang:nif_error(nif_not_loaded).

-spec hafnian_repeated(pq_core:packed_cmatrix(), [non_neg_integer()], opts()) ->
          {ok, {pq_core:nif_complex(), haf_info()}} | kernel_error().
hafnian_repeated(_M, _Mult, _Opts) -> erlang:nif_error(nif_not_loaded).

%% the hafnian of M with its diagonal replaced by Diag, indices repeated by Mult (FS A.9)
-spec loop_hafnian(pq_core:packed_cmatrix(), pq_core:packed_cvector(), [non_neg_integer()], opts()) ->
          {ok, {pq_core:nif_complex(), haf_info()}} | kernel_error().
loop_hafnian(_M, _Diag, _Mult, _Opts) -> erlang:nif_error(nif_not_loaded).

%% tor(O) for O = I - Q^-1 restricted to the clicked modes, laid out (modes, conjugate modes) (FS A.10)
-spec torontonian(pq_core:packed_cmatrix(), opts()) ->
          {ok, {float(), #{scale := float(), terms := non_neg_integer()}}} | kernel_error().
torontonian(_M, _Opts) -> erlang:nif_error(nif_not_loaded).
