%% The only module that calls Rust (FS 2.1). The functions below are replaced by the NIF library at load time.
-module(pq_nif).
-export([nif_version/0, check_version/0, cancel_token_new/0, cancel_token_set/1, rng_u64/3]).
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
          {ok, [non_neg_integer()]} | {error, {bad_count, non_neg_integer()}}.
rng_u64(_Seed, _Stream, _Count) -> erlang:nif_error(nif_not_loaded).
