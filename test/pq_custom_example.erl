%% A custom instruction for the tests: a 50:50 beamsplitter written as a macro (FS 4.3).
-module(pq_custom_example).
-behaviour(pq_custom_instruction).
-export([name/0, validate/2, expand/2]).

name() -> half_mirror.
validate(#{modes := [A, B]}, N) when A < N, B < N, A =/= B -> ok;
validate(_, _) -> {error, bad_modes}.
expand(#{modes := [A, B]}, _Binding) -> {ok, [{beamsplitter, A, B, {mul, [{rat, 1, 4}, {const, pi}]}, 0}]}.
