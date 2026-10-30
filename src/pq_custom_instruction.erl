%% Behaviour of a custom instruction (FS 4.3): a macro that lowering replaces by built-in instructions. Modules
%% are registered in the application environment of `pq`, key `custom_instructions`.
-module(pq_custom_instruction).

-callback name() -> atom().
-callback validate(Args :: map(), NModes :: pos_integer()) -> ok | {error, term()}.
-callback expand(Args :: map(), Binding :: #{atom() => term()}) -> {ok, [tuple()]} | {error, term()}.
