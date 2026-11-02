%% An optional bridge to FORM (J.A.M. Vermaseren et al.; version 5, arXiv:2601.19982) for large polynomial work:
%% expansion and symbolic permanents beyond what pq_symbolic is meant to hold (docs/decisions/SY-FORM-comparison.md).
%%
%% FORM is GPL-3 and PQCF Apache-2.0: FORM runs as a separate program (the executable `form`, or the option
%% form_executable), never linked or bundled. Without it, available/0 is false and every call returns
%% {error, {symbolic_unsupported, #{form => not_available}}}.
%%
%% The fragment: symbols, exact integers, rationals and Gaussian rationals, i, sums, products and integer powers
%% (negative ones included). Anything else -- exponentials, radicals, pi, floats, matrices where a scalar is
%% expected -- is {symbolic_unsupported, ...}: those stay with pq_symbolic. Symbols are renamed s1, s2, ... for FORM
%% (it reserves names with '_') and back. FORM writes the result to a file; it is read by a parser of FORM's
%% sum-of-terms output (megabytes, beyond pq_sym_parse's 64 KB) and returned in pq_symbolic's canonical form, so a
%% result from FORM and one from pq_symbolic compare with =:=.
-module(pq_form).
-export([available/0, available/1, version/0, expand/1, expand/2, permanent/1, permanent/2]).

-define(TIMEOUT, 600000).

-spec available() -> boolean().
available() -> available(#{}).

-spec available(map()) -> boolean().
available(Opts) -> executable(Opts) =/= false.

%% the first line FORM prints, or not_available
-spec version() -> binary() | not_available.
version() ->
    case executable(#{}) of
        false -> not_available;
        Exe ->
            Dir = tmpdir(),
            F = filename:join(Dir, "v.frm"),
            ok = file:write_file(F, <<".end\n">>),
            Out = case run(Exe, F, #{quiet => false}) of {ok, O} -> O; {error, _} -> <<>> end,
            _ = file:del_dir_r(Dir),
            hd(binary:split(string:trim(Out), <<"\n">>) ++ [<<>>])
    end.

%% E expanded by FORM, in canonical form. Opts: form_executable, form_timeout_ms (default 600 s), and the
%% pq_symbolic options for the canonicalisation of the result.
-spec expand(pq_symbolic:sym_expr()) -> {ok, pq_symbolic:sym_expr()} | {error, pq_core:reason()}.
expand(E) -> expand(E, #{}).

-spec expand(pq_symbolic:sym_expr(), map()) -> {ok, pq_symbolic:sym_expr()} | {error, pq_core:reason()}.
expand(E, Opts) ->
    with_form(Opts, fun(Names) -> "Local F = " ++ to_form(E, Names) ++ ";\n" end, [E]).

%% the permanent of a square matrix ({matrix, Rows}) by FORM, row by row with the terms that use a column twice
%% dropped as they appear (the expansion of the comparison benchmark)
-spec permanent(pq_symbolic:sym_expr()) -> {ok, pq_symbolic:sym_expr()} | {error, pq_core:reason()}.
permanent(M) -> permanent(M, #{}).

-spec permanent(pq_symbolic:sym_expr(), map()) -> {ok, pq_symbolic:sym_expr()} | {error, pq_core:reason()}.
permanent({matrix, Rows} = M, Opts) when is_list(Rows), Rows =/= [] ->
    K = length(Rows),
    case lists:all(fun(R) -> is_list(R) andalso length(R) =:= K end, Rows) of
        false -> {error, {bad_shape, permanent}};
        true ->
            with_form(Opts, fun(Names) ->
                                    Xs = [io_lib:format("x~b", [J]) || J <- lists:seq(1, K)],
                                    "Symbols " ++ lists:join(",", Xs) ++ ";\nLocal F = 1;\n.sort\n" ++
                                        lists:append(
                                          [["Multiply (", lists:join("+", [["(", to_form(Ejk, Names), ")*", X]
                                                                            || {Ejk, X} <- lists:zip(R, Xs)]), ");\n",
                                            [["id ", X, "^2 = 0;\n"] || X <- Xs], ".sort\n"] || R <- Rows]) ++
                                        "id " ++ lists:join("*", Xs) ++ " = 1;\n"
                            end, [M])
    end;
permanent(_, _Opts) -> {error, {bad_shape, permanent}}.

%%% The run

with_form(Opts, Body, Exprs) ->
    case executable(Opts) of
        false -> {error, {symbolic_unsupported, #{form => not_available}}};
        Exe ->
            Dir = tmpdir(),
            try
                %% the names as written (free_symbols would simplify first: x - x has none, but FORM reads x)
                Syms = lists:usort(lists:append([syms(X) || X <- Exprs])),
                Names = maps:from_list([{S, "s" ++ integer_to_list(K)} || {K, S} <- lists:enumerate(Syms)]),
                Script = iolist_to_binary(
                           [[["Symbols ", lists:join(",", maps:values(Names)), ";\n"] || map_size(Names) > 0],
                            Body(Names), ".sort\n",
                            "#write <", filename:join(Dir, "out.txt"), "> \"%E\", F\n.end\n"]),
                File = filename:join(Dir, "in.frm"),
                ok = file:write_file(File, Script),
                case run(Exe, File, Opts) of
                    {ok, _Log} ->
                        {ok, Text} = file:read_file(filename:join(Dir, "out.txt")),
                        Back = maps:from_list([{list_to_binary(N), S} || S := N <- Names]),
                        pq_symbolic:canonicalise(from_form(Text, Back), Opts);
                    {error, _} = Err -> Err
                end
            catch
                throw:{unsupported, Node} -> {error, {symbolic_unsupported, #{form => Node}}};
                throw:{form_output, What} -> {error, {symbolic_unsupported, #{form_output => What}}}
            after
                _ = file:del_dir_r(Dir)
            end
    end.

syms({sym, S}) -> [S];
syms(T) when is_tuple(T) -> syms(tuple_to_list(T));
syms(L) when is_list(L) -> lists:append([syms(X) || X <- L]);
syms(_) -> [].

executable(Opts) ->
    case maps:get(form_executable, Opts, pq_core_get(form_executable)) of
        undefined -> os:find_executable("form");
        Path when is_list(Path) -> case filelib:is_regular(Path) of true -> Path; false -> false end;
        _ -> false
    end.

pq_core_get(K) ->
    case application:get_env(pq, K) of {ok, V} -> V; undefined -> undefined end.

tmpdir() ->
    Base = os:getenv("TMPDIR", "/tmp"),
    Dir = filename:join(Base, "pq_form_" ++ integer_to_list(erlang:unique_integer([positive])) ++ "_" ++
                              os:getpid()),
    ok = file:make_dir(Dir),
    Dir.

%% FORM as a separate OS process; killed at the timeout
run(Exe, File, Opts) ->
    Timeout = maps:get(form_timeout_ms, Opts, ?TIMEOUT),
    Args = case maps:get(quiet, Opts, true) of true -> ["-q", File]; false -> [File] end,
    Port = open_port({spawn_executable, Exe}, [{args, Args}, exit_status, stderr_to_stdout, binary,
                                               {cd, filename:dirname(File)}]),
    collect(Port, Timeout, erlang:monotonic_time(millisecond) + Timeout, []).

collect(Port, Timeout, Deadline, Acc) ->
    Left = max(0, Deadline - erlang:monotonic_time(millisecond)),
    receive
        {Port, {data, D}} -> collect(Port, Timeout, Deadline, [D | Acc]);
        {Port, {exit_status, 0}} -> {ok, iolist_to_binary(lists:reverse(Acc))};
        {Port, {exit_status, S}} ->
            Log = iolist_to_binary(lists:reverse(Acc)),
            {error, {symbolic_unsupported, #{form_exit => S, log => binary:part(Log, max(0, byte_size(Log) - 400),
                                                                                  min(400, byte_size(Log)))}}}
    after Left ->
            _ = case erlang:port_info(Port, os_pid) of
                    {os_pid, Pid} -> os:cmd("kill -9 " ++ integer_to_list(Pid));
                    undefined -> ok
                end,
            catch port_close(Port),
            {error, {timeout, #{form => Timeout}}}
    end.

%%% pq_symbolic expressions to FORM's language: the polynomial fragment only

to_form({sym, S}, Names) -> maps:get(S, Names);
to_form(N, _) when is_integer(N) -> "(" ++ integer_to_list(N) ++ ")";
to_form({rat, N, D}, _) when is_integer(N), is_integer(D), D > 0 -> io_lib:format("(~b/~b)", [N, D]);
to_form({complex, Re, Im}, Names) -> "(" ++ to_form(Re, Names) ++ "+" ++ to_form(Im, Names) ++ "*i_)";
to_form({add, L}, Names) when is_list(L), L =/= [] -> "(" ++ lists:join("+", [to_form(X, Names) || X <- L]) ++ ")";
to_form({add, []}, _) -> "(0)";
to_form({mul, L}, Names) when is_list(L), L =/= [] -> "(" ++ lists:join("*", [to_form(X, Names) || X <- L]) ++ ")";
to_form({mul, []}, _) -> "(1)";
to_form({pow, B, N}, Names) when is_integer(N) -> "(" ++ to_form(B, Names) ++ ")^(" ++ integer_to_list(N) ++ ")";
to_form(Other, _) -> throw({unsupported, kind(Other)}).

kind(T) when is_tuple(T) -> element(1, T);
kind(T) -> T.

%%% FORM's output back: a sum of terms, each [sign] [coefficient] factors joined by '*', a factor `sK`, `sK^E`,
%%% `i_` or an integer or fraction; no brackets (the scripts never ask for them)

from_form(Text, Back) ->
    Flat = << <<C>> || <<C>> <= Text, C =/= $\s, C =/= $\n, C =/= $\r, C =/= $\t, C =/= $; >>,
    case Flat of
        <<>> -> 0;
        <<"0">> -> 0;
        _ -> {add, [term(T, Back) || T <- split_terms(Flat)]}
    end.

%% the terms with their signs: split before every '+' or '-' that is not an exponent's sign
split_terms(B) -> split_terms(B, <<>>, []).

split_terms(<<>>, Cur, Acc) -> lists:reverse([Cur | Acc]);
split_terms(<<$^, $-, Rest/binary>>, Cur, Acc) -> split_terms(Rest, <<Cur/binary, "^-">>, Acc);
split_terms(<<S, Rest/binary>>, Cur, Acc) when (S =:= $+ orelse S =:= $-), Cur =/= <<>> ->
    split_terms(Rest, <<S>>, [Cur | Acc]);
split_terms(<<C, Rest/binary>>, Cur, Acc) -> split_terms(Rest, <<Cur/binary, C>>, Acc).

term(<<$-, T/binary>>, Back) -> {mul, [-1, product(T, Back)]};
term(<<$+, T/binary>>, Back) -> product(T, Back);
term(T, Back) -> product(T, Back).

product(T, Back) -> {mul, [factor(F, Back) || F <- binary:split(T, <<"*">>, [global])]}.

factor(<<"i_">>, _) -> {complex, 0, 1};
factor(F, Back) ->
    case binary:split(F, <<"^">>) of
        [Base, Exp] -> {pow, base(Base, Back), int(Exp)};
        [Base] -> base(Base, Back)
    end.

base(<<C, _/binary>> = B, _Back) when C >= $0, C =< $9 ->
    case binary:split(B, <<"/">>) of
        [N, D] -> {rat, int(N), int(D)};
        [N] -> int(N)
    end;
base(<<"i_">>, _) -> {complex, 0, 1};
base(B, Back) ->
    case maps:find(B, Back) of
        {ok, S} -> {sym, S};
        error -> throw({form_output, B})
    end.

int(B) ->
    try binary_to_integer(B)
    catch error:badarg -> throw({form_output, B})
    end.
