%% Recursive-descent parser for the parameter grammar of FS B.3, with amendment B-9 (exponents ^(p/q)). At most
%% 65,536 bytes and nesting depth 200. Errors: {sym_parse, Position} (1-based byte offset). Atom safety: symbol
%% names become atoms only after every name has been checked (allowed_symbols, or at most max_symbols names).
-module(pq_sym_parse).
-export([parse/2]).

-define(MAX_BYTES, 65536).
-define(MAX_DEPTH, 200).
-define(MAX_NAME, 64).

parse(Bin, Opts) when is_binary(Bin) ->
    try
        byte_size(Bin) > ?MAX_BYTES andalso throw({sym_parse, ?MAX_BYTES + 1}),
        Tokens = tokens(Bin, 1, []),
        {E, Rest} = expr(Tokens, 0),
        case Rest of
            [{eof, _, _}] -> ok;
            [{_, _, Pos} | _] -> throw({sym_parse, Pos})
        end,
        {ok, atoms(E, Opts)}
    catch
        throw:{sym_parse, _} = R -> {error, R};
        throw:{symbolic_limit, _} = R -> {error, R}
    end;
parse(_, _) -> {error, {sym_parse, 1}}.

%%% Tokens: {num, Value, Pos} | {id, Binary, Pos} | {op, Char, Pos} | {eof, none, Pos}

tokens(<<>>, Pos, Acc) -> lists:reverse([{eof, none, Pos} | Acc]);
tokens(<<C, Rest/binary>>, Pos, Acc) when C =:= $\s; C =:= $\t; C =:= $\n; C =:= $\r -> tokens(Rest, Pos + 1, Acc);
tokens(<<C, _/binary>> = B, Pos, Acc) when C >= $0, C =< $9 ->
    {Tok, Len} = number(B),
    tokens(binary:part(B, Len, byte_size(B) - Len), Pos + Len, [{num, Tok, Pos} | Acc]);
tokens(<<C, _/binary>> = B, Pos, Acc) when (C >= $a andalso C =< $z) orelse (C >= $A andalso C =< $Z) ->
    Len = ident_len(B, 0),
    Len > ?MAX_NAME andalso throw({sym_parse, Pos}),
    <<Name:Len/binary, Rest/binary>> = B,
    tokens(Rest, Pos + Len, [{id, Name, Pos} | Acc]);
tokens(<<C, Rest/binary>>, Pos, Acc) when C =:= $+; C =:= $-; C =:= $*; C =:= $/; C =:= $^; C =:= $(; C =:= $) ->
    tokens(Rest, Pos + 1, [{op, C, Pos} | Acc]);
tokens(_, Pos, _) -> throw({sym_parse, Pos}).

ident_len(<<C, Rest/binary>>, N) when (C >= $a andalso C =< $z) orelse (C >= $A andalso C =< $Z)
                                      orelse (C >= $0 andalso C =< $9) orelse C =:= $_ -> ident_len(Rest, N + 1);
ident_len(_, N) -> N.

%% the longest match of digits [. digits] [(e|E) [+|-] digits]: an integer, or a float
number(B) ->
    I = digits(B, 0),
    <<Int:I/binary, R1/binary>> = B,
    {Frac, R2} = case R1 of
                     <<$., C, _/binary>> when C >= $0, C =< $9 ->
                         F = digits(binary:part(R1, 1, byte_size(R1) - 1), 0),
                         {binary:part(R1, 1, F), binary:part(R1, 1 + F, byte_size(R1) - 1 - F)};
                     _ -> {none, R1}
                 end,
    Exp = case R2 of
              <<E, S, C2, _/binary>> when (E =:= $e orelse E =:= $E), (S =:= $+ orelse S =:= $-), C2 >= $0, C2 =< $9 ->
                  D = digits(binary:part(R2, 2, byte_size(R2) - 2), 0),
                  {<<S>>, binary:part(R2, 2, D), 2 + D};
              <<E, C2, _/binary>> when (E =:= $e orelse E =:= $E), C2 >= $0, C2 =< $9 ->
                  D = digits(binary:part(R2, 1, byte_size(R2) - 1), 0),
                  {<<>>, binary:part(R2, 1, D), 1 + D};
              _ -> none
          end,
    FracLen = case Frac of none -> 0; _ -> 1 + byte_size(Frac) end,
    case {Frac, Exp} of
        {none, none} -> {binary_to_integer(Int), I};
        _ ->
            F1 = case Frac of none -> <<"0">>; _ -> Frac end,
            {ES, ED, EL} = case Exp of none -> {<<>>, <<"0">>, 0}; _ -> Exp end,
            {binary_to_float(<<Int/binary, ".", F1/binary, "e", ES/binary, ED/binary>>), I + FracLen + EL}
    end.

digits(<<C, Rest/binary>>, N) when C >= $0, C =< $9 -> digits(Rest, N + 1);
digits(_, N) -> N.

%%% Grammar

depth(D, [{_, _, Pos} | _]) when D > ?MAX_DEPTH -> throw({sym_parse, Pos});
depth(_, _) -> ok.

%% expr = term { ("+" | "-") term }
expr(Ts, D) ->
    depth(D, Ts),
    {T, Rest} = term(Ts, D + 1),
    expr_tail(Rest, D, [T]).

expr_tail([{op, $+, _} | Ts], D, Acc) -> {T, R} = term(Ts, D + 1), expr_tail(R, D, [T | Acc]);
expr_tail([{op, $-, _} | Ts], D, Acc) -> {T, R} = term(Ts, D + 1), expr_tail(R, D, [{mul, [-1, T]} | Acc]);
expr_tail(Ts, _D, [T]) -> {T, Ts};
expr_tail(Ts, _D, Acc) -> {{add, lists:reverse(Acc)}, Ts}.

%% term = unary { ("*" | "/") unary }
term(Ts, D) ->
    {U, Rest} = unary(Ts, D + 1),
    term_tail(Rest, D, [U]).

term_tail([{op, $*, _} | Ts], D, Acc) -> {U, R} = unary(Ts, D + 1), term_tail(R, D, [U | Acc]);
term_tail([{op, $/, _} | Ts], D, Acc) -> {U, R} = unary(Ts, D + 1), term_tail(R, D, [{pow, U, -1} | Acc]);
term_tail(Ts, _D, [U]) -> {U, Ts};
term_tail(Ts, _D, Acc) -> {{mul, lists:reverse(Acc)}, Ts}.

%% unary = "-" unary | power
unary([{op, $-, _} | Ts], D) ->
    depth(D, Ts),
    {U, R} = unary(Ts, D + 1),
    {{mul, [-1, U]}, R};
unary(Ts, D) -> power(Ts, D).

%% power = primary [ "^" exponent ];  exponent = ["-"] integer | "(" ["-"] integer [ "/" integer ] ")"
power(Ts, D) ->
    {P, Rest} = primary(Ts, D + 1),
    case Rest of
        [{op, $^, _} | R1] -> {E, R2} = exponent(R1), {{pow, P, E}, R2};
        _ -> {P, Rest}
    end.

exponent([{op, $-, _}, {num, N, _} | R]) when is_integer(N) -> {-N, R};
exponent([{num, N, _} | R]) when is_integer(N) -> {N, R};
exponent([{op, $(, _} | R]) ->
    {Sign, R1} = case R of [{op, $-, _} | X] -> {-1, X}; _ -> {1, R} end,
    case R1 of
        [{num, N, _}, {op, $/, _}, {num, M, Pos}, {op, $), _} | R2] when is_integer(N), is_integer(M) ->
            M > 0 orelse throw({sym_parse, Pos}),
            G = pq_sym_exact:gcd(N, M),
            case M div G of
                1 -> {Sign * N div G, R2};
                Den -> {{rat, Sign * N div G, Den}, R2}
            end;
        [{num, N, _}, {op, $), _} | R2] when is_integer(N) -> {Sign * N, R2};
        [{_, _, Pos} | _] -> throw({sym_parse, Pos})
    end;
exponent([{_, _, Pos} | _]) -> throw({sym_parse, Pos}).

primary([{num, N, _} | R], _D) -> {N, R};
primary([{id, <<"pi">>, _} | R], _D) -> {{const, pi}, R};
primary([{id, <<"e">>, _} | R], _D) -> {{const, e}, R};
primary([{id, <<"i">>, _} | R], _D) -> {{complex, 0, 1}, R};
primary([{id, F, _}, {op, $(, _} | R], D) when F =:= <<"sin">>; F =:= <<"cos">>; F =:= <<"exp">>;
                                               F =:= <<"conj">>; F =:= <<"sqrt">> ->
    {E, R1} = expr(R, D + 1),
    R2 = close(R1),
    case F of
        <<"sqrt">> -> {{pow, E, {rat, 1, 2}}, R2};
        _ -> {{binary_to_atom(F), E}, R2}
    end;
primary([{id, F, Pos} | _], _D) when F =:= <<"sin">>; F =:= <<"cos">>; F =:= <<"exp">>; F =:= <<"conj">>;
                                     F =:= <<"sqrt">> -> throw({sym_parse, Pos + byte_size(F)});
primary([{id, Name, Pos} | R], _D) -> {{sym_name, Name, Pos}, R};
primary([{op, $(, _} | R], D) ->
    {E, R1} = expr(R, D + 1),
    {E, close(R1)};
primary([{_, _, Pos} | _], _D) -> throw({sym_parse, Pos}).

close([{op, $), _} | R]) -> R;
close([{_, _, Pos} | _]) -> throw({sym_parse, Pos}).

%%% Atom safety

atoms(E, Opts) ->
    Names = lists:ukeysort(1, names(E, [])),
    Max = maps:get(max_symbols, Opts, 1024),
    length(Names) > Max andalso throw({symbolic_limit, #{reason => max_symbols, max => Max}}),
    Map = case maps:find(allowed_symbols, Opts) of
              {ok, Allowed} ->
                  maps:from_list([{N, allowed(N, Pos, Allowed)} || {N, Pos} <- Names]);
              error ->
                  maps:from_list([{N, binary_to_atom(N)} || {N, _} <- Names])
          end,
    replace(E, Map).

allowed(Name, Pos, Allowed) ->
    try binary_to_existing_atom(Name) of
        A ->
            lists:member(A, Allowed) orelse throw({sym_parse, Pos}),
            A
    catch
        error:badarg -> throw({sym_parse, Pos})
    end.

names({sym_name, N, Pos}, Acc) -> [{N, Pos} | Acc];
names(T, Acc) when is_tuple(T) -> names(tuple_to_list(T), Acc);
names(L, Acc) when is_list(L) -> lists:foldl(fun names/2, Acc, L);
names(_, Acc) -> Acc.

replace({sym_name, N, _}, Map) -> {sym, maps:get(N, Map)};
replace({add, L}, Map) -> {add, [replace(X, Map) || X <- L]};
replace({mul, L}, Map) -> {mul, [replace(X, Map) || X <- L]};
replace({pow, B, E}, Map) -> {pow, replace(B, Map), E};
replace({F, A}, Map) when F =:= sin; F =:= cos; F =:= exp; F =:= conj -> {F, replace(A, Map)};
replace(X, _Map) -> X.
