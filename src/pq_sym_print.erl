%% Printing (SY 5.6, FS B.3 with amendment B-9): ASCII that pq_sym_parse reads back, and LaTeX. Floats print in
%% their shortest round-trip form. {trig, fold} is not implemented yet: exponentials print as they are.
-module(pq_sym_print).
-export([ascii/1, latex/1, in_subset/1]).

-define(ADD, 1).
-define(MUL, 2).
-define(NEG, 3).
-define(POW, 4).
-define(ATOM, 5).

ascii(E) -> iolist_to_binary(element(2, p(E))).

%% {Precedence, IoList}
p(N) when is_integer(N), N >= 0 -> {?ATOM, integer_to_list(N)};
p(N) when is_integer(N) -> {?NEG, ["-", integer_to_list(-N)]};
p(X) when is_float(X), X >= 0 -> {?ATOM, float_to_list(X, [short])};
p(X) when is_float(X) -> {?NEG, ["-", float_to_list(-X, [short])]};
p({rat, N, D}) when N >= 0 -> {?MUL, [integer_to_list(N), "/", integer_to_list(D)]};
p({rat, N, D}) -> {?MUL, ["-", integer_to_list(-N), "/", integer_to_list(D)]};
p({complex, Re, Im}) ->
    case is_zero(Re) of
        true -> imag(Im);
        false ->
            {_, R} = p(Re),
            case negative(Im) of
                true -> {?ATOM, ["(", R, " - ", element(2, imag(negate(Im))), ")"]};
                false -> {?ATOM, ["(", R, " + ", element(2, imag(Im)), ")"]}
            end
    end;
p({const, pi}) -> {?ATOM, "pi"};
p({const, e}) -> {?ATOM, "e"};
p({const, i}) -> {?ATOM, "i"};
p({sym, X}) -> {?ATOM, atom_to_list(X)};
p({add, []}) -> {?ATOM, "0"};
p({add, [T | Ts]}) ->
    {?ADD, [wrap(T, ?ADD) | [case negative_term(X) of
                                 true -> [" - ", wrap(negate_term(X), ?MUL)];
                                 false -> [" + ", wrap(X, ?MUL)]
                             end || X <- Ts]]};
p({mul, []}) -> {?ATOM, "1"};
p({mul, [C | Fs]}) when Fs =/= [], (C =:= 1 orelse C =:= 1.0) -> p({mul, Fs});
p({mul, [C | Fs]}) when Fs =/= [], (C =:= -1 orelse C =:= -1.0) -> {?NEG, ["-", wrap({mul, Fs}, ?NEG)]};
p({mul, [F]}) -> p(F);
p({mul, [{complex, Z, Im} | Fs]}) when Fs =/= [], (Z == 0 orelse Z =:= {rat, 0, 1}) -> p({mul, [Im, {const, i} | Fs]});
p({mul, [F | Fs]}) ->
    First = case F of
                {rat, _, _} -> element(2, p(F));     % (n/d)*x reads back the same: * and / associate left
                _ -> wrap(F, ?NEG)
            end,
    {?MUL, [First | [["*", wrap(X, ?POW)] || X <- Fs]]};
p({pow, B, {rat, 1, 2}}) -> {?ATOM, ["sqrt(", element(2, p(B)), ")"]};
p({pow, B, E}) -> {?POW, [wrap(B, ?ATOM), "^", exponent(E)]};
p({F, A}) when F =:= exp; F =:= sin; F =:= cos; F =:= conj ->
    {?ATOM, [atom_to_list(F), "(", element(2, p(A)), ")"]};
p({matrix, Rows}) ->
    {?ATOM, ["[", lists:join(", ", [["[", lists:join(", ", [element(2, p(E)) || E <- R]), "]"] || R <- Rows]), "]"]};
p({identity, N}) -> {?ATOM, ["identity(", integer_to_list(N), ")"]};
p({Op, A, B}) when Op =:= matmul; Op =:= kron; Op =:= smul ->
    {?ATOM, [atom_to_list(Op), "(", element(2, p(A)), ", ", element(2, p(B)), ")"]};
p({Op, A}) when Op =:= transpose; Op =:= dagger -> {?ATOM, [atom_to_list(Op), "(", element(2, p(A)), ")"]};
p({madd, Ms}) -> {?ATOM, ["madd(", lists:join(", ", [element(2, p(M)) || M <- Ms]), ")"]};
p(Other) -> {?ATOM, io_lib:format("~p", [Other])}.

wrap(E, Min) ->
    case p(E) of
        {P, S} when P >= Min -> S;
        {_, S} -> ["(", S, ")"]
    end.

exponent(N) when is_integer(N), N >= 0 -> integer_to_list(N);
exponent(N) when is_integer(N) -> ["-", integer_to_list(-N)];
exponent({rat, N, D}) -> ["(", integer_to_list(N), "/", integer_to_list(D), ")"];
exponent(X) when is_float(X) -> ["(", float_to_list(X, [short]), ")"].

imag(1) -> {?ATOM, "i"};
imag(-1) -> {?NEG, "-i"};
imag(Im) -> p({mul, [Im, {const, i}]}).

is_zero(X) -> X == 0 orelse X =:= {rat, 0, 1}.

negative(N) when is_number(N) -> N < 0;
negative({rat, N, _}) -> N < 0;
negative(_) -> false.

negate(N) when is_number(N) -> -N;
negate({rat, N, D}) -> {rat, -N, D}.

negative_term(T) when is_number(T) -> T < 0;
negative_term({rat, N, _}) -> N < 0;
negative_term({mul, [C | _]}) -> negative(C);
negative_term(_) -> false.

negate_term({mul, [C | Fs]}) -> {mul, [negate(C) | Fs]};
negate_term(T) -> negate(T).

%% Whether the expression is inside the parameter subset of B.3 (with B-9): no matrices, exponents integer or
%% rational, no float exponent.
in_subset(N) when is_number(N) -> true;
in_subset({rat, _, _}) -> true;
in_subset({complex, _, _}) -> true;
in_subset({const, _}) -> true;
in_subset({sym, _}) -> true;
in_subset({add, Ts}) -> lists:all(fun in_subset/1, Ts);
in_subset({mul, Ts}) -> lists:all(fun in_subset/1, Ts);
in_subset({pow, B, E}) -> (is_integer(E) orelse element(1, E) =:= rat) andalso in_subset(B);
in_subset({F, A}) when F =:= exp; F =:= sin; F =:= cos; F =:= conj -> in_subset(A);
in_subset(_) -> false.

%%% LaTeX

latex(E) -> iolist_to_binary(element(2, l(E))).

l(N) when is_integer(N), N >= 0 -> {?ATOM, integer_to_list(N)};
l(N) when is_integer(N) -> {?NEG, ["-", integer_to_list(-N)]};
l(X) when is_float(X), X >= 0 -> {?ATOM, float_to_list(X, [short])};
l(X) when is_float(X) -> {?NEG, ["-", float_to_list(-X, [short])]};
l({rat, N, D}) when N >= 0 -> {?ATOM, ["\\frac{", integer_to_list(N), "}{", integer_to_list(D), "}"]};
l({rat, N, D}) -> {?NEG, ["-\\frac{", integer_to_list(-N), "}{", integer_to_list(D), "}"]};
l({complex, Re, Im}) ->
    case is_zero(Re) of
        true -> {?MUL, [element(2, l(Im)), " i"]};
        false -> {?ATOM, ["\\left(", element(2, l(Re)), " + ", element(2, l(Im)), " i\\right)"]}
    end;
l({const, pi}) -> {?ATOM, "\\pi"};
l({const, C}) -> {?ATOM, atom_to_list(C)};
l({sym, X}) -> {?ATOM, atom_to_list(X)};
l({add, [T | Ts]}) ->
    {?ADD, [lwrap(T, ?ADD) | [case negative_term(X) of
                                  true -> [" - ", lwrap(negate_term(X), ?MUL)];
                                  false -> [" + ", lwrap(X, ?MUL)]
                              end || X <- Ts]]};
l({mul, [C | Fs]}) when Fs =/= [], (C =:= 1 orelse C =:= 1.0) -> l({mul, Fs});
l({mul, [C | Fs]}) when Fs =/= [], (C =:= -1 orelse C =:= -1.0) -> {?NEG, ["-", lwrap({mul, Fs}, ?NEG)]};
l({mul, [F]}) -> l(F);
l({mul, Fs}) -> {?MUL, lists:join(" \\cdot ", [lwrap(F, ?NEG) || F <- Fs])};
l({pow, B, {rat, 1, 2}}) -> {?ATOM, ["\\sqrt{", element(2, l(B)), "}"]};
l({pow, B, E}) -> {?POW, [lwrap(B, ?ATOM), "^{", element(2, l(E)), "}"]};
l({exp, A}) -> {?POW, ["e^{", element(2, l(A)), "}"]};
l({conj, A}) -> {?ATOM, ["\\overline{", element(2, l(A)), "}"]};
l({F, A}) when F =:= sin; F =:= cos -> {?ATOM, ["\\", atom_to_list(F), "\\left(", element(2, l(A)), "\\right)"]};
l({matrix, Rows}) ->
    {?ATOM, ["\\begin{pmatrix}", lists:join(" \\\\ ", [lists:join(" & ", [element(2, l(E)) || E <- R]) || R <- Rows]),
             "\\end{pmatrix}"]};
l(Other) -> {?ATOM, ["\\text{", element(2, p(Other)), "}"]}.

lwrap(E, Min) ->
    case l(E) of
        {P, S} when P >= Min -> S;
        {_, S} -> ["\\left(", S, "\\right)"]
    end.
