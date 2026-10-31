%% Evaluation of a canonical value under a binding (SY 5.2): complex doubles term by term in canonical key order,
%% principal branch for powers; or, with exact => true, exact Gaussian rationals for the algebraic fragment.
-module(pq_sym_eval).
-export([eval/3, cpow/2, rpow/2, free_symbols/1, check_binding/2, value/1]).

-define(ZERO, {0.0, 0.0}).
-define(ONEC, {1.0, 0.0}).

eval(V, Binding, Opts) ->
    Syms = free_symbols(V),
    case [S || S <- Syms, not maps:is_key(S, Binding)] of
        [] -> ok;
        Missing -> throw({unbound_symbols, Missing})
    end,
    check_binding(maps:with(Syms, Binding), maps:get(domains, Opts, #{})),
    case maps:get(exact, Opts, false) of
        true -> exact_value(V, Binding);
        false -> float_value(V, Binding)
    end.

%%% Symbols

free_symbols({scalar, P}) -> lists:usort(poly_symbols(P));
free_symbols({identity, _}) -> [];
free_symbols({mat, _, _, Rows}) -> lists:usort(lists:append([poly_symbols(E) || Row <- Rows, E <- Row])).

poly_symbols(P) -> lists:append([mono_symbols(M) || M <- maps:keys(P)]).

mono_symbols({Alg, Exp}) ->
    lists:append([atom_symbols(K) || {K, _} <- Alg]) ++ lists:append([mono_symbols(B) || {B, _} <- Exp]).

atom_symbols({sym, X}) -> [X];
atom_symbols({conj, X}) -> [X];
atom_symbols({base, B}) -> poly_symbols(B);
atom_symbols({cbase, B}) -> poly_symbols(B);
atom_symbols(_) -> [].

%% domain_violation unless every value fits its symbol's domain (real symbols by default)
check_binding(Binding, Domains) ->
    maps:foreach(
      fun(X, V) ->
              D = maps:get(X, Domains, real),
              fits(value(V), D) orelse throw({domain_violation, {X, D}})
      end, Binding).

fits(_, complex) -> true;
fits({_, Im}, real) -> Im == 0;
fits({Re, Im}, nonnegative) -> Im == 0 andalso Re >= 0;
fits({Re, Im}, positive) -> Im == 0 andalso Re > 0.

%% a binding's value as a complex pair
value(N) when is_integer(N) -> {float(N), 0.0};
value(X) when is_float(X) -> {X, 0.0};
value({rat, N, D}) when is_integer(N), is_integer(D), D > 0 -> {N / D, 0.0};
value({complex, Re, Im}) -> {element(1, value(Re)), element(1, value(Im))};
value(Other) -> throw({bad_expr, {binding, Other}}).

%%% Floats

float_value({scalar, P}, B) -> result(poly(P, B));
float_value({identity, N}, _B) ->
    float_value({mat, N, N, [[case I =:= J of true -> pq_sym_poly:one(); false -> #{} end
                              || J <- lists:seq(1, N)] || I <- lists:seq(1, N)]}, #{});
float_value({mat, _, _, Rows}, B) ->
    {ok, Packed} = pq_core:pack({cmat, [[begin {Re, Im} = poly(E, B), {complex, Re, Im} end || E <- Row]
                                        || Row <- Rows]}),
    Packed.

result({Re, Im}) when Im == 0 -> Re + 0.0;
result({Re, Im}) -> {complex, Re + 0.0, Im + 0.0}.

poly(P, B) ->
    lists:foldl(fun({M, C}, Acc) -> add(Acc, mul(pq_sym_exact:to_complex(C), mono(M, B))) end,
                ?ZERO, pq_sym_poly:sorted(P)).

mono({Alg, Exp}, B) ->
    V1 = lists:foldl(fun({K, E}, Acc) -> mul(Acc, atom(K, E, B)) end, ?ONEC, Alg),
    lists:foldl(fun({Base, C}, Acc) -> mul(Acc, cexp(mul(pq_sym_exact:to_complex(C), mono(Base, B)))) end,
                V1, Exp).

atom({sym, X}, E, B) -> rpow(value(maps:get(X, B)), E);
atom({conj, X}, E, B) -> {Re, Im} = value(maps:get(X, B)), rpow({Re, -Im}, E);
atom(pi, E, _B) -> rpow({math:pi(), 0.0}, E);
atom({rad, P}, {1, 1}, _B) -> {math:sqrt(P), 0.0};
atom({base, P}, E, B) -> rpow(poly(P, B), E);
atom({cbase, P}, E, B) -> {Re, Im} = rpow(poly(P, B), E), {Re, -Im}.

%% z^(N/D): integer powers by repeated squaring; a positive real base by math:pow; otherwise exp(w Log z)
rpow(Z, {N, 1}) -> ipow(Z, N);
rpow(Z, {N, 2}) -> ipow(csqrt(Z), N);
rpow({Re, Im}, {N, D}) when Im == 0, Re > 0 -> {math:pow(Re, N / D), 0.0};
rpow(Z, {N, D}) -> cpow(Z, {N / D, 0.0}).

%% the principal square root, by the stable formula: sqrt(-1) is exactly i
csqrt({Re, Im}) ->
    R = math:sqrt(Re * Re + Im * Im),
    A = math:sqrt((R + Re) / 2),
    B = math:sqrt((R - Re) / 2),
    case Im + 0.0 < 0 of
        true -> {A, -B};
        false -> {A, B}
    end.

ipow(_Z, 0) -> ?ONEC;
ipow(Z, N) when N < 0 -> inv(ipow(Z, -N));
ipow(Z, 1) -> Z;
ipow(Z, N) ->
    H = ipow(Z, N div 2),
    H2 = mul(H, H),
    case N rem 2 of 0 -> H2; 1 -> mul(H2, Z) end.

cpow({Re, Im}, {WRe, _}) when Re == 0, Im == 0 ->
    case WRe > 0 of
        true -> ?ZERO;
        false -> throw({non_finite, division_by_zero})
    end;
cpow(Z, W) -> cexp(mul(W, clog(Z))).

%% principal branch, arg in (-pi, pi]: a zero imaginary part counts as +0.0
clog({Re, Im}) -> {math:log(math:sqrt(Re * Re + Im * Im)), math:atan2(Im + 0.0, Re)}.

cexp({Re, Im}) ->
    M = math:exp(Re),
    {M * math:cos(Im), M * math:sin(Im)}.

add({A, B}, {C, D}) -> {A + C, B + D}.
mul({A, B}, {C, D}) -> {A * C - B * D, A * D + B * C}.

inv({Re, Im}) ->
    N = Re * Re + Im * Im,
    N == 0 andalso throw({non_finite, division_by_zero}),
    {Re / N, -Im / N}.

%%% Exact (the test oracle of SY 9.2): integer powers of symbols, exact bindings, no pi, radicals or exponentials

exact_value({scalar, P}, B) -> exact_result(exact_poly(P, B));
exact_value(_, _B) -> throw({symbolic_unsupported, exact}).

exact_poly(P, B) ->
    lists:foldl(fun({M, C}, Acc) -> pq_sym_exact:cadd(Acc, pq_sym_exact:cmul(C, exact_mono(M, B))) end,
                pq_sym_exact:c(0), pq_sym_poly:sorted(P)).

exact_mono({Alg, []}, B) ->
    lists:foldl(fun({K, E}, Acc) -> pq_sym_exact:cmul(Acc, exact_atom(K, E, B)) end, pq_sym_exact:c(1), Alg);
exact_mono(_, _B) -> throw({symbolic_unsupported, exact}).

exact_atom({sym, X}, {N, 1}, B) -> pq_sym_exact:cpow_int(exact(maps:get(X, B)), N);
exact_atom({conj, X}, {N, 1}, B) -> pq_sym_exact:cpow_int(pq_sym_exact:cconj(exact(maps:get(X, B))), N);
exact_atom({base, P}, {N, 1}, B) -> pq_sym_exact:cpow_int(exact_poly(P, B), N);
exact_atom(_, _, _B) -> throw({symbolic_unsupported, exact}).

exact(N) when is_integer(N) -> pq_sym_exact:c(N);
exact({rat, N, D}) -> pq_sym_exact:c(pq_sym_exact:rat(N, D));
exact({complex, Re, Im}) -> pq_sym_exact:cq(exact_part(Re), exact_part(Im));
exact(_) -> throw({symbolic_unsupported, exact}).

exact_part(N) when is_integer(N) -> {N, 1};
exact_part({rat, N, D}) -> pq_sym_exact:rat(N, D);
exact_part(_) -> throw({symbolic_unsupported, exact}).

exact_result(C) -> pq_sym_norm:coeff_to_ast(C).
