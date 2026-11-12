%% Laurent series in epsilon (d = 4 - 2 epsilon) with coefficients polynomial in zeta values, for the comparison with
%% FORM's Forcer (comparison/forcer/README.md): the expansion of a reduction's coefficients times its master
%% integrals, in Forcer's G-scheme and in MS-bar (B. Ruijl, T. Ueda, J.A.M. Vermaseren,
%% https://arxiv.org/abs/1704.06650, sections "One-loop sub integrals" and "The master integrals").
%%
%% A series is {Coeffs, Prec}: Coeffs #{Order => Poly}, known exactly below the order Prec (infinity: exact). A Poly
%% is #{Monomial => {N, D}}, a monomial a sorted list of {Symbol, Exponent} (z3, z5, z5z3, ...); [] is 1.
-module(epsexp).
-export([eval/3, gc/3, msbar/1, series_of_form/2, truncate/2, to_string/1, sub/2, equal/2, is_zero/1, leading/1,
         sadd/2, smul/2, spow/3, drop/2, sinv/2, exact/1, scale/2, coeffs/1]).

-define(INF, 1000000).

%%% Polynomials in the zeta values

pconst({0, _}) -> #{};
pconst(Q) -> #{[] => Q}.
psym(S) -> #{[{S, 1}] => {1, 1}}.
padd(A, B) ->
    maps:fold(fun(K, V, Acc) ->
                      case qadd(maps:get(K, Acc, {0, 1}), V) of
                          {0, _} -> maps:remove(K, Acc);
                          S -> Acc#{K => S}
                      end
              end, A, B).
pmul(A, B) ->
    lists:foldl(fun(T, Acc) -> padd(Acc, T) end, #{}, [#{mmul(Ka, Kb) => qmul(Va, Vb)} || Ka := Va <- A, Kb := Vb <- B]).
pscale(_, {0, _}) -> #{};
pscale(P, Q) -> maps:map(fun(_, V) -> qmul(V, Q) end, P).
mmul(A, B) ->
    M = lists:foldl(fun({S, E}, Acc) -> Acc#{S => maps:get(S, Acc, 0) + E} end, maps:from_list(A), B),
    lists:sort([{S, E} || S := E <- M, E =/= 0]).

%%% Series

exact(Coeffs) -> {maps:filter(fun(_, P) -> map_size(P) > 0 end, Coeffs), ?INF}.
sconst(Q) -> exact(#{0 => pconst(Q)}).
seps() -> exact(#{1 => pconst({1, 1})}).

valuation({C, Prec}) ->
    case maps:keys(C) of
        [] -> Prec;
        Ks -> lists:min(Ks)
    end.

sadd({A, Pa}, {B, Pb}) ->
    P = min(Pa, Pb),
    {maps:filter(fun(K, V) -> K < P andalso map_size(V) > 0 end,
                 maps:fold(fun(K, V, Acc) -> Acc#{K => padd(maps:get(K, Acc, #{}), V)} end, A, B)), P}.

smul(X = {A, Pa}, Y = {B, Pb}) ->
    P = min(sat(Pa + valuation(Y)), sat(Pb + valuation(X))),
    C = lists:foldl(fun({K, V}, Acc) ->
                            case K < P of
                                true -> Acc#{K => padd(maps:get(K, Acc, #{}), V)};
                                false -> Acc
                            end
                    end, #{}, [{Ka + Kb, pmul(Va, Vb)} || Ka := Va <- A, Kb := Vb <- B]),
    {maps:filter(fun(_, V) -> map_size(V) > 0 end, C), P}.

sat(X) when X > ?INF -> ?INF;
sat(X) -> X.

%% 1 / S to N orders beyond its leading one (the leading coefficient a rational number)
sinv({C, ?INF}, _) when map_size(C) =:= 1 ->
    %% an exact monomial a eps^k: exactly 1/a eps^-k
    [{K, #{[] := A} = P}] = maps:to_list(C),
    1 = map_size(P),
    exact(#{-K => pconst(qinv(A))});
sinv(S = {C, Prec}, N) ->
    V = valuation(S),
    case maps:get(V, C, undefined) of
        #{[] := Lead} = L when map_size(L) =:= 1 ->
            Inv = qinv(Lead),
            P = case Prec of ?INF -> -V + N; _ -> min(-V + N, Prec - 2 * V) end,
            %% b_0 = 1/a_0, b_k = -(1/a_0) sum_{j=1..k} a_j b_{k-j}
            Bs = lists:foldl(fun(K, Acc) ->
                                     Sum = lists:foldl(fun(J, S2) ->
                                                               padd(S2, pmul(maps:get(V + J, C, #{}), maps:get(K - J, Acc)))
                                                       end, #{}, lists:seq(1, K)),
                                     Acc#{K => pscale(Sum, qneg(Inv))}
                             end, #{0 => pconst(Inv)}, lists:seq(1, max(0, P + V - 1))),
            {maps:filter(fun(K, X) -> K < P andalso map_size(X) > 0 end, maps:from_list([{K - V, B} || K := B <- Bs])), P};
        _ -> error({no_inverse, S})
    end.

spow(S, K, _N) when K >= 0 -> lists:foldl(fun(_, Acc) -> smul(Acc, S) end, sconst({1, 1}), lists:seq(1, K));
spow(S, K, N) -> spow(sinv(S, N), -K, N).

%% exp of a series with no constant or negative orders, to order N
sexp(S = {C, _}, N) ->
    0 = length([K || K := _ <- C, K =< 0]),
    Trunc = truncate(S, N + 1),
    {Sum, _} = lists:foldl(fun(K, {Acc, Term}) ->
                                   T2 = scale(smul(Term, Trunc), {1, K}),
                                   {sadd(Acc, T2), T2}
                           end, {sconst({1, 1}), sconst({1, 1})}, lists:seq(1, N)),
    truncate(Sum, N + 1).

scale({C, P}, Q) -> {maps:filter(fun(_, V) -> map_size(V) > 0 end, maps:map(fun(_, V) -> pscale(V, Q) end, C)), P}.

truncate({C, P}, N) -> {maps:filter(fun(K, _) -> K < N end, C), min(P, N)}.
%% the terms of order N and above set to zero (Forcer's msbarexpand drops them): the precision is kept
drop({C, P}, N) -> {maps:filter(fun(K, _) -> K < N end, C), P}.

sub(A, B) -> sadd(A, scale(B, {-1, 1})).
is_zero({C, _}) -> map_size(C) =:= 0.
equal(A, B) -> is_zero(sub(A, B)).
leading(S) -> valuation(S).

%%% Gamma functions

%% sum_{k>=2} (-1)^k zeta_k x^k / k for x = c epsilon, to order N (log Gamma(1 + c eps) without the gamma_E term,
%% which cancels in every ratio used here); zeta_2 is dropped, as Forcer drops it (it never appears in massless
%% propagators' physical results: https://arxiv.org/abs/1704.06650)
lngamma1(C, N) ->
    exact(maps:from_list([{K, pscale(psym(zeta(K)), qnorm({sign(K) * pow(C, K), K}))} || K <- lists:seq(3, N)]
                         ++ [{2, #{}}])).

sign(K) -> case K rem 2 of 0 -> 1; 1 -> -1 end.
pow(_, 0) -> 1;
pow(C, K) -> C * pow(C, K - 1).
zeta(K) -> list_to_atom("z" ++ integer_to_list(K)).

%% Forcer's GschemeConstants(a, b) = (1 + a + b) G(1 + a eps, 1 + b eps) / (eps G(1, 1)) (forcer-aux.h; the
%% G-scheme sets GschemeConstants(0, 0) = 1/eps), G the one-loop function
%%   G(al, be) = Gamma(al + be - D/2) Gamma(D/2 - al) Gamma(D/2 - be) / (Gamma(al) Gamma(be) Gamma(D - al - be))
%% = (1 - 2 eps) / (eps (1 - (2 + a + b) eps)) exp(L(a, b) - L(0, 0)), L the sum of the log Gammas
-spec gc(integer(), integer(), pos_integer()) -> tuple().
gc(A, B, N) ->
    L = fun(X, Y) ->
                lists:foldl(fun({S, C}, Acc) -> sadd(Acc, scale(lngamma1(C, N + 2), {S, 1})) end, exact(#{}),
                            [{1, X + Y + 1}, {1, -(1 + X)}, {1, -(1 + Y)}, {-1, X}, {-1, Y}, {-1, -(2 + X + Y)}])
        end,
    E = sexp(sub(L(A, B), L(0, 0)), N + 2),
    Pre = smul(exact(#{0 => pconst({1, 1}), 1 => pconst({-2, 1})}),
               sinv(smul(seps(), exact(#{0 => pconst({1, 1}), 1 => pconst({-(2 + A + B), 1})})), N + 2)),
    truncate(smul(Pre, E), N).

%% the MS-bar factor per loop: eps G(1, 1) e^(gamma_E eps) = exp(sum_k zeta_k/k ((-1)^k + 2 - 2^k) eps^k) / (1 - 2 eps)
-spec msbar(pos_integer()) -> tuple().
msbar(N) ->
    Log = exact(maps:from_list([{K, pscale(psym(zeta(K)), qnorm({sign(K) + 2 - pow(2, K), K}))} || K <- lists:seq(3, N)])),
    truncate(smul(sexp(Log, N), sinv(exact(#{0 => pconst({1, 1}), 1 => pconst({-2, 1})}), N)), N).

%%% Evaluating expressions

%% An AST in d, ep, zeta symbols, gc_A_B (GschemeConstants(A, B)) and master_X (Forcer's Master(X), from Masters)
%% as a series to N orders: #{Name => Series} for the masters.
-spec eval(term(), map(), pos_integer()) -> tuple().
eval(N, _, _) when is_integer(N) -> sconst({N, 1});
eval({rat, A, B}, _, _) -> sconst(qnorm({A, B}));
eval({sym, ep}, _, _) -> seps();
eval({sym, d}, _, _) -> exact(#{0 => pconst({4, 1}), 1 => pconst({-2, 1})});
eval({sym, S}, Masters, N) ->
    case atom_to_list(S) of
        "gc_" ++ AB ->
            [A, B] = [list_to_integer(X) || X <- string:split(AB, "_")],
            gc(A, B, N);
        "master_" ++ Name -> maps:get(Name, Masters);
        %% zeta values, colour factors, flavours: constants in epsilon
        _ -> exact(#{0 => psym(S)})
    end;
eval({add, Es}, M, N) -> lists:foldl(fun(E, Acc) -> sadd(Acc, eval(E, M, N)) end, exact(#{}), Es);
eval({mul, Es}, M, N) -> lists:foldl(fun(E, Acc) -> smul(Acc, eval(E, M, N)) end, sconst({1, 1}), Es);
eval({pow, B, K}, M, N) when is_integer(K) -> spow(eval(B, M, N), K, N).

%% FORM's printed expression (Forcer's output, or its master tables) as a series: function calls renamed to
%% symbols, then parsed by pq_symbolic
-spec series_of_form(binary() | string(), map()) -> {term(), fun((pos_integer()) -> tuple())}.
series_of_form(Text, Masters) ->
    T0 = re:replace(Text, "\\s+", "", [global, {return, list}]),
    T1 = re:replace(T0, "GschemeConstants\\(([0-9]+),([0-9]+)\\)", "gc_\\1_\\2", [global, {return, list}]),
    T2 = re:replace(T1, "Master\\(([a-z0-9]+)\\)", "master_\\1", [global, {return, list}]),
    T3a = re:replace(T2, "rat\\(", "ratf(", [global, {return, list}]),
    %% a sign right after an opening bracket (FORM's bracketed output): a 0 in front, for the parser
    T3 = re:replace(T3a, "\\(([-+])", "(0\\1", [global, {return, list}]),
    %% FORM writes a leading unary plus; a 0 in front keeps the parser's grammar
    {ok, Ast} = pq_symbolic:parse(list_to_binary("0" ++ unrat(T3))),
    {Ast, fun(N) -> eval(Ast, Masters, N) end}.

%% ratf(A,B) -> ((A)/(B)), the arguments split at the top-level comma
unrat(S) ->
    case string:find(S, "ratf(") of
        nomatch -> S;
        Rest0 ->
            Pre = string:slice(S, 0, length(S) - length(Rest0)),
            Inner = string:slice(Rest0, 5),
            {A, B, After} = split_args(Inner, 0, [], []),
            unrat(Pre ++ "((" ++ A ++ ")/(" ++ B ++ "))" ++ After)
    end.

split_args([$, | T], 0, Acc, []) -> split_args(T, 0, [], lists:reverse(Acc));
split_args([$) | T], 0, Acc, A) -> {A, lists:reverse(Acc), T};
split_args([$( | T], D, Acc, A) -> split_args(T, D + 1, [$( | Acc], A);
split_args([$) | T], D, Acc, A) -> split_args(T, D - 1, [$) | Acc], A);
split_args([C | T], D, Acc, A) -> split_args(T, D, [C | Acc], A).

%% a series as FORM prints it: + ep^-3 * ( 1/24 ) ... (for the report)
-spec to_string(tuple()) -> iolist().
to_string({C, P}) ->
    Terms = [io_lib:format("  ep^~b: ~s~n", [K, poly_string(V)]) || {K, V} <- lists:sort(maps:to_list(C))],
    [Terms, case P > ?INF div 2 of true -> ""; false -> io_lib:format("  + O(ep^~b)~n", [P]) end].

poly_string(P) ->
    lists:join(" + ", [[case Q of {N, 1} -> integer_to_list(N); {N, D} -> io_lib:format("~b/~b", [N, D]) end,
                        [["*", atom_to_list(S), case E of 1 -> ""; _ -> "^" ++ integer_to_list(E) end] || {S, E} <- M]]
                       || {M, Q} <- lists:sort(maps:to_list(P))]).

%%% Rationals

qnorm({N, D}) when D < 0 -> qnorm({-N, -D});
qnorm({0, _}) -> {0, 1};
qnorm({N, D}) -> G = gcd(abs(N), D), {N div G, D div G}.
qadd({A, B}, {C, D}) -> qnorm({A * D + C * B, B * D}).
qmul({A, B}, {C, D}) -> qnorm({A * C, B * D}).
qinv({A, B}) -> qnorm({B, A}).
qneg({A, B}) -> {-A, B}.
gcd(A, 0) -> A;
gcd(A, B) -> gcd(B, A rem B).

%% the coefficients: #{Order => Poly}, a Poly #{Monomial => {N, D}}
coeffs({C, _}) -> C.
