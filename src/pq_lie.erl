%% Lie algebras for physics (plan "theoretical physics instruments", F1): angular-momentum coupling coefficients
%% exactly, and SU(N) colour algebra with N symbolic.
%%
%% Wigner symbols. Racah's formulas, the angular momenta given doubled (2j, 2m as integers) so that half-integers stay
%% exact. A 3j or 6j symbol is s sqrt(q) with s = -1, 0, 1 and q rational: it is returned as {S, Q} (wigner_3j/6,
%% wigner_6j/6, clebsch_gordan/6) or as an expression of the engine (to_expr/1). A 9j symbol is a sum over a product
%% of three 6j symbols (wigner_9j/9, an expression).
%%   3j:  (j1 j2 j3; m1 m2 m3) = delta(m1 + m2 + m3) (-1)^(j1 - j2 - m3) sqrt(Delta(j1 j2 j3))
%%        sqrt((j1+m1)! (j1-m1)! (j2+m2)! (j2-m2)! (j3+m3)! (j3-m3)!)
%%        sum_k (-1)^k / (k! (j3-j2+k+m1)! (j3-j1+k-m2)! (j1+j2-j3-k)! (j1-k-m1)! (j2-k+m2)!),
%%        Delta(a b c) = (a+b-c)! (a-b+c)! (-a+b+c)! / (a+b+c+1)!
%%   CG:  <j1 m1 j2 m2 | J M> = (-1)^(j1 - j2 + M) sqrt(2J + 1) (j1 j2 J; m1 m2 -M)
%%   6j:  {j1 j2 j3; j4 j5 j6} = sqrt(Delta(j1 j2 j3) Delta(j1 j5 j6) Delta(j4 j2 j6) Delta(j4 j5 j3))
%%        sum_t (-1)^t (t+1)! / ((t-a1)! (t-a2)! (t-a3)! (t-a4)! (b1-t)! (b2-t)! (b3-t)!)
%%   9j:  sum_x (-1)^(2x) (2x+1) {j1 j4 j7; j8 j9 x} {j2 j5 j8; j4 x j6} {j3 j6 j9; x j1 j2}
%%
%% Colour algebra. A colour expression is a sum of products of traces of SU(N) generators in the fundamental
%% representation, Tr(T^a T^b ...), normalised Tr(T^a T^b) = delta^ab/2, the adjoint labels summed when they appear
%% twice. The Fierz identity sum_a T^a_ij T^a_kl = (1/2)(delta_il delta_kj - delta_ij delta_kl/N) removes them one by one:
%%   Tr(A T^a B T^a) = (1/2) (Tr(A) Tr(B) - Tr(A B)/N),   Tr(A T^a) Tr(B T^a) = (1/2) (Tr(A B) - Tr(A) Tr(B)/N),
%% with Tr() = N and Tr(T^a) = 0, so every fully contracted expression is a Laurent polynomial in N, exactly.
%% f^abc = -2i Tr([T^a, T^b] T^c) and d^abc = 2 Tr({T^a, T^b} T^c) are written as traces.
-module(pq_lie).
-export([wigner_3j/6, wigner_6j/6, wigner_9j/9, clebsch_gordan/6, to_expr/1,
         contract/1, colour/1, invariant/1]).

%%% Wigner symbols (arguments doubled)

wigner_3j(J1, J2, J3, M1, M2, M3) ->
    case M1 + M2 + M3 =:= 0 andalso triangle(J1, J2, J3) andalso abs(M1) =< J1 andalso abs(M2) =< J2 andalso abs(M3) =< J3
         andalso even(J1 + M1) andalso even(J2 + M2) andalso even(J3 + M3) of
        false -> {0, {0, 1}};
        true ->
            H = fun(X) -> X div 2 end,                       % halves of sums that are even
            Sq = q_mul(delta(J1, J2, J3),
                       {fact(H(J1 + M1)) * fact(H(J1 - M1)) * fact(H(J2 + M2)) * fact(H(J2 - M2)) * fact(H(J3 + M3)) * fact(H(J3 - M3)), 1}),
            Kmin = max(0, max(H(J2 - J3 - M1), H(J1 - J3 + M2))),
            Kmax = min(H(J1 + J2 - J3), min(H(J1 - M1), H(J2 + M2))),
            Sum = lists:foldl(fun(K, Acc) ->
                                      D = fact(K) * fact(H(J3 - J2 + M1) + K) * fact(H(J3 - J1 - M2) + K) * fact(H(J1 + J2 - J3) - K)
                                          * fact(H(J1 - M1) - K) * fact(H(J2 + M2) - K),
                                      q_add(Acc, {sign(K) , D})
                              end, {0, 1}, lists:seq(Kmin, Kmax)),
            signed(H(J1 - J2 - M3), Sum, Sq)
    end.

clebsch_gordan(J1, M1, J2, M2, J, M) ->
    {S, Q} = wigner_3j(J1, J2, J, M1, M2, -M),
    {S * phase((J1 - J2 + M) div 2), q_mul(Q, {J + 1, 1})}.

wigner_6j(J1, J2, J3, J4, J5, J6) ->
    case triangle(J1, J2, J3) andalso triangle(J1, J5, J6) andalso triangle(J4, J2, J6) andalso triangle(J4, J5, J3) of
        false -> {0, {0, 1}};
        true ->
            H = fun(X) -> X div 2 end,
            Sq = lists:foldl(fun q_mul/2, {1, 1}, [delta(J1, J2, J3), delta(J1, J5, J6), delta(J4, J2, J6), delta(J4, J5, J3)]),
            A = [H(J1 + J2 + J3), H(J1 + J5 + J6), H(J4 + J2 + J6), H(J4 + J5 + J3)],
            B = [H(J1 + J2 + J4 + J5), H(J2 + J3 + J5 + J6), H(J3 + J1 + J6 + J4)],
            Sum = lists:foldl(fun(T, Acc) ->
                                      D = lists:foldl(fun(X, P) -> P * fact(T - X) end, 1, A) * lists:foldl(fun(X, P) -> P * fact(X - T) end, 1, B),
                                      q_add(Acc, {sign(T) * fact(T + 1), D})
                              end, {0, 1}, lists:seq(lists:max(A), lists:min(B))),
            signed(0, Sum, Sq)
    end.

%% a sum of signed square roots: an expression of the engine
wigner_9j(J1, J2, J3, J4, J5, J6, J7, J8, J9) ->
    Lo = lists:max([abs(J1 - J9), abs(J4 - J8), abs(J2 - J6)]),
    Hi = lists:min([J1 + J9, J4 + J8, J2 + J6]),
    Terms = [{mul, [phase(X) * (X + 1), to_expr(wigner_6j(J1, J4, J7, J8, J9, X)), to_expr(wigner_6j(J2, J5, J8, J4, X, J6)),
                    to_expr(wigner_6j(J3, J6, J9, X, J1, J2))]}
             || X <- lists:seq(Lo, Hi, 2)],
    {ok, E} = pq_symbolic:canonicalise({add, Terms}),
    E.

%% s sqrt(q) as an expression
to_expr({0, _}) -> 0;
to_expr({S, {N, D}}) -> {mul, [S, {pow, {rat, N, D}, {rat, 1, 2}}]}.

%% S * Sum * sqrt(Sq) as {sign, square}, Phase the exponent of -1
signed(Phase, {Sn, _} = Sum, Sq) ->
    case Sn of
        0 -> {0, {0, 1}};
        _ -> {phase(Phase) * sgn(Sn), q_mul(q_mul(Sum, Sum), Sq)}
    end.

delta(A, B, C) -> H = fun(X) -> X div 2 end,
                  {fact(H(A + B - C)) * fact(H(A - B + C)) * fact(H(-A + B + C)), fact(H(A + B + C) + 1)}.

triangle(A, B, C) -> A >= 0 andalso B >= 0 andalso C >= 0 andalso C =< A + B andalso C >= abs(A - B) andalso even(A + B + C).

even(X) -> X rem 2 =:= 0.
phase(K) -> case K rem 2 of 0 -> 1; _ -> -1 end.
sign(K) -> phase(K).
sgn(X) when X > 0 -> 1;
sgn(X) when X < 0 -> -1.

fact(N) when N =< 1 -> 1;
fact(N) -> N * fact(N - 1).

q_add({A, B}, {C, D}) -> norm(A * D + C * B, B * D).
q_mul({A, B}, {C, D}) -> norm(A * C, B * D).
norm(N, D) when D < 0 -> norm(-N, -D);
norm(N, D) -> G = gcd(abs(N), D), {N div G, D div G}.
gcd(A, 0) -> A;
gcd(A, B) -> gcd(B, A rem B).

%%% Colour algebra
%%%
%%% A term is {Coef, [Trace]}: Coef a Laurent polynomial in N (a map Power => {Num, Den}, possibly with a factor i^k
%%% kept apart as {I, Poly}), Trace a list of adjoint labels (any terms), in order. The result of a full contraction is
%%% a Laurent polynomial.

%% sum of {Coefficient, [Trace]}: the full contraction as a Laurent polynomial #{Power => {Num, Den}}; the
%% coefficients are rationals {N, D} or integers
contract(Terms) ->
    lists:foldl(fun({C, Traces}, Acc) -> padd(Acc, pscale(reduce(Traces), q(C))) end, #{}, Terms).

q({N, D}) -> norm(N, D);
q(N) when is_integer(N) -> {N, 1}.

%% a product of traces, fully contracted
reduce(Traces) ->
    Key = canon(Traces),
    case get({?MODULE, Key}) of
        undefined -> V = reduce1(Key), put({?MODULE, Key}, V), V;
        V -> V
    end.

reduce1(Traces) ->
    case lists:keytake([], 1, [{T, x} || T <- Traces]) of
        {value, _, Rest} -> pmul(#{1 => {1, 1}}, reduce([T || {T, _} <- Rest]));    % Tr() = N
        false ->
            case [T || T <- Traces, length(T) =:= 1] of
                [_ | _] -> #{};                                                      % Tr(T^a) = 0
                [] ->
                    case pair(Traces) of
                        none -> case Traces of [] -> #{0 => {1, 1}}; _ -> error({pq_lie, {open_labels, Traces}}) end;
                        {same, Others, A, B} ->
                            % Tr(A T B T) = (1/2)(Tr A Tr B - Tr(A B)/N)
                            padd(pscale(reduce([A, B | Others]), {1, 2}),
                                 pscale(pmul(#{-1 => {1, 1}}, reduce([A ++ B | Others])), {-1, 2}));
                        {apart, Others, A, B} ->
                            % Tr(A T) Tr(B T) = (1/2)(Tr(A B) - Tr A Tr B / N)
                            padd(pscale(reduce([A ++ B | Others]), {1, 2}),
                                 pscale(pmul(#{-1 => {1, 1}}, reduce([A, B | Others])), {-1, 2}))
                    end
            end
    end.

%% a label that appears twice: in one trace (the trace rotated to A T B T, returning A and B) or in two (each rotated
%% to end with it, returning their fronts)
pair(Traces) ->
    Counts = lists:foldl(fun(L, M) -> M#{L => maps:get(L, M, 0) + 1} end, #{}, lists:append(Traces)),
    case [L || L := 2 <- Counts] of
        [] -> none;
        Ls ->
            L = lists:min(Ls),
            case [T || T <- Traces, length([X || X <- T, X =:= L]) =:= 2] of
                [T | _] ->
                    Others = Traces -- [T],
                    R = rotate_to(T, L),                        % [L | A ++ [L] ++ B]
                    [L | Rest] = R,
                    {A, [L | B]} = lists:splitwith(fun(X) -> X =/= L end, Rest),
                    {same, Others, B, A};
                [] ->
                    [T1, T2] = [T || T <- Traces, lists:member(L, T)],
                    Others = (Traces -- [T1]) -- [T2],
                    [L | A] = rotate_to(T1, L), [L | B] = rotate_to(T2, L),
                    {apart, Others, A, B}
            end
    end.

rotate_to(T, L) ->
    {Before, After} = lists:splitwith(fun(X) -> X =/= L end, T),
    After ++ Before.

%% a canonical form of a product of traces: each trace rotated to its least rotation, the traces sorted (labels are
%% dummies, but renaming them is not attempted: the cache only needs to be correct)
canon(Traces) -> lists:sort([least_rotation(T) || T <- Traces]).
least_rotation([]) -> [];
least_rotation(T) -> lists:min([begin {A, B} = lists:split(I, T), B ++ A end || I <- lists:seq(0, length(T) - 1)]).

%% Laurent polynomials
padd(A, B) -> maps:filter(fun(_, {N, _}) -> N =/= 0 end,
                          maps:fold(fun(K, V, M) -> M#{K => case M of #{K := W} -> q_add(V, W); _ -> V end} end, A, B)).
pscale(A, C) -> maps:filter(fun(_, {N, _}) -> N =/= 0 end, maps:map(fun(_, V) -> q_mul(V, C) end, A)).
pmul(A, B) -> maps:fold(fun(K1, V1, M) -> maps:fold(fun(K2, V2, M2) ->
                                                            K = K1 + K2, V = q_mul(V1, V2),
                                                            M2#{K => case M2 of #{K := W} -> q_add(V, W); _ -> V end}
                                                    end, M, B) end, #{}, A).

%% a Laurent polynomial as an expression of the engine in {sym, n}, cancelled
colour(Poly) ->
    E = {add, [{mul, [{rat, N, D}, {pow, {sym, n}, K}]} || K := {N, D} <- Poly]},
    {ok, C} = pq_symbolic:cancel(E),
    C.

%%% The invariants, from their definitions by traces

%% f^abc = -2i (Tr(abc) - Tr(bac)) and d^abc = 2 (Tr(abc) + Tr(bac)), as [{Coef, Trace}] with Coef = {Re, Im} rationals
f(A, B, C) -> [{{0, -2}, [A, B, C]}, {{0, 2}, [B, A, C]}].
d(A, B, C) -> [{{2, 0}, [A, B, C]}, {{2, 0}, [B, A, C]}].

%% the product of factors (each a list of {Complex, Trace}), expanded: [{Complex, [Trace]}]
expand(Factors) ->
    lists:foldl(fun(F, Acc) -> [{cmul(C1, C2), Ts ++ [T]} || {C1, Ts} <- Acc, {C2, T} <- F] end, [{{1, 0}, []}], Factors).

cmul({A, B}, {C, D}) -> {A * C - B * D, A * D + B * C}.

%% the contraction of complex-weighted products; the result must be real
contract_c(Terms) ->
    Re = contract([{C, Ts} || {{C, _}, Ts} <- Terms, C =/= 0]),
    Im = contract([{C, Ts} || {{_, C}, Ts} <- Terms, C =/= 0]),
    case map_size(Im) of 0 -> Re; _ -> error({pq_lie, {not_real, Im}}) end.

%% the standard invariants as expressions in n, from their definitions by traces (Tr(T^a T^b) = delta^ab/2):
%%   tf      Tr(T^a T^b) = TF delta^ab, TF = 1/2;           na  N_A = N^2 - 1, the adjoint's dimension
%%   cf      T^a T^a = CF 1                                ca  f^acd f^bcd = CA delta^ab
%%   d33     d^abc d^abc
%%   d44_ff, d44_fa, d44_aa   d_R1^abcd d_R2^abcd / N_A for the fundamental (F) and adjoint (A) representations,
%%           d_R^abcd the trace of four generators symmetrised over their order (the quartic Casimirs of the
%%           four-loop beta function); the adjoint's generators (F^a)_xy = -i f^axy
invariant(na) -> colour(na());
invariant(tf) -> per_na(contract([{1, [[a, a]]}]));
invariant(cf) -> colour(pmul(contract([{1, [[a, a]]}]), #{-1 => {1, 1}}));
invariant(ca) -> per_na(contract_c(expand([f(a, c, d), f(a, c, d)])));
invariant(d33) -> colour(contract_c(expand([d(a, b, c), d(a, b, c)])));
invariant(d44_ff) -> per_na(sym_contract(fun(P) -> [fund(P)] end, fun(P) -> [fund(P)] end));
invariant(d44_fa) -> per_na(sym_contract(fun(P) -> [fund(P)] end, fun(P) -> adj(P, 1) end));
invariant(d44_aa) -> per_na(sym_contract(fun(P) -> adj(P, 1) end, fun(P) -> adj(P, 2) end)).

na() -> #{2 => {1, 1}, 0 => {-1, 1}}.

%% P / (N^2 - 1) as an expression, cancelled
per_na(P) ->
    {ok, C} = pq_symbolic:cancel({mul, [colour(P), {pow, {add, [{pow, {sym, n}, 2}, -1]}, -1}]}),
    C.

%% sum over the six orders p of (b, c, d) for each factor, weight 1/36: d_R1^abcd d_R2^abcd
sym_contract(F1, F2) ->
    Orders = [[a, B, C, D] || [B, C, D] <- perms([b, c, d])],
    Terms = lists:append([[{cmul({1, 0}, C), Ts} || {C, Ts} <- expand(F1(P1) ++ F2(P2))] || P1 <- Orders, P2 <- Orders]),
    pscale(contract_c(Terms), {1, 36}).

%% Tr(T^a T^b T^c T^d) as a factor
fund(Labels) -> [{{1, 0}, Labels}].

%% Tr(F^a F^b F^c F^d) = (-i)^4 f^(a x1 x2) f^(b x2 x3) f^(c x3 x4) f^(d x4 x1), dummies tagged by K: four factors
adj([A, B, C, D], K) ->
    X = fun(I) -> {x, K, I} end,
    [f(A, X(1), X(2)), f(B, X(2), X(3)), f(C, X(3), X(4)), f(D, X(4), X(1))].

perms([]) -> [[]];
perms(L) -> [[H | T] || H <- L, T <- perms(L -- [H])].
