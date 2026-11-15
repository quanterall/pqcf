%% Exact steady states and spectra of master equations (quantum optics plan, Q3): the Liouvillian of pq_qop operators
%% with rational parameters, in Gaussian rationals (pq_sym_exact's {q, Re, Im}), solved by exact elimination.
%%
%% In the Fock basis the matrix elements of a and a^dagger are square roots, so even a rational model has an
%% irrational Liouvillian. In the basis f_n = sqrt(n!) |n> they are integers: a^dagger f_n = f_(n+1), a f_n = n f_(n-1).
%% With S = diag(sqrt(n!)) (the atoms untouched) and rho' = S^-1 rho S^-1, a term A rho B of the master equation becomes
%% A' rho' B'' with A' = S^-1 A S (the "left" form) and B'' = S B S^-1 (the "right" form), both rational; a trace
%% Tr(O rho) = Tr(S O S rho') (the "trace" form, rational too). The jump terms L rho L^dagger take the products of the
%% coefficients c conj(c') symbolically, under the declared domains, so that sqrt(gamma) conj(sqrt(gamma)) = gamma.
%%
%% Vectors are maps Index => coefficient (the zeros left out); rho' is vectorised row by row (rho'_ij at i D + j);
%% matrices are maps {Row, Col} => coefficient.
-module(pq_qo_exact).
-export([liouvillian/5, steady_state/2, operator/5, expect/3, apply_left/3, spectrum/4, real_value/1, kerr_moment/3, kerr2_moment/4, steady_state_symbolic/6]).

-define(E, pq_sym_exact).

%%% Operators in the rescaled basis

%% the matrix of a pq_qop expression in one of the forms left (S^-1 O S), right (S O S^-1) or trace (S O S), the
%% coefficients exact at the rational Params (#{omega => 2, gamma => {1, 2}}); Basis from pq_lindblad:basis/2
operator(Basis, Spaces, Expr, Params, Form) ->
    pq_qop:with_spaces(Spaces, fun() ->
        lists:foldl(fun({{R, C}, V}, Acc) -> madd(Acc, R, C, V) end, #{},
                    lists:append([word(W, coef(Cf, Params, Spaces), Basis, Form) || W := Cf <- Expr]))
    end).

%% the elements {{Row, Col}, Value} of one word: the tensor product over the subsystems, the last name fastest
word(W, Cv, Basis, Form) ->
    lists:foldl(fun({N, S, Kind}, Acc) ->
                        Local = case proplists:get_value(N, W, none) of
                                    % the identity: S 1 S = diag(n!) in the trace form, 1 otherwise
                                    none -> [{I, I, case {Kind, Form} of {boson, trace} -> fact(I); _ -> 1 end}
                                             || I <- lists:seq(0, S - 1)];
                                    F -> factor(Kind, F, S, Form)
                                end,
                        [{{R * S + R1, C * S + C1}, ?E:cmul(V, ?E:c(V1))} || {{R, C}, V} <- Acc, {R1, C1, V1} <- Local]
                end, [{{0, 0}, Cv}], Basis).

%% (a^dagger)^P a^Q from f_n to f_m, m = n - Q + P: left n!/(n-Q)!, right m!/(n-Q)!, trace n! m!/(n-Q)!
factor(boson, {P, Q}, S, Form) ->
    [{M, N, case Form of
                left -> fact(N) div fact(N - Q);
                right -> fact(M) div fact(N - Q);
                trace -> fact(N) * fact(M) div fact(N - Q)
            end} || N <- lists:seq(Q, S - 1), M <- [N - Q + P], M =< S - 1];
factor(atom, {I, J}, _, _) -> [{I - 1, J - 1, 1}].

fact(0) -> 1;
fact(N) -> N * fact(N - 1).

%%% The Liouvillian

%% -i [H, rho] + sum_k (L_k rho L_k^dagger - (1/2) {L_k^dagger L_k, rho}) on rho', as {liouvillian, D, Rows}: Rows a
%% map Row => [{Col, Value}]. H and Ls are pq_qop expressions, Spaces as for pq_qop (domains included).
liouvillian(Basis, Spaces, H, Ls, Params) ->
    D = pq_lindblad:dim(Basis),
    Op = fun(E, Form) -> operator(Basis, Spaces, E, Params, Form) end,
    Id = maps:from_list([{{I, I}, ?E:c(1)} || I <- lists:seq(0, D - 1)]),
    MinusI = ?E:cq({0, 1}, {-1, 1}), PlusI = ?E:cq({0, 1}, {1, 1}), Half = ?E:c({-1, 2}),
    Terms0 = [{Op(H, left), Id, MinusI}, {Id, Op(H, right), PlusI}],
    Cf = fun(E) -> coef(E, Params, Spaces) end,
    Jumps = lists:append([jump_terms(L, Op, Cf, Id, Half, Spaces) || L <- Ls]),
    Elems = lists:foldl(fun({A, B, K}, Acc) -> super(A, B, K, D, Acc) end, #{}, Terms0 ++ Jumps),
    Rows = maps:fold(fun({R, C}, V, Acc) -> Acc#{R => [{C, V} | maps:get(R, Acc, [])]} end, #{}, Elems),
    {liouvillian, D, Rows}.

%% L rho L^dagger word by word: the coefficient c_w conj(c_w') of each pair taken symbolically (the domains apply);
%% then -(1/2) (L^dagger L rho + rho L^dagger L)
jump_terms(L, Op, Cf, Id, Half, Spaces) ->
    Pairs = [{#{W1 => 1}, #{dag_word(W2) => 1}, {mul, [C1, {conj, C2}]}} || W1 := C1 <- L, W2 := C2 <- L],
    K = pq_qop:with_spaces(Spaces, fun() -> pq_qop:mul(pq_qop:dag(L), L) end),
    [{Op(A, left), Op(B, right), Cf(C)} || {A, B, C} <- Pairs]
        ++ [{Op(K, left), Id, Half}, {Id, Op(K, right), Half}].

dag_word(W) -> [{N, {Y, X}} || {N, {X, Y}} <- W].

%% A rho B into the superoperator: the element ((i,j),(k,l)) = K A_ik B_lj
super(A, B, Kv, D, Acc) ->
    ByRowB = maps:fold(fun({L, J}, V, M) -> M#{L => [{J, V} | maps:get(L, M, [])]} end, #{}, B),
    maps:fold(fun({I, Kk}, Va, Acc1) ->
                      lists:foldl(fun({L, J}, Acc2) ->
                                          Vb = maps:get({L, J}, B),
                                          madd(Acc2, I * D + J, Kk * D + L, ?E:cmul(Kv, ?E:cmul(Va, Vb)))
                                  end, Acc1, [{L, J} || L <- lists:seq(0, D - 1), {J, _} <- maps:get(L, ByRowB, [])])
              end, Acc, A).

madd(M, R, C, V) ->
    case ?E:is_zero(V) of
        true -> M;
        false ->
            case maps:find({R, C}, M) of
                {ok, V0} ->
                    S = ?E:cadd(V0, V),
                    case ?E:is_zero(S) of true -> maps:remove({R, C}, M); false -> M#{{R, C} => S} end;
                error -> M#{{R, C} => V}
            end
    end.

%%% Coefficients: an engine expression at rational parameters, exactly; the domains of the spaces apply

coef(Expr, Params, Spaces) ->
    Domains = maps:get(domains, Spaces, #{}),
    Sub = pq_symbolic:substitute(Expr, maps:map(fun(_, V) -> num(V) end, Params)),
    case pq_symbolic:canonicalise(Sub, #{domains => Domains}) of
        {ok, C} -> exact(C);
        {error, E} -> error({pq_qo_exact, E})
    end.

num({N, D}) when is_integer(N), is_integer(D) -> {rat, N, D};
num({rat, _, _} = R) -> R;
num(N) when is_integer(N) -> N.

exact(N) when is_integer(N) -> ?E:c(N);
exact({rat, N, D}) -> ?E:c({N, D});
exact({complex, A, B}) -> ?E:cadd(exact(A), ?E:cmul(?E:cq({0, 1}, {1, 1}), exact(B)));
exact({add, L}) -> lists:foldl(fun(X, S) -> ?E:cadd(S, exact(X)) end, ?E:c(0), L);
exact({mul, L}) -> lists:foldl(fun(X, S) -> ?E:cmul(S, exact(X)) end, ?E:c(1), L);
exact({pow, B, K}) when is_integer(K) -> ?E:cpow_int(exact(B), K);
exact({conj, X}) -> ?E:cconj(exact(X));
exact(E) -> error({pq_qo_exact, {not_rational, E}}).

%%% The steady state: L rho' = 0 with Tr rho = 1, by exact elimination

%% rho' as a vector map; the first equation replaced by the trace, sum_i w_i rho'_ii = 1, w_i the product of n! over
%% the bosons of state i (Tr rho = Tr(S rho' S))
steady_state({liouvillian, D, Rows}, Basis) ->
    W = weights(Basis),
    Trace = maps:from_list([{I * D + I, ?E:c(maps:get(I, W))} || I <- lists:seq(0, D - 1)]),
    Eqs = [{Trace, ?E:c(1)} | [{maps:from_list(maps:get(R, Rows, [])), ?E:c(0)} || R <- lists:seq(1, D * D - 1)]],
    solve(Eqs, D * D).

weights(Basis) ->
    {_, W} = lists:foldl(fun({_, S, Kind}, {Sz, Acc}) ->
                                 Loc = [case Kind of boson -> fact(N); atom -> 1 end || N <- lists:seq(0, S - 1)],
                                 {Sz * S, maps:from_list([{I * S + J, maps:get(I, Acc) * lists:nth(J + 1, Loc)}
                                                          || I <- lists:seq(0, Sz - 1), J <- lists:seq(0, S - 1)])}
                         end, {1, #{0 => 1}}, Basis),
    W.

%% Gaussian elimination on sparse rows [{RowMap, Rhs}], column by column; the pivot the shortest row (exact
%% arithmetic: any non-zero pivot will do, the shortest keeps the fill-in small)
solve(Eqs, N) -> back(elim(Eqs, 0, N, []), #{}).

elim(_, Col, N, Done) when Col >= N -> Done;
elim(Eqs, Col, N, Done) ->
    {With, Without} = lists:partition(fun({Row, _}) -> is_map_key(Col, Row) end, Eqs),
    case With of
        [] -> elim(Without, Col + 1, N, Done);
        _ ->
            [{PRow, PRhs} = Pivot | _] = lists:sort(fun({A, _}, {B, _}) -> map_size(A) =< map_size(B) end, With),
            Inv = ?E:cinv(maps:get(Col, PRow)),
            Rest = [begin
                        F = ?E:cmul(maps:get(Col, Row), Inv),
                        {row_sub(Row, PRow, F, Col), ?E:cadd(Rhs, ?E:cneg(?E:cmul(F, PRhs)))}
                    end || {Row, Rhs} = E <- With, E =/= Pivot],
            elim(Rest ++ Without, Col + 1, N, [{Col, PRow, PRhs} | Done])
    end.

row_sub(Row, PRow, F, Col) ->
    maps:remove(Col, maps:fold(fun(C, V, Acc) ->
                                       New = ?E:cadd(maps:get(C, Acc, ?E:c(0)), ?E:cneg(?E:cmul(F, V))),
                                       case ?E:is_zero(New) of true -> maps:remove(C, Acc); false -> Acc#{C => New} end
                               end, Row, PRow)).

back([], Sol) -> Sol;
back([{Col, Row, Rhs} | T], Sol) ->
    S = maps:fold(fun(C, V, Acc) when C =/= Col -> ?E:cadd(Acc, ?E:cneg(?E:cmul(V, maps:get(C, Sol, ?E:c(0)))));
                     (_, _, Acc) -> Acc
                  end, Rhs, Row),
    X = ?E:cmul(S, ?E:cinv(maps:get(Col, Row))),
    back(T, case ?E:is_zero(X) of true -> Sol; false -> Sol#{Col => X} end).

%%% Expectation values and spectra

%% Tr(O rho) from the trace form of O (operator(..., trace)) and rho': sum_ij O_ij rho'_ji
expect(Otr, Rho, D) ->
    maps:fold(fun({I, J}, V, S) -> ?E:cadd(S, ?E:cmul(V, maps:get(J * D + I, Rho, ?E:c(0)))) end, ?E:c(0), Otr).

%% B rho' from the left form of B: (B rho')_ij = sum_k B_ik rho'_kj
apply_left(Bl, Rho, D) ->
    maps:fold(fun({I, K}, V, Acc) ->
                      lists:foldl(fun(J, A2) ->
                                          case maps:find(K * D + J, Rho) of
                                              {ok, R} -> vadd(A2, I * D + J, ?E:cmul(V, R));
                                              error -> A2
                                          end
                                  end, Acc, lists:seq(0, D - 1))
              end, #{}, Bl).

vadd(M, K, V) ->
    S = ?E:cadd(maps:get(K, M, ?E:c(0)), V),
    case ?E:is_zero(S) of true -> maps:remove(K, M); false -> M#{K => S} end.

%% The incoherent spectrum by the quantum regression theorem: S(w) = Re Tr[A (-(L + i w))^-1 (B rho - <B> rho)] =
%% Re of the integral over tau > 0 of e^(i w tau) <dA(tau) dB(0)>, at the rational frequencies Ws (non-zero), exactly.
%% Spec = {Atrace, Bleft, Btrace, Rho, Basis}: A's trace form, B's left and trace forms. [{W, S(W)}], S as a rational.
spectrum({liouvillian, D, Rows}, {Atr, Bl, Btr, Rho, _Basis}, Ws, _Opts) ->
    Bmean = expect(Btr, Rho, D),
    X0 = apply_left(Bl, Rho, D),
    X = maps:fold(fun(K, V, Acc) -> vadd(Acc, K, ?E:cneg(?E:cmul(Bmean, V))) end, X0, Rho),
    [begin
         Wq = ?E:cq({0, 1}, to_rat(W)),
         % (L + i w) y = -X
         Eqs = [{row_plus(maps:get(R, Rows, []), R, Wq), ?E:cneg(maps:get(R, X, ?E:c(0)))} || R <- lists:seq(0, D * D - 1)],
         Y = solve(Eqs, D * D),
         {W, real_value(expect(Atr, Y, D))}
     end || W <- Ws].

row_plus(Row, R, Wq) ->
    M = maps:from_list(Row),
    S = ?E:cadd(maps:get(R, M, ?E:c(0)), Wq),
    case ?E:is_zero(S) of true -> maps:remove(R, M); false -> M#{R => S} end.

to_rat(N) when is_integer(N) -> {N, 1};
to_rat({N, D}) -> ?E:rat(N, D).

real_value({q, Re, _}) -> Re.

%%% The driven Kerr cavity, exactly (Drummond and Walls, J. Phys. A 13, 725 (1980), by the complex-P representation)
%%%
%%% For H = Delta a^dagger a + chi a^dagger^2 a^2 + E (a + a^dagger) and L = sqrt(kappa) a, the steady state's normally
%%% ordered moments are
%%%   <a^dagger^j a^j> = |x|^(2j) / prod_(k<j) |c + k|^2 * S(j) / S(0),
%%%   S(j) = sum_n z^n / (n! prod_(k<n) |c + j + k|^2),   x = E/chi, c = (Delta - i kappa/2)/chi, z = 2 |x|^2.
%%% The Gamma functions of the published form, |Gamma(c)|^2 / |Gamma(c + j)|^2, are these finite products, so every term is a
%%% real rational for rational parameters: the series (0F2-type, converging like z^n/n!^3) is summed to the precision P.

%% <a^dagger^j a^j> at Params #{delta, chi, e, kappa} (rationals {N, D} or integers), as an mp at P bits
kerr_moment(#{delta := Dl, chi := Chi, e := E, kappa := K}, J, P) ->
    Q = fun(V) -> pq_mp:from(case V of {_, _} -> V; N -> {N, 1} end, P) end,
    Cr = pq_mp:divide(Q(Dl), Q(Chi), P),
    Ci = pq_mp:divide(pq_mp:mul(Q(K), pq_mp:from({1, 2}, P), P), Q(Chi), P),
    Abs2 = fun(Shift) -> R = pq_mp:add(Cr, pq_mp:from(Shift, P), P),
                         pq_mp:add(pq_mp:mul(R, R, P), pq_mp:mul(Ci, Ci, P), P) end,
    X = pq_mp:divide(Q(E), Q(Chi), P),
    X2 = pq_mp:mul(X, X, P),
    Z = pq_mp:mul(pq_mp:from(2, P), X2, P),
    Pre = lists:foldl(fun(Kk, Acc) -> pq_mp:divide(pq_mp:mul(Acc, X2, P), Abs2(Kk), P) end, pq_mp:from(1, P), lists:seq(0, J - 1)),
    pq_mp:divide(pq_mp:mul(Pre, kerr_series(Abs2, Z, J, P), P), kerr_series(Abs2, Z, 0, P), P).

%% S(j), summed until a term is below 2^-(P + 8) of the sum
kerr_series(Abs2, Z, J, P) -> kerr_series(Abs2, Z, J, 0, pq_mp:from(1, P), pq_mp:from(0, P), P).

kerr_series(Abs2, Z, J, N, T, S, P) ->
    S2 = pq_mp:add(S, T, P),
    T2 = pq_mp:divide(pq_mp:mul(T, Z, P), pq_mp:mul(pq_mp:from(N + 1, P), Abs2(J + N), P), P),
    case small(T2, S2, P) of
        true -> S2;
        false -> kerr_series(Abs2, Z, J, N + 1, T2, S2, P)
    end.

small({mp, 0, _}, _, _) -> true;
small(T, S, P) -> mag(T) < mag(S) - P - 8.

mag({mp, 0, _}) -> -(1 bsl 40);
mag({mp, M, E}) -> length(integer_to_list(abs(M), 2)) + E.

%%% Steady states as rational functions of the parameters (finite-field reconstruction, decision record R1)
%%%
%%% rho' is Hermitian, so its real coordinates are x_ii, and x_ij, y_ij for i < j (rho'_ij = x_ij + i y_ij): D^2 real
%%% unknowns. The equations are the real parts of (L rho')_ii and the real and imaginary parts of (L rho')_ij, i < j; the
%%% one dependent equation (row (0,0): the trace of L rho' vanishes) is replaced by Tr rho = 1. The entries are
%%% expressions of the engine in the parameters, and pq_recon:linear/4 reconstructs each unknown as a reduced rational
%%% function from eliminations modulo primes at sample points: no symbolic elimination.

%% Tr(O rho) for each O of Os (pq_qop expressions, Hermitian) as a reduced rational function of Vars (a list of the
%% parameters' symbols, declared real or positive in the spaces' domains): [{O, Expr}]
steady_state_symbolic(Basis, Spaces, H, Ls, Vars, Os) ->
    D = pq_lindblad:dim(Basis),
    Doms = maps:get(domains, Spaces, #{}),
    Can = fun(E) -> {ok, C} = pq_symbolic:canonicalise(E, #{domains => Doms}), C end,
    Op = fun(E, Form) -> operator_sym(Basis, Spaces, E, Form) end,
    Id = maps:from_list([{{I, I}, 1} || I <- lists:seq(0, D - 1)]),
    MinusI = {complex, 0, -1}, PlusI = {complex, 0, 1}, Half = {rat, -1, 2},
    Jumps = lists:append([begin
                              Pairs = [{#{W1 => 1}, #{dag_word(W2) => 1}, Can({mul, [C1, {conj, C2}]})} || W1 := C1 <- L, W2 := C2 <- L],
                              K = pq_qop:with_spaces(Spaces, fun() -> pq_qop:mul(pq_qop:dag(L), L) end),
                              [{Op(A, left), Op(B, right), Cf} || {A, B, Cf} <- Pairs]
                                  ++ [{Op(K, left), Id, Half}, {Id, Op(K, right), Half}]
                          end || L <- Ls]),
    Terms = [{Op(H, left), Id, MinusI}, {Id, Op(H, right), PlusI} | Jumps],
    %% the superoperator: {(i,j), (k,l)} => [products], summed and canonicalised once
    Acc = lists:foldl(fun({A, B, K}, M) ->
                              maps:fold(fun({I, Kk}, Va, M1) ->
                                                maps:fold(fun({L, J}, Vb, M2) ->
                                                                  Key = {I * D + J, Kk * D + L},
                                                                  M2#{Key => [{mul, [K, Va, Vb]} | maps:get(Key, M2, [])]}
                                                          end, M1, B)
                                        end, M, A)
                      end, #{}, Terms),
    Sup = maps:filtermap(fun(_, Es) -> case Can({add, Es}) of 0 -> false; E -> {true, E} end end, Acc),
    %% the real unknowns, numbered: x_ij (i =< j), then y_ij (i < j)
    Xs = [{x, I, J} || I <- lists:seq(0, D - 1), J <- lists:seq(I, D - 1)],
    Ys = [{y, I, J} || I <- lists:seq(0, D - 1), J <- lists:seq(I + 1, D - 1)],
    Unknowns = Xs ++ Ys,
    Index = maps:from_list(lists:zip(Unknowns, lists:seq(0, length(Unknowns) - 1))),
    %% rho'_kl in the unknowns: [{Unknown, Factor}]
    Coord = fun(K, L) when K =:= L -> [{{x, K, K}, 1}];
                (K, L) when K < L -> [{{x, K, L}, 1}, {{y, K, L}, {complex, 0, 1}}];
                (K, L) -> [{{x, L, K}, 1}, {{y, L, K}, {complex, 0, -1}}]
            end,
    ByRow = maps:fold(fun({R, C}, E, M) -> M#{R => [{C, E} | maps:get(R, M, [])]} end, #{}, Sup),
    Row = fun(I, J) ->
                  Cs = lists:foldl(fun({C, E}, M) ->
                                           lists:foldl(fun({U, F}, M2) -> M2#{U => [{mul, [F, E]} | maps:get(U, M2, [])]} end,
                                                       M, Coord(C div D, C rem D))
                                   end, #{}, maps:get(I * D + J, ByRow, [])),
                  maps:map(fun(_, Es) -> Can({add, Es}) end, Cs)
          end,
    Re = fun(E) -> Can({mul, [{rat, 1, 2}, {add, [E, {conj, E}]}]}) end,
    Im = fun(E) -> Can({mul, [{complex, 0, -1}, {rat, 1, 2}, {add, [E, {mul, [-1, {conj, E}]}]}]}) end,
    Dense = fun(M) -> Full = maps:from_list([{maps:get(U, Index), V} || U := V <- M, V =/= 0]),
                      [maps:get(C, Full, 0) || C <- lists:seq(0, length(Unknowns) - 1)] end,
    W = weights(Basis),
    Trace = Dense(maps:from_list([{{x, I, I}, maps:get(I, W)} || I <- lists:seq(0, D - 1)])),
    Rows = [Trace | [Dense(maps:map(fun(_, E) -> Re(E) end, Row(I, I))) || I <- lists:seq(1, D - 1)]]
        ++ lists:append([[Dense(maps:map(fun(_, E) -> Re(E) end, Row(I, J))),
                          Dense(maps:map(fun(_, E) -> Im(E) end, Row(I, J)))]
                         || I <- lists:seq(0, D - 1), J <- lists:seq(I + 1, D - 1)]),
    Rhs = [1 | lists:duplicate(length(Unknowns) - 1, 0)],
    {ok, Sol} = pq_recon:linear(Rows, Rhs, Vars, #{}),
    Val = maps:from_list([{lists:nth(C + 1, Unknowns),
                           case map_size(N) of 0 -> 0; _ -> {mul, [pq_sym_norm:to_ast(N), {pow, pq_sym_norm:to_ast(Dn), -1}]} end}
                          || {C, {N, Dn}} <- Sol]),
    Rho = fun(K, L) -> {add, [{mul, [F, maps:get(U, Val)]} || {U, F} <- Coord(K, L)]} end,
    [{O, begin
             Otr = operator_sym(Basis, Spaces, O, trace),
             {ok, E} = pq_symbolic:cancel(Can({add, [{mul, [V, Rho(J, I)]} || {I, J} := V <- Otr]})),
             E
         end} || O <- Os].

%% operator/5 with the coefficients kept as expressions: {Row, Col} => Expr
operator_sym(Basis, Spaces, Expr, Form) ->
    pq_qop:with_spaces(Spaces, fun() ->
        lists:foldl(fun({{R, C}, V}, Acc) -> Acc#{{R, C} => {add, [V | case Acc of #{{R, C} := V0} -> [V0]; _ -> [] end]}} end, #{},
                    lists:append([[{RC, {mul, [Cf, N]}} || {RC, {q, {N0, 1}, {0, 1}}} <- word(W, ?E:c(1), Basis, Form), N <- [N0]]
                                  || W := Cf <- Expr]))
    end).

%%% The Kerr resonator with one- and two-photon driving and dissipation, exactly (Bartolo, Minganti, Casteels, Ciuti,
%%% PRA 94, 033841 (2016), arXiv:1607.06739; Eqs. 5-11 there). For
%%%   H = -Delta ad a + (U/2) ad^2 a^2 + F (ad + a) + (G/2) (ad^2 + a^2),   L1 = sqrt(gamma) a,   L2 = sqrt(eta) a^2
%%% (F and G real here), with c = (Delta + i gamma/2)/(U - i eta), f = F/(U - i eta), g = G/(U - i eta):
%%%   <ad^i a^j> = (1/N) sum_m (2^m/m!) F_(m+j) conj(F_(m+i)),   N = sum_m (2^m/m!) |F_m|^2,
%%%   F_m = (i sqrt(g))^m 2F1(-m, -c - i f/sqrt(g); -2c; 2),
%%% the hypergeometric function terminating (a polynomial of degree m in 2). With U = 0 and F = 0 it is the degenerate
%%% parametric oscillator with its pump adiabatically eliminated (two-photon driving and two-photon loss). G must not be
%%% 0 (the limit g -> 0 is the Kerr cavity of kerr_moment/3). The sums are taken until a term falls below 2^-(P + 8) of
%%% the normalisation, at precision P + 64 bits (the terms of N alternate in size before they fall).

%% <ad^i a^j> as {Re, Im} mps at P bits; Params #{delta, u, f, g, gamma, eta} (rationals {N, D} or integers)
kerr2_moment(Params, I, J, P) ->
    Pw = P + 64,
    R = fun(K) -> V = maps:get(K, Params), pq_mp:from(case V of {_, _} -> V; N -> {N, 1} end, Pw) end,
    Den = {R(u), pq_mp:neg(R(eta))},
    Over = fun(Z) -> cdiv(Z, Den, Pw) end,
    C = Over({R(delta), pq_mp:mul(R(gamma), pq_mp:from({1, 2}, Pw), Pw)}),
    Fp = Over({R(f), zero_mp()}),
    G = Over({R(g), zero_mp()}),
    Sg = csqrt(G, Pw),
    Ip = {zero_mp(), pq_mp:from(1, Pw)},
    B = csub2(cneg(C), cmul2(Ip, cdiv(Fp, Sg, Pw), Pw), Pw),            % -c - i f/sqrt(g)
    Cc = cneg(cscale2(C, 2, Pw)),                                        % -2c
    Pre = cmul2(Ip, Sg, Pw),
    Fm = fun(M) -> cmul2(cpow2(Pre, M, Pw), hyp(M, B, Cc, Pw), Pw) end,
    Sum = fun(Term) -> sum_series(Term, 0, {zero_mp(), zero_mp()}, zero_mp(), Pw, 0) end,
    Norm = Sum(fun(M) -> cscale_r(cmul2(Fm(M), cconj2(Fm(M)), Pw), w(M, Pw), Pw) end),
    Corr = Sum(fun(M) -> cscale_r(cmul2(Fm(M + J), cconj2(Fm(M + I)), Pw), w(M, Pw), Pw) end),
    {Re, Im} = cdiv(Corr, Norm, Pw),
    {pq_mp:from(Re, P), pq_mp:from(Im, P)}.

%% 2^m / m!
w(M, P) -> pq_mp:divide(pq_mp:from(1 bsl M, P), pq_mp:from(fact(M), P), P).

%% 2F1(-m, b; c; 2) = sum_(k=0..m) (-m)_k (b)_k / ((c)_k k!) 2^k
hyp(M, B, C, P) ->
    {S, _} = lists:foldl(fun(K, {Acc, T}) ->
                                 T2 = case K of
                                          0 -> T;
                                          _ -> cdiv(cmul2(cscale2(T, 2 * (K - 1 - M), P), cadd2(B, {pq_mp:from(K - 1, P), zero_mp()}, P), P),
                                                    cscale2(cadd2(C, {pq_mp:from(K - 1, P), zero_mp()}, P), K, P), P)
                                      end,
                                 {cadd2(Acc, T2, P), T2}
                         end, {{zero_mp(), zero_mp()}, {pq_mp:from(1, P), zero_mp()}}, lists:seq(0, M)),
    S.

%% sum_m Term(m) until the terms (after the first few) stay below 2^-(P + 8) of the sum for 8 consecutive m
sum_series(_, M, Acc, _, _, Quiet) when Quiet >= 8, M > 16 -> Acc;
sum_series(Term, M, Acc, Max, P, Quiet) ->
    T = Term(M),
    Acc2 = cadd2(Acc, T, P),
    Mag = max(cm(T), -(1 bsl 40)),
    Ref = max(cm(Acc2), Mag),
    sum_series(Term, M + 1, Acc2, Max, P, case Mag < Ref - P - 8 of true -> Quiet + 1; false -> 0 end).

cm({A, B}) -> max(mag(A), mag(B)).

cadd2({A, B}, {C, D}, P) -> {pq_mp:add(A, C, P), pq_mp:add(B, D, P)}.
csub2({A, B}, {C, D}, P) -> {pq_mp:sub(A, C, P), pq_mp:sub(B, D, P)}.
cmul2({A, B}, {C, D}, P) -> {pq_mp:sub(pq_mp:mul(A, C, P), pq_mp:mul(B, D, P), P), pq_mp:add(pq_mp:mul(A, D, P), pq_mp:mul(B, C, P), P)}.
cneg({A, B}) -> {pq_mp:neg(A), pq_mp:neg(B)}.
cconj2({A, B}) -> {A, pq_mp:neg(B)}.
cscale2({A, B}, K, P) -> {pq_mp:mul(A, pq_mp:from(K, P), P), pq_mp:mul(B, pq_mp:from(K, P), P)}.
cscale_r({A, B}, R, P) -> {pq_mp:mul(A, R, P), pq_mp:mul(B, R, P)}.
cdiv({A, B}, {C, D}, P) ->
    N = pq_mp:add(pq_mp:mul(C, C, P), pq_mp:mul(D, D, P), P),
    {pq_mp:divide(pq_mp:add(pq_mp:mul(A, C, P), pq_mp:mul(B, D, P), P), N, P),
     pq_mp:divide(pq_mp:sub(pq_mp:mul(B, C, P), pq_mp:mul(A, D, P), P), N, P)}.
cpow2(_, 0, P) -> {pq_mp:from(1, P), zero_mp()};
cpow2(Z, K, P) -> cmul2(Z, cpow2(Z, K - 1, P), P).
%% the principal square root
csqrt({X, Y}, P) ->
    R = pq_mp:sqrt(pq_mp:add(pq_mp:mul(X, X, P), pq_mp:mul(Y, Y, P), P), P),
    Re = pq_mp:sqrt(pq_mp:mul(pq_mp:add(R, X, P), pq_mp:from({1, 2}, P), P), P),
    Im0 = pq_mp:sqrt(pq_mp:mul(pq_mp:sub(R, X, P), pq_mp:from({1, 2}, P), P), P),
    {Re, case pq_mp:cmp(Y, zero_mp()) of -1 -> pq_mp:neg(Im0); _ -> Im0 end}.
zero_mp() -> {mp, 0, 0}.

