%% Equations of motion and spectra (decision record SY2, phase 3).
%%
%% - schrodinger/5: from a Lagrangian L(z, psi, psi') of one mode, its Euler-Lagrange equation a psi'' + b psi' + c psi
%%   = 0 (linear in psi), and the Liouville transformation psi = exp(-int h) phi, h = b / (2a), to Schrödinger form
%%   -phi'' + V phi = lambda phi, V = lambda - (c/a - h^2 - h'). Symbolic and exact (pq_symbolic).
%% - eigenvalues/3: the eigenvalues lambda of -phi'' + V phi = lambda phi, V a rational function of z, to any number
%%   of digits (pq_mp). The equation is written D phi'' + (lambda D - N) phi = 0 with V = N / D (polynomials).
%%   - The left end is a regular singular point at 0 (Frobenius series, the larger exponent) or a Dirichlet point.
%%   - The integration is the Taylor-series method: at each point the solution's Taylor coefficients by the
%%     recurrence of the polynomial equation, each step half the distance to the nearest zero of D.
%%   - At the right end: phi = 0 (Dirichlet), or alpha phi + beta phi' = 0 (Robin).
%%   - The eigenvalues are bracketed on a grid at low precision, then refined by the Illinois method at full
%%     precision.
%% - eigenvalues_collocation/3: the same eigenvalues by an independent method, for the left end at a regular singular
%%   point 0. phi = z^s u, u analytic, is discretised by Chebyshev-Lobatto collocation on [0, B]; the equation divided
%%   by z, (D/z) u'' + 2 s d u' + ((s(s-1) d - N)/z) u + lambda z d u = 0 (D = z^2 d), makes z = 0 its natural regularity
%%   condition; the generalised eigenproblem is solved by shifted inverse iteration in arbitrary precision.
-module(pq_ode).
-export([euler_lagrange/5, schrodinger/5, eigenvalues/3, eigenvalues_collocation/3, frobenius_solution/4, eigenvalues_sl/3]).

%%% Symbolic: Euler-Lagrange and Schrödinger form

%% d/dz (dL/dpsi') - dL/dpsi, with psi'' the symbol D2
-spec euler_lagrange(term(), atom(), atom(), atom(), atom()) -> {ok, term()} | {error, term()}.
euler_lagrange(L, Z, Psi, DPsi, D2) ->
    try
        Lp = d(L, DPsi),
        Lq = d(L, Psi),
        Total = can({add, [d(Lp, Z), {mul, [{sym, DPsi}, d(Lp, Psi)]}, {mul, [{sym, D2}, d(Lp, DPsi)]}]}),
        {ok, can({add, [Total, {mul, [-1, Lq]}]})}
    catch throw:{ode, Why} -> {error, Why}
    end.

%% #{equation => EL, potential => V, h => H}: V free of Lambda when c/a = Lambda + (free of Lambda)
-spec schrodinger(term(), atom(), atom(), atom(), atom()) -> {ok, map()} | {error, term()}.
schrodinger(L, Z, Psi, DPsi, Lambda) ->
    D2 = '$psi_zz',
    case euler_lagrange(L, Z, Psi, DPsi, D2) of
        {ok, EL} ->
            try
                A = d(EL, D2), B = d(EL, DPsi), C = d(EL, Psi),
                % linear and homogeneous in psi: nothing left of EL after the three
                Rest = can({add, [EL, {mul, [-1, {sym, D2}, A]}, {mul, [-1, {sym, DPsi}, B]}, {mul, [-1, {sym, Psi}, C]}]}),
                Rest =:= 0 orelse throw({ode, not_linear_homogeneous}),
                H = can({mul, [{rat, 1, 2}, B, {pow, A, -1}]}),
                Q = can({add, [{mul, [C, {pow, A, -1}]}, {mul, [-1, {pow, H, 2}]}, {mul, [-1, d(H, Z)]}]}),
                V = simplify(can({add, [{sym, Lambda}, {mul, [-1, Q]}]})),
                lists:member(Lambda, pq_symbolic:free_symbols(V)) andalso throw({ode, {weight_not_one, V}}),
                {ok, #{equation => EL, potential => V, h => H}}
            catch throw:{ode, Why} -> {error, Why}
            end;
        Err -> Err
    end.

d(E, X) ->
    case pq_symbolic:differentiate(E, X) of
        {ok, R} -> R;
        {error, Why} -> throw({ode, Why})
    end.

can(E) ->
    case pq_symbolic:canonicalise(E) of
        {ok, R} -> R;
        {error, Why} -> throw({ode, Why})
    end.

%% a rational function in normal form when it is one
simplify(V) ->
    case pq_symbolic:cancel(V) of
        {ok, R} -> R;
        {error, _} -> V
    end.

%%% Numeric: eigenvalues

%% Opts:
%%   params => #{Symbol => number}   the other symbols of V
%%   left => {frobenius, 0} | {dirichlet, A}
%%   right => {dirichlet, B} | {robin, B, Alpha, Beta}
%%   range => {Lo, Hi}, grid => K (bracketing steps, default 200), count => N (at most N eigenvalues)
%%   digits => D (default 30)
%% {ok, [Digits]} (strings, ascending)
-spec eigenvalues(term(), atom(), map()) -> {ok, [string()]} | {error, term()}.
eigenvalues(V0, Z, Opts) ->
    try
        Params = maps:get(params, Opts, #{}),
        V = case map_size(Params) of
                0 -> V0;
                _ -> {ok, S} = pq_symbolic:subs(V0, maps:map(fun(_, X) -> num_ast(X) end, Params)), S
            end,
        {N, Dd} = case pq_symbolic:numer_denom(V) of
                      {ok, ND} -> ND;
                      {error, _} -> throw({ode, potential_not_rational})
                  end,
        Np = coeffs(N, Z), Dp = coeffs(Dd, Z),
        Digits = maps:get(digits, Opts, 30),
        P = ceil(Digits * 3.3219280948873622) + 24,
        Prob = #{n => Np, d => Dp, roots => roots(Dp), left => maps:get(left, Opts, {frobenius, 0}),
                 right => maps:get(right, Opts)},
        {Lo, Hi} = maps:get(range, Opts),
        K = maps:get(grid, Opts, 200),
        Count = maps:get(count, Opts, 1000),
        Brackets = brackets(Prob, Lo, Hi, K, 80),
        Es = [refine(Prob, A, B, P) || {A, B} <- lists:sublist(Brackets, Count)],
        {ok, [pq_mp:to_string(E, Digits) || E <- Es]}
    catch throw:{ode, Why} -> {error, {symbolic_unsupported, #{eigenvalues => Why}}}
    end.

num_ast(N) when is_integer(N) -> N;
num_ast({N, D}) -> {rat, N, D};
num_ast({rat, _, _} = R) -> R.

%% a polynomial in Z (canonical expression) as exact coefficients [c0, c1, ...] ({N, D} each)
coeffs(E, Z) ->
    P = pq_sym_gcd:from_ast(E, [Z]),
    Deg = lists:max([K || [K] <- maps:keys(P)]),
    [q(maps:get([K], P, 0)) || K <- lists:seq(0, Deg)].

q(N) when is_integer(N) -> {N, 1};
q({N, D}) -> {N, D}.

%% the complex roots of D (floats, Durand-Kerner), for the step sizes
roots(Dp) ->
    % exact zeros at 0 first (a pole of V at the left end), then the others
    {Zeros, Rest} = lists:splitwith(fun({N, _}) -> N =:= 0 end, Dp),
    [{0.0, 0.0} || _ <- Zeros] ++ roots_rest(Rest).

roots_rest(Dp) ->
    Fs = [N / Dn || {N, Dn} <- Dp],
    Trim = lists:reverse(lists:dropwhile(fun(X) -> X == 0.0 end, lists:reverse(Fs))),
    case length(Trim) - 1 of
        Deg when Deg < 1 -> [];
        Deg ->
            Lead = lists:last(Trim),
            Mon = [X / Lead || X <- Trim],
            Init = [cpow({0.4, 0.9}, I) || I <- lists:seq(0, Deg - 1)],
            dk(Mon, Init, 500)
    end.

dk(_, Rs, 0) -> Rs;
dk(Mon, Rs, It) ->
    New = [csub(R, cdiv(peval(Mon, R), lists:foldl(fun(S, Acc) -> cmul(Acc, csub(R, S)) end, {1.0, 0.0},
                                                    [S || S <- Rs, S =/= R])))
           || R <- Rs],
    Delta = lists:max([cabs(csub(A, B)) || {A, B} <- lists:zip(New, Rs)]),
    case Delta < 1.0e-14 of
        true -> New;
        false -> dk(Mon, New, It - 1)
    end.

peval(C, X) -> lists:foldl(fun(Ci, Acc) -> cadd(cmul(Acc, X), {Ci, 0.0}) end, {0.0, 0.0}, lists:reverse(C)).
cadd({A, B}, {C, D}) -> {A + C, B + D}.
csub({A, B}, {C, D}) -> {A - C, B - D}.
cmul({A, B}, {C, D}) -> {A * C - B * D, A * D + B * C}.
cdiv({A, B}, {C, D}) -> M = C * C + D * D, {(A * C + B * D) / M, (B * C - A * D) / M}.
cabs({A, B}) -> math:sqrt(A * A + B * B).
cpow(_, 0) -> {1.0, 0.0};
cpow(Z, K) -> cmul(Z, cpow(Z, K - 1)).

%% sign changes of the shooting function on a grid of K steps, at precision P0
brackets(Prob, Lo, Hi, K, P0) ->
    Es = [Lo + (Hi - Lo) * I / K || I <- lists:seq(0, K)],
    Fs = [{E, pq_mp:to_float(shoot(Prob, pq_mp:from(E, P0), P0))} || E <- Es],
    [{A, B} || {{A, Fa}, {B, Fb}} <- lists:zip(lists:droplast(Fs), tl(Fs)), Fa * Fb < 0].

%% the Illinois method on [A, B] at precision P
refine(Prob, A0, B0, P) ->
    A = pq_mp:from(A0, P), B = pq_mp:from(B0, P),
    illinois(Prob, A, shoot(Prob, A, P), B, shoot(Prob, B, P), P, 0, none).

illinois(_, A, _, B, _, P, 200, _) -> pq_mp:ldexp(pq_mp:add(A, B, P), -1);
illinois(Prob, A, Fa, B, Fb, P, I, Side) ->
    % the secant point
    C = pq_mp:sub(B, pq_mp:divide(pq_mp:mul(Fb, pq_mp:sub(B, A, P), P), pq_mp:sub(Fb, Fa, P), P), P),
    Width = pq_mp:sub(B, A, P),
    case small(Width, C, P) of
        true -> C;
        false ->
            Fc = shoot(Prob, C, P),
            case sgn(Fc) of
                0 -> C;
                S ->
                    case S =:= sgn(Fb) of
                        true ->
                            % C replaces B; A kept: halve Fa if A was kept last time too
                            Fa2 = case Side of a -> pq_mp:ldexp(Fa, -1); _ -> Fa end,
                            illinois(Prob, A, Fa2, C, Fc, P, I + 1, a);
                        false ->
                            Fb2 = case Side of b -> pq_mp:ldexp(Fb, -1); _ -> Fb end,
                            illinois(Prob, C, Fc, B, Fb2, P, I + 1, b)
                    end
            end
    end.

small(W, C, P) ->
    {mp, Mw, Ew} = W, {mp, Mc, Ec} = C,
    Mw =:= 0 orelse (Mc =/= 0 andalso Ew + bitlen(abs(Mw)) < Ec + bitlen(abs(Mc)) - P + 8).

bitlen(0) -> 0;
bitlen(N) -> length(integer_to_list(N, 2)).

sgn({mp, M, _}) when M > 0 -> 1;
sgn({mp, M, _}) when M < 0 -> -1;
sgn(_) -> 0.

%% the boundary function at the right end for lambda = E (an mp), at precision P
shoot(#{right := Right} = Prob, E, P) ->
    W = P + 32,
    {Z0, Phi0, Dphi0} = start(Prob, E, W),
    Zb = case Right of {dirichlet, B} -> B; {robin, B, _, _} -> B end,
    {Phi, Dphi} = integrate(Prob, E, Z0, pq_mp:from(Zb, W), Phi0, Dphi0, W),
    case Right of
        {dirichlet, _} -> pq_mp:from(Phi, P);
        {robin, _, Al, Be} -> pq_mp:from(pq_mp:add(pq_mp:mul(pq_mp:from(Al, W), Phi, W), pq_mp:mul(pq_mp:from(Be, W), Dphi, W), W), P)
    end.

%% the solution and its derivative at the first point
start(Prob, E, W) ->
    {Z0, Phi, Dphi, _} = start(Prob, E, W, larger),
    {Z0, Phi, Dphi}.

%% {Z0, phi(Z0), phi'(Z0), int_0^Z0 phi^2} for the Frobenius solution of the larger or smaller exponent (c_0 = 1)
start(#{left := {dirichlet, A}}, _, W, _) -> {pq_mp:from(A, W), {mp, 0, 0}, {mp, 1, 0}, {mp, 0, 0}};
start(#{left := {frobenius, 0}, n := Np, d := Dp, roots := Rs}, E, W, Root) ->
    % D = z^2 d(z), d(0) /= 0; M = E D - N; phi = z^s sum c_k z^k
    case Dp of
        [{0, _}, {0, _}, {D2n, _} | _] when D2n =/= 0 -> ok;
        _ -> throw({ode, {not_a_double_pole_at_0, Dp}})
    end,
    Dc = [pq_mp:from(X, W) || X <- lists:nthtail(2, Dp)],          % d_j = D_(j+2)
    Mc = mcoeffs(Np, Dp, E, W),
    [M0 | _] = Mc ++ [{mp, 0, 0}],
    D0 = hd(Dc),
    % the indicial equation d_0 s (s - 1) + M_0 = 0: the larger root
    Disc = pq_mp:sub({mp, 1, 0}, pq_mp:ldexp(pq_mp:divide(M0, D0, W), 2), W),
    pq_mp:to_float(Disc) >= 0 orelse throw({ode, complex_indicial_exponent}),
    Sq = pq_mp:sqrt(Disc, W),
    S = pq_mp:ldexp(case Root of larger -> pq_mp:add({mp, 1, 0}, Sq, W); smaller -> pq_mp:sub({mp, 1, 0}, Sq, W) end, -1),
    Rho = lists:min([cabs(R) || R <- Rs, cabs(R) > 1.0e-12] ++ [1.0e9]),
    Z0f = min(0.25 * Rho, 0.25),
    Z0 = pq_mp:from(Z0f, W),
    Cs = frob(Dc, Mc, S, W, Z0),
    % phi = z^s sum c_k z^k, phi' = z^(s-1) sum (k + s) c_k z^k
    Zs = pq_mp:exp(pq_mp:mul(S, pq_mp:log(Z0, W), W), W),
    {Sum, DSum} = lists:foldl(fun({K, C}, {A, B}) ->
                                      Zk = ipow(Z0, K, W),
                                      T = pq_mp:mul(C, Zk, W),
                                      {pq_mp:add(A, T, W), pq_mp:add(B, pq_mp:mul(pq_mp:add(S, pq_mp:from(K, W), W), T, W), W)}
                              end, {{mp, 0, 0}, {mp, 0, 0}}, lists:enumerate(0, Cs)),
    % int_0^Z0 phi^2 = sum_m (sum_(i+j=m) c_i c_j) Z0^(2s+m+1) / (2s+m+1)
    Ct = list_to_tuple(Cs),
    Nc = tuple_size(Ct),
    Int = lists:foldl(fun(_, undefined) -> undefined;
                         (M, Acc) ->
                              Cm = lists:foldl(fun(I, A) -> pq_mp:add(A, pq_mp:mul(element(I + 1, Ct), element(M - I + 1, Ct), W), W) end,
                                               {mp, 0, 0}, lists:seq(max(0, M - Nc + 1), min(M, Nc - 1))),
                              Ex = pq_mp:add(pq_mp:ldexp(S, 1), pq_mp:from(M + 1, W), W),
                              % phi^2 ~ z^(2s+m) with 2s + m + 1 =< 0 (the smaller exponent): not integrable at 0
                              case pq_mp:to_float(Ex) =< 1.0e-30 of
                                  true -> undefined;
                                  false ->
                                      Zp = pq_mp:exp(pq_mp:mul(Ex, pq_mp:log(Z0, W), W), W),
                                      pq_mp:add(Acc, pq_mp:divide(pq_mp:mul(Cm, Zp, W), Ex, W), W)
                              end
                      end, {mp, 0, 0}, lists:seq(0, 2 * Nc - 2)),
    {Z0, pq_mp:mul(Zs, Sum, W), pq_mp:divide(pq_mp:mul(Zs, DSum, W), Z0, W), Int}.

%% the coefficients of M = E D - N
mcoeffs(Np, Dp, E, W) ->
    L = max(length(Np), length(Dp)),
    [pq_mp:sub(pq_mp:mul(E, pq_mp:from(at(Dp, J), W), W), pq_mp:from(at(Np, J), W), W) || J <- lists:seq(0, L - 1)].

at(L, J) when J < length(L) -> lists:nth(J + 1, L);
at(_, _) -> {0, 1}.

%% the Frobenius coefficients c_0 = 1, c_m [d_0 (m+s)(m+s-1) + M_0] = -sum_(j>=1) [d_j (m-j+s)(m-j+s-1) + M_j] c_(m-j),
%% until the terms at Z0 are negligible
frob(Dc, Mc, S, W, Z0) -> frob(Dc, Mc, S, W, Z0, 1, [{mp, 1, 0}], {mp, 1, 0}).

frob(Dc, Mc, S, W, Z0, M, Acc, Zm) ->
    Ms = pq_mp:add(S, pq_mp:from(M, W), W),
    Den = pq_mp:add(pq_mp:mul(hd(Dc), pq_mp:mul(Ms, pq_mp:sub(Ms, {mp, 1, 0}, W), W), W), hd(Mc ++ [{mp, 0, 0}]), W),
    Rev = Acc,   % c_(m-1), c_(m-2), ..., c_0
    Num = lists:foldl(fun({J, C}, Sum) ->
                              K = pq_mp:add(S, pq_mp:from(M - J, W), W),
                              Dj = idx(Dc, J), Mj = idx(Mc, J),
                              F = pq_mp:add(pq_mp:mul(Dj, pq_mp:mul(K, pq_mp:sub(K, {mp, 1, 0}, W), W), W), Mj, W),
                              pq_mp:add(Sum, pq_mp:mul(F, C, W), W)
                      end, {mp, 0, 0}, lists:zip(lists:seq(1, M), Rev)),
    Cm = case mag_rel(Den, Ms, W) of
             true ->
                 % a resonance (the exponents differ by the integer M): no logarithm only if the right side vanishes
                 % too; the free coefficient is then 0
                 case negligible(Num, W - 40) of
                     true -> {mp, 0, 0};
                     false -> throw({ode, logarithmic_frobenius_solution})
                 end;
             false -> pq_mp:neg(pq_mp:divide(Num, Den, W))
         end,
    Zm2 = pq_mp:mul(Zm, Z0, W),
    T = pq_mp:mul(Cm, Zm2, W),
    Acc2 = [Cm | Acc],
    case M > 8 andalso negligible(T, W) andalso negligible(pq_mp:mul(hd(tl(Acc2)), Zm, W), W) of
        true -> lists:reverse(Acc2);
        false when M > 4000 -> throw({ode, frobenius_not_converging});
        false -> frob(Dc, Mc, S, W, Z0, M + 1, Acc2, Zm2)
    end.

%% Den negligible against the size of its terms
mag_rel({mp, 0, _}, _, _) -> true;
mag_rel(Den, Ms, W) -> mag(Den) < mag(pq_mp:mul(Ms, Ms, W)) - W + 40.

idx(L, J) when J < length(L) -> lists:nth(J + 1, L);
idx(_, _) -> {mp, 0, 0}.

negligible({mp, 0, _}, _) -> true;
negligible({mp, M, E}, W) -> E + bitlen(abs(M)) < -W - 4.

ipow(_, 0, _) -> {mp, 1, 0};
ipow(X, K, W) -> pq_mp:mul(X, ipow(X, K - 1, W), W).

%% the Taylor-series method from Z0 to Zb
integrate(Prob, E, Z0, Zb, Phi, Dphi, W) ->
    case pq_mp:cmp(Z0, Zb) >= 0 of
        true -> {Phi, Dphi};
        false ->
            {H, Last} = step_size(Prob, Z0, Zb, W),
            {Phi2, Dphi2} = taylor_step(Prob, E, Z0, H, Phi, Dphi, W),
            Z1 = case Last of true -> Zb; false -> pq_mp:add(Z0, H, W) end,
            integrate(Prob, E, Z1, Zb, Phi2, Dphi2, W)
    end.

%% the same with int phi^2 accumulated: each step's Taylor polynomial squared, integrated over [0, h]
integrate_int(Prob, E, Z0, Zb, Phi, Dphi, undefined, W) ->
    {P2, D2} = integrate(Prob, E, Z0, Zb, Phi, Dphi, W),
    {P2, D2, undefined};
integrate_int(Prob, E, Z0, Zb, Phi, Dphi, Int, W) ->
    case pq_mp:cmp(Z0, Zb) >= 0 of
        true -> {Phi, Dphi, Int};
        false ->
            {H, Last} = step_size(Prob, Z0, Zb, W),
            {Phi2, Dphi2} = taylor_step(Prob, E, Z0, H, Phi, Dphi, W),
            Coef = get({?MODULE, last_coefs}),
            K = maps:size(Coef),
            Int2 = lists:foldl(fun(M, Acc) ->
                                       Cm = lists:foldl(fun(I, A) -> pq_mp:add(A, pq_mp:mul(maps:get(I, Coef), maps:get(M - I, Coef), W), W) end,
                                                        {mp, 0, 0}, lists:seq(max(0, M - K + 1), min(M, K - 1))),
                                       pq_mp:add(Acc, pq_mp:divide(pq_mp:mul(Cm, ipow(H, M + 1, W), W), pq_mp:from(M + 1, W), W), W)
                               end, Int, lists:seq(0, 2 * K - 2)),
            Z1 = case Last of true -> Zb; false -> pq_mp:add(Z0, H, W) end,
            integrate_int(Prob, E, Z1, Zb, Phi2, Dphi2, Int2, W)
    end.

%% half the distance to the nearest zero of D, at most the rest of the way
step_size(#{roots := Rs}, Z0, Zb, W) ->
    Zf = pq_mp:to_float(Z0),
    Rho = lists:min([cabs(csub({Zf, 0.0}, R)) || R <- Rs] ++ [1.0e9]),
    H0 = min(0.5 * Rho, 1.0),
    Rest = pq_mp:sub(Zb, Z0, W),
    case pq_mp:to_float(Rest) =< H0 of
        true -> {Rest, true};
        false -> {pq_mp:from(H0, W), false}
    end.

%% phi(a + h), phi'(a + h) from phi(a), phi'(a): the Taylor coefficients by
%%   phi_(m+2) = -[sum_(j=1..m) D_j (m-j+2)(m-j+1) phi_(m-j+2) + sum_(j=0..m) M_j phi_(m-j)] / (D_0 (m+2)(m+1)),
%% D_j, M_j the coefficients of D(a + t), M(a + t)
taylor_step(#{n := Np, d := Dp}, E, A, H, Phi, Dphi, W) ->
    Ds = shift_poly([pq_mp:from(X, W) || X <- Dp], A, W),
    Ms = shift_poly(mcoeffs(Np, Dp, E, W), A, W),
    D0 = hd(Ds),
    Dt = list_to_tuple(Ds), Mt = list_to_tuple(Ms),
    Coef = taylor_coefs(Dt, Mt, D0, W, H, #{0 => Phi, 1 => Dphi}, 0),
    {Sum, DSum} = maps:fold(fun(K, C, {S, Ds1}) ->
                                    Hk = ipow(H, K, W),
                                    S2 = pq_mp:add(S, pq_mp:mul(C, Hk, W), W),
                                    Ds2 = case K of
                                              0 -> Ds1;
                                              _ -> pq_mp:add(Ds1, pq_mp:mul(pq_mp:mul(pq_mp:from(K, W), C, W), ipow(H, K - 1, W), W), W)
                                          end,
                                    {S2, Ds2}
                            end, {{mp, 0, 0}, {mp, 0, 0}}, Coef),
    put({?MODULE, last_coefs}, Coef),
    {Sum, DSum}.

taylor_coefs(Dt, Mt, D0, W, H, Acc, M) ->
    Num = lists:foldl(fun(J, S) ->
                              case J >= 1 andalso J =< tuple_size(Dt) - 1 of
                                  true ->
                                      K = M - J + 2,
                                      pq_mp:add(S, pq_mp:mul(pq_mp:mul(element(J + 1, Dt), pq_mp:from(K * (K - 1), W), W), maps:get(K, Acc), W), W);
                                  false -> S
                              end
                      end, {mp, 0, 0}, lists:seq(1, M)),
    Num2 = lists:foldl(fun(J, S) ->
                               case J =< tuple_size(Mt) - 1 of
                                   true -> pq_mp:add(S, pq_mp:mul(element(J + 1, Mt), maps:get(M - J, Acc), W), W);
                                   false -> S
                               end
                       end, Num, lists:seq(0, M)),
    C = pq_mp:neg(pq_mp:divide(Num2, pq_mp:mul(D0, pq_mp:from((M + 2) * (M + 1), W), W), W)),
    Acc2 = Acc#{M + 2 => C},
    % stop when two successive terms, times h^k, are negligible relative to the sum's scale
    Scale = max(mag(maps:get(0, Acc)), mag(pq_mp:mul(maps:get(1, Acc), H, W))),
    T1 = mag(pq_mp:mul(C, ipow(H, M + 2, W), W)),
    T0 = mag(pq_mp:mul(maps:get(M + 1, Acc2), ipow(H, M + 1, W), W)),
    case M > 6 andalso T1 < Scale - W - 4 andalso T0 < Scale - W - 4 of
        true -> Acc2;
        false when M > 3000 -> throw({ode, taylor_not_converging});
        false -> taylor_coefs(Dt, Mt, D0, W, H, Acc2, M + 1)
    end.

mag({mp, 0, _}) -> -1 bsl 30;
mag({mp, M, E}) -> E + bitlen(abs(M)).

%% the coefficients of p(a + t) from those of p(z): Horner's synthetic division repeated
shift_poly(Cs, A, W) ->
    N = length(Cs),
    shift_poly(lists:reverse(Cs), A, W, N, []).

shift_poly(_, _, _, 0, Acc) -> lists:reverse(Acc);
shift_poly(Rev, A, W, N, Acc) ->
    % divide by (t - a) style: evaluate and deflate; Rev from the highest coefficient
    {Quot, Rem} = synth(Rev, A, W),
    shift_poly(Quot, A, W, N - 1, [Rem | Acc]).

synth([H | T], A, W) ->
    {Q, Last} = lists:foldl(fun(C, {Qs, Prev}) ->
                                    V = pq_mp:add(C, pq_mp:mul(Prev, A, W), W),
                                    {[Prev | Qs], V}
                            end, {[], H}, T),
    {lists:reverse(Q), Last}.

%%% Numeric: eigenvalues by Chebyshev collocation and inverse iteration (independent of the shooting method)

%% Opts as for eigenvalues/3 (left must be {frobenius, 0}), plus points => N (default 64); the starting shifts are
%% the midpoints of the low-precision brackets (only where to look: the values are the collocation's)
-spec eigenvalues_collocation(term(), atom(), map()) -> {ok, [string()]} | {error, term()}.
eigenvalues_collocation(V0, Z, Opts) ->
    try
        Params = maps:get(params, Opts, #{}),
        V = case map_size(Params) of
                0 -> V0;
                _ -> {ok, S0} = pq_symbolic:subs(V0, maps:map(fun(_, X) -> num_ast(X) end, Params)), S0
            end,
        {ok, {N, Dd}} = pq_symbolic:numer_denom(V),
        Np = coeffs(N, Z), Dp = coeffs(Dd, Z),
        Digits = maps:get(digits, Opts, 30),
        P = ceil(Digits * 3.3219280948873622) + 24,
        W = P + 32,
        Prob = #{n => Np, d => Dp, roots => roots(Dp), left => {frobenius, 0}, right => maps:get(right, Opts)},
        {Lo, Hi} = maps:get(range, Opts),
        Count = maps:get(count, Opts, 1000),
        Shifts = [(A + B) / 2 || {A, B} <- lists:sublist(brackets(Prob, Lo, Hi, maps:get(grid, Opts, 200), 80), Count)],
        {A0, A1} = collocation(Np, Dp, maps:get(right, Opts), maps:get(points, Opts, 64), W),
        Es = [inverse_iteration(A0, A1, pq_mp:from(S, W), W, P) || S <- Shifts],
        {ok, [pq_mp:to_string(E, Digits) || E <- Es]}
    catch throw:{ode, Why} -> {error, {symbolic_unsupported, #{eigenvalues => Why}}}
    end.

%% the matrices A0, A1 (lists of rows of mp) of the collocation: A0 u + lambda A1 u = 0
collocation(Np, Dp, Right, N, W) ->
    case Dp of
        [{0, _}, {0, _}, {D2n, _} | _] when D2n =/= 0 -> ok;
        _ -> throw({ode, {not_a_double_pole_at_0, Dp}})
    end,
    B = case Right of {dirichlet, Bb} -> Bb; {robin, Bb, _, _} -> Bb end,
    L = pq_mp:from(B, W),
    Pi = pq_mp:pi(W),
    % Chebyshev-Lobatto points x_j = cos(pi j / N) on [-1, 1]; z_j = L (1 - x_j) / 2, from 0 to L
    Xs = [pq_mp:cos(pq_mp:divide(pq_mp:mul(pq_mp:from(J, W), Pi, W), pq_mp:from(N, W), W), W) || J <- lists:seq(0, N)],
    Zs = [pq_mp:ldexp(pq_mp:mul(L, pq_mp:sub({mp, 1, 0}, X, W), W), -1) || X <- Xs],
    Dx = cheb_matrix(Xs, N, W),
    % d/dz = -(2/L) d/dx
    F = pq_mp:neg(pq_mp:divide({mp, 2, 0}, L, W)),
    D1 = [[pq_mp:mul(F, E, W) || E <- Row] || Row <- Dx],
    D2 = matmul(D1, D1, W),
    % the indicial exponent: d0 s (s - 1) - N0 = 0, the larger root
    Dc = [pq_mp:from(X, W) || X <- lists:nthtail(2, Dp)],      % d(z) = D / z^2
    Nc = [pq_mp:from(X, W) || X <- Np],
    D0 = hd(Dc), N0 = hd(Nc),
    Disc = pq_mp:add({mp, 1, 0}, pq_mp:ldexp(pq_mp:divide(N0, D0, W), 2), W),
    S = pq_mp:ldexp(pq_mp:add({mp, 1, 0}, pq_mp:sqrt(Disc, W), W), -1),
    % polynomials: D/z = z d, 2 s d, (s(s-1) d - N)/z, z d
    Ss1 = pq_mp:mul(S, pq_mp:sub(S, {mp, 1, 0}, W), W),
    Q0 = padd([pq_mp:mul(Ss1, X, W) || X <- Dc], [pq_mp:neg(X) || X <- Nc], W),   % s(s-1) d - N, Q0(0) = 0
    Qz = tl(Q0),                                                                 % divided by z
    Zd = [{mp, 0, 0} | Dc],                                                      % z d
    TwoSd = [pq_mp:mul(pq_mp:ldexp(S, 1), X, W) || X <- Dc],
    Ev = fun(Poly, Zj) -> horner(Poly, Zj, W) end,
    Rows = lists:zip3(lists:seq(0, N), Zs, lists:zip(D1, D2)),
    A0 = [case J of
              N -> bc_row(Right, N, S, L, lists:nth(N + 1, D1), W);
              _ ->
                  C2 = Ev(Zd, Zj), C1 = Ev(TwoSd, Zj), C0 = Ev(Qz, Zj),
                  [pq_mp:add(pq_mp:add(pq_mp:mul(C2, E2, W), pq_mp:mul(C1, E1, W), W),
                             case K of J -> C0; _ -> {mp, 0, 0} end, W)
                   || {K, E1, E2} <- lists:zip3(lists:seq(0, N), R1, R2)]
          end || {J, Zj, {R1, R2}} <- Rows],
    A1 = [case J of
              N -> [{mp, 0, 0} || _ <- lists:seq(0, N)];
              _ -> C = Ev(Zd, Zj), [case K of J -> C; _ -> {mp, 0, 0} end || K <- lists:seq(0, N)]
          end || {J, Zj, _} <- Rows],
    {A0, A1}.

bc_row({dirichlet, _}, N, _, _, _, _) -> [case K of N -> {mp, 1, 0}; _ -> {mp, 0, 0} end || K <- lists:seq(0, N)];
bc_row({robin, _, Al, Be}, N, S, L, D1row, W) ->
    % alpha phi + beta phi' = z^s [ (alpha + beta s / z) u + beta u' ] at z = L
    A = pq_mp:add(pq_mp:from(Al, W), pq_mp:divide(pq_mp:mul(pq_mp:from(Be, W), S, W), L, W), W),
    Bm = pq_mp:from(Be, W),
    [pq_mp:add(pq_mp:mul(Bm, E, W), case K of N -> A; _ -> {mp, 0, 0} end, W) || {K, E} <- lists:enumerate(0, D1row)].

%% the Chebyshev differentiation matrix on the points Xs (Trefethen, Spectral Methods in MATLAB, cheb.m), the
%% diagonal as minus the row sums
cheb_matrix(Xs, N, W) ->
    C = fun(I) -> (case I of 0 -> 2; N -> 2; _ -> 1 end) * (case I rem 2 of 0 -> 1; 1 -> -1 end) end,
    Xt = list_to_tuple(Xs),
    Off = [[case I of
                J -> {mp, 0, 0};
                _ -> pq_mp:divide(pq_mp:from({C(I) * sign(C(J)), abs(C(J))}, W),
                                  pq_mp:sub(element(I + 1, Xt), element(J + 1, Xt), W), W)
            end || J <- lists:seq(0, N)] || I <- lists:seq(0, N)],
    [begin
         Sum = lists:foldl(fun(E, A) -> pq_mp:add(A, E, W) end, {mp, 0, 0}, Row),
         setnth(I + 1, Row, pq_mp:neg(Sum))
     end || {I, Row} <- lists:enumerate(0, Off)].

sign(X) when X < 0 -> -1;
sign(_) -> 1.

setnth(1, [_ | T], X) -> [X | T];
setnth(K, [H | T], X) -> [H | setnth(K - 1, T, X)].

matmul(A, B, W) ->
    Bt = transpose(B),
    [[dot(Row, Col, W) || Col <- Bt] || Row <- A].

transpose([[] | _]) -> [];
transpose(M) -> [[hd(R) || R <- M] | transpose([tl(R) || R <- M])].

dot(A, B, W) -> lists:foldl(fun({X, Y}, S) -> pq_mp:add(S, pq_mp:mul(X, Y, W), W) end, {mp, 0, 0}, lists:zip(A, B)).

padd(A, B, W) ->
    N = max(length(A), length(B)),
    Pad = fun(L) -> L ++ [{mp, 0, 0} || _ <- lists:seq(1, N - length(L))] end,
    [pq_mp:add(X, Y, W) || {X, Y} <- lists:zip(Pad(A), Pad(B))].

horner(Poly, X, W) -> lists:foldl(fun(C, Acc) -> pq_mp:add(pq_mp:mul(Acc, X, W), C, W) end, {mp, 0, 0}, lists:reverse(Poly)).

%% (A0 + sigma A1) y = A1 x; lambda = sigma - x_k / y_k (k the largest component); sigma := lambda, repeated
inverse_iteration(A0, A1, Sigma0, W, P) ->
    N = length(A0),
    X0 = [{mp, 1, 0} || _ <- lists:seq(1, N)],
    inv_it(A0, A1, Sigma0, X0, W, P, 0).

inv_it(_, _, Sigma, _, _, _, 60) -> Sigma;
inv_it(A0, A1, Sigma, X, W, P, I) ->
    M = [[pq_mp:add(E0, pq_mp:mul(Sigma, E1, W), W) || {E0, E1} <- lists:zip(R0, R1)] || {R0, R1} <- lists:zip(A0, A1)],
    Rhs = [dot(R1, X, W) || R1 <- A1],
    Y = solve(M, Rhs, W),
    {K, _} = lists:foldl(fun({J, Yj}, {Kb, Mb}) -> case abs(pq_mp:to_float(Yj)) > Mb of true -> {J, abs(pq_mp:to_float(Yj))}; false -> {Kb, Mb} end end,
                         {1, -1.0}, lists:enumerate(Y)),
    Lambda = pq_mp:sub(Sigma, pq_mp:divide(lists:nth(K, X), lists:nth(K, Y), W), W),
    Yk = lists:nth(K, Y),
    X2 = [pq_mp:divide(Yj, Yk, W) || Yj <- Y],
    Change = pq_mp:sub(Lambda, Sigma, W),
    case small(Change, Lambda, P) andalso I > 1 of
        true -> Lambda;
        false -> inv_it(A0, A1, Lambda, X2, W, P, I + 1)
    end.

%% Gaussian elimination with partial pivoting
solve(M, B, W) ->
    Aug = [Row ++ [Bi] || {Row, Bi} <- lists:zip(M, B)],
    back(elim(Aug, W), W).

elim([], _) -> [];
elim(Rows, W) ->
    {Piv, Rest} = pivot(Rows),
    [P0 | Pt] = Piv,
    Red = [begin
               [R0 | Rt] = R,
               F = pq_mp:divide(R0, P0, W),
               [pq_mp:sub(E, pq_mp:mul(F, Pe, W), W) || {E, Pe} <- lists:zip(Rt, Pt)]
           end || R <- Rest],
    [Piv | elim(Red, W)].

pivot(Rows) ->
    {_, I} = lists:max([{abs(pq_mp:to_float(hd(R))), I} || {I, R} <- lists:enumerate(Rows)]),
    Piv = lists:nth(I, Rows),
    {Piv, lists:sublist(Rows, I - 1) ++ lists:nthtail(I, Rows)}.

%% the triangular rows (each one shorter, the right-hand side last) back-substituted
back(Tri, W) ->
    lists:foldr(fun(Row, Xs) ->
                        [Diag | Rest] = Row,
                        {Coefs, [Rhs]} = lists:split(length(Rest) - 1, Rest),
                        S = lists:foldl(fun({C, X}, A) -> pq_mp:add(A, pq_mp:mul(C, X, W), W) end, {mp, 0, 0}, lists:zip(Coefs, Xs)),
                        [pq_mp:divide(pq_mp:sub(Rhs, S, W), Diag, W) | Xs]
                end, [], Tri).

%%% The Frobenius solution itself

%% the solution phi = z^s (1 + ...) of -phi'' + V phi = E phi from the regular singular point 0, s the larger or the
%% smaller exponent (Opts exponent => larger | smaller; a smaller one with a logarithm is refused), at the right end B
%% (Opts b => B): {ok, #{phi, dphi, norm}} as decimal strings, norm = int_0^B phi^2 dz (undefined when phi^2 is not
%% integrable at 0). E an
%% integer, a fraction {N, D}, or a decimal string.
-spec frobenius_solution(term(), atom(), term(), map()) -> {ok, map()} | {error, term()}.
frobenius_solution(V0, Z, E0, Opts) ->
    try
        Params = maps:get(params, Opts, #{}),
        V = case map_size(Params) of
                0 -> V0;
                _ -> {ok, S0} = pq_symbolic:subs(V0, maps:map(fun(_, X) -> num_ast(X) end, Params)), S0
            end,
        {ok, {N, Dd}} = pq_symbolic:numer_denom(V),
        Np = coeffs(N, Z), Dp = coeffs(Dd, Z),
        Digits = maps:get(digits, Opts, 30),
        P = ceil(Digits * 3.3219280948873622) + 24,
        W = P + 32,
        E = case E0 of
                S when is_list(S) -> pq_mp:from_string(S, W);
                _ -> pq_mp:from(E0, W)
            end,
        B = maps:get(b, Opts),
        Prob = #{n => Np, d => Dp, roots => roots(Dp), left => {frobenius, 0}},
        {Z0, Phi0, Dphi0, Int0} = start(Prob, E, W, maps:get(exponent, Opts, larger)),
        {Phi, Dphi, Int} = integrate_int(Prob, E, Z0, pq_mp:from(B, W), Phi0, Dphi0, Int0, W),
        erase({?MODULE, last_coefs}),
        {ok, #{phi => pq_mp:to_string(Phi, Digits), dphi => pq_mp:to_string(Dphi, Digits),
               norm => case Int of undefined -> undefined; _ -> pq_mp:to_string(Int, Digits) end}}
    catch throw:{ode, Why} -> {error, {symbolic_unsupported, #{frobenius => Why}}}
    end.

%%% Sturm-Liouville problems with two regular singular ends (one of them at infinity)

%% The eigenvalues lambda of p2 psi'' + p1 psi' + (p0 + lambda q) psi = 0, p2, p1, p0, q polynomials in X (Coefs:
%% #{p2, p1, p0, q} expressions), on [A, infinity):
%%   left => {A, larger | smaller}: A a simple zero of p2 (a regular singular point), the Frobenius exponent there;
%%   right => larger | smaller: the exponent at infinity, in w = 1/x (where the equation is brought to the same form);
%%   match => M (A < M): the two Frobenius solutions, integrated by the Taylor-series method to x = M, and their
%%   Wronskian, which vanishes at the eigenvalues;
%%   range, grid, count, digits, params as for eigenvalues/3.
%% {ok, [Digits]}.
-spec eigenvalues_sl(map(), atom(), map()) -> {ok, [string()]} | {error, term()}.
eigenvalues_sl(Coefs0, X, Opts) ->
    try
        Params = maps:get(params, Opts, #{}),
        Sub = fun(E) -> case map_size(Params) of
                            0 -> E;
                            _ -> {ok, S1} = pq_symbolic:subs(E, maps:map(fun(_, V) -> num_ast(V) end, Params)), S1
                        end end,
        Coefs = maps:map(fun(_, E) -> Sub(E) end, Coefs0),
        Digits = maps:get(digits, Opts, 25),
        P = ceil(Digits * 3.3219280948873622) + 24,
        {A, RootL} = maps:get(left, Opts),
        M = maps:get(match, Opts),
        Left = sl_problem(Coefs, X, A),
        Right = sl_problem(at_infinity(Coefs, X), '$w', 0),
        Prob = #{left => Left, right => Right, root_l => RootL, root_r => maps:get(right, Opts), a => A, m => M},
        {Lo, Hi} = maps:get(range, Opts),
        K = maps:get(grid, Opts, 200),
        Es = [Lo + (Hi - Lo) * I / K || I <- lists:seq(0, K)],
        Fs = [{E, pq_mp:to_float(sl_wronskian(Prob, pq_mp:from(E, 90), 90))} || E <- Es],
        Brackets = [{Ea, Eb} || {{Ea, Fa}, {Eb, Fb}} <- lists:zip(lists:droplast(Fs), tl(Fs)), Fa * Fb < 0],
        Shoot = fun(Prb, E, Pp) -> sl_wronskian(Prb, E, Pp) end,
        Out = [refine_with(Shoot, Prob, Ea, Eb, P) || {Ea, Eb} <- lists:sublist(Brackets, maps:get(count, Opts, 1000))],
        {ok, [pq_mp:to_string(E, Digits) || E <- Out]}
    catch throw:{ode, Why} -> {error, {symbolic_unsupported, #{eigenvalues_sl => Why}}}
    end.

%% the coefficients as exact polynomials in X: #{p2, p1, p0, q} => [{N, D}]
sl_problem(Coefs, X, A) ->
    C = fun(K) -> coeffs_or_zero(maps:get(K, Coefs), X) end,
    #{p2 => C(p2), p1 => C(p1), p0 => C(p0), q => C(q), a => A}.

coeffs_or_zero(E, X) ->
    {ok, Ce} = pq_symbolic:canonicalise(E),
    case Ce of
        0 -> [{0, 1}];
        _ -> coeffs(Ce, X)
    end.

%% the equation in w = 1/x: psi_x = -w^2 psi_w, psi_xx = w^4 psi_ww + 2 w^3 psi_w; the coefficients multiplied by the
%% power of w that makes them polynomials
at_infinity(Coefs, X) ->
    W = {sym, '$w'},
    At = fun(K) -> {ok, S} = pq_symbolic:subs(maps:get(K, Coefs), #{X => {pow, W, -1}}), S end,
    E2 = {mul, [At(p2), {pow, W, 4}]},
    E1 = {add, [{mul, [2, At(p2), {pow, W, 3}]}, {mul, [-1, At(p1), {pow, W, 2}]}]},
    Es = #{p2 => E2, p1 => E1, p0 => At(p0), q => At(q)},
    % the largest power of 1/w in any of them
    Deg = lists:max([neg_degree(E) || E <- maps:values(Es)]),
    maps:map(fun(_, E) -> {ok, R} = pq_symbolic:cancel({mul, [E, {pow, W, Deg}]}), R end, Es).

neg_degree(E) ->
    {ok, {_, D}} = pq_symbolic:numer_denom(E),
    case D of
        1 -> 0;
        {sym, _} -> 1;
        {pow, {sym, _}, K} -> K;
        {mul, [_, {pow, {sym, _}, K}]} -> K;
        {mul, [_, {sym, _}]} -> 1;
        _ -> throw({ode, {denominator_not_a_power, D}})
    end.

%% the Wronskian at the match point x = M of the left and right Frobenius solutions, for lambda = E
sl_wronskian(#{left := L, right := R, root_l := Rl, root_r := Rr, m := M}, E, P) ->
    W = P + 32,
    {PhiL, DPhiL} = sl_solve(L, Rl, E, pq_mp:from(M, W), W),
    {PhiR, DPhiRw} = sl_solve(R, Rr, E, pq_mp:divide({mp, 1, 0}, pq_mp:from(M, W), W), W),
    % d/dx = -w^2 d/dw
    Wm = pq_mp:divide({mp, 1, 0}, pq_mp:from(M, W), W),
    DPhiR = pq_mp:neg(pq_mp:mul(pq_mp:mul(Wm, Wm, W), DPhiRw, W)),
    % normalised so that the sign does not depend on the solutions' sizes
    Wr = pq_mp:sub(pq_mp:mul(PhiL, DPhiR, W), pq_mp:mul(PhiR, DPhiL, W), W),
    Scale = pq_mp:sqrt(pq_mp:add(pq_mp:mul(PhiL, PhiL, W), pq_mp:mul(DPhiL, DPhiL, W), W), W),
    Scale2 = pq_mp:sqrt(pq_mp:add(pq_mp:mul(PhiR, PhiR, W), pq_mp:mul(DPhiR, DPhiR, W), W), W),
    pq_mp:from(pq_mp:divide(Wr, pq_mp:mul(Scale, Scale2, W), W), P).

%% the Frobenius solution at the problem's point a (a simple zero of p2), integrated to Xe (> a): {psi, psi'}
sl_solve(#{p2 := P2, p1 := P1, p0 := P0, q := Q, a := A}, Root, E, Xe, W) ->
    Am = pq_mp:from(A, W),
    R = padd([pq_mp:from(C, W) || C <- P0], [pq_mp:mul(E, pq_mp:from(C, W), W) || C <- Q], W),
    P2m = [pq_mp:from(C, W) || C <- P2], P1m = [pq_mp:from(C, W) || C <- P1],
    % shifted to t = x - a
    S2 = shift_poly(P2m, Am, W), S1 = shift_poly(P1m, Am, W), S0 = shift_poly(R, Am, W),
    is_negligible(hd(S2), W) orelse throw({ode, {not_a_singular_point, A}}),
    p2_ok = case tl(S2) of [] -> throw({ode, degenerate_p2}); [L | _] -> case is_negligible(L, W) of true -> throw({ode, not_a_simple_zero}); false -> p2_ok end end,
    D2 = tl(S2),                                   % p2 = t d2(t)
    % the indicial equation d2_0 s (s - 1) + p1_0 s = 0: s = 0 or 1 - p1_0 / d2_0
    S1r = pq_mp:sub({mp, 1, 0}, pq_mp:divide(hd(S1 ++ [{mp, 0, 0}]), hd(D2), W), W),
    {Sl, Ss} = case pq_mp:to_float(S1r) >= 0 of true -> {S1r, {mp, 0, 0}}; false -> {{mp, 0, 0}, S1r} end,
    Sx = case Root of larger -> Sl; smaller -> Ss end,
    % the other zeros of p2 limit the series' radius
    Rs = [Z || Z <- roots([float_q(C) || C <- P2]), cabs(csub(Z, {float(A), 0.0})) > 1.0e-9],
    Rho = lists:min([cabs(csub(Z, {float(A), 0.0})) || Z <- Rs] ++ [1.0e9]),
    H0 = pq_mp:from(min(0.25 * Rho, 0.25), W),
    Cs = sl_frob(D2, S1, S0, Sx, W, H0),
    Ts = pq_mp:exp(pq_mp:mul(Sx, pq_mp:log(H0, W), W), W),
    {Sum, DSum} = lists:foldl(fun({K, C}, {Ac, Bc}) ->
                                      T = pq_mp:mul(C, ipow(H0, K, W), W),
                                      {pq_mp:add(Ac, T, W), pq_mp:add(Bc, pq_mp:mul(pq_mp:add(Sx, pq_mp:from(K, W), W), T, W), W)}
                              end, {{mp, 0, 0}, {mp, 0, 0}}, lists:enumerate(0, Cs)),
    Psi0 = pq_mp:mul(Ts, Sum, W),
    DPsi0 = pq_mp:divide(pq_mp:mul(Ts, DSum, W), H0, W),
    Prob = #{p2 => P2m, p1 => P1m, r => R, roots => [{float(A), 0.0} | Rs]},
    sl_integrate(Prob, pq_mp:add(Am, H0, W), Xe, Psi0, DPsi0, W).

float_q({N, D}) -> {N, D}.

is_negligible(X, W) -> negligible(X, W - 40).

%% c_0 = 1; c_m [d_0 (m+s)(m+s-1) + p1_0 (m+s)] = -[sum_(j>=1) (d_j (m-j+s)(m-j+s-1) + p1_j (m-j+s)) c_(m-j) +
%% sum_(j>=0) r_j c_(m-1-j)]
sl_frob(D2, S1, S0, S, W, H0) -> sl_frob(D2, S1, S0, S, W, H0, 1, [{mp, 1, 0}], {mp, 1, 0}).

sl_frob(D2, S1, S0, S, W, H0, M, Acc, Hm) ->
    F = fun(J, K) -> % the coefficient of c at offset J with K = m - j + s
                pq_mp:add(pq_mp:mul(idx(D2, J), pq_mp:mul(K, pq_mp:sub(K, {mp, 1, 0}, W), W), W), pq_mp:mul(idx(S1, J), K, W), W)
        end,
    Ms = pq_mp:add(S, pq_mp:from(M, W), W),
    Den = F(0, Ms),
    Rev = Acc,   % c_(m-1), ..., c_0
    Num1 = lists:foldl(fun({J, C}, Sm) -> pq_mp:add(Sm, pq_mp:mul(F(J, pq_mp:add(S, pq_mp:from(M - J, W), W)), C, W), W) end,
                       {mp, 0, 0}, lists:zip(lists:seq(1, M), Rev)),
    Num = lists:foldl(fun({J, C}, Sm) -> pq_mp:add(Sm, pq_mp:mul(idx(S0, J), C, W), W) end,
                      Num1, lists:zip(lists:seq(0, M - 1), Rev)),
    Cm = case mag_rel(Den, Ms, W) of
             true -> case negligible(Num, W - 40) of true -> {mp, 0, 0}; false -> throw({ode, logarithmic_frobenius_solution}) end;
             false -> pq_mp:neg(pq_mp:divide(Num, Den, W))
         end,
    Hm2 = pq_mp:mul(Hm, H0, W),
    Acc2 = [Cm | Acc],
    case M > 8 andalso negligible(pq_mp:mul(Cm, Hm2, W), W) andalso negligible(pq_mp:mul(hd(tl(Acc2)), Hm, W), W) of
        true -> lists:reverse(Acc2);
        false when M > 4000 -> throw({ode, frobenius_not_converging});
        false -> sl_frob(D2, S1, S0, S, W, H0, M + 1, Acc2, Hm2)
    end.

%% the Taylor-series method for p2 psi'' + p1 psi' + r psi = 0
sl_integrate(Prob, X0, Xe, Psi, DPsi, W) ->
    case pq_mp:cmp(X0, Xe) >= 0 of
        true -> {Psi, DPsi};
        false ->
            {H, Last} = step_size(Prob, X0, Xe, W),
            {Psi2, DPsi2} = sl_step(Prob, X0, H, Psi, DPsi, W),
            X1 = case Last of true -> Xe; false -> pq_mp:add(X0, H, W) end,
            sl_integrate(Prob, X1, Xe, Psi2, DPsi2, W)
    end.

%% psi_(m+2) = -[sum_(j>=1) P2_j (m-j+2)(m-j+1) psi_(m-j+2) + sum_(j>=0) P1_j (m-j+1) psi_(m-j+1) + sum_(j>=0) R_j psi_(m-j)]
%%             / (P2_0 (m+2)(m+1))
sl_step(#{p2 := P2, p1 := P1, r := R}, A, H, Psi, DPsi, W) ->
    S2 = list_to_tuple(shift_poly(P2, A, W)), S1 = list_to_tuple(shift_poly(P1, A, W)), S0 = list_to_tuple(shift_poly(R, A, W)),
    El = fun(T, J) -> case J < tuple_size(T) of true -> element(J + 1, T); false -> {mp, 0, 0} end end,
    D0 = element(1, S2),
    Coef = sl_coefs(S2, S1, S0, El, D0, W, H, #{0 => Psi, 1 => DPsi}, 0),
    maps:fold(fun(K, C, {Sm, Dm}) ->
                      Sm2 = pq_mp:add(Sm, pq_mp:mul(C, ipow(H, K, W), W), W),
                      Dm2 = case K of 0 -> Dm; _ -> pq_mp:add(Dm, pq_mp:mul(pq_mp:mul(pq_mp:from(K, W), C, W), ipow(H, K - 1, W), W), W) end,
                      {Sm2, Dm2}
              end, {{mp, 0, 0}, {mp, 0, 0}}, Coef).

sl_coefs(S2, S1, S0, El, D0, W, H, Acc, M) ->
    N2 = lists:foldl(fun(J, Sm) -> K = M - J + 2, pq_mp:add(Sm, pq_mp:mul(pq_mp:mul(El(S2, J), pq_mp:from(K * (K - 1), W), W), maps:get(K, Acc), W), W) end,
                     {mp, 0, 0}, lists:seq(1, M)),
    N1 = lists:foldl(fun(J, Sm) -> K = M - J + 1, pq_mp:add(Sm, pq_mp:mul(pq_mp:mul(El(S1, J), pq_mp:from(K, W), W), maps:get(K, Acc), W), W) end,
                     N2, lists:seq(0, M)),
    N0 = lists:foldl(fun(J, Sm) -> pq_mp:add(Sm, pq_mp:mul(El(S0, J), maps:get(M - J, Acc), W), W) end, N1, lists:seq(0, M)),
    C = pq_mp:neg(pq_mp:divide(N0, pq_mp:mul(D0, pq_mp:from((M + 2) * (M + 1), W), W), W)),
    Acc2 = Acc#{M + 2 => C},
    Scale = max(mag(maps:get(0, Acc)), mag(pq_mp:mul(maps:get(1, Acc), H, W))),
    T1 = mag(pq_mp:mul(C, ipow(H, M + 2, W), W)),
    T0 = mag(pq_mp:mul(maps:get(M + 1, Acc2), ipow(H, M + 1, W), W)),
    case M > 6 andalso T1 < Scale - W - 4 andalso T0 < Scale - W - 4 of
        true -> Acc2;
        false when M > 3000 -> throw({ode, taylor_not_converging});
        false -> sl_coefs(S2, S1, S0, El, D0, W, H, Acc2, M + 1)
    end.

%% the Illinois method on [A, B] with the shooting function Shoot(Prob, E, P)
refine_with(Shoot, Prob, A0, B0, P) ->
    A = pq_mp:from(A0, P), B = pq_mp:from(B0, P),
    illinois_f(Shoot, Prob, A, Shoot(Prob, A, P), B, Shoot(Prob, B, P), P, 0, none).

illinois_f(_, _, A, _, B, _, P, 200, _) -> pq_mp:ldexp(pq_mp:add(A, B, P), -1);
illinois_f(Shoot, Prob, A, Fa, B, Fb, P, I, Side) ->
    C = pq_mp:sub(B, pq_mp:divide(pq_mp:mul(Fb, pq_mp:sub(B, A, P), P), pq_mp:sub(Fb, Fa, P), P), P),
    case small(pq_mp:sub(B, A, P), C, P) of
        true -> C;
        false ->
            Fc = Shoot(Prob, C, P),
            case sgn(Fc) of
                0 -> C;
                S ->
                    case S =:= sgn(Fb) of
                        true ->
                            Fa2 = case Side of a -> pq_mp:ldexp(Fa, -1); _ -> Fa end,
                            illinois_f(Shoot, Prob, A, Fa2, C, Fc, P, I + 1, a);
                        false ->
                            Fb2 = case Side of b -> pq_mp:ldexp(Fb, -1); _ -> Fb end,
                            illinois_f(Shoot, Prob, C, Fc, B, Fb2, P, I + 1, b)
                    end
            end
    end.
