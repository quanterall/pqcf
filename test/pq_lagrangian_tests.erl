%% pq_lagrangian (plan "theoretical physics instruments", F4): Feynman rules from Lagrangians against Peskin and
%% Schroeder, An Introduction to Quantum Field Theory (1995): phi^4, Yukawa, QED, scalar QED and QCD (quark-gluon,
%% three- and four-gluon, ghost-gluon vertices), all momenta incoming. Leg indices m1.. (Lorentz) and a1.. (colour).
-module(pq_lagrangian_tests).
-include_lib("eunit/include/eunit.hrl").

-define(I(X), {{0, 1}, X}).          % i X
-define(R(X), {X, {0, 1}}).          % X

kinds() -> #{phi => [], phistar => [], psi => [], psibar => [], a => [lorentz], gl => [lorentz, internal],
             q => [], qbar => [], c => [internal], cbar => [internal]}.

v(L, Legs) -> pq_lagrangian:vertex(L, Legs, kinds()).
eq(Want, Got) -> ?assertEqual(pq_lagrangian:canonical(Want), Got).

%% L = -lambda/4! phi^4: -i lambda
phi4_test() ->
    L = [{{-1, 24}, [{c, lambda}], [{phi, [], []}, {phi, [], []}, {phi, [], []}, {phi, [], []}]}],
    eq([{?I({-1, 1}), [{c, lambda}]}], v(L, [{phi, []}, {phi, []}, {phi, []}, {phi, []}])),
    % phi^3 with -g/3!: -i g
    L3 = [{{-1, 6}, [{c, g}], [{phi, [], []}, {phi, [], []}, {phi, [], []}]}],
    eq([{?I({-1, 1}), [{c, g}]}], v(L3, [{phi, []}, {phi, []}, {phi, []}])).

%% Yukawa, L = -g phi psibar psi: -i g
yukawa_test() ->
    L = [{-1, [{c, g}], [{phi, [], []}, {psibar, [], []}, {psi, [], []}]}],
    eq([{?I({-1, 1}), [{c, g}]}], v(L, [{psibar, []}, {psi, []}, {phi, []}])).

%% QED, L = -e psibar gamma^mu psi A_mu: -i e gamma^mu
qed_test() ->
    L = [{-1, [{c, e}, {gamma, mu}], [{psibar, [], []}, {psi, [], []}, {a, [mu], []}]}],
    eq([{?I({-1, 1}), [{c, e}, {gamma, m1}]}], v(L, [{psibar, []}, {psi, []}, {a, [m1]}])).

%% scalar QED, |D phi|^2 with D = d + i e A:
%%   L3 = i e A_mu (phi d^mu phi* - phi* d^mu phi),   L4 = e^2 g^(mu nu) A_mu A_nu phi* phi
%% P&S: -i e (p + p')^mu (p the incoming phi, p' the outgoing one: the phi* leg's incoming momentum is -p'), and
%% 2 i e^2 g^(mu nu)
scalar_qed_test() ->
    L3 = [{{complex, 0, 1}, [{c, e}], [{a, [mu], []}, {phi, [], []}, {phistar, [], [mu]}]},
          {{complex, 0, -1}, [{c, e}], [{a, [mu], []}, {phistar, [], []}, {phi, [], [mu]}]}],
    eq([{?I({-1, 1}), [{c, e}, {p, 1, m1}]}, {?I({1, 1}), [{c, e}, {p, 2, m1}]}],
       v(L3, [{phi, []}, {phistar, []}, {a, [m1]}])),
    L4 = [{1, [{c, e}, {c, e}, {g, mu, nu}], [{a, [mu], []}, {a, [nu], []}, {phistar, [], []}, {phi, [], []}]}],
    eq([{?I({2, 1}), [{c, e}, {c, e}, {g, m1, m2}]}], v(L4, [{a, [m1]}, {a, [m2]}, {phistar, []}, {phi, []}])).

%% QCD (P&S chapter 16), the gauge field A^a_mu = gl with indices [mu, a]
qcd() ->
    [% quark: qbar (i Dslash) q with D = d - i g A^a t^a: + g qbar gamma^mu t^a q A^a_mu
     {1, [{c, g}, {gamma, mu}, {t, a}], [{qbar, [], []}, {q, [], []}, {gl, [mu, a], []}]},
     % -1/4 F^a F^a, F^a_(mu nu) = d_mu A^a_nu - d_nu A^a_mu + g f^abc A^b_mu A^c_nu:
     % cubic -g f^abc (d_mu A^a_nu) A^(b mu) A^(c nu), quartic -g^2/4 f^abe f^cde A^a_mu A^b_nu A^(c mu) A^(d nu)
     {-1, [{c, g}, {f, a, b, c}], [{gl, [nu, a], [mu]}, {gl, [mu, b], []}, {gl, [nu, c], []}]},
     {{-1, 4}, [{c, g}, {c, g}, {f, a, b, e}, {f, c, d, e}],
      [{gl, [mu, a], []}, {gl, [nu, b], []}, {gl, [mu, c], []}, {gl, [nu, d], []}]},
     % ghosts: (d^mu cbar^a) D^ac_mu c^c, D^ac_mu = d_mu delta^ac + g f^abc A^b_mu: cubic g f^abc (d^mu cbar^a) A^b_mu c^c
     {1, [{c, g}, {f, a, b, c}], [{cbar, [a], [mu]}, {gl, [mu, b], []}, {c, [c], []}]}].

%% quark-gluon: i g gamma^mu t^a
quark_gluon_test() ->
    eq([{?I({1, 1}), [{c, g}, {gamma, m1}, {t, a1}]}], v(qcd(), [{qbar, []}, {q, []}, {gl, [m1, a1]}])).

%% three gluons, (a1, m1, p1), (a2, m2, p2), (a3, m3, p3) incoming:
%%   g f^(a1 a2 a3) [g^(m1 m2) (p1 - p2)^m3 + g^(m2 m3) (p2 - p3)^m1 + g^(m3 m1) (p3 - p1)^m2]
three_gluon_test() ->
    F = {f, a1, a2, a3},
    Want = [{?R({S, 1}), [{c, g}, F, {g, X, Y}, {p, K, Z}]}
            || {{X, Y, Z}, Ks} <- [{{m1, m2, m3}, [{1, 1}, {2, -1}]}, {{m2, m3, m1}, [{2, 1}, {3, -1}]},
                                   {{m3, m1, m2}, [{3, 1}, {1, -1}]}],
               {K, S} <- Ks],
    eq(Want, v(qcd(), [{gl, [m1, a1]}, {gl, [m2, a2]}, {gl, [m3, a3]}])).

%% four gluons: -i g^2 [f^(abe) f^(cde) (g^mu rho g^nu sigma - g^mu sigma g^nu rho)
%%                    + f^(ace) f^(bde) (g^mu nu g^rho sigma - g^mu sigma g^nu rho)
%%                    + f^(ade) f^(bce) (g^mu nu g^rho sigma - g^mu rho g^nu sigma)], (a, mu) .. (d, sigma) the legs
four_gluon_test() ->
    G = fun(X, Y) -> {g, X, Y} end,
    Ff = fun(A, B, C, D) -> [{f, A, B, e}, {f, C, D, e}] end,
    Want = [{?I({S, 1}), [{c, g}, {c, g}] ++ Ff(A, B, C, D) ++ [G(X1, Y1), G(X2, Y2)]}
            || {{A, B, C, D}, Pairs} <- [{{a1, a2, a3, a4}, [{-1, {m1, m3}, {m2, m4}}, {1, {m1, m4}, {m2, m3}}]},
                                          {{a1, a3, a2, a4}, [{-1, {m1, m2}, {m3, m4}}, {1, {m1, m4}, {m2, m3}}]},
                                          {{a1, a4, a2, a3}, [{-1, {m1, m2}, {m3, m4}}, {1, {m1, m3}, {m2, m4}}]}],
               {S, {X1, Y1}, {X2, Y2}} <- Pairs],
    eq(Want, v(qcd(), [{gl, [m1, a1]}, {gl, [m2, a2]}, {gl, [m3, a3]}, {gl, [m4, a4]}])).

%% ghost-gluon: -g f^abc p^mu with p the outgoing ghost momentum; with the antighost leg (a1) incoming momentum p1
%% (= -p): g f^(a1 a2 a3) p1^m2
ghost_test() ->
    eq([{?R({1, 1}), [{c, g}, {f, a1, a2, a3}, {p, 1, m2}]}], v(qcd(), [{cbar, [a1]}, {gl, [m2, a2]}, {c, [a3]}])).

%% the quadratic parts give the inverse propagators: i (p^2 - m^2) for 1/2 (d phi)^2 - 1/2 m^2 phi^2, with the two
%% legs' momenta p1 and p2 = -p1 the result is in p1.p2 = -p^2
quadratic_test() ->
    L = [{{1, 2}, [{g, mu, nu}], [{phi, [], [mu]}, {phi, [], [nu]}]}, {{-1, 2}, [{c, m}, {c, m}], [{phi, [], []}, {phi, [], []}]}],
    % i [(-i p1)(-i p2) g] (two assignments, each 1/2) - i m^2 = -i p1.p2 - i m^2; with p2 = -p1: i (p^2 - m^2)
    eq([{?I({-1, 1}), [{dot, 1, 2}]}, {?I({-1, 1}), [{c, m}, {c, m}]}], v(L, [{phi, []}, {phi, []}])).

%% canonical forms: f is antisymmetric, the metric symmetric, dummy names do not matter
canonical_test() ->
    ?assertEqual(pq_lagrangian:canonical([{?R({1, 1}), [{f, a, b, c}]}]),
                 pq_lagrangian:canonical([{?R({-1, 1}), [{f, b, a, c}]}])),
    ?assertEqual([], pq_lagrangian:canonical([{?R({1, 1}), [{f, a, a, c}]}])),
    ?assertEqual(pq_lagrangian:canonical([{?R({1, 1}), [{f, a, b, x}, {f, c, d, x}]}]),
                 pq_lagrangian:canonical([{?R({1, 1}), [{f, a, b, y}, {f, c, d, y}]}])),
    ?assertEqual(pq_lagrangian:canonical([{?R({1, 1}), [{g, m1, nu}, {p, 1, nu}]}]),
                 pq_lagrangian:canonical([{?R({1, 1}), [{p, 1, m1}]}])).

%%% Propagators: Gamma_2 D = -1, numerically

%% Gamma_2 from the quadratic Lagrangian (two legs), then evaluated at a momentum p with p_2 = -p_1, the Minkowski
%% metric diag(1, -1, -1, -1) and the Dirac matrices in the Dirac representation; complex numbers as {Re, Im} floats
%% of exact integers (so the comparisons are exact)
-define(P, [3, 1, -2, 5]).
eta(0) -> 1;
eta(_) -> -1.
pv(K) -> [X * case K of 1 -> 1; 2 -> -1 end || X <- ?P].
psq() -> lists:sum([eta(M) * X * X || {M, X} <- lists:zip([0, 1, 2, 3], ?P)]).

%% 4x4 complex matrices as maps {R, C} => {Re, Im}
gamma(0) -> diag([1, 1, -1, -1]);
gamma(K) ->
    S = case K of 1 -> [[{0, 0}, {1, 0}], [{1, 0}, {0, 0}]]; 2 -> [[{0, 0}, {0, -1}], [{0, 1}, {0, 0}]];
                  3 -> [[{1, 0}, {0, 0}], [{0, 0}, {-1, 0}]] end,
    maps:from_list([{{R, C + 2}, lists:nth(C + 1, lists:nth(R + 1, S))} || R <- [0, 1], C <- [0, 1]] ++
                   [{{R + 2, C}, cneg(lists:nth(C + 1, lists:nth(R + 1, S)))} || R <- [0, 1], C <- [0, 1]]).
diag(L) -> maps:from_list([{{I, I}, {X, 0}} || {I, X} <- lists:zip([0, 1, 2, 3], L)]).
mget(M, K) -> maps:get(K, M, {0, 0}).
cneg({A, B}) -> {-A, -B}.
cm({A, B}, {C, D}) -> {A * C - B * D, A * D + B * C}.
ca({A, B}, {C, D}) -> {A + C, B + D}.
madd(A, B) -> maps:from_list([{{R, C}, ca(mget(A, {R, C}), mget(B, {R, C}))} || R <- [0, 1, 2, 3], C <- [0, 1, 2, 3]]).
mscale(A, S) -> maps:map(fun(_, V) -> cm(V, S) end, A).
mmul(A, B) -> maps:from_list([{{R, C}, lists:foldl(fun(K, Acc) -> ca(Acc, cm(mget(A, {R, K}), mget(B, {K, C}))) end, {0, 0}, [0, 1, 2, 3])}
                              || R <- [0, 1, 2, 3], C <- [0, 1, 2, 3]]).
ident() -> diag([1, 1, 1, 1]).

%% the value of a term list with free indices bound by Env, a repeated index summed with its metric factor; the
%% result a complex number (scalar terms) or a matrix (one gamma)
eval(Terms, Env, Params) ->
    lists:foldl(fun({{Re, Im}, Fs}, Acc) -> add_v(Acc, cm_v(eval_term(Fs, Env, Params), {fl(Re), fl(Im)})) end, {0, 0}, Terms).
fl({N, 1}) -> N.
eval_term(Fs, Env, Params) ->
    All = lists:append([idx(F) || F <- Fs]),
    Free = [I || I <- lists:usort(All), maps:is_key(I, Env)],
    Dummies = lists:usort([I || I <- All, not lists:member(I, Free)]),
    lists:foldl(fun(Vals, Acc) ->
                        E = maps:merge(Env, maps:from_list(lists:zip(Dummies, Vals))),
                        Metric = lists:foldl(fun(V, P) -> P * eta(V) end, 1, Vals),
                        add_v(Acc, cm_v(lists:foldl(fun(F, X) -> mul_v(X, factor(F, E, Params)) end, {1, 0}, Fs), {Metric, 0}))
                end, {0, 0}, tuples(length(Dummies))).
tuples(0) -> [[]];
tuples(N) -> [[V | T] || V <- [0, 1, 2, 3], T <- tuples(N - 1)].
idx({g, A, B}) -> [A, B];
idx({p, _, A}) -> [A];
idx({gamma, A}) -> [A];
idx(_) -> [].
factor({c, X}, _, Params) -> {maps:get(X, Params), 0};
factor({g, A, B}, E, _) -> case maps:get(A, E) =:= maps:get(B, E) of true -> {eta(maps:get(A, E)), 0}; false -> {0, 0} end;
factor({p, K, A}, E, _) -> {lists:nth(maps:get(A, E) + 1, pv(K)), 0};
factor({dot, 1, 2}, _, _) -> {-psq(), 0};
factor({dot, K, K}, _, _) -> {psq(), 0};
factor({gamma, A}, E, _) -> gamma(maps:get(A, E)).
%% values: complex {Re, Im} or matrices (maps)
mul_v(X, Y) when is_map(X), is_map(Y) -> mmul(X, Y);
mul_v(X, Y) when is_map(X) -> mscale(X, Y);
mul_v(X, Y) when is_map(Y) -> mscale(Y, X);
mul_v(X, Y) -> cm(X, Y).
cm_v(X, S) when is_map(X) -> mscale(X, S);
cm_v(X, S) -> cm(X, S).
add_v(X, Y) when is_map(X), is_map(Y) -> madd(X, Y);
add_v({0, 0}, Y) -> Y;
add_v(X, {0, 0}) -> X;
add_v(X, Y) when is_map(X) -> madd(X, mscale(ident(), Y));
add_v(X, Y) when is_map(Y) -> madd(mscale(ident(), X), Y);
add_v(X, Y) -> ca(X, Y).

%% scalar: Gamma_2 = i (p^2 - m^2), D = i / (p^2 - m^2): Gamma_2 D = -1
scalar_propagator_test() ->
    L = [{{1, 2}, [{g, mu, nu}], [{phi, [], [mu]}, {phi, [], [nu]}]}, {{-1, 2}, [{c, m}, {c, m}], [{phi, [], []}, {phi, [], []}]}],
    Prm = #{m => 7},
    G2 = eval(v(L, [{phi, []}, {phi, []}]), #{}, Prm),
    {Num, Den} = pq_lagrangian:propagator(scalar, m),
    % D = Num / Den: Gamma_2 Num = -Den
    ?assertEqual(cneg(eval(Den, #{}, Prm)), cm(G2, eval(Num, #{}, Prm))).

%% Dirac: Gamma_2 = i (pslash - m) with p the psi leg's incoming momentum (leg 2 here, = -p_1)
dirac_propagator_test() ->
    L = [{{complex, 0, 1}, [{gamma, mu}], [{psibar, [], []}, {psi, [], [mu]}]}, {-1, [{c, m}], [{psibar, [], []}, {psi, [], []}]}],
    Prm = #{m => 7},
    % legs: psi first, so its momentum is p_1 = p
    G2 = eval(v(L, [{psi, []}, {psibar, []}]), #{}, Prm),
    {Num, Den} = pq_lagrangian:propagator(dirac, m),
    ?assertEqual(nz(mscale(ident(), cneg(eval(Den, #{}, Prm)))), nz(mul_v(G2, eval(Num, #{}, Prm)))).

%% a matrix without its zero entries
nz(M) -> maps:filter(fun(_, V) -> V =/= {0, 0} end, M).

%% vector in R_xi gauge (xi = 3): Gamma_2^(mu nu) D_(nu rho) = -delta^mu_rho, the indices lowered with eta
%% (Lagrangian: -1/4 F^2 - 1/(2 xi) (d.A)^2 written out: -1/2 (d_mu A_nu)^2 + 1/2 (d_mu A_nu)(d^nu A^mu) - (d.A)^2/(2 xi))
vector_propagator_test() ->
    Xi = 3,
    % the Lagrangian times 2 xi = 6, to keep the arithmetic in integers: then Gamma_2 D = -6
    L = [{-Xi, [{g, mu, al}, {g, nu, be}], [{a, [nu], [mu]}, {a, [be], [al]}]},
         {Xi, [{g, mu, be}, {g, nu, al}], [{a, [nu], [mu]}, {a, [be], [al]}]},
         {-1, [], [{a, [mu], [mu]}, {a, [nu], [nu]}]}],
    G2 = v(L, [{a, [m1]}, {a, [m2]}]),
    {Num, Den} = pq_lagrangian:propagator({vector, n1, n2}, xi),
    Prm = #{xi => Xi, {one_minus, xi} => 1 - Xi},
    P2 = psq(),
    DenV = eval(Den, #{}, Prm),
    [begin
         Sum = lists:foldl(fun(Nu, Acc) ->
                                   G = eval(G2, #{m1 => Mu, m2 => Nu}, Prm),
                                   D = eval([inv_fix(T, P2) || T <- Num], #{n1 => Nu, n2 => Rho}, Prm),
                                   ca(Acc, cm(cm(G, D), {eta(Nu), 0}))
                           end, {0, 0}, [0, 1, 2, 3]),
         % Sum / (p^2 * p^2 for the 1/p^2 carried in the numerator) = -delta
         % all indices upper: Gamma eta D = -eta, here times 2 xi (the Lagrangian's factor) and p^4 (both 1/p^2 cleared)
         Want = case Mu =:= Rho of true -> cneg(cm(DenV, {2 * Xi * P2 * eta(Mu), 0})); false -> {0, 0} end,
         ?assertEqual(Want, Sum)
     end || Mu <- [0, 1, 2, 3], Rho <- [0, 1, 2, 3]].

%% the numerator's 1/p^2: everything multiplied by p^2 instead (so the check stays in integers)
inv_fix({C, Fs}, P2) ->
    case lists:member({inv, dot11}, Fs) of
        true -> {C, Fs -- [{inv, dot11}]};
        false -> {mulc(C, P2), Fs}
    end.
mulc({{A, 1}, {B, 1}}, K) -> {{A * K, 1}, {B * K, 1}}.
