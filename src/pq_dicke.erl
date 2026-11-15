%% N identical two-level atoms in the permutation-invariant (Dicke) representation (quantum optics plan, Q5), as PIQS
%% does (Shammah, Ahmed, Lambert, De Liberato, Nori, PRA 98, 063815 (2018), arXiv:1805.05129; QuTiP's qutip.piqs).
%%
%% For permutation-symmetric dynamics the density matrix needs only the elements rho_(j, m, m') within each Dicke
%% ladder j (j = N/2, N/2 - 1, ..., 0 or 1/2; m, m' = -j .. j), the degenerate copies of a ladder carried together.
%% Collective processes (rates gCE, gCD, gCP for J_-, J_z, J_+) act within a ladder; local ones (each atom's sigma_-,
%% sigma_z, sigma_+ at gE, gD, gP) also move between ladders. The master equation is
%%   d rho_(j,m,m')/dt = sum over the nine neighbours (j + dj, m + dm, m' + dm), dj, dm in {-1, 0, 1}, of
%%                       Gamma_k(neighbour) rho_(neighbour),
%% with PIQS's coefficients gamma1 .. gamma9 (Shammah et al., Eqs. 10-12), here in arbitrary precision; plus
%% -i [H, rho] for H = omega J_z + omega_x J_x within each ladder. Every process shifts m and m' together, so m - m' is
%% conserved; without the J_x drive a diagonal initial state needs only the populations (about N^2/4 numbers).
%%
%% j and m are carried doubled (J = 2j, M = 2m) so that they stay integers. The vector of rho's elements is indexed by
%% the triples {J, M, M1} kept; pq_lindblad:evolve/4 integrates it.
-module(pq_dicke).
-export([triples/2, liouvillian/4, state/3, expect/4, j_ladders/1]).

%%% The basis

%% the doubled j of the ladders, from N/2 down
j_ladders(N) -> lists:seq(N, N rem 2, -2).

%% the triples {J, M, M1} kept: all (Mode = all) or the populations only (Mode = diagonal)
triples(N, Mode) ->
    [{J, M, M1} || J <- j_ladders(N), M <- lists:seq(J, -J, -2), M1 <- lists:seq(J, -J, -2),
                   Mode =:= all orelse M =:= M1].

%%% The Liouvillian

%% Rates #{emission, dephasing, pumping, collective_emission, collective_dephasing, collective_pumping, omega, omega_x}
%% (rationals {N, D} or integers; missing ones 0); Triples from triples/2. A sparse matrix in pq_lindblad's form.
liouvillian(N, Rates, Triples, P) ->
    R = fun(K) -> case maps:get(K, Rates, 0) of {A, B} -> pq_mp:from({A, B}, P); X -> pq_mp:from(X, P) end end,
    Y = #{e => R(emission), d => R(dephasing), p => R(pumping), ce => R(collective_emission),
          cd => R(collective_dephasing), cp => R(collective_pumping), w => R(omega), wx => R(omega_x)},
    Index = maps:from_list(lists:zip(Triples, lists:seq(0, length(Triples) - 1))),
    Nm = pq_mp:from(N, P),
    Rows = [row(T, Index, Y, Nm, P) || T <- Triples],
    {sparse, length(Triples), list_to_tuple(Rows)}.

row({J, M, M1} = T, Index, Y, N, P) ->
    Neigh = [{1, T}, {2, {J, M + 2, M1 + 2}}, {3, {J + 2, M + 2, M1 + 2}}, {4, {J - 2, M + 2, M1 + 2}},
             {5, {J + 2, M, M1}}, {6, {J - 2, M, M1}}, {7, {J + 2, M - 2, M1 - 2}}, {8, {J, M - 2, M1 - 2}},
             {9, {J - 2, M - 2, M1 - 2}}],
    Diss = [{maps:get(C, Index), {gamma(K, C, Y, N, P), zero_mp()}} || {K, C} <- Neigh, is_map_key(C, Index)],
    Ham = hamiltonian(T, Index, Y, P),
    merge(Diss ++ Ham, P).

merge(Entries, P) ->
    M = lists:foldl(fun({C, V}, Acc) -> Acc#{C => case Acc of #{C := V0} -> cadd(V0, V, P); _ -> V end} end, #{}, Entries),
    lists:sort([{C, V} || C := V <- M, not czero(V)]).

%% PIQS's gamma1 .. gamma9 at the triple {J, M, M1} (doubled), as an mp
gamma(K, {J2, M2, Mp2}, Y, N, P) ->
    J = half(J2, P), M = half(M2, P), M1 = half(Mp2, P),
    #{e := E, d := D, p := Pu, ce := CE, cd := CD, cp := CP} = Y,
    Half = pq_mp:from({1, 2}, P),
    Nh = pq_mp:mul(N, Half, P),
    Add = fun(L) -> lists:foldl(fun(X, S) -> pq_mp:add(S, X, P) end, zero_mp(), L) end,
    Mul = fun(L) -> lists:foldl(fun(X, S) -> pq_mp:mul(S, X, P) end, pq_mp:from(1, P), L) end,
    Sub = fun(A, B) -> pq_mp:sub(A, B, P) end,
    I = fun(X) -> pq_mp:from(X, P) end,
    Div = fun(A, B) -> pq_mp:divide(A, B, P) end,
    Sq = fun(X) -> pq_mp:sqrt(X, P) end,
    Jp = Mul([J, Add([J, I(1)])]),                                     % j (j + 1)
    Pos = fun(X) -> pq_mp:cmp(X, zero_mp()) =:= 1 end,
    case K of
        1 ->
            Spont = Mul([CE, Half, Sub(Sub(Mul([I(2), Jp]), Mul([M, Sub(M, I(1))])), Mul([M1, Sub(M1, I(1))]))]),
            Loss = Mul([E, Half, Add([N, M, M1])]),
            Pump = Mul([Pu, Half, Sub(Sub(N, M), M1)]),
            CPump = Mul([CP, Half, Sub(Sub(Mul([I(2), Jp]), Mul([M, Add([M, I(1)])])), Mul([M1, Add([M1, I(1)])]))]),
            CDeph = Mul([CD, Half, Mul([Sub(M, M1), Sub(M, M1)])]),
            Deph = case Pos(J) of
                       false -> Mul([D, N, I({1, 4})]);
                       true -> Mul([D, Half, Sub(Nh, Div(Mul([M, M1, Add([Nh, I(1)])]), Jp))])
                   end,
            pq_mp:neg(Add([Spont, Loss, Pump, CPump, CDeph, Deph]));
        2 ->
            Root = Sq(Mul([Add([J, M]), Add([Sub(J, M), I(1)]), Add([J, M1]), Add([Sub(J, M1), I(1)])])),
            Loss = case Pos(J) of true -> Div(Mul([E, Half, Root, Add([Nh, I(1)])]), Jp); false -> zero_mp() end,
            Add([Mul([CE, Root]), Loss]);
        3 ->
            nz(Pos(J), fun() -> Div(Mul([E, Half, Sq(Mul([Add([J, M]), Sub(Add([J, M]), I(1)), Add([J, M1]), Sub(Add([J, M1]), I(1))])),
                                         Add([Nh, J, I(1)])]), Mul([J, Add([Mul([I(2), J]), I(1)])])) end);
        4 ->
            nz(Pos(Add([J, I(1)])), fun() -> Div(Mul([E, Half, Sq(Mul([Add([Sub(J, M), I(1)]), Add([Sub(J, M), I(2)]),
                                                                       Add([Sub(J, M1), I(1)]), Add([Sub(J, M1), I(2)])])),
                                                     Sub(Nh, J)]), Mul([Add([J, I(1)]), Add([Mul([I(2), J]), I(1)])])) end);
        5 ->
            nz(Pos(J), fun() -> Div(Mul([D, Half, Sq(Mul([Sub(Mul([J, J]), Mul([M, M])), Sub(Mul([J, J]), Mul([M1, M1]))])),
                                         Add([Nh, J, I(1)])]), Mul([J, Add([Mul([I(2), J]), I(1)])])) end);
        6 ->
            J1 = Add([J, I(1)]),
            Div(Mul([D, Half, Sq(Mul([Sub(Mul([J1, J1]), Mul([M, M])), Sub(Mul([J1, J1]), Mul([M1, M1]))])), Sub(Nh, J)]),
                Mul([J1, Add([Mul([I(2), J]), I(1)])]));
        7 ->
            nz(Pos(J), fun() -> Div(Mul([Pu, Half, Sq(Mul([Sub(Sub(J, M), I(1)), Sub(J, M), Sub(Sub(J, M1), I(1)), Sub(J, M1)])),
                                         Add([Nh, J, I(1)])]), Mul([J, Add([Mul([I(2), J]), I(1)])])) end);
        8 ->
            Root = Sq(Mul([Add([J, M, I(1)]), Sub(J, M), Add([J, M1, I(1)]), Sub(J, M1)])),
            Pump = case Pos(J) of true -> Div(Mul([Pu, Half, Root, Add([Nh, I(1)])]), Jp); false -> zero_mp() end,
            Add([Pump, Mul([CP, Root])]);
        9 ->
            J1 = Add([J, I(1)]),
            Div(Mul([Pu, Half, Sq(Mul([Add([J, M, I(1)]), Add([J, M, I(2)]), Add([J, M1, I(1)]), Add([J, M1, I(2)])])), Sub(Nh, J)]),
                Mul([J1, Add([Mul([I(2), J]), I(1)])]))
    end.

nz(true, F) -> F();
nz(false, _) -> zero_mp().

half(X2, P) -> pq_mp:from({X2, 2}, P).

%% -i (H rho - rho H) for H = omega J_z + omega_x J_x, within the ladder of T; J_x = (J_+ + J_-)/2,
%% <j, m +- 1| J_+- |j, m> = sqrt(j(j+1) - m(m +- 1))
hamiltonian({J, M, M1}, Index, #{w := W, wx := Wx}, P) ->
    Mi = {zero_mp(), pq_mp:from(-1, P)},
    Wz = case W of {mp, 0, _} -> []; _ ->
             % (H rho - rho H)_(m,m1) = omega (m - m1) rho_(m,m1)
             [{maps:get({J, M, M1}, Index), cmul(Mi, {pq_mp:mul(W, half(M - M1, P), P), zero_mp()}, P)}]
         end,
    Wxs = case Wx of {mp, 0, _} -> []; _ ->
              Half = pq_mp:from({1, 2}, P),
              Me = fun(Ma, Mb) -> % <m_a| J_x |m_b>, m_a = m_b +- 1 (doubled indices)
                           Jh = half(J, P), Mbh = half(Mb, P),
                           S = case Ma - Mb of 2 -> 1; -2 -> -1 end,
                           pq_mp:mul(Half, pq_mp:sqrt(pq_mp:sub(pq_mp:mul(Jh, pq_mp:add(Jh, pq_mp:from(1, P), P), P),
                                                                pq_mp:mul(Mbh, pq_mp:add(Mbh, pq_mp:from(S, P), P), P), P), P), P)
                   end,
              % (J_x rho)_(m, m1) = sum_k <m|J_x|k> rho_(k, m1); (rho J_x)_(m, m1) = sum_k rho_(m, k) <k|J_x|m1>
              Left = [{{J, K, M1}, Me(M, K)} || K <- [M - 2, M + 2], abs(K) =< J],
              Right = [{{J, M, K}, Me(K, M1)} || K <- [M1 - 2, M1 + 2], abs(K) =< J],
              [{maps:get(T, Index), cmul(Mi, {pq_mp:mul(Wx, V, P), zero_mp()}, P)} || {T, V} <- Left, is_map_key(T, Index)]
                  ++ [{maps:get(T, Index), cmul(Mi, {pq_mp:neg(pq_mp:mul(Wx, V, P)), zero_mp()}, P)} || {T, V} <- Right, is_map_key(T, Index)]
          end,
    Wz ++ Wxs.

%%% States and averages

%% the vector of an initial state: excited (all atoms up: j = m = N/2), ground (j = N/2, m = -N/2), or {dicke, J2, M2}
state(Triples, Which, P) ->
    {J, M} = case Which of
                 {excited, N} -> {N, N};
                 {ground, N} -> {N, -N};
                 {dicke, J2, M2} -> {J2, M2}
             end,
    list_to_tuple([case T of {J, M, M} -> {pq_mp:from(1, P), zero_mp()}; _ -> zero() end || T <- Triples]).

%% <O> for O = jz | jp_jm (J_+ J_-) | jz2 (J_z^2) | trace, from the populations (the operators are diagonal in m)
expect(O, Triples, Rho, P) ->
    lists:foldl(fun({{J, M, M1}, {Re, _}}, S) when M =:= M1 ->
                        Jh = half(J, P), Mh = half(M, P),
                        F = case O of
                                trace -> pq_mp:from(1, P);
                                jz -> Mh;
                                jz2 -> pq_mp:mul(Mh, Mh, P);
                                jp_jm -> pq_mp:mul(pq_mp:add(Jh, Mh, P), pq_mp:add(pq_mp:sub(Jh, Mh, P), pq_mp:from(1, P), P), P)
                            end,
                        pq_mp:add(S, pq_mp:mul(F, Re, P), P);
                   (_, S) -> S
                end, zero_mp(), lists:zip(Triples, tuple_to_list(Rho))).

%%% Complex helpers ({Re, Im} of mps)

zero_mp() -> {mp, 0, 0}.
zero() -> {zero_mp(), zero_mp()}.
cadd({A, B}, {C, D}, P) -> {pq_mp:add(A, C, P), pq_mp:add(B, D, P)}.
cmul({A, B}, {C, D}, P) -> {pq_mp:sub(pq_mp:mul(A, C, P), pq_mp:mul(B, D, P), P), pq_mp:add(pq_mp:mul(A, D, P), pq_mp:mul(B, C, P), P)}.
czero({{mp, 0, _}, {mp, 0, _}}) -> true;
czero(_) -> false.
