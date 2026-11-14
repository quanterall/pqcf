%% Master equations in a truncated basis, in arbitrary precision (quantum optics plan, Q2).
%%
%% The operators of pq_qop become sparse matrices in the product basis of the subsystems (Fock states 0 .. N of each
%% bosonic mode, levels 1 .. L of each atom), their coefficients evaluated at the parameters given. A word is projected
%% on the truncated basis as it stands (normal ordered), so a a^dagger = a^dagger a + 1 is n + 1 up to n = N.
%%
%% The Schrödinger equation psi' = -i H psi and the master equation
%%   rho' = L rho = -i [H, rho] + sum_k (L_k rho L_k^dagger - (1/2) {L_k^dagger L_k, rho})
%% are linear with constant coefficients, y' = M y: they are integrated by the Taylor series, term_k = (h/k) M term_(k-1),
%% summed until the terms fall below the working precision, with steps h =< 1/|M| so that the series converges like
%% 1/k!. The result is exact to the precision P (bits) apart from the rounding of the arithmetic.
%%
%% Numbers are {Re, Im} pairs of pq_mp values; vectors are tuples; rho is vectorised row by row (rho_ij at i D + j);
%% sparse matrices are {sparse, D, Rows}, Rows a tuple of [{Col, Value}] (0-based columns).
-module(pq_lindblad).
-export([basis/2, dim/1, op/4, state/3, fock/1, coherent/3, dm/2, hamiltonian_generator/2, liouvillian/3,
         evolve/4, sesolve/4, mesolve/5, expect/3, expect_dm/3, trace/2, steady_state/2, trajectories/7,
         to_complex/1, cmag2/2]).

%%% Bases: [{Name, Size, Kind}] in name order, Kind boson | atom; the state index in mixed radix, the last name fastest

%% Spaces: the pq_qop spaces (atoms {atom, L}); Dims: #{Mode => N}, the largest photon number of each bosonic mode
basis(Spaces, Dims) ->
    Atoms = [N || N := {atom, _} <- Spaces],
    Names = lists:usort(maps:keys(Dims) ++ Atoms),
    [case maps:get(N, Spaces, boson) of
         {atom, L} -> {N, L, atom};
         _ -> {N, maps:get(N, Dims) + 1, boson}
     end || N <- Names].

dim(Basis) -> lists:foldl(fun({_, S, _}, Acc) -> Acc * S end, 1, Basis).

%%% Operators

%% the matrix of a pq_qop expression; Params binds the coefficients' symbols (#{g => 1, kappa => {1, 10}})
op(Basis, Spaces, Expr, {Params, P}) ->
    D = dim(Basis),
    Elems = pq_qop:with_spaces(Spaces, fun() ->
        lists:append([word_elements(W, ceval(C, Params, P), Basis, P) || W := C <- Expr])
    end),
    sparse(D, Elems, P).

%% the elements {Row, Col, Value} of one word: the tensor product over the subsystems
word_elements(W, Cv, Basis, P) ->
    lists:foldl(fun({N, S, Kind}, Acc) ->
                        Local = case proplists:get_value(N, W, none) of
                                    none -> [{I, I, one()} || I <- lists:seq(0, S - 1)];
                                    F -> factor_elements(Kind, F, S, P)
                                end,
                        [{R * S + R1, C * S + C1, cmul(V, V1, P)} || {R, C, V} <- Acc, {R1, C1, V1} <- Local]
                end, [{0, 0, Cv}], Basis).

%% (a^dagger)^Pp a^Q |n> = sqrt(n!/(n-Q)!) sqrt(m!/(n-Q)!) |m>, m = n - Q + Pp (dropped if m > N); |I><J| of an atom
factor_elements(boson, {Pp, Q}, S, P) ->
    [{N - Q + Pp, N, {pq_mp:sqrt(pq_mp:from(ffact(N, Q) * ffact(N - Q + Pp, Pp), P), P), zero_mp()}}
     || N <- lists:seq(Q, S - 1), N - Q + Pp =< S - 1];
factor_elements(atom, {I, J}, _, _) -> [{I - 1, J - 1, one()}].

%% n (n-1) ... (n-k+1)
ffact(N, K) -> lists:foldl(fun(X, A) -> A * X end, 1, lists:seq(N - K + 1, N)).

%% a sparse matrix from elements, duplicates summed, zeros dropped
sparse(D, Elems, P) ->
    M = lists:foldl(fun({R, C, V}, Acc) ->
                            Acc#{{R, C} => case Acc of #{{R, C} := V0} -> cadd(V0, V, P); _ -> V end}
                    end, #{}, Elems),
    ByRow = maps:fold(fun({R, C}, V, Acc) ->
                              case czero(V) of
                                  true -> Acc;
                                  false -> Acc#{R => [{C, V} | maps:get(R, Acc, [])]}
                              end
                      end, #{}, M),
    {sparse, D, list_to_tuple([lists:sort(maps:get(R, ByRow, [])) || R <- lists:seq(0, D - 1)])}.

elements({sparse, D, Rows}) ->
    [{R, C, V} || R <- lists:seq(0, D - 1), {C, V} <- element(R + 1, Rows)].

%%% States

%% a basis state of one subsystem: {Label, Amplitude} lists; a photon number (from 0) or an atomic level (from 1)
fock(N) -> [{N, one()}].

%% |alpha> of a mode truncated at N photons, renormalised; alpha a number or {complex, Re, Im}
coherent(N, Alpha, P) ->
    A = cnum(Alpha, P),
    {Amps, _} = lists:mapfoldl(fun(K, Pow) ->
                                       V = cscale(Pow, pq_mp:divide(pq_mp:from(1, P), pq_mp:sqrt(pq_mp:from(ffact(K, K), P), P), P), P),
                                       {{K, V}, cmul(Pow, A, P)}
                               end, one(), lists:seq(0, N)),
    Norm = pq_mp:sqrt(lists:foldl(fun({_, V}, S) -> pq_mp:add(S, cmag2(V, P), P) end, zero_mp(), Amps), P),
    [{K, cscale(V, pq_mp:divide(pq_mp:from(1, P), Norm, P), P)} || {K, V} <- Amps].

%% the product state, States #{Name => Amps}; a subsystem not given is in its ground state
state(Basis, States, P) ->
    D = dim(Basis),
    Combos = lists:foldl(fun({N, S, Kind}, Acc) ->
                                 Amps = maps:get(N, States, [{case Kind of atom -> 1; boson -> 0 end, one()}]),
                                 [{Idx * S + local(Kind, L), cmul(V, Va, P)} || {Idx, V} <- Acc, {L, Va} <- Amps]
                         end, [{0, one()}], Basis),
    Map = maps:from_list(Combos),
    list_to_tuple([maps:get(I, Map, zero()) || I <- lists:seq(0, D - 1)]).

local(atom, L) -> L - 1;
local(boson, N) -> N.

%% |psi><psi|, vectorised
dm(Psi, P) ->
    L = tuple_to_list(Psi),
    list_to_tuple([cmul(X, cconj(Y), P) || X <- L, Y <- L]).

%%% Generators

%% -i H, the generator of the Schrödinger equation
hamiltonian_generator(H, P) -> scale_sparse(H, {zero_mp(), pq_mp:from(-1, P)}, P).

%% the Liouvillian, a sparse D^2 x D^2 matrix: (A rho B)_ij = sum A_ik rho_kl B_lj, the element ((i,j),(k,l)) A_ik B_lj
liouvillian(H, Ls, P) ->
    {sparse, D, _} = H,
    MI = {zero_mp(), pq_mp:from(-1, P)},
    PI = {zero_mp(), pq_mp:from(1, P)},
    Half = {pq_mp:from({-1, 2}, P), zero_mp()},
    Id = [{I, I, one()} || I <- lists:seq(0, D - 1)],
    Hs = elements(H),
    Terms = [{Hs, Id, MI}, {Id, Hs, PI}
             | lists:append([begin
                                 Ld = [{C, R, cconj(V)} || {R, C, V} <- elements(L)],
                                 K = elements(matmul(sparse(D, Ld, P), L, P)),
                                 [{elements(L), Ld, one()}, {K, Id, Half}, {Id, K, Half}]
                             end || L <- Ls])],
    Elems = [{I * D + J, K * D + L, cmul(Coef, cmul(Va, Vb, P), P)}
             || {As, Bs, Coef} <- Terms, {I, K, Va} <- As, {L, J, Vb} <- Bs],
    sparse(D * D, Elems, P).

matmul(A, {sparse, D, RowsB}, P) ->
    sparse(D, [{I, J, cmul(X, Y, P)} || {I, K, X} <- elements(A), {J, Y} <- element(K + 1, RowsB)], P).

scale_sparse({sparse, D, Rows}, K, P) ->
    {sparse, D, list_to_tuple([[{C, cmul(V, K, P)} || {C, V} <- R] || R <- tuple_to_list(Rows)])}.

%%% Time evolution

%% y(t) of y' = M y for each t of Times (ascending, from 0; integers, fractions {N, D} or floats): [{T, Y}]
evolve(M, Y0, Times, P) ->
    H0 = step_bound(M, P),
    {Out, _, _} = lists:foldl(fun(T, {Acc, Y, T0}) ->
                                      Tm = pq_mp:from(T, P),
                                      Y2 = advance(M, Y, pq_mp:sub(Tm, T0, P), H0, P),
                                      {[{T, Y2} | Acc], Y2, Tm}
                              end, {[], Y0, zero_mp()}, Times),
    lists:reverse(Out).

sesolve(H, Psi0, Times, P) -> evolve(hamiltonian_generator(H, P), Psi0, Times, P).

%% from a state vector or a vectorised rho
mesolve(H, Ls, Psi0OrRho0, Times, P) ->
    {sparse, D, _} = H,
    Rho0 = case tuple_size(Psi0OrRho0) of D -> dm(Psi0OrRho0, P); _ -> Psi0OrRho0 end,
    evolve(liouvillian(H, Ls, P), Rho0, Times, P).

%% steps of at most 1/|M|, |M| the largest row sum of magnitudes (a power of two below it, as an mp)
step_bound({sparse, _, Rows}, P) ->
    Norm = lists:max([0.0 | [lists:sum([math:sqrt(pq_mp:to_float(cmag2(V, P))) || {_, V} <- R]) || R <- tuple_to_list(Rows)]]),
    case Norm > 0 of
        true -> pq_mp:ldexp(pq_mp:from(1, P), -ceil(math:log2(Norm)));
        false -> pq_mp:from(1, P)
    end.

advance(M, Y, Dt, H0, P) ->
    case pq_mp:cmp(Dt, zero_mp()) of
        1 ->
            H = case pq_mp:cmp(Dt, H0) of 1 -> H0; _ -> Dt end,
            advance(M, taylor_step(M, Y, H, P), pq_mp:sub(Dt, H, P), H0, P);
        _ -> Y
    end.

taylor_step(M, Y, H, P) ->
    Floor = vmag(Y) - P - 4,
    taylor_step(M, Y, Y, H, 1, Floor, P).

taylor_step(M, Sum, Term, H, K, Floor, P) ->
    T2 = vscale(matvec(M, Term, P), {pq_mp:divide(H, pq_mp:from(K, P), P), zero_mp()}, P),
    Sum2 = vadd(Sum, T2, P),
    case vmag(T2) < Floor andalso K > 2 of
        true -> Sum2;
        false -> taylor_step(M, Sum2, T2, H, K + 1, Floor, P)
    end.

%%% Expectation values

%% <psi|O|psi>
expect(O, Psi, P) -> vdot(Psi, matvec(O, Psi, P), P).

%% Tr(O rho), rho vectorised: sum_ij O_ij rho_ji
expect_dm({sparse, D, _} = O, Rho, P) ->
    lists:foldl(fun({I, J, V}, S) -> cadd(S, cmul(V, element(J * D + I + 1, Rho), P), P) end, zero(), elements(O)).

trace(Rho, P) ->
    D = isqrt(tuple_size(Rho)),
    lists:foldl(fun(I, S) -> cadd(S, element(I * D + I + 1, Rho), P) end, zero(), lists:seq(0, D - 1)).

isqrt(N) -> round(math:sqrt(N)).

%%% Steady states: L rho = 0 with Tr rho = 1, by Gaussian elimination on the sparse rows (the first equation replaced
%%% by the trace: L has the one dependent row Tr(L rho) = 0)

steady_state({sparse, N, Rows}, P) ->
    D = isqrt(N),
    TraceRow = maps:from_list([{I * D + I, one()} || I <- lists:seq(0, D - 1)]),
    Eqs = [{TraceRow, one()} | [{maps:from_list(element(R + 1, Rows)), zero()} || R <- lists:seq(1, N - 1)]],
    Sol = eliminate(Eqs, 0, N, [], P),
    list_to_tuple([maps:get(I, Sol, zero()) || I <- lists:seq(0, N - 1)]).

%% forward elimination column by column, partial pivoting on the magnitude; then back substitution
eliminate(_, Col, N, Done, P) when Col >= N -> back(Done, #{}, P);
eliminate(Eqs, Col, N, Done, P) ->
    {With, Without} = lists:partition(fun({Row, _}) -> is_map_key(Col, Row) end, Eqs),
    case With of
        [] -> eliminate(Without, Col + 1, N, Done, P);
        [First | More] ->
            Mag = fun({Row, _}) -> mag2(maps:get(Col, Row), P) end,
            {Pivot, Others} = lists:foldl(fun(E, {Best, Rest}) ->
                                                  case Mag(E) > Mag(Best) of
                                                      true -> {E, [Best | Rest]};
                                                      false -> {Best, [E | Rest]}
                                                  end
                                          end, {First, []}, More),
            {PRow, PRhs} = Pivot,
            Inv = cinv(maps:get(Col, PRow), P),
            Reduced = [begin
                           F = cmul(maps:get(Col, Row), Inv, P),
                           {row_sub(Row, PRow, F, Col, P), csub(Rhs, cmul(F, PRhs, P), P)}
                       end || {Row, Rhs} <- Others],
            eliminate(Reduced ++ Without, Col + 1, N, [{Col, PRow, PRhs} | Done], P)
    end.

%% a magnitude for pivoting (a float of the binary exponent, no overflow)
mag2(V, P) -> mag(cmag2(V, P)).

row_sub(Row, PRow, F, Col, P) ->
    R = maps:fold(fun(C, V, Acc) ->
                          New = csub(maps:get(C, Acc, zero()), cmul(F, V, P), P),
                          case czero(New) of true -> maps:remove(C, Acc); false -> Acc#{C => New} end
                  end, Row, PRow),
    maps:remove(Col, R).

back([], Sol, _) -> Sol;
back([{Col, Row, Rhs} | T], Sol, P) ->
    S = maps:fold(fun(C, V, Acc) when C =/= Col -> csub(Acc, cmul(V, maps:get(C, Sol, zero()), P), P);
                     (_, _, Acc) -> Acc
                  end, Rhs, Row),
    back(T, Sol#{Col => cmul(S, cinv(maps:get(Col, Row), P), P)}, P).

%%% Quantum trajectories (Dalibard, Castin and Mølmer, PRL 68, 580 (1992); Plenio and Knight, RMP 70, 101 (1998)):
%%% psi evolves under H_eff = H - (i/2) sum L_k^dagger L_k until |psi|^2 falls to a uniform random r, then jumps to
%%% L_k psi with probability proportional to |L_k psi|^2. The jump time is found by bisection to 2^-40 of a step.

%% the averages of <O_j> over Ntraj trajectories and their standard errors: [{T, [{Mean, StdErr}]}] (floats, real parts)
trajectories(H, Ls, Psi0, Times, Os, {Ntraj, Seed}, P) ->
    {sparse, D, _} = H,
    K = sparse(D, lists:append([elements(matmul(dagger(L, P), L, P)) || L <- Ls]), P),
    Heff = sparse(D, elements(H) ++ [{R, C, cmul({zero_mp(), pq_mp:from({-1, 2}, P)}, V, P)} || {R, C, V} <- elements(K)], P),
    M = hamiltonian_generator(Heff, P),
    H0 = step_bound(M, P),
    S0 = rand:seed_s(exsss, Seed),
    {Sums, _} = lists:foldl(fun(_, {Acc, S}) ->
                                    {Vals, S2} = trajectory(M, Ls, Psi0, Times, Os, H0, S, P),
                                    {[[{A + X, B + X * X} || {{A, B}, X} <- lists:zip(Ra, Rv)] || {Ra, Rv} <- lists:zip(Acc, Vals)], S2}
                            end, {[[{0.0, 0.0} || _ <- Os] || _ <- Times], S0}, lists:seq(1, Ntraj)),
    [{T, [{A / Ntraj, math:sqrt(max(B / Ntraj - (A / Ntraj) * (A / Ntraj), 0.0) / Ntraj)} || {A, B} <- Row]}
     || {T, Row} <- lists:zip(Times, Sums)].

dagger({sparse, D, _} = L, P) -> sparse(D, [{C, R, cconj(V)} || {R, C, V} <- elements(L)], P).

trajectory(M, Ls, Psi0, Times, Os, H0, S, P) ->
    {R, S1} = rand:uniform_s(S),
    {Vals, _, _, _, S2} =
        lists:foldl(fun(T, {Acc, Psi, T0, Rv, St}) ->
                            Tm = pq_mp:from(T, P),
                            {Psi2, Rv2, St2} = run(M, Ls, Psi, pq_mp:sub(Tm, T0, P), H0, Rv, St, P),
                            Nrm = pq_mp:to_float(vnorm2(Psi2, P)),
                            V = [pq_mp:to_float(element(1, expect(O, Psi2, P))) / Nrm || O <- Os],
                            {[V | Acc], Psi2, Tm, Rv2, St2}
                    end, {[], Psi0, zero_mp(), R, S1}, Times),
    {lists:reverse(Vals), S2}.

%% evolve over Dt with jumps; Rv the current threshold for |psi|^2
run(M, Ls, Psi, Dt, H0, Rv, S, P) ->
    case pq_mp:cmp(Dt, zero_mp()) of
        1 ->
            H = case pq_mp:cmp(Dt, H0) of 1 -> H0; _ -> Dt end,
            Psi2 = taylor_step(M, Psi, H, P),
            case pq_mp:to_float(vnorm2(Psi2, P)) > Rv of
                true -> run(M, Ls, Psi2, pq_mp:sub(Dt, H, P), H0, Rv, S, P);
                false ->
                    {Tj, PsiJ} = bisect(M, Psi, zero_mp(), H, Rv, 40, P),
                    {PsiN, S2} = jump(Ls, PsiJ, S, P),
                    {Rv2, S3} = rand:uniform_s(S2),
                    run(M, Ls, PsiN, pq_mp:sub(Dt, Tj, P), H0, Rv2, S3, P)
            end;
        _ -> {Psi, Rv, S}
    end.

%% the time in [Lo, Hi] where |psi|^2 falls to Rv, and the state there
bisect(M, Psi, _, Hi, _, 0, P) -> {Hi, taylor_step(M, Psi, Hi, P)};
bisect(M, Psi, Lo, Hi, Rv, N, P) ->
    Mid = pq_mp:ldexp(pq_mp:add(Lo, Hi, P), -1),
    case pq_mp:to_float(vnorm2(taylor_step(M, Psi, Mid, P), P)) > Rv of
        true -> bisect(M, Psi, Mid, Hi, Rv, N - 1, P);
        false -> bisect(M, Psi, Lo, Mid, Rv, N - 1, P)
    end.

jump(Ls, Psi, S, P) ->
    Cands = [begin X = matvec(L, Psi, P), {X, pq_mp:to_float(vnorm2(X, P))} end || L <- Ls],
    Total = lists:sum([W || {_, W} <- Cands]),
    {U, S2} = rand:uniform_s(S),
    X = pick(Cands, U * Total),
    Inv = pq_mp:divide(pq_mp:from(1, P), pq_mp:sqrt(vnorm2(X, P), P), P),
    {vscale(X, {Inv, zero_mp()}, P), S2}.

pick([{X, _}], _) -> X;
pick([{X, W} | T], U) -> case U =< W of true -> X; false -> pick(T, U - W) end.

%%% Vectors

matvec({sparse, _, Rows}, V, P) ->
    list_to_tuple([lists:foldl(fun({C, X}, S) -> cadd(S, cmul(X, element(C + 1, V), P), P) end, zero(), R)
                   || R <- tuple_to_list(Rows)]).

vadd(A, B, P) -> list_to_tuple(lists:zipwith(fun(X, Y) -> cadd(X, Y, P) end, tuple_to_list(A), tuple_to_list(B))).
vscale(A, K, P) -> list_to_tuple([cmul(X, K, P) || X <- tuple_to_list(A)]).
vdot(A, B, P) -> lists:foldl(fun({X, Y}, S) -> cadd(S, cmul(cconj(X), Y, P), P) end, zero(), lists:zip(tuple_to_list(A), tuple_to_list(B))).
vnorm2(A, P) -> lists:foldl(fun(X, S) -> pq_mp:add(S, cmag2(X, P), P) end, zero_mp(), tuple_to_list(A)).

%% the binary order of magnitude of the largest component (a very small number for the zero vector)
vmag(V) -> lists:max([-1 bsl 40 | [max(mag(A), mag(B)) || {A, B} <- tuple_to_list(V)]]).
mag({mp, 0, _}) -> -1 bsl 40;
mag({mp, M, E}) -> bitlen(abs(M)) + E.

bitlen(N) -> bitlen(N, 0).
bitlen(0, Acc) -> Acc;
bitlen(N, Acc) when N >= 1 bsl 64 -> bitlen(N bsr 64, Acc + 64);
bitlen(N, Acc) -> bitlen(N bsr 1, Acc + 1).

%%% Complex numbers {Re, Im}

zero_mp() -> {mp, 0, 0}.
zero() -> {zero_mp(), zero_mp()}.
one() -> {{mp, 1, 0}, zero_mp()}.
cadd({A, B}, {C, D}, P) -> {pq_mp:add(A, C, P), pq_mp:add(B, D, P)}.
csub({A, B}, {C, D}, P) -> {pq_mp:sub(A, C, P), pq_mp:sub(B, D, P)}.
cmul({A, B}, {C, D}, P) ->
    {pq_mp:sub(pq_mp:mul(A, C, P), pq_mp:mul(B, D, P), P), pq_mp:add(pq_mp:mul(A, D, P), pq_mp:mul(B, C, P), P)}.
cscale({A, B}, K, P) -> {pq_mp:mul(A, K, P), pq_mp:mul(B, K, P)}.
cconj({A, B}) -> {A, pq_mp:neg(B)}.
cmag2({A, B}, P) -> pq_mp:add(pq_mp:mul(A, A, P), pq_mp:mul(B, B, P), P).
cinv({A, B} = Z, P) -> N = cmag2(Z, P), {pq_mp:divide(A, N, P), pq_mp:neg(pq_mp:divide(B, N, P))}.
czero({{mp, 0, _}, {mp, 0, _}}) -> true;
czero(_) -> false.

%% a value as floats {Re, Im}
to_complex({A, B}) -> {pq_mp:to_float(A), pq_mp:to_float(B)}.

cnum({complex, A, B}, P) -> {pq_mp:from(A, P), pq_mp:from(B, P)};
cnum(A, P) -> {pq_mp:from(A, P), zero_mp()}.

%% a coefficient (an expression of the engine) at the parameters
ceval(N, _, P) when is_integer(N); is_float(N) -> {pq_mp:from(N, P), zero_mp()};
ceval({rat, A, B}, _, P) -> {pq_mp:from({A, B}, P), zero_mp()};
ceval({complex, A, B}, Ps, P) ->
    {Ar, Ai} = ceval(A, Ps, P), {Br, Bi} = ceval(B, Ps, P),
    {pq_mp:sub(Ar, Bi, P), pq_mp:add(Ai, Br, P)};
ceval({sym, X}, Ps, P) ->
    case maps:find(X, Ps) of
        {ok, V} -> cnum(V, P);
        error -> error({pq_lindblad, {unbound, X}})
    end;
ceval({add, L}, Ps, P) -> lists:foldl(fun(T, S) -> cadd(S, ceval(T, Ps, P), P) end, zero(), L);
ceval({mul, L}, Ps, P) -> lists:foldl(fun(T, S) -> cmul(S, ceval(T, Ps, P), P) end, one(), L);
ceval({pow, B, K}, Ps, P) when is_integer(K), K >= 0 -> cpow(ceval(B, Ps, P), K, P);
ceval({pow, B, K}, Ps, P) when is_integer(K) -> cinv(cpow(ceval(B, Ps, P), -K, P), P);
ceval({pow, B, E}, Ps, P) ->
    % a real power of a positive real base
    {Br, _} = ceval(B, Ps, P), {Er, _} = ceval(E, Ps, P),
    {pq_mp:exp(pq_mp:mul(Er, pq_mp:log(Br, P), P), P), zero_mp()};
ceval({conj, X}, Ps, P) -> cconj(ceval(X, Ps, P));
ceval({exp, X}, Ps, P) ->
    {A, B} = ceval(X, Ps, P),
    R = pq_mp:exp(A, P),
    {pq_mp:mul(R, pq_mp:cos(B, P), P), pq_mp:mul(R, pq_mp:sin(B, P), P)};
ceval({F, X}, Ps, P) when F =:= sin; F =:= cos ->
    {A, _} = ceval(X, Ps, P),
    {pq_mp:F(A, P), zero_mp()};
ceval(E, _, _) -> error({pq_lindblad, {unsupported_coefficient, E}}).

cpow(_, 0, _) -> one();
cpow(Z, K, P) -> cmul(Z, cpow(Z, K - 1, P), P).
