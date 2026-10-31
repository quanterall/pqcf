%% FR-F06: a Gaussian state as a Fock density matrix within a cutoff. Closed forms (coherent state, squeezed vacuum)
%% and cross-validation against the density-matrix backend evolving the same circuit (GT-N01's pairing).
-module(fr_f06_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

c(N, Is) -> {ok, C} = pq_circuit:new(N), C#circuit{instrs = Is}.

rho(S, Cut) ->
    {ok, #fock_dm{rho = R}, Info} = pq_state:gaussian_to_fock_dm(S, Cut),
    {ok, {cmat, Rows}} = pq_core:unpack(R),
    {list_to_tuple([list_to_tuple(Row) || Row <- Rows]), Info}.

at(T, I, J) -> element(J + 1, element(I + 1, T)).

cmul({A, B}, {C, D}) -> {A * C - B * D, A * D + B * C}.
cpow(_, 0) -> {1.0, 0.0};
cpow(Z, N) -> cmul(Z, cpow(Z, N - 1)).
fact(N) -> lists:foldl(fun erlang:'*'/2, 1, lists:seq(1, N)).
close({complex, A, B}, {C, D}, Tol) -> abs(A - C) =< Tol andalso abs(B - D) =< Tol.

%% rho(m, n) = <m|alpha><alpha|n> = e^(-|alpha|^2) alpha^m conj(alpha)^n / sqrt(m! n!), alpha = r e^(i phi)
coherent_test() ->
    {R, Phi} = {0.7, 0.9},
    {ok, S, _} = pq_api:evolve(c(1, [{displacement, 0, R, Phi}]), #{}),
    Cut = 6,
    {T, #{truncated_mass := M}} = rho(S, Cut),
    Al = {R * math:cos(Phi), R * math:sin(Phi)},
    Ac = {R * math:cos(Phi), -R * math:sin(Phi)},
    [begin
         {Re, Im} = cmul(cpow(Al, Mi), cpow(Ac, Ni)),
         F = math:exp(-R * R) / math:sqrt(fact(Mi) * fact(Ni)),
         ?assert(close(at(T, Mi, Ni), {Re * F, Im * F}, 1.0e-14))
     end || Mi <- lists:seq(0, Cut), Ni <- lists:seq(0, Cut)],
    Kept = lists:sum([math:exp(-R * R) * math:pow(R * R, K) / fact(K) || K <- lists:seq(0, Cut)]),
    ?assert(abs(M - (1 - Kept)) < 1.0e-14).

%% S(z)|0>, z = r e^(i phi): c(2n) = (-e^(i phi) tanh r)^n sqrt((2n)!) / (2^n n! sqrt(cosh r))
squeezed_test() ->
    {R, Phi} = {0.6, -0.4},
    {ok, S, _} = pq_api:evolve(c(1, [{squeezing, 0, R, Phi}]), #{}),
    Cut = 8,
    {T, _} = rho(S, Cut),
    Amp = fun(K) when K rem 2 =:= 1 -> {0.0, 0.0};
             (K) ->
                  N = K div 2,
                  {Zr, Zi} = cpow({-math:tanh(R) * math:cos(Phi), -math:tanh(R) * math:sin(Phi)}, N),
                  F = math:sqrt(fact(2 * N)) / (math:pow(2, N) * fact(N) * math:sqrt(math:cosh(R))),
                  {Zr * F, Zi * F}
          end,
    [begin
         {Ar, Ai} = Amp(Mi), {Br, Bi} = Amp(Ni),
         ?assert(close(at(T, Mi, Ni), cmul({Ar, Ai}, {Br, -Bi}), 1.0e-14))
     end || Mi <- lists:seq(0, Cut), Ni <- lists:seq(0, Cut)].

%% the same mixed two-mode state through the Gaussian backend and converted, and through the density-matrix
%% backend at a high cutoff: the elements up to 3 photons agree
cross_validation_test() ->
    Is = [{squeezing, 0, 0.3, 0.2}, {displacement, 1, 0.25, -0.6}, {beamsplitter, 0, 1, 0.7, 0.3}, {loss, 0, 0.8},
          {squeezing, 1, 0.15, 1.1}],
    {ok, G, _} = pq_api:evolve(c(2, Is), #{}),
    {T, #{truncated_mass := M}} = rho(G, 3),
    ?assert(M > 0 andalso M < 0.05),
    {ok, #fock_dm{rho = Rf}, _} = pq_api:evolve(c(2, Is), #{backend => fock_dm, cutoff => 14}),
    {ok, {cmat, Rows}} = pq_core:unpack(Rf),
    F = list_to_tuple([list_to_tuple(R) || R <- Rows]),
    D = pq_state:basis_size(2, 3),
    Err = lists:max([begin
                         {ok, Pi} = pq_state:pattern(I, 2, 3), {ok, Pj} = pq_state:pattern(J, 2, 3),
                         {ok, Ki} = pq_state:index(Pi, 2, 14), {ok, Kj} = pq_state:index(Pj, 2, 14),
                         {complex, A, B} = at(T, I, J), {complex, C, E} = at(F, Ki, Kj),
                         abs(A - C) + abs(B - E)
                     end || I <- lists:seq(0, D - 1), J <- lists:seq(0, D - 1)]),
    io:format(user, "~nFR-F06: Gaussian -> Fock against the density-matrix backend: max error ~.3g~n", [Err]),
    ?assert(Err < 1.0e-9).

limits_test() ->
    {ok, S, _} = pq_api:evolve(c(1, [{squeezing, 0, 0.3, 0}]), #{}),
    ?assertMatch({error, {bad_shape, _}}, pq_state:gaussian_to_fock_dm(S, -1)),
    {ok, S8, _} = pq_api:evolve(c(8, []), #{}),
    ?assertMatch({error, {dimension_limit, 6435, 4096}}, pq_state:gaussian_to_fock_dm(S8, 7)).
