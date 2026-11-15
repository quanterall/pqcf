%% pq_ode (decision record SY2, phase 3): the soft-wall and hard-wall models of holographic QCD, from the action to
%% the meson spectrum.
%% - Soft wall (Karch, Katz, Son, Stephanov, hep-ph/0602229, https://arxiv.org/abs/hep-ph/0602229): vector mesons,
%%   dilaton z^2 (c = 1): m_n^2 = 4 (n + 1).
%% - Hard wall (Erlich, Katz, Son, Stephanov, hep-ph/0501128, https://arxiv.org/abs/hep-ph/0501128): vector mesons
%%   on 0 < z =< z_m with psi'(z_m) = 0: m_n z_m = the zeros of J_0.
-module(pq_ode_tests).
-include_lib("eunit/include/eunit.hrl").

-define(Z, {sym, z}).
-define(PSI, {sym, psi}).
-define(DPSI, {sym, dpsi}).

%% the radial Lagrangian of a vector mode: (1/2) e^(-Phi) (1/z) (psi'^2 - m^2 psi^2), Phi the dilaton
lagrangian(Dilaton) ->
    {mul, [{rat, 1, 2}, {exp, {mul, [-1, Dilaton]}}, {pow, ?Z, -1},
           {add, [{pow, ?DPSI, 2}, {mul, [-1, {sym, m2}, {pow, ?PSI, 2}]}]}]}.

soft_wall_potential_test() ->
    {ok, #{potential := V}} = pq_ode:schrodinger(lagrangian({pow, ?Z, 2}), z, psi, dpsi, m2),
    % V = z^2 + 3/(4 z^2)
    ?assertEqual({ok, true}, pq_symbolic:equiv(V, {add, [{pow, ?Z, 2}, {mul, [{rat, 3, 4}, {pow, ?Z, -2}]}]}, #{})).

hard_wall_potential_test() ->
    {ok, #{potential := V}} = pq_ode:schrodinger(lagrangian(0), z, psi, dpsi, m2),
    ?assertEqual({ok, true}, pq_symbolic:equiv(V, {mul, [{rat, 3, 4}, {pow, ?Z, -2}]}, #{})).

%% the decimal strings agree to D significant digits
agree(A, B, D) ->
    F = fun(S) -> [C || C <- S, C >= $0, C =< $9] end,
    {Ma, Ea} = mant_exp(A), {Mb, Eb} = mant_exp(B),
    Ea =:= Eb andalso lists:sublist(F(Ma), D) =:= lists:sublist(F(Mb), D).

mant_exp(S) -> [M, E] = string:split(S, "e"), {M, list_to_integer(E)}.

soft_wall_spectrum_test_() ->
    {timeout, 600, fun() ->
        {ok, #{potential := V}} = pq_ode:schrodinger(lagrangian({pow, ?Z, 2}), z, psi, dpsi, m2),
        % the decaying solution at z = 10: e^(-100) below the 30 digits asked for
        {ok, Es} = pq_ode:eigenvalues(V, z, #{right => {dirichlet, 10}, range => {1, 17}, grid => 64, digits => 30}),
        ?assertEqual(4, length(Es)),
        [?assert(agree(E, pq_mp:to_string(pq_mp:from(4 * (N + 1), 120), 30), 28), {N, E})
         || {N, E} <- lists:enumerate(0, Es)]
    end}.

hard_wall_spectrum_test_() ->
    {timeout, 600, fun() ->
        {ok, #{potential := V}} = pq_ode:schrodinger(lagrangian(0), z, psi, dpsi, m2),
        % z_m = 1; psi'(1) = 0 is phi'(1) + phi(1)/2 = 0 (psi = z^(1/2) phi)
        {ok, Es} = pq_ode:eigenvalues(V, z, #{right => {robin, 1, {1, 2}, 1}, range => {1, 80}, grid => 80, digits => 30}),
        ?assertEqual(3, length(Es)),
        [begin
             J = pq_mp:besselj_zero(0, K, 130),
             ?assert(agree(E, pq_mp:to_string(pq_mp:mul(J, J, 130), 30), 28), {K, E})
         end || {K, E} <- lists:enumerate(1, Es)]
    end}.

%% the second, independent solver (Chebyshev collocation, inverse iteration) agrees with the closed forms to 25 digits
collocation_test_() ->
    {timeout, 600, fun() ->
        Vs = {add, [{pow, ?Z, 2}, {mul, [{rat, 3, 4}, {pow, ?Z, -2}]}]},
        {ok, Es} = pq_ode:eigenvalues_collocation(Vs, z, #{right => {dirichlet, 10}, range => {1, 17}, grid => 64, digits => 25, points => 64}),
        [?assert(agree(E, pq_mp:to_string(pq_mp:from(4 * (N + 1), 120), 25), 25), {soft, N, E}) || {N, E} <- lists:enumerate(0, Es)],
        Vh = {mul, [{rat, 3, 4}, {pow, ?Z, -2}]},
        {ok, Eh} = pq_ode:eigenvalues_collocation(Vh, z, #{right => {robin, 1, {1, 2}, 1}, range => {1, 80}, grid => 80, digits => 25, points => 32}),
        [begin
             J = pq_mp:besselj_zero(0, K, 130),
             ?assert(agree(E, pq_mp:to_string(pq_mp:mul(J, J, 130), 25), 24), {hard, K, E})
         end || {K, E} <- lists:enumerate(1, Eh)]
    end}.

%% the hard wall's axial and pion sector (comparison/holographic/hardwall.erl), against an independent floating-point
%% solution (SciPy solve_ivp and brentq, double precision): m_a1 = 1358.2431775 (A), 1219.5033847 (B) MeV; f_pi in the
%% chiral limit 90.903 (A), 82.255 (B) MeV; and the paper's F_rho^(1/2) = 329 (A), 353 (B), F_a1^(1/2) = 486, 440
hardwall_test_() ->
    {timeout, 900, fun() ->
        Dir = filename:dirname(code:which(?MODULE)),
        Src = hd([F || F <- [filename:join([Dir, "..", "comparison", "holographic", "hardwall.erl"]),
                             filename:join([Dir, "..", "..", "..", "..", "comparison", "holographic", "hardwall.erl"]),
                             "comparison/holographic/hardwall.erl"], filelib:is_regular(F)]),
        {ok, hardwall, Bin} = compile:file(Src, [binary]),
        {module, hardwall} = code:load_binary(hardwall, Src, Bin),
        [A, B] = hardwall:models(),
        F = fun(M, K) -> pq_mp:to_float(maps:get(K, M)) end,
        Near = fun(X, Y, Tol) -> abs(X - Y) =< Tol * abs(Y) end,
        ?assert(Near(F(A, m_a1), 1358.2431775, 1.0e-9)),
        ?assert(Near(F(B, m_a1), 1219.5033847, 1.0e-9)),
        % the chiral-limit f_pi: the exact power series of A(z) (comparison/holographic/scipy_hardwall.py agrees to 1e-9)
        ?assert(Near(F(A, f_pi), 90.9033372812686, 1.0e-12)),
        ?assert(Near(F(B, f_pi), 82.2545523007217, 1.0e-12)),
        [?assert(Near(F(M, K), V, 0.005)) || {M, K, V} <- [{A, f_rho, 329}, {B, f_rho, 353}, {A, f_a1, 486}, {B, f_a1, 440},
                                                           {B, m_rho, 832}, {A, m_pi, 139.6}, {B, m_pi, 141}]]
    end}.

%% Sakai-Sugimoto (hep-th/0412141, https://arxiv.org/abs/hep-th/0412141): the vector and axial mesons,
%% -K^(1/3) d_z (K d_z psi) = lambda psi, K = 1 + z^2, in u = K^(1/3): u (u^3 - 1) psi'' + (5u^3/2 - 1) psi' + (9u/4)
%% lambda psi = 0 on [1, infinity), exponent 0 (even, vectors) or 1/2 (odd, axials) at u = 1, decaying u^(-3/2) at
%% infinity. The paper: lambda = 0.67, 1.6, 2.9, 4.5. An independent SciPy shooting in z, extrapolated in its cutoff
%% (comparison/holographic/scipy_sakai.py), gives 0.669314492244 and 1.568765751391; the result must not depend on the
%% match point.
sakai_sugimoto_test_() ->
    {timeout, 900, fun() ->
        U = {sym, u},
        C = #{p2 => {add, [{pow, U, 4}, {mul, [-1, U]}]}, p1 => {add, [{mul, [{rat, 5, 2}, {pow, U, 3}]}, -1]}, p0 => 0,
              q => {mul, [{rat, 9, 4}, U]}},
        Opts = #{right => larger, range => {0.1, 5}, grid => 40, digits => 20},
        {ok, V2} = pq_ode:eigenvalues_sl(C, u, Opts#{left => {1, smaller}, match => 2}),
        {ok, V3} = pq_ode:eigenvalues_sl(C, u, Opts#{left => {1, smaller}, match => 3}),
        {ok, Ax} = pq_ode:eigenvalues_sl(C, u, Opts#{left => {1, larger}, match => 2}),
        ?assertEqual(V2, V3),
        [L1, L3] = [list_to_float(hd(string:split(S, "e"))) * math:pow(10, list_to_integer(lists:last(string:split(S, "e")))) || S <- V2],
        [L2, L4] = [list_to_float(hd(string:split(S, "e"))) * math:pow(10, list_to_integer(lists:last(string:split(S, "e")))) || S <- Ax],
        ?assert(abs(L1 - 0.669314492244) < 1.0e-11),
        ?assert(abs(L2 - 1.568765751391) < 1.0e-11),
        [?assert(abs(X - Y) < 0.05 * Y) || {X, Y} <- [{L1, 0.67}, {L2, 1.6}, {L3, 2.9}, {L4, 4.5}]]
    end}.

%% A Frobenius series with gaps: -phi'' + (3/(4 u^2) + u^4) phi = 0 has phi = u^(-1/2) A(u), A = sum a_n u^n with
%% a_0 = 1, a_n = a_(n-6) / (n (n-2)): every coefficient but each sixth is zero. The series must not stop at the first
%% two zeros in a row (it once stopped after a_6, an error of 1e-6 at any precision: the hard wall's f_pi).
frobenius_gaps_test() ->
    U = {sym, u},
    V = {add, [{mul, [{rat, 3, 4}, {pow, U, -2}]}, {pow, U, 4}]},
    {ok, #{phi := Phi, dphi := Dphi}} = pq_ode:frobenius_solution(V, u, 0, #{b => 1, digits => 30, exponent => smaller}),
    % A(1) and A'(1) exactly, as fractions: phi(1) = A(1), phi'(1) = A'(1) - A(1)/2
    Terms = lists:foldl(fun(N, [{Pn, Pd} | _] = Acc) -> [{Pn, Pd * N * (N - 2)} | Acc] end, [{1, 1}], lists:seq(6, 120, 6)),
    Ns = lists:seq(120, 0, -6),
    Sum = fun(F) -> lists:foldl(fun({{Tn, Td}, N}, {An, Ad}) -> {An * Td + F(N) * Tn * Ad, Ad * Td} end, {0, 1}, lists:zip(Terms, Ns)) end,
    A = Sum(fun(_) -> 1 end),
    Ap = Sum(fun(N) -> N end),
    P = 120,
    Mp = fun(S) -> pq_mp:from_string(S, P) end,
    Rel = fun(X, Y) -> pq_mp:to_float(pq_mp:divide(pq_mp:sub(X, Y, P), Y, P)) end,
    Aexact = pq_mp:from(A, P),
    Dexact = pq_mp:sub(pq_mp:from(Ap, P), pq_mp:ldexp(Aexact, -1), P),
    ?assert(abs(Rel(Mp(Phi), Aexact)) < 1.0e-28),
    ?assert(abs(Rel(Mp(Dphi), Dexact)) < 1.0e-28).
