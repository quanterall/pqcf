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
        ?assert(Near(F(A, f_pi), 90.903, 2.0e-4)),
        ?assert(Near(F(B, f_pi), 82.255, 2.0e-4)),
        [?assert(Near(F(M, K), V, 0.005)) || {M, K, V} <- [{A, f_rho, 329}, {B, f_rho, 353}, {A, f_a1, 486}, {B, f_a1, 440},
                                                           {B, m_rho, 832}, {A, m_pi, 139.6}, {B, m_pi, 141}]]
    end}.

%% Sakai-Sugimoto (hep-th/0412141, https://arxiv.org/abs/hep-th/0412141): the vector and axial mesons,
%% -K^(1/3) d_z (K d_z psi) = lambda psi, K = 1 + z^2, in u = K^(1/3): u (u^3 - 1) psi'' + (5u^3/2 - 1) psi' + (9u/4)
%% lambda psi = 0 on [1, infinity), exponent 0 (even, vectors) or 1/2 (odd, axials) at u = 1, decaying u^(-3/2) at
%% infinity. The paper: lambda = 0.67, 1.6, 2.9, 4.5. An independent SciPy shooting in z, extrapolated in its cutoff,
%% gives 0.6693144 and 1.5687644 (to its 1e-7); the result must not depend on the match point.
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
        ?assert(abs(L1 - 0.6693144) < 2.0e-7),
        ?assert(abs(L2 - 1.5687644) < 2.0e-6),
        [?assert(abs(X - Y) < 0.05 * Y) || {X, Y} <- [{L1, 0.67}, {L2, 1.6}, {L3, 2.9}, {L4, 4.5}]]
    end}.
