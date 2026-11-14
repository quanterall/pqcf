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
