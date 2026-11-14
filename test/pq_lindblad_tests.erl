%% pq_lindblad (quantum optics plan, Q2): master equations against closed forms. Jaynes-Cummings Rabi oscillations and
%% collapse-revival (Eberly, Narozhny and Sanchez-Mondragon, PRL 44, 1323 (1980)), a damped coherent state, atomic decay,
%% the steady state of resonance fluorescence, and quantum trajectories against the master equation.
-module(pq_lindblad_tests).
-include_lib("eunit/include/eunit.hrl").

-define(L, pq_lindblad).
-define(Q, pq_qop).

re(Z) -> element(1, ?L:to_complex(Z)).

%% at P bits, |x - y| < 2^-(P - 20) relative to 1
mp_close(Z, Want, P) ->
    {A, _} = Z,
    D = pq_mp:sub(A, Want, P),
    case D of
        {mp, 0, _} -> true;
        {mp, M, E} -> bitlen(abs(M)) + E < -(P - 24)
    end.

bitlen(0) -> 0;
bitlen(N) -> 1 + bitlen(N bsr 1).

jc() ->
    % H = g (ad sigma- + a sigma+), on resonance in the rotating frame
    Sp = #{q => {atom, 2}},
    H = ?Q:with_spaces(Sp, fun() -> ?Q:scale({sym, g}, ?Q:add(?Q:mul(?Q:ad(c), ?Q:sm(q)), ?Q:mul(?Q:a(c), ?Q:sp(q)))) end),
    Pe = ?Q:with_spaces(Sp, fun() -> ?Q:proj(q, 2, 2) end),
    {Sp, H, Pe}.

%% vacuum Rabi oscillations, |e, 0> -> P_e(t) = cos^2(g t), to 200 bits
rabi_test() ->
    P = 200,
    {Sp, H, Pe} = jc(),
    B = ?L:basis(Sp, #{c => 3}),
    Hm = ?L:op(B, Sp, H, {#{g => 1}, P}),
    Pem = ?L:op(B, Sp, Pe, {#{}, P}),
    Psi0 = ?L:state(B, #{q => ?L:fock(2)}, P),
    [begin
         Want = pq_mp:mul(pq_mp:cos(pq_mp:from(T, P), P), pq_mp:cos(pq_mp:from(T, P), P), P),
         ?assert(mp_close(?L:expect(Pem, Psi, P), Want, P))
     end || {T, Psi} <- ?L:sesolve(Hm, Psi0, [{1, 2}, 1, 3], P)].

%% collapse and revival: |e> with the field coherent, P_e(t) = sum_n p_n cos^2(g sqrt(n+1) t), p_n Poisson of the
%% truncated coherent state (to 112 bits), at t = 0 .. 30 (the revival near 2 pi |alpha| = 12.6 for alpha = 2)
collapse_revival_test_() ->
    {timeout, 120, fun() ->
        P = 120, N = 20,
        {Sp, H, Pe} = jc(),
        B = ?L:basis(Sp, #{c => N}),
        Hm = ?L:op(B, Sp, H, {#{g => 1}, P}),
        Pem = ?L:op(B, Sp, Pe, {#{}, P}),
        % the field truncated one below the basis, so that |e, n> always has its partner |g, n + 1>: the closed form is then
        % exact for the truncated state
        Amps = ?L:coherent(N - 1, 2, P),
        Psi0 = ?L:state(B, #{q => ?L:fock(2), c => Amps}, P),
        Times = [0, 3, 6, 10, 13, 20, 30],
        [begin
             Want = lists:foldl(fun({K, V}, S) ->
                                        C = pq_mp:cos(pq_mp:mul(pq_mp:sqrt(pq_mp:from(K + 1, P), P), pq_mp:from(T, P), P), P),
                                        pq_mp:add(S, pq_mp:mul(?L:cmag2(V, P), pq_mp:mul(C, C, P), P), P)
                                end, {mp, 0, 0}, Amps),
             ?assert(mp_close(?L:expect(Pem, Psi, P), Want, P - 8))
         end || {T, Psi} <- ?L:sesolve(Hm, Psi0, Times, P)]
    end}.

%% a damped cavity, L = sqrt(kappa) a: a coherent state stays coherent, <a> = alpha e^(-kappa t/2) and
%% <ad a> = |alpha|^2 e^(-kappa t). The truncation (N = 12, alpha = 1) is exact for the decay: the photon number never grows.
damped_cavity_test_() ->
    {timeout, 120, fun() ->
        P = 100, N = 12,
        Sp = #{},
        B = ?L:basis(Sp, #{c => N}),
        Hm = ?L:op(B, Sp, ?Q:scale({sym, w}, ?Q:mul(?Q:ad(c), ?Q:a(c))), {#{w => 0}, P}),
        Lm = ?L:op(B, Sp, ?Q:scale({pow, {sym, kappa}, {rat, 1, 2}}, ?Q:a(c)), {#{kappa => {1, 2}}, P}),
        Am = ?L:op(B, Sp, ?Q:a(c), {#{}, P}),
        Nm = ?L:op(B, Sp, ?Q:mul(?Q:ad(c), ?Q:a(c)), {#{}, P}),
        Amps = ?L:coherent(N, 1, P),
        Psi0 = ?L:state(B, #{c => Amps}, P),
        % the truncated state's own <a> and <n> at t = 0, then the exact decay
        A0 = element(1, ?L:expect(Am, Psi0, P)), N0 = element(1, ?L:expect(Nm, Psi0, P)),
        [begin
             E = pq_mp:exp(pq_mp:neg(pq_mp:from({T, 2}, P)), P),
             ?assert(mp_close(?L:expect_dm(Nm, Rho, P), pq_mp:mul(N0, E, P), P)),
             ?assert(mp_close(?L:expect_dm(Am, Rho, P), pq_mp:mul(A0, pq_mp:sqrt(E, P), P), P - 4)),
             ?assert(mp_close(?L:trace(Rho, P), pq_mp:from(1, P), P))
         end || {T, Rho} <- ?L:mesolve(Hm, [Lm], Psi0, [1, 2, 4], P)]
    end}.

%% spontaneous decay, L = sqrt(gamma) sigma-: P_e(t) = e^(-gamma t)
decay_test() ->
    P = 150,
    Sp = #{q => {atom, 2}},
    B = ?L:basis(Sp, #{}),
    Hm = ?L:op(B, Sp, ?Q:with_spaces(Sp, fun() -> ?Q:num(0) end), {#{}, P}),
    Lm = ?L:op(B, Sp, ?Q:with_spaces(Sp, fun() -> ?Q:sm(q) end), {#{}, P}),
    Pem = ?L:op(B, Sp, ?Q:with_spaces(Sp, fun() -> ?Q:proj(q, 2, 2) end), {#{}, P}),
    Psi0 = ?L:state(B, #{q => ?L:fock(2)}, P),
    [?assert(mp_close(?L:expect_dm(Pem, Rho, P), pq_mp:exp(pq_mp:from(-T, P), P), P))
     || {T, Rho} <- ?L:mesolve(Hm, [Lm], Psi0, [1, 3, 5], P)].

%% resonance fluorescence: H = (Omega/2)(sigma+ + sigma-) + Delta sigma+ sigma-... here on resonance with the drive,
%% L = sqrt(gamma) sigma-: the steady state rho_ee = (Omega^2/4) / (gamma^2/4 + Omega^2/2) (Mollow; Loudon, The
%% Quantum Theory of Light, ch. 7); with Delta: (Omega^2/4) / (Delta^2 + gamma^2/4 + Omega^2/2)
steady_state_test() ->
    P = 150,
    Sp = #{q => {atom, 2}},
    B = ?L:basis(Sp, #{}),
    [begin
         H = ?Q:with_spaces(Sp, fun() ->
                 ?Q:add(?Q:scale({mul, [{rat, 1, 2}, {sym, omega}]}, ?Q:add(?Q:sp(q), ?Q:sm(q))),
                        ?Q:scale({mul, [-1, {sym, delta}]}, ?Q:proj(q, 2, 2)))
             end),
         Hm = ?L:op(B, Sp, H, {#{omega => Om, delta => Dl}, P}),
         Lm = ?L:op(B, Sp, ?Q:with_spaces(Sp, fun() -> ?Q:sm(q) end), {#{}, P}),
         Pem = ?L:op(B, Sp, ?Q:with_spaces(Sp, fun() -> ?Q:proj(q, 2, 2) end), {#{}, P}),
         Rho = ?L:steady_state(?L:liouvillian(Hm, [Lm], P), P),
         % gamma = 1
         Want = pq_mp:from({Om * Om * 4, 4 * (4 * Dl * Dl + 1 + 2 * Om * Om)}, P),
         ?assert(mp_close(?L:expect_dm(Pem, Rho, P), Want, P - 8))
     end || {Om, Dl} <- [{1, 0}, {3, 0}, {2, 1}]].

%% the driven damped cavity: H = E (ad + a), L = sqrt(kappa) a, steady state coherent with alpha = -2 i E / kappa
driven_cavity_test_() ->
    {timeout, 120, fun() ->
        P = 100, N = 14,
        Sp = #{},
        B = ?L:basis(Sp, #{c => N}),
        Hm = ?L:op(B, Sp, ?Q:scale({sym, e}, ?Q:add(?Q:ad(c), ?Q:a(c))), {#{e => {1, 4}}, P}),
        Lm = ?L:op(B, Sp, ?Q:a(c), {#{}, P}),
        Nm = ?L:op(B, Sp, ?Q:mul(?Q:ad(c), ?Q:a(c)), {#{}, P}),
        Am = ?L:op(B, Sp, ?Q:a(c), {#{}, P}),
        Rho = ?L:steady_state(?L:liouvillian(Hm, [Lm], P), P),
        % alpha = -i/2, |alpha|^2 = 1/4; the truncation at 14 photons errs by ~ |alpha|^30 / 15! ~ 1e-21
        {Ar, Ai} = ?L:to_complex(?L:expect_dm(Am, Rho, P)),
        ?assert(abs(Ar) < 1.0e-18),
        ?assert(abs(Ai + 0.5) < 1.0e-18),
        ?assert(abs(re(?L:expect_dm(Nm, Rho, P)) - 0.25) < 1.0e-18)
    end}.

%% trajectories: the mean of P_e over 400 trajectories of a driven, decaying atom against the master equation, within
%% four standard errors
trajectories_test_() ->
    {timeout, 300, fun() ->
        P = 53,
        Sp = #{q => {atom, 2}},
        B = ?L:basis(Sp, #{}),
        H = ?Q:with_spaces(Sp, fun() -> ?Q:scale({rat, 1, 1}, ?Q:add(?Q:sp(q), ?Q:sm(q))) end),
        Hm = ?L:op(B, Sp, H, {#{}, P}),
        Lm = ?L:op(B, Sp, ?Q:with_spaces(Sp, fun() -> ?Q:sm(q) end), {#{}, P}),
        Pem = ?L:op(B, Sp, ?Q:with_spaces(Sp, fun() -> ?Q:proj(q, 2, 2) end), {#{}, P}),
        Psi0 = ?L:state(B, #{}, P),
        Times = [{1, 2}, 1, 2, 4],
        Me = [re(?L:expect_dm(Pem, Rho, P)) || {_, Rho} <- ?L:mesolve(Hm, [Lm], Psi0, Times, P)],
        Tr = ?L:trajectories(Hm, [Lm], Psi0, Times, [Pem], {400, {1, 2, 3}}, P),
        [?assert(abs(M - X) < 4 * max(Se, 0.005)) || {M, {_, [{X, Se}]}} <- lists:zip(Me, Tr)]
    end}.
