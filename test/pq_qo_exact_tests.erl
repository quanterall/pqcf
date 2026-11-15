%% pq_qo_exact (quantum optics plan, Q3): exact steady states and spectra against SymPy's dense, exact solutions in the
%% Fock basis (test/oracle/gt_qo_exact.py, test/fixtures/gt_qo_exact.json), equal as rationals; the Mollow triplet's
%% strong-field limit; resonance fluorescence's g2(tau) against its closed form (Kimble and Mandel, Phys. Rev. A
%% 13, 2123 (1976)); and the driven Kerr cavity against Drummond and Walls's exact solution.
-module(pq_qo_exact_tests).
-include_lib("eunit/include/eunit.hrl").

-define(X, pq_qo_exact).
-define(Q, pq_qop).
-define(L, pq_lindblad).

fixture() ->
    Dir = filename:dirname(code:which(?MODULE)),
    Paths = [filename:join([Dir, "..", "test", "fixtures", "gt_qo_exact.json"]),
             filename:join([Dir, "..", "..", "..", "..", "test", "fixtures", "gt_qo_exact.json"]),
             "test/fixtures/gt_qo_exact.json"],
    [P | _] = [P || P <- Paths, filelib:is_regular(P)],
    {ok, B} = file:read_file(P),
    json:decode(B).

%% "p/q" or "p" as a rational {P, Q}
q(B) when is_binary(B) ->
    case binary:split(B, <<"/">>) of
        [N, D] -> pq_sym_exact:rat(binary_to_integer(N), binary_to_integer(D));
        [N] -> {binary_to_integer(N), 1}
    end.
c([Re, Im]) -> {q, q(Re), q(Im)}.
par(B) -> q(B).

sym(X) -> {sym, X}.

two_level() ->
    Sp = #{q => {atom, 2}, domains => #{gamma => positive}},
    H = ?Q:with_spaces(Sp, fun() ->
            ?Q:add(?Q:scale({mul, [{rat, 1, 2}, sym(omega)]}, ?Q:add(?Q:sp(q), ?Q:sm(q))),
                   ?Q:scale({mul, [-1, sym(delta)]}, ?Q:proj(q, 2, 2)))
        end),
    L = ?Q:with_spaces(Sp, fun() -> ?Q:scale({pow, sym(gamma), {rat, 1, 2}}, ?Q:sm(q)) end),
    {Sp, H, L}.

op(B, Sp, E, P, F) -> ?X:operator(B, Sp, ?Q:with_spaces(Sp, fun() -> E() end), P, F).

two_level_test_() ->
    {timeout, 120, fun() ->
        #{<<"cases">> := Cases} = fixture(),
        {Sp, H, L} = two_level(),
        B = ?L:basis(Sp, #{}),
        [begin
             P = #{omega => par(Om), delta => par(Dl), gamma => par(G)},
             Lv = ?X:liouvillian(B, Sp, H, [L], P),
             Rho = ?X:steady_state(Lv, B),
             Pe = op(B, Sp, fun() -> ?Q:proj(q, 2, 2) end, P, trace),
             Sm = op(B, Sp, fun() -> ?Q:sm(q) end, P, trace),
             ?assertEqual(c(WantPe), ?X:expect(Pe, Rho, 2)),
             ?assertEqual(c(WantSm), ?X:expect(Sm, Rho, 2)),
             Spec = ?X:spectrum(Lv, {op(B, Sp, fun() -> ?Q:sp(q) end, P, trace), op(B, Sp, fun() -> ?Q:sm(q) end, P, left),
                                     Sm, Rho, B}, [q(W) || [W, _] <- WantSpec], #{}),
             ?assertEqual([q(S) || [_, S] <- WantSpec], [S || {_, S} <- Spec])
         end || #{<<"model">> := <<"two_level">>, <<"params">> := #{<<"omega">> := Om, <<"delta">> := Dl, <<"gamma">> := G},
                  <<"pe">> := WantPe, <<"sm">> := WantSm, <<"spectrum">> := WantSpec} <- Cases]
    end}.

kerr_test_() ->
    {timeout, 120, fun() ->
        #{<<"cases">> := Cases} = fixture(),
        Sp = #{domains => #{kappa => positive}},
        A = ?Q:a(c), Ad = ?Q:ad(c),
        H = ?Q:add([?Q:scale(sym(delta), ?Q:mul(Ad, A)), ?Q:scale(sym(chi), ?Q:mul([Ad, Ad, A, A])), ?Q:scale(sym(e), ?Q:add(A, Ad))]),
        L = ?Q:scale({pow, sym(kappa), {rat, 1, 2}}, A),
        [begin
             B = ?L:basis(Sp, #{c => N}),
             P = #{e => par(E), chi => par(Chi), kappa => par(K), delta => par(Dl)},
             Rho = ?X:steady_state(?X:liouvillian(B, Sp, H, [L], P), B),
             D = N + 1,
             ?assertEqual(c(WantN), ?X:expect(op(B, Sp, fun() -> ?Q:mul(Ad, A) end, P, trace), Rho, D)),
             ?assertEqual(c(WantA), ?X:expect(op(B, Sp, fun() -> A end, P, trace), Rho, D)),
             ?assertEqual(c(WantN2), ?X:expect(op(B, Sp, fun() -> ?Q:mul([Ad, Ad, A, A]) end, P, trace), Rho, D))
         end || #{<<"model">> := <<"kerr">>, <<"params">> := #{<<"n">> := N, <<"e">> := E, <<"chi">> := Chi,
                                                               <<"kappa">> := K, <<"delta">> := Dl},
                  <<"n">> := WantN, <<"a">> := WantA, <<"n2">> := WantN2} <- Cases]
    end}.

jaynes_cummings_test_() ->
    {timeout, 120, fun() ->
        #{<<"cases">> := Cases} = fixture(),
        [#{<<"params">> := #{<<"n">> := N, <<"g">> := G, <<"e">> := E, <<"kappa">> := K, <<"gamma">> := Gm},
           <<"n">> := WantN, <<"pe">> := WantPe, <<"a">> := WantA}] = [C || #{<<"model">> := <<"jc">>} = C <- Cases],
        Sp = #{q => {atom, 2}, domains => #{kappa => positive, gamma => positive}},
        {H, Ls} = ?Q:with_spaces(Sp, fun() ->
            {?Q:add(?Q:scale(sym(g), ?Q:add(?Q:mul(?Q:ad(c), ?Q:sm(q)), ?Q:mul(?Q:a(c), ?Q:sp(q)))),
                    ?Q:scale(sym(e), ?Q:add(?Q:a(c), ?Q:ad(c)))),
             [?Q:scale({pow, sym(kappa), {rat, 1, 2}}, ?Q:a(c)), ?Q:scale({pow, sym(gamma), {rat, 1, 2}}, ?Q:sm(q))]}
        end),
        B = ?L:basis(Sp, #{c => N}),
        P = #{g => par(G), e => par(E), kappa => par(K), gamma => par(Gm)},
        Rho = ?X:steady_state(?X:liouvillian(B, Sp, H, Ls, P), B),
        D = 2 * (N + 1),
        ?assertEqual(c(WantN), ?X:expect(op(B, Sp, fun() -> ?Q:mul(?Q:ad(c), ?Q:a(c)) end, P, trace), Rho, D)),
        ?assertEqual(c(WantPe), ?X:expect(op(B, Sp, fun() -> ?Q:proj(q, 2, 2) end, P, trace), Rho, D)),
        ?assertEqual(c(WantA), ?X:expect(op(B, Sp, fun() -> ?Q:a(c) end, P, trace), Rho, D))
    end}.

%% the Mollow triplet: for Omega >> gamma the incoherent spectrum has a central peak and two side peaks at about
%% +-Omega, the central one three times as high (widths gamma/2 and 3 gamma/4, weights 1/2 and 1/4)
mollow_triplet_test() ->
    {Sp, H, L} = two_level(),
    B = ?L:basis(Sp, #{}),
    P = #{omega => {200, 1}, delta => {0, 1}, gamma => {1, 1}},
    Lv = ?X:liouvillian(B, Sp, H, [L], P),
    Rho = ?X:steady_state(Lv, B),
    Sm = op(B, Sp, fun() -> ?Q:sm(q) end, P, trace),
    [{_, S0}, {_, Sside}, {_, Sbelow}, {_, Sabove}, {_, Smid}] =
        ?X:spectrum(Lv, {op(B, Sp, fun() -> ?Q:sp(q) end, P, trace), op(B, Sp, fun() -> ?Q:sm(q) end, P, left), Sm, Rho, B},
                    [{1, 1000}, {200, 1}, {195, 1}, {205, 1}, {100, 1}], #{}),
    F = fun pq_sym_exact:to_float/1,
    ?assert(abs(F(S0) / F(Sside) - 3) < 0.01),
    ?assert(F(Sside) > F(Sbelow) andalso F(Sside) > F(Sabove)),
    ?assert(F(Smid) < F(Sside) / 100).

%% g2(tau) of resonance fluorescence, H = (Omega/2)(sigma+ + sigma-), L = sqrt(gamma) sigma-:
%%   g2(tau) = 1 - e^(-3 gamma tau/4) (cos(mu tau) + (3 gamma/(4 mu)) sin(mu tau)), mu = sqrt(Omega^2 - gamma^2/16),
%% g2(0) = 0 (antibunching). Tr[sigma+ sigma- e^(L tau)(sigma- rho sigma+)] / <sigma+ sigma->^2, evolved at 150 bits.
g2_test() ->
    Pb = 150,
    {Sp, H, L} = two_level(),
    B = ?L:basis(Sp, #{}),
    P = #{omega => 2, delta => 0, gamma => 1},
    Hm = ?L:op(B, Sp, H, {P, Pb}), Lm = ?L:op(B, Sp, L, {P, Pb}),
    Pe = ?L:op(B, Sp, ?Q:with_spaces(Sp, fun() -> ?Q:proj(q, 2, 2) end), {#{}, Pb}),
    Lv = ?L:liouvillian(Hm, [Lm], Pb),
    Rho = ?L:steady_state(Lv, Pb),
    {N, _} = ?L:expect_dm(Pe, Rho, Pb),
    % sigma- rho sigma+ = rho_ee |1><1|: the vector with rho_ee at index 0
    Z = {{mp, 0, 0}, {mp, 0, 0}},
    X0 = {element(4, Rho), Z, Z, Z},
    Mu = pq_mp:sqrt(pq_mp:sub(pq_mp:from(4, Pb), pq_mp:from({1, 16}, Pb), Pb), Pb),
    [begin
         {Num, _} = ?L:expect_dm(Pe, X, Pb),
         G2 = pq_mp:divide(Num, pq_mp:mul(N, N, Pb), Pb),
         Tm = pq_mp:from(T, Pb),
         MT = pq_mp:mul(Mu, Tm, Pb),
         Want = pq_mp:sub(pq_mp:from(1, Pb),
                          pq_mp:mul(pq_mp:exp(pq_mp:neg(pq_mp:mul(pq_mp:from({3, 4}, Pb), Tm, Pb)), Pb),
                                    pq_mp:add(pq_mp:cos(MT, Pb),
                                              pq_mp:mul(pq_mp:divide(pq_mp:from({3, 4}, Pb), Mu, Pb), pq_mp:sin(MT, Pb), Pb), Pb), Pb), Pb),
         ?assert(abs(pq_mp:to_float(pq_mp:sub(G2, Want, Pb))) < 1.0e-35)
     end || {T, X} <- ?L:evolve(Lv, X0, [{1, 2}, 1, 2, 5], Pb)],
    % g2(0) = 0 exactly: sigma- sigma- = 0
    {N0, _} = ?L:expect_dm(Pe, X0, Pb),
    ?assertEqual(0.0, pq_mp:to_float(N0)).


%% the driven Kerr cavity, H = Delta ad a + chi ad^2 a^2 + E (a + ad), L = sqrt(kappa) a, at Delta = -2, chi = 1/2,
%% kappa = 1 (semiclassically bistable: Delta < -sqrt(3) kappa/2), the drive scanned across the bistable region:
%% <ad a> and <ad^2 a^2> (so g2(0)) of Drummond and Walls's exact solution (pq_qo_exact:kerr_moment/3) against the
%% steady state of the master equation truncated at 40 photons (pq_lindblad, 200 bits), to 30 digits
kerr_exact_test_() ->
    {timeout, 600, fun() ->
        P = 200, Q = pq_qop, L = pq_lindblad,
        A = Q:a(c), Ad = Q:ad(c),
        H = Q:add([Q:scale(sym(delta), Q:mul(Ad, A)), Q:scale(sym(chi), Q:mul([Ad, Ad, A, A])), Q:scale(sym(e), Q:add(A, Ad))]),
        Lop = Q:scale({pow, sym(kappa), {rat, 1, 2}}, A),
        B = L:basis(#{}, #{c => 40}),
        Nop = L:op(B, #{}, Q:mul(Ad, A), {#{}, P}),
        N2op = L:op(B, #{}, Q:mul([Ad, Ad, A, A]), {#{}, P}),
        Close = fun(X, Y) -> D = pq_mp:sub(X, Y, P), element(2, D) =:= 0 orelse
                             abs(pq_mp:to_float(D)) < 1.0e-30 * max(1.0, abs(pq_mp:to_float(Y))) end,
        [begin
             Par = #{delta => {-2, 1}, chi => {1, 2}, e => E, kappa => {1, 1}},
             Num = maps:map(fun(_, {X, Y}) -> {rat, X, Y} end, Par),
             Rho = L:steady_state(L:liouvillian(L:op(B, #{}, H, {Num, P}), [L:op(B, #{}, Lop, {Num, P})], P), P),
             {N1, _} = L:expect_dm(Nop, Rho, P),
             {N2, _} = L:expect_dm(N2op, Rho, P),
             ?assert(Close(N1, pq_qo_exact:kerr_moment(Par, 1, P))),
             ?assert(Close(N2, pq_qo_exact:kerr_moment(Par, 2, P)))
         end || E <- [{1, 4}, {1, 1}, {3, 2}, {2, 1}]]
    end}.

%% the steady state as rational functions of the parameters (finite-field reconstruction): resonance fluorescence's
%% optical Bloch steady state with Omega, Delta and gamma symbolic, exactly
symbolic_two_level_test_() ->
    {timeout, 120, fun() ->
        {Sp0, H, L} = two_level(),
        Sp = Sp0#{domains => #{omega => real, delta => real, gamma => positive}},
        Os = ?Q:with_spaces(Sp, fun() -> [?Q:proj(q, 2, 2), ?Q:add(?Q:sp(q), ?Q:sm(q)),
                                          ?Q:scale({complex, 0, 1}, ?Q:sub(?Q:sm(q), ?Q:sp(q)))] end),
        [Pe, Sx, Sy] = [E || {_, E} <- ?X:steady_state_symbolic(?L:basis(Sp, #{}), Sp, H, [L], [omega, delta, gamma], Os)],
        Den = {add, [{mul, [4, {pow, sym(delta), 2}]}, {pow, sym(gamma), 2}, {mul, [2, {pow, sym(omega), 2}]}]},
        Is = fun(E, Num) -> {ok, Z} = pq_symbolic:cancel({add, [E, {mul, [-1, Num, {pow, Den, -1}]}]}), Z =:= 0 end,
        ?assert(Is(Pe, {pow, sym(omega), 2})),
        ?assert(Is(Sx, {mul, [4, sym(delta), sym(omega)]})),
        ?assert(Is(Sy, {mul, [2, sym(gamma), sym(omega)]}))
    end}.

%% the driven Kerr cavity truncated at 2 photons, Delta, chi, E and kappa symbolic: <ad a> and <ad^2 a^2> as rational
%% functions, equal to the exact steady state (pq_qo_exact:steady_state/2) at three rational points
symbolic_kerr_test_() ->
    {timeout, 300, fun() ->
        A = ?Q:a(c), Ad = ?Q:ad(c),
        Sp = #{domains => #{delta => real, chi => real, e => real, kappa => positive}},
        H = ?Q:add([?Q:scale(sym(delta), ?Q:mul(Ad, A)), ?Q:scale(sym(chi), ?Q:mul([Ad, Ad, A, A])), ?Q:scale(sym(e), ?Q:add(A, Ad))]),
        L = ?Q:scale({pow, sym(kappa), {rat, 1, 2}}, A),
        B = ?L:basis(Sp, #{c => 2}),
        Os = [?Q:mul(Ad, A), ?Q:mul([Ad, Ad, A, A])],
        Sym = [E || {_, E} <- ?X:steady_state_symbolic(B, Sp, H, [L], [delta, chi, e, kappa], Os)],
        [begin
             Rho = ?X:steady_state(?X:liouvillian(B, Sp, H, [L], P), B),
             Bind = maps:map(fun(_, {X, Y}) -> {rat, X, Y} end, P),
             [begin
                  {ok, V} = pq_symbolic:evaluate(E, Bind, #{exact => true}),
                  {q, Want, {0, 1}} = ?X:expect(op(B, Sp, fun() -> O end, P, trace), Rho, 3),
                  ?assertEqual(Want, case V of {rat, X, Y} -> {X, Y}; N when is_integer(N) -> {N, 1} end)
              end || {O, E} <- lists:zip(Os, Sym)]
         end || P <- [#{delta => {0, 1}, chi => {1, 1}, e => {1, 1}, kappa => {1, 1}},
                      #{delta => {-2, 1}, chi => {1, 3}, e => {1, 2}, kappa => {3, 2}},
                      #{delta => {5, 7}, chi => {2, 1}, e => {3, 1}, kappa => {1, 4}}]]
    end}.
