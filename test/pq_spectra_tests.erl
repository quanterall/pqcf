%% Spectra to high precision (plan "theoretical physics instruments", QM2): the quartic oscillator H = p^2 + x^4 by
%% pq_ode's shooting on the half-line, even states with psi'(0) = 0 (the Neumann left end) and odd ones with
%% psi(0) = 0, against mpmath's independent shooting (test/oracle/gt_quartic.py, test/fixtures/gt_quartic.json), to 38
%% digits; and hydrogen's radial equation, the s-states included (a simple pole at r = 0), against -1/n^2.
-module(pq_spectra_tests).
-include_lib("eunit/include/eunit.hrl").

fixture() ->
    Dir = filename:dirname(code:which(?MODULE)),
    Paths = [filename:join([Dir, "..", "test", "fixtures", "gt_quartic.json"]),
             filename:join([Dir, "..", "..", "..", "..", "test", "fixtures", "gt_quartic.json"]), "test/fixtures/gt_quartic.json"],
    [P | _] = [P || P <- Paths, filelib:is_regular(P)],
    {ok, B} = file:read_file(P),
    json:decode(B).

quartic_test_() ->
    {timeout, 300, fun() ->
        #{<<"B">> := B, <<"states">> := States} = fixture(),
        Want = [binary_to_list(E) || #{<<"E">> := E} <- States],
        Opts = fun(L) -> #{left => L, right => {dirichlet, B}, range => {0.5, 13}, grid => 60, count => 2, digits => 40} end,
        {ok, [E0, E2]} = pq_ode:eigenvalues({pow, {sym, x}, 4}, x, Opts({neumann, 0})),
        {ok, [E1, E3]} = pq_ode:eigenvalues({pow, {sym, x}, 4}, x, Opts({dirichlet, 0})),
        [begin
             D = pq_mp:sub(pq_mp:from_string(Got, 160), pq_mp:from_string(W, 160), 160),
             ?assert(element(2, D) =:= 0 orelse abs(pq_mp:to_float(D)) < 1.0e-37)
         end || {Got, W} <- lists:zip([E0, E1, E2, E3], Want)]
    end}.

%% hydrogen: -u'' + (l(l+1)/r^2 - 2/r) u = E u has E = -1/n^2 (n = l + 1, l + 2, ...). The wall at r = 200 costs
%% e^(-2 r/n): below 1e-50 for n =< 3. l = 0 has only a simple pole at r = 0 (a regular singular point all the same)
hydrogen_test_() ->
    {timeout, 600, fun() ->
        R = {sym, r},
        [begin
             V = {add, [{mul, [L * (L + 1), {pow, R, -2}]}, {mul, [-2, {pow, R, -1}]}]},
             {ok, Es} = pq_ode:eigenvalues(V, r, #{left => {frobenius, 0}, right => {dirichlet, 200}, range => {-1.2, -0.09},
                                                   grid => 240, count => 3 - L, digits => 30}),
             [begin
                  D = pq_mp:add(pq_mp:from_string(E, 140), pq_mp:from({1, N * N}, 140), 140),
                  ?assert(element(2, D) =:= 0 orelse abs(pq_mp:to_float(D)) < 1.0e-29)
              end || {E, N} <- lists:zip(Es, lists:seq(L + 1, 3))]
         end || L <- [0, 1, 2]]
    end}.
