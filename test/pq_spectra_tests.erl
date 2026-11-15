%% Spectra to high precision (plan "theoretical physics instruments", QM2): the quartic oscillator H = p^2 + x^4 by
%% pq_ode's shooting on the half-line, even states with psi'(0) = 0 (the Neumann left end) and odd ones with
%% psi(0) = 0, against mpmath's independent shooting (test/oracle/gt_quartic.py, test/fixtures/gt_quartic.json), to 38
%% digits.
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
