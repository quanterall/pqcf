%% pq_qop (quantum optics plan, Q1): the operator algebra against SymPy's normal ordering (test/oracle/gt_qop.py,
%% test/fixtures/gt_qop.json), textbook identities, and the Heisenberg-Lindblad equations of the Jaynes-Cummings model.
-module(pq_qop_tests).
-include_lib("eunit/include/eunit.hrl").

-define(Q, pq_qop).

eq(A, B) -> ?Q:is_zero(?Q:sub(A, B)).

fixture() ->
    Dir = filename:dirname(code:which(?MODULE)),
    Paths = [filename:join([Dir, "..", "test", "fixtures", "gt_qop.json"]),
             filename:join([Dir, "..", "..", "..", "..", "test", "fixtures", "gt_qop.json"]),
             "test/fixtures/gt_qop.json"],
    [P | _] = [P || P <- Paths, filelib:is_regular(P)],
    {ok, B} = file:read_file(P),
    json:decode(B).

%% every word's normal-ordered form equals SymPy's, term for term
sympy_test() ->
    #{<<"cases">> := Cases} = fixture(),
    Letter = fun(<<"a">>) -> ?Q:a(a); (<<"ad">>) -> ?Q:ad(a); (<<"b">>) -> ?Q:a(b); (<<"bd">>) -> ?Q:ad(b) end,
    Mono = fun(Pa, Qa, Pb, Qb) -> ?Q:mul([?Q:pow(?Q:ad(a), Pa), ?Q:pow(?Q:a(a), Qa), ?Q:pow(?Q:ad(b), Pb), ?Q:pow(?Q:a(b), Qb)]) end,
    [begin
         Got = ?Q:mul([Letter(L) || L <- Word]),
         Want = ?Q:add([?Q:scale(C, Mono(Pa, Qa, Pb, Qb)) || [C, Pa, Qa, Pb, Qb] <- Terms]),
         ?assert(eq(Got, Want), Word)
     end || #{<<"word">> := Word, <<"terms">> := Terms} <- Cases].

boson_identities_test() ->
    A = ?Q:a(c), Ad = ?Q:ad(c),
    % a a^dagger = a^dagger a + 1
    ?assert(eq(?Q:mul(A, Ad), ?Q:add(?Q:mul(Ad, A), ?Q:one()))),
    % (a + a^dagger)^4 = ad^4 + 4 ad^3 a + 6 ad^2 a^2 + 4 ad a^3 + a^4 + 6 ad^2 + 12 ad a + 6 a^2 + 3
    X = ?Q:add(A, Ad),
    M = fun(P, Q) -> ?Q:mul(?Q:pow(Ad, P), ?Q:pow(A, Q)) end,
    Want = ?Q:add([M(4, 0), ?Q:scale(4, M(3, 1)), ?Q:scale(6, M(2, 2)), ?Q:scale(4, M(1, 3)), M(0, 4),
                   ?Q:scale(6, M(2, 0)), ?Q:scale(12, M(1, 1)), ?Q:scale(6, M(0, 2)), ?Q:num(3)]),
    ?assert(eq(?Q:pow(X, 4), Want)),
    % [a, ad^n] = n ad^(n-1)
    [?assert(eq(?Q:comm(A, ?Q:pow(Ad, N)), ?Q:scale(N, ?Q:pow(Ad, N - 1)))) || N <- lists:seq(1, 6)],
    % e^(-alpha ad) a e^(alpha ad) = a + alpha: the series sum_n (1/n!) [X, [X, ... a]] with X = -alpha ad
    Xa = ?Q:scale({mul, [-1, {sym, alpha}]}, Ad),
    Nested = lists:foldl(fun(N, {Acc, T}) -> T2 = ?Q:comm(Xa, T), {?Q:add(Acc, ?Q:scale({rat, 1, fact(N)}, T2)), T2} end,
                         {A, A}, lists:seq(1, 4)),
    ?assert(eq(element(1, Nested), ?Q:add(A, ?Q:num({sym, alpha})))),
    % different modes commute; the adjoint of a ad a is ad a ad = ad^2 a + ad
    ?assert(eq(?Q:mul(?Q:a(c), ?Q:ad(d)), ?Q:mul(?Q:ad(d), ?Q:a(c)))),
    ?assert(eq(?Q:dag(?Q:mul([A, Ad, A])), ?Q:add(M(2, 1), Ad))).

fact(0) -> 1;
fact(N) -> N * fact(N - 1).

atom_identities_test() ->
    ?Q:with_spaces(#{q => {atom, 2}, r => {atom, 3}}, fun() ->
        Sp = ?Q:sp(q), Sm = ?Q:sm(q), Sz = ?Q:sz(q),
        % [sigma+, sigma-] = sigma_z; sigma_z^2 = 1; sigma+ sigma- = (1 + sigma_z)/2; sigma-^2 = 0
        ?assert(eq(?Q:comm(Sp, Sm), Sz)),
        ?assert(eq(?Q:mul(Sz, Sz), ?Q:one())),
        ?assert(eq(?Q:mul(Sp, Sm), ?Q:scale({rat, 1, 2}, ?Q:add(?Q:one(), Sz)))),
        ?assert(?Q:is_zero(?Q:mul(Sm, Sm))),
        % three levels: |1><2| |2><1| = |1><1| = 1 - |2><2| - |3><3|
        ?assert(eq(?Q:mul(?Q:proj(r, 1, 2), ?Q:proj(r, 2, 1)),
                   ?Q:sub(?Q:one(), ?Q:add(?Q:proj(r, 2, 2), ?Q:proj(r, 3, 3)))))
    end).

%% the Jaynes-Cummings model, H = w ad a + (w0/2) sigma_z + g (ad sigma- + a sigma+), a damped cavity L = sqrt(kappa) a:
%%   d a/dt = -i w a - i g sigma- - (kappa/2) a
%%   d sigma-/dt = -i w0 sigma- + i g a sigma_z
jaynes_cummings_test() ->
    Sp = #{q => {atom, 2}, domains => #{kappa => positive, w => real, w0 => real, g => real}},
    S = fun(X) -> {sym, X} end,
    I = {complex, 0, 1},
    H = fun() -> ?Q:add([?Q:scale(S(w), ?Q:mul(?Q:ad(c), ?Q:a(c))), ?Q:scale({mul, [{rat, 1, 2}, S(w0)]}, ?Q:sz(q)),
                         ?Q:scale(S(g), ?Q:add(?Q:mul(?Q:ad(c), ?Q:sm(q)), ?Q:mul(?Q:a(c), ?Q:sp(q))))]) end,
    L = ?Q:scale({pow, S(kappa), {rat, 1, 2}}, ?Q:a(c)),
    Da = ?Q:heisenberg(Sp, H(), [L], ?Q:a(c)),
    ?assert(eq(Da, ?Q:add([?Q:scale({mul, [-1, I, S(w)]}, ?Q:a(c)), ?Q:scale({mul, [-1, I, S(g)]}, ?Q:sm(q)),
                           ?Q:scale({mul, [{rat, -1, 2}, S(kappa)]}, ?Q:a(c))]))),
    Ds = ?Q:heisenberg(Sp, H(), [], ?Q:sm(q)),
    Want = ?Q:with_spaces(Sp, fun() ->
                                      ?Q:add(?Q:scale({mul, [-1, I, S(w0)]}, ?Q:sm(q)), ?Q:scale({mul, [I, S(g)]}, ?Q:mul(?Q:a(c), ?Q:sz(q))))
                              end),
    ?assert(eq(Ds, Want)).
