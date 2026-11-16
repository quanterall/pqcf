%% pq_dirac (plan "theoretical physics instruments", F2): d-dimensional Dirac traces against FORM's tracen
%% (test/oracle/gt_dirac.sh, test/fixtures/gt_dirac.json) exactly, and against textbook identities.
-module(pq_dirac_tests).
-include_lib("eunit/include/eunit.hrl").

fixture() ->
    Dir = filename:dirname(code:which(?MODULE)),
    Paths = [filename:join([Dir, "..", "test", "fixtures", "gt_dirac.json"]),
             filename:join([Dir, "..", "..", "..", "..", "test", "fixtures", "gt_dirac.json"]), "test/fixtures/gt_dirac.json"],
    [P | _] = [P || P <- Paths, filelib:is_regular(P)],
    {ok, B} = file:read_file(P),
    json:decode(B).

item(<<"i", _/binary>> = I) -> {idx, binary_to_atom(I)};
item(<<"p", _/binary>> = P) -> {vec, binary_to_atom(P)}.

zero(E) -> {ok, Z} = pq_symbolic:cancel(E), Z =:= 0.

%% 60 fully contracted traces of up to 12 gamma matrices, equal to FORM's
form_test_() ->
    {timeout, 120, fun() ->
        [begin
             {ok, Want} = pq_symbolic:parse(F),
             Got = pq_dirac:to_expr(pq_dirac:trace([item(G) || G <- Gs])),
             ?assert(zero({add, [Got, {mul, [-1, Want]}]}))
         end || #{<<"gammas">> := Gs, <<"form">> := F} <- fixture()]
    end}.

%% textbook identities in d dimensions
identities_test() ->
    T = fun(L) -> pq_dirac:to_expr(pq_dirac:trace(L)) end,
    D = {sym, d}, Dot = fun(P, Q) -> {sym, pq_dirac:dot_symbol(P, Q)} end,
    I = fun(M) -> {idx, M} end, V = fun(P) -> {vec, P} end,
    ?assert(zero({add, [T([I(mu), I(mu)]), {mul, [-4, D]}]})),
    ?assert(zero({add, [T([V(p), V(q)]), {mul, [-4, Dot(p, q)]}]})),
    % gamma^mu p gamma_mu = (2 - d) p
    ?assert(zero({add, [T([I(mu), V(p), I(mu), V(q)]), {mul, [-4, {add, [2, {mul, [-1, D]}]}, Dot(p, q)]}]})),
    % Tr(p q r s) = 4 (p.q r.s - p.r q.s + p.s q.r)
    ?assert(zero({add, [T([V(p), V(q), V(r), V(s)]),
                        {mul, [-4, {add, [{mul, [Dot(p, q), Dot(r, s)]}, {mul, [-1, Dot(p, r), Dot(q, s)]}, {mul, [Dot(p, s), Dot(q, r)]}]}]}]})),
    % an odd number of gamma matrices is refused (its trace vanishes; asked for, it is a mistake)
    ?assertError({pq_dirac, odd_number_of_gammas}, pq_dirac:trace([V(p)])).

%% a free index: Tr(gamma^mu p) = 4 p^mu, kept as a factor
free_index_test() ->
    ?assertEqual(#{{0, [], [{p, p, mu}]} => 4}, pq_dirac:trace([{idx, mu}, {vec, p}])).
