%% pq_cumulant (quantum optics plan, Q4): the equations for averages and their cumulant closure against an exact
%% solution, against QuantumCumulants.jl (test/fixtures/qc_cumulants.json, written by comparison/cumulants/run.sh), and
%% against the exact master equation.
-module(pq_cumulant_tests).
-include_lib("eunit/include/eunit.hrl").

-define(C, pq_cumulant).
-define(Q, pq_qop).

fixture() ->
    Dir = filename:dirname(code:which(?MODULE)),
    Paths = [filename:join([Dir, "..", "test", "fixtures", "qc_cumulants.json"]),
             filename:join([Dir, "..", "..", "..", "..", "test", "fixtures", "qc_cumulants.json"]),
             "test/fixtures/qc_cumulants.json"],
    [P | _] = [P || P <- Paths, filelib:is_regular(P)],
    {ok, B} = file:read_file(P),
    json:decode(B).

sym(X) -> {sym, X}.

%% a linear system: the driven, damped cavity. From vacuum <a>(t) = -(2 i E/kappa)(1 - e^(-kappa t/2)) and
%% <ad a> = |<a>|^2; order 2 is exact
driven_cavity_test() ->
    A = ?Q:a(c), Ad = ?Q:ad(c),
    Sp = #{domains => #{e => real, kappa => positive}},
    Eqs = ?C:equations(Sp, ?Q:scale(sym(e), ?Q:add(A, Ad)), [?Q:scale({pow, sym(kappa), {rat, 1, 2}}, A)], [?Q:mul(Ad, A)], 2),
    ?assertEqual(2, length(Eqs)),
    [begin
         Ex = 1 - math:exp(-T / 2),
         {Ar, Ai} = maps:get(avg__c_0_1, V), {N, _} = maps:get(avg__c_1_1, V),
         ?assert(abs(Ar) < 1.0e-15), ?assert(abs(Ai + Ex) < 1.0e-15), ?assert(abs(N - Ex * Ex) < 1.0e-15)
     end || {T, V} <- ?C:solve(Eqs, #{e => {1, 2}, kappa => 1}, #{}, [1.0, 3.0], 20)].

laser() ->
    Sp = #{q => {atom, 2}, domains => #{delta => real, g => real, kappa => positive, gamma => positive, nu => positive}},
    {H, Ls, Ops} = ?Q:with_spaces(Sp, fun() ->
        {?Q:add(?Q:scale(sym(delta), ?Q:mul(?Q:ad(c), ?Q:a(c))),
                ?Q:scale(sym(g), ?Q:add(?Q:mul(?Q:ad(c), ?Q:sm(q)), ?Q:mul(?Q:a(c), ?Q:sp(q))))),
         [?Q:scale({pow, sym(kappa), {rat, 1, 2}}, ?Q:a(c)), ?Q:scale({pow, sym(gamma), {rat, 1, 2}}, ?Q:sm(q)),
          ?Q:scale({pow, sym(nu), {rat, 1, 2}}, ?Q:sp(q))],
         [?Q:mul(?Q:ad(c), ?Q:a(c)), ?Q:proj(q, 2, 2)]}
    end),
    {Sp, H, Ls, Ops}.

%% QuantumCumulants.jl's single-atom laser (its documentation's first example), order 2: the same 8 equations and the
%% same <ad a> and P_e at t = 1, 2, 5, 10 (to its solver's tolerance, 1e-13)
laser_test_() ->
    {timeout, 120, fun() ->
        #{<<"laser">> := #{<<"equations">> := Ne, <<"t">> := Ts, <<"n">> := Ns, <<"pe">> := Pes}} = fixture(),
        {Sp, H, Ls, Ops} = laser(),
        Eqs = ?C:equations(Sp, H, Ls, Ops, 2),
        ?assertEqual(Ne, length(Eqs)),
        Sol = ?C:solve(Eqs, #{delta => {1, 2}, g => {3, 2}, kappa => 1, gamma => {1, 4}, nu => 4}, #{}, Ts, 20),
        [begin
             ?assert(abs(element(1, maps:get(avg__c_1_1, V)) - N) < 1.0e-11),
             ?assert(abs(element(1, maps:get(avg__q_2_2, V)) - Pe) < 1.0e-11)
         end || {{_, V}, N, Pe} <- lists:zip3(Sol, Ns, Pes)]
    end}.

jc() ->
    Sp = #{q => {atom, 2}, domains => #{g => real, e => real, kappa => positive, gamma => positive}},
    {H, Ls} = ?Q:with_spaces(Sp, fun() ->
        {?Q:add(?Q:scale(sym(g), ?Q:add(?Q:mul(?Q:ad(c), ?Q:sm(q)), ?Q:mul(?Q:a(c), ?Q:sp(q)))), ?Q:scale(sym(e), ?Q:add(?Q:a(c), ?Q:ad(c)))),
         [?Q:scale({pow, sym(kappa), {rat, 1, 2}}, ?Q:a(c)), ?Q:scale({pow, sym(gamma), {rat, 1, 2}}, ?Q:sm(q))]}
    end),
    {Sp, H, Ls, #{g => 1, e => {1, 2}, kappa => {1, 2}, gamma => {1, 4}}}.

%% the driven, damped Jaynes-Cummings model at orders 2, 3, 4: as many equations as QuantumCumulants.jl and the same
%% <ad a>(2); and against the exact master equation (pq_lindblad, 25 photons, 120 bits): order 4 within 1e-5, order 2
%% within 1e-4 (order 3 is further than order 2 here, in QuantumCumulants.jl too: the expansion is not monotone)
jc_orders_test_() ->
    {timeout, 300, fun() ->
        #{<<"jc_orders">> := Qc} = fixture(),
        {Sp, H, Ls, Par} = jc(),
        P = 120, L = pq_lindblad,
        B = L:basis(Sp, #{c => 25}),
        Num = maps:map(fun(_, {X, Y}) -> {rat, X, Y}; (_, X) -> X end, Par),
        [{_, Rho}] = L:mesolve(L:op(B, Sp, H, {Num, P}), [L:op(B, Sp, X, {Num, P}) || X <- Ls], L:state(B, #{}, P), [2], P),
        Exact = pq_mp:to_float(element(1, L:expect_dm(L:op(B, Sp, ?Q:mul(?Q:ad(c), ?Q:a(c)), {#{}, P}), Rho, P))),
        Errors = [begin
                      Eqs = ?C:equations(Sp, H, Ls, [?Q:mul(?Q:ad(c), ?Q:a(c))], K),
                      #{<<"equations">> := Ne, <<"n">> := Want} = maps:get(integer_to_binary(K), Qc),
                      ?assertEqual(Ne, length(Eqs)),
                      [{_, V}] = ?C:solve(Eqs, Par, #{}, [2.0], 20),
                      {N, _} = maps:get(avg__c_1_1, V),
                      ?assert(abs(N - Want) < 1.0e-11),
                      abs(N - Exact)
                  end || K <- [2, 3, 4]],
        [E2, _, E4] = Errors,
        ?assert(E4 < 1.0e-5), ?assert(E2 < 1.0e-4)
    end}.

%% an average of more than K operators has no equation at order K
order_too_low_test() ->
    {Sp, H, Ls, _} = jc(),
    ?assertError({pq_cumulant, {order_too_low, 1, _}}, ?C:equations(Sp, H, Ls, [?Q:mul(?Q:ad(c), ?Q:a(c))], 1)).
