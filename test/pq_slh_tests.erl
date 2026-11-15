%% pq_slh (quantum optics plan, Q6): SLH networks reduced symbolically, against the cascaded-system results (Gardiner,
%% PRL 70, 2269 (1993); Carmichael, PRL 70, 2273 (1993)), a feedback loop worked by hand from Gough and James's
%% reduction, PQCF's own symbolic circuit unitary, and a driven cascade's steady state.
-module(pq_slh_tests).
-include_lib("eunit/include/eunit.hrl").

-define(S, pq_slh).
-define(Q, pq_qop).

sym(X) -> {sym, X}.
sqrtk(K) -> {pow, sym(K), {rat, 1, 2}}.
eq(A, B) -> ?Q:is_zero(?Q:sub(A, B)).

cavity(Mode, Kappa, Delta) ->
    ?S:new([[?Q:num(1)]], [?Q:scale(sqrtk(Kappa), ?Q:a(Mode))], ?Q:scale(sym(Delta), ?Q:mul(?Q:ad(Mode), ?Q:a(Mode)))).

%% two cavities in cascade: L = sqrt(ka) a + sqrt(kb) b and H = Da a+a + Db b+b + (sqrt(ka kb)/2i)(b+ a - a+ b)
cascade_test() ->
    Sp = #{domains => #{ka => positive, kb => positive, da => real, db => real}},
    ?Q:with_spaces(Sp, fun() ->
        G = ?S:series(cavity(b, kb, db), cavity(a, ka, da)),
        {H, [L]} = ?S:master(G),
        ?assert(eq(L, ?Q:add(?Q:scale(sqrtk(ka), ?Q:a(a)), ?Q:scale(sqrtk(kb), ?Q:a(b))))),
        K = {mul, [sqrtk(ka), sqrtk(kb), {rat, 1, 2}, {complex, 0, -1}]},
        Want = ?Q:add([?Q:scale(sym(da), ?Q:mul(?Q:ad(a), ?Q:a(a))), ?Q:scale(sym(db), ?Q:mul(?Q:ad(b), ?Q:a(b))),
                       ?Q:scale(K, ?Q:sub(?Q:mul(?Q:ad(b), ?Q:a(a)), ?Q:mul(?Q:ad(a), ?Q:a(b))))]),
        ?assert(eq(H, Want))
    end).

%% a cavity with two mirrors (L = (sqrt(k1) a, sqrt(k2) a)), its output 1 sent through a phase z = e^(i phi) back into
%% input 2: by Gough and James's reduction, L' = (sqrt(k2) + z sqrt(k1)) a, S' = z, H' = (Delta + sqrt(k1 k2) sin phi) a+a,
%% so the decay rate is k1 + k2 + 2 sqrt(k1 k2) cos phi. Phases with exact cos and sin: z = i, -1, (3 + 4i)/5
feedback_test() ->
    Sp = #{domains => #{k1 => positive, k2 => positive, delta => real}},
    ?Q:with_spaces(Sp, fun() ->
        N = ?Q:mul(?Q:ad(c), ?Q:a(c)),
        G = ?S:new([[?Q:num(1), ?Q:num(0)], [?Q:num(0), ?Q:num(1)]],
                   [?Q:scale(sqrtk(k1), ?Q:a(c)), ?Q:scale(sqrtk(k2), ?Q:a(c))], ?Q:scale(sym(delta), N)),
        [begin
             Phase = ?S:static([[Z, 0], [0, 1]]),
             R = ?S:feedback(?S:series(Phase, G), 1, 2),
             {H, [L]} = ?S:master(R),
             ?assert(eq(L, ?Q:scale({add, [sqrtk(k2), {mul, [Z, sqrtk(k1)]}]}, ?Q:a(c)))),
             ?assert(eq(H, ?Q:scale({add, [sym(delta), {mul, [sqrtk(k1), sqrtk(k2), Sin]}]}, N))),
             ?assert(eq(?Q:mul(?Q:dag(L), L),
                        ?Q:scale({add, [sym(k1), sym(k2), {mul, [2, Cos, sqrtk(k1), sqrtk(k2)]}]}, N))),
             #{s := [[Sr]]} = R,
             ?assert(eq(Sr, ?Q:num(Z)))
         end || {Z, Cos, Sin} <- [{{complex, 0, 1}, 0, 1}, {-1, -1, 0}, {{complex, {rat, 3, 5}, {rat, 4, 5}}, {rat, 3, 5}, {rat, 4, 5}}]]
    end).

%% a three-mode mesh of beamsplitters and a phase: the SLH network of static elements has, as its S, exactly PQCF's
%% symbolic circuit unitary (pq_symbolic:circuit_unitary/1), entry by entry
circuit_test() ->
    I = {complex, 0, 1},
    Bs = fun(A, B, T, P) ->
                 G = [[{cos, sym(T)}, {mul, [-1, {exp, {mul, [-1, I, sym(P)]}}, {sin, sym(T)}]}],
                      [{mul, [{exp, {mul, [I, sym(P)]}}, {sin, sym(T)}]}, {cos, sym(T)}]],
                 embed(3, [A, B], G)
         end,
    Ph = fun(M, P) -> embed(3, [M], [[{exp, {mul, [I, sym(P)]}}]]) end,
    Steps = [{{beamsplitter, 0, 1, sym(t1), sym(p1)}, Bs(0, 1, t1, p1)},
             {{phaseshift, 1, sym(f)}, Ph(1, f)},
             {{beamsplitter, 1, 2, sym(t2), sym(p2)}, Bs(1, 2, t2, p2)},
             {{beamsplitter, 0, 1, sym(t3), sym(p3)}, Bs(0, 1, t3, p3)}],
    {ok, C0} = pq_circuit:new(3),
    C = lists:foldl(fun({Ins, _}, A) -> pq_circuit:add(A, Ins) end, C0, [S || S <- Steps]),
    {ok, {matrix, U}} = pq_symbolic:circuit_unitary(C),
    Net = lists:foldl(fun({_, M}, Acc) -> ?S:series(?S:static(M), Acc) end, ?S:identity(3), Steps),
    #{s := S} = Net,
    [begin
         {ok, D} = pq_symbolic:cancel({add, [Want, {mul, [-1, case maps:to_list(Got) of [] -> 0; [{[], X}] -> X end]}]}),
         ?assertEqual(0, D)
     end || {RowS, RowU} <- lists:zip(S, U), {Got, Want} <- lists:zip(RowS, RowU)].

%% the scalar matrix G on the modes Ms (0-based) of N, the identity elsewhere
embed(N, Ms, G) ->
    Pos = maps:from_list(lists:zip(Ms, lists:seq(1, length(Ms)))),
    [[case {maps:find(I, Pos), maps:find(J, Pos)} of
          {{ok, A}, {ok, B}} -> lists:nth(B, lists:nth(A, G));
          _ -> case I =:= J of true -> 1; false -> 0 end
      end || J <- lists:seq(0, N - 1)] || I <- lists:seq(0, N - 1)].

%% a coherent field alpha driving cavity a, whose output drives cavity b (both on resonance): the steady state is
%% <a> = -2 alpha/sqrt(ka) and, the field leaving a being -alpha, <b> = 2 alpha/sqrt(kb). The reduced network's
%% equations (pq_cumulant at order 1, exact for a linear system) relax to them
driven_cascade_test_() ->
    {timeout, 120, fun() ->
        Sp = #{domains => #{ka => positive, kb => positive, da => real, db => real}},
        {H, Ls} = ?Q:with_spaces(Sp, fun() ->
            ?S:master(?S:series(cavity(b, kb, db), ?S:series(cavity(a, ka, da), ?S:coherent({rat, 1, 2}))))
        end),
        Eqs = pq_cumulant:equations(Sp, H, Ls, [?Q:a(a), ?Q:a(b)], 1),
        [{_, V}] = pq_cumulant:solve(Eqs, #{ka => 1, kb => 4, da => 0, db => 0}, #{}, [80.0], 20),
        {Ar, Ai} = maps:get(avg__a_0_1, V), {Br, Bi} = maps:get(avg__b_0_1, V),
        ?assert(abs(Ar + 1.0) < 1.0e-12), ?assert(abs(Ai) < 1.0e-12),
        ?assert(abs(Br - 0.5) < 1.0e-12), ?assert(abs(Bi) < 1.0e-12)
    end}.
