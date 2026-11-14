%% PQCF's side of the quantum optics head-to-head (README.md): every task of tasks.json through pq_qop and pq_lindblad
%% at 128 bits, timed; writes pqcf.json, the values as decimal strings of 36 digits.
-module(qo_pqcf).
-export([main/2]).

-define(P, 128).
-define(BUDGET, 120.0).
-define(L, pq_lindblad).
-define(Q, pq_qop).

main(TasksFile, Dir) ->
    {ok, Bin} = file:read_file(TasksFile),
    T = json:decode(Bin),
    Tasks = #{<<"rabi">> => task(fun() -> rabi(maps:get(<<"rabi">>, T)) end),
              <<"revival">> => task(fun() -> revival(maps:get(<<"revival">>, T)) end),
              <<"cavity_decay">> => task(fun() -> cavity_decay(maps:get(<<"cavity_decay">>, T)) end),
              <<"fluorescence">> => task(fun() -> fluorescence(maps:get(<<"fluorescence">>, T)) end),
              <<"damped_jc">> => damped_jc(maps:get(<<"damped_jc">>, T))},
    Out = #{<<"tool">> => <<"PQCF (128 bits)">>, <<"version">> => version(), <<"tasks">> => Tasks},
    ok = file:write_file(filename:join(Dir, "pqcf.json"), json:encode(Out)).

version() ->
    {ok, V} = application:get_key(pq, vsn),
    list_to_binary(V).

task(F) ->
    {S, V} = timed(F),
    #{<<"values">> => [str(X) || X <- V], <<"seconds">> => S}.

str({A, _}) -> list_to_binary(pq_mp:to_string(A, 36)).

%% the inputs: integers and floats (binary fractions, exact)
num(X) -> X.

%%% The tasks, as in qutip_side.py

jc(N, G) ->
    Sp = #{q => {atom, 2}},
    B = ?L:basis(Sp, #{c => N}),
    H = ?Q:with_spaces(Sp, fun() -> ?Q:scale({sym, g}, ?Q:add(?Q:mul(?Q:ad(c), ?Q:sm(q)), ?Q:mul(?Q:a(c), ?Q:sp(q)))) end),
    Pe = ?L:op(B, Sp, ?Q:with_spaces(Sp, fun() -> ?Q:proj(q, 2, 2) end), {#{}, ?P}),
    {Sp, B, H, ?L:op(B, Sp, H, {#{g => num(G)}, ?P}), Pe}.

rabi(#{<<"n">> := N, <<"g">> := G, <<"times">> := Ts}) ->
    {_, B, _, Hm, Pe} = jc(N, G),
    Psi0 = ?L:state(B, #{q => ?L:fock(2)}, ?P),
    [?L:expect(Pe, Psi, ?P) || {_, Psi} <- ?L:sesolve(Hm, Psi0, Ts, ?P)].

revival(#{<<"n">> := N, <<"g">> := G, <<"alpha">> := A, <<"times">> := Ts}) ->
    {_, B, _, Hm, Pe} = jc(N, G),
    Psi0 = ?L:state(B, #{q => ?L:fock(2), c => ?L:coherent(N - 1, A, ?P)}, ?P),
    [?L:expect(Pe, Psi, ?P) || {_, Psi} <- ?L:sesolve(Hm, Psi0, Ts, ?P)].

cavity_decay(#{<<"n">> := N, <<"alpha">> := A, <<"kappa">> := K, <<"times">> := Ts}) ->
    Sp = #{},
    B = ?L:basis(Sp, #{c => N}),
    H = ?L:op(B, Sp, ?Q:num(0), {#{}, ?P}),
    L = ?L:op(B, Sp, ?Q:scale({pow, {sym, kappa}, {rat, 1, 2}}, ?Q:a(c)), {#{kappa => num(K)}, ?P}),
    Nm = ?L:op(B, Sp, ?Q:mul(?Q:ad(c), ?Q:a(c)), {#{}, ?P}),
    Psi0 = ?L:state(B, #{c => ?L:coherent(N, A, ?P)}, ?P),
    [?L:expect_dm(Nm, Rho, ?P) || {_, Rho} <- ?L:mesolve(H, [L], Psi0, Ts, ?P)].

fluorescence(#{<<"gamma">> := Gm, <<"cases">> := Cases}) ->
    Sp = #{q => {atom, 2}},
    B = ?L:basis(Sp, #{}),
    H = ?Q:with_spaces(Sp, fun() ->
            ?Q:add(?Q:scale({mul, [{rat, 1, 2}, {sym, omega}]}, ?Q:add(?Q:sp(q), ?Q:sm(q))),
                   ?Q:scale({mul, [-1, {sym, delta}]}, ?Q:proj(q, 2, 2)))
        end),
    L = ?L:op(B, Sp, ?Q:with_spaces(Sp, fun() -> ?Q:scale({pow, {sym, gamma}, {rat, 1, 2}}, ?Q:sm(q)) end),
              {#{gamma => num(Gm)}, ?P}),
    Pe = ?L:op(B, Sp, ?Q:with_spaces(Sp, fun() -> ?Q:proj(q, 2, 2) end), {#{}, ?P}),
    [begin
         Hm = ?L:op(B, Sp, H, {#{omega => num(Om), delta => num(Dl)}, ?P}),
         ?L:expect_dm(Pe, ?L:steady_state(?L:liouvillian(Hm, [L], ?P), ?P), ?P)
     end || [Om, Dl] <- Cases].

%% the driven, damped Jaynes-Cummings model, <ad a> and P_e at t, sizes in order until one exceeds the budget
damped_jc(#{<<"sizes">> := Sizes} = T) -> damped_jc(Sizes, T, []).

damped_jc([], _, Acc) -> lists:reverse(Acc);
damped_jc([N | Rest], #{<<"g">> := G, <<"kappa">> := K, <<"gamma">> := Gm, <<"drive">> := E, <<"t">> := Tf} = T, Acc) ->
    F = fun() ->
                {Sp, B, H, _, Pe} = jc(N, G),
                Hd = ?Q:add(H, ?Q:scale({sym, e}, ?Q:add(?Q:a(c), ?Q:ad(c)))),
                Hm = ?L:op(B, Sp, Hd, {#{g => num(G), e => num(E)}, ?P}),
                La = ?L:op(B, Sp, ?Q:scale({pow, {sym, kappa}, {rat, 1, 2}}, ?Q:a(c)), {#{kappa => num(K)}, ?P}),
                Ls = ?L:op(B, Sp, ?Q:with_spaces(Sp, fun() -> ?Q:scale({pow, {sym, gamma}, {rat, 1, 2}}, ?Q:sm(q)) end),
                           {#{gamma => num(Gm)}, ?P}),
                Nm = ?L:op(B, Sp, ?Q:mul(?Q:ad(c), ?Q:a(c)), {#{}, ?P}),
                Psi0 = ?L:state(B, #{}, ?P),
                [{_, Rho}] = ?L:mesolve(Hm, [La, Ls], Psi0, [Tf], ?P),
                [?L:expect_dm(Nm, Rho, ?P), ?L:expect_dm(Pe, Rho, ?P)]
        end,
    {S, V} = once(F),
    Row = #{<<"size">> => N, <<"values">> => [str(X) || X <- V], <<"seconds">> => S},
    case S > ?BUDGET of
        true -> lists:reverse([Row | Acc]);
        false -> damped_jc(Rest, T, [Row | Acc])
    end.

%% the median of repeated runs: once if it takes a second or more, else up to 7 times within about two seconds
timed(F) ->
    {T1, V} = once(F),
    Ts = repeat(F, T1, [T1]),
    {median(Ts), V}.
repeat(_, _, Ts) when length(Ts) >= 7 -> Ts;
repeat(F, T1, Ts) when T1 * (length(Ts) + 1) < 2.0, T1 < 1.0 -> {T, _} = once(F), repeat(F, T1, [T | Ts]);
repeat(_, _, Ts) -> Ts.
once(F) ->
    T0 = erlang:monotonic_time(),
    V = F(),
    {erlang:convert_time_unit(erlang:monotonic_time() - T0, native, microsecond) / 1.0e6, V}.
median(Ts) -> lists:nth((length(Ts) + 1) div 2, lists:sort(Ts)).
