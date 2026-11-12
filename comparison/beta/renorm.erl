%% The QCD beta function from the gluon and ghost self-energies in Landau gauge (comparison/beta/README.md): the
%% bare self-energies Gamma/Q^2 = 1 + sum_L a_B^L P_L (MS-bar, a = g^2/(16 pi^2)) made finite by Z_3 and Z~_3,
%% order by order, with the bare coupling a_B = Z_a a and Z_a = 1/(Z~_3^2 Z_3) -- the ghost-gluon vertex is finite
%% in Landau gauge (Taylor), and Landau gauge needs no gauge-parameter renormalisation. Then
%% beta_(L-1) = -L [a^L eps^-1] Z_a, the higher poles checked against the lower coefficients.
-module(renorm).
-export([run/2, beta/1, write/3]).

-define(L, 3).

%% Files: fun(Field, Loops) -> FORM output text of P (the gluon's times D - 1)
run(Dir, MaxL) ->
    {ok, _} = application:ensure_all_started(pq),
    Read = fun(F, L) ->
                   {ok, Bin} = file:read_file(filename:join(Dir, io_lib:format("res_~s_~b.log", [F, L]))),
                   {match, [E]} = re:run(Bin, "\\n   F =\\n(.*?);", [dotall, {capture, all_but_first, list}]),
                   {_, S} = epsexp:series_of_form(E, #{}),
                   S(8)
           end,
    %% the gluon's projection carried D - 1 = 3 - 2 eps
    Dm1 = epsexp:sinv(epsexp:exact(#{0 => #{[] => {3, 1}}, 1 => #{[] => {-2, 1}}}), 10),
    G = maps:from_list([{L, epsexp:smul(Read("glu", L), Dm1)} || L <- lists:seq(1, MaxL)]),
    H = maps:from_list([{L, Read("gho", L)} || L <- lists:seq(1, MaxL)]),
    beta(#{gluon => G, ghost => H, loops => MaxL}).

beta(#{gluon := G, ghost := H, loops := N}) ->
    One = #{0 => one()},
    {Z3, Zt} = lists:foldl(fun(L, {Z3a, Zta}) ->
                                   Za = za(Z3a, Zta, N),
                                   Ab = amul(#{1 => one()}, Za, N),
                                   Fg = amul(Z3a, gamma(G, Ab, N), N),
                                   Fh = amul(Zta, gamma(H, Ab, N), N),
                                   {Z3a#{L => neg(pole(maps:get(L, Fg)))}, Zta#{L => neg(pole(maps:get(L, Fh)))}}
                           end, {One, One}, lists:seq(1, N)),
    Za = za(Z3, Zt, N),
    %% finiteness, now with every order's Z
    Ab = amul(#{1 => one()}, Za, N),
    Finite = [epsexp:is_zero(pole(maps:get(L, amul(Z, gamma(P, Ab, N), N)))) || {Z, P} <- [{Z3, G}, {Zt, H}], L <- lists:seq(1, N)],
    Betas = [{L - 1, scale(coef(maps:get(L, Za), -1), {-L, 1})} || L <- lists:seq(1, N)],
    #{z3 => Z3, zghost => Zt, za => Za, finite => lists:all(fun(X) -> X end, Finite), betas => Betas, loops => N}.

%% Gamma(a_B) = 1 + sum_L P_L a_B^L
gamma(P, Ab, N) ->
    lists:foldl(fun(L, Acc) -> aadd(Acc, ascale(apow(Ab, L, N), maps:get(L, P), N)) end, #{0 => one()}, lists:seq(1, N)).

za(Z3, Zt, N) -> ainv(amul(amul(Zt, Zt, N), Z3, N), N).

%%% power series in a with epsilon-series coefficients (#{Order => Series}), truncated at order N

one() -> epsexp:exact(#{0 => #{[] => {1, 1}}}).
zero() -> epsexp:exact(#{}).
neg(S) -> epsexp:scale(S, {-1, 1}).
scale(P, Q) -> maps:map(fun(_, V) -> qmul(V, Q) end, P).
pole(S) -> epsexp:exact(maps:filter(fun(K, _) -> K < 0 end, epsexp:coeffs(S))).
coef(S, K) -> maps:get(K, epsexp:coeffs(S), #{}).

aadd(A, B) -> maps:fold(fun(K, V, Acc) -> Acc#{K => epsexp:sadd(maps:get(K, Acc, zero()), V)} end, A, B).
amul(A, B, N) ->
    lists:foldl(fun({K, V}, Acc) -> Acc#{K => epsexp:sadd(maps:get(K, Acc, zero()), V)} end, #{},
                [{Ka + Kb, epsexp:smul(Va, Vb)} || Ka := Va <- A, Kb := Vb <- B, Ka + Kb =< N]).
ascale(A, S, N) -> amul(A, #{0 => S}, N).
apow(_, 0, _) -> #{0 => one()};
apow(A, K, N) -> amul(A, apow(A, K - 1, N), N).
%% 1/A for A = 1 + O(a)
ainv(A, N) ->
    lists:foldl(fun(K, B) ->
                        S = lists:foldl(fun(J, Acc) -> epsexp:sadd(Acc, epsexp:smul(maps:get(J, A, zero()), maps:get(K - J, B))) end,
                                        zero(), lists:seq(1, K)),
                        B#{K => neg(S)}
                end, #{0 => one()}, lists:seq(1, N)).

qmul({A, B}, {C, D}) -> qn({A * C, B * D}).
qn({N, D}) when D < 0 -> qn({-N, -D});
qn({0, _}) -> {0, 1};
qn({N, D}) -> G = gcd(abs(N), D), {N div G, D div G}.
gcd(A, 0) -> A;
gcd(A, B) -> gcd(B, A rem B).

%% the renormalisation constants and the beta function as text, into File: Z_3, Z~_3 and Z_a = Z_g^2 to a^N,
%% each a^L coefficient's poles in epsilon, and beta_0 .. beta_(N-1)
write(File, R = #{loops := N}, Header) ->
    Z = fun(Name, Key) -> [io_lib:format("~s, coefficient of a^~b:~n~s", [Name, L, epsexp:to_string(pole(maps:get(L, maps:get(Key, R))))])
                           || L <- lists:seq(1, N)] end,
    Betas = [io_lib:format("beta_~b = ~s~n", [K, poly(B)]) || {K, B} <- maps:get(betas, R)],
    ok = file:write_file(File, [Header, "\n# beta(a) = -sum_k beta_k a^(k+2), a = g^2/(16 pi^2); cA = C_A, cR = C_F, I2R = T_F, nf = n_f\n",
                                Betas, "\n", Z("Z_3 (gluon)", z3), "\n", Z("Z~_3 (ghost)", zghost), "\n", Z("Z_a = Z_g^2 = 1/(Z~_3^2 Z_3)", za)]).

poly(P) -> lists:join(" + ", [[case Q of {Nu, 1} -> integer_to_list(Nu); {Nu, D} -> io_lib:format("~b/~b", [Nu, D]) end,
                               [["*", atom_to_list(S), case E of 1 -> ""; _ -> "^" ++ integer_to_list(E) end] || {S, E} <- M]]
                              || {M, Q} <- lists:sort(maps:to_list(P))]).
