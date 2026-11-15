%% Equations of motion for averages, closed by a cumulant expansion (quantum optics plan, Q4), as QuantumCumulants.jl
%% does (Plankensteiner, Hotter, Ritsch, Quantum 6, 617 (2022), arXiv:2105.01657).
%%
%% For an operator O, d<O>/dt = <i [H, O] + sum_k (L_k^dagger O L_k - (1/2) {L_k^dagger L_k, O})>, exactly
%% (pq_qop:heisenberg/4): a sum of averages of normal-ordered words. A word is a product of elementary operators (each
%% a^dagger, a, or an atomic |i><j|); its order is their number. At order K, an average of more than K elementary
%% operators is replaced by the moment whose cumulants above order K vanish:
%%   <X_1 ... X_n> = sum over partitions pi of {1..n} with every block of size =< K of prod_(B in pi) kappa(B),
%%   kappa(B) = sum over partitions sigma of B of (|sigma| - 1)! (-1)^(|sigma| - 1) prod_(C in sigma) <C>,
%% each block keeping the operators' order (a subsequence of a normal-ordered word is normal ordered). Starting from
%% the operators asked for, every new average that appears gets its own equation, until the set closes. An average
%% and that of the adjoint word are complex conjugates: one of the two is a variable, the other its conjugate.
%%
%% The averages are symbols of the engine (avg_...), so that the equations are expressions it can canonicalise,
%% differentiate, evaluate and integrate (pq_jet).
-module(pq_cumulant).
-export([equations/5, average_symbol/1, word_of/1, order/2, solve/5]).

-define(Q, pq_qop).

%%% The equations

%% [{Word, Symbol, Rhs}] for the averages of the words of Ops (pq_qop expressions; each of their words is followed)
%% and every average they need, at order K. Spaces as for pq_qop (atoms, and domains => the parameters' domains).
%% Symbol is the variable (an atom); Rhs an engine expression in the variables ({sym, avg_...}, {conj, {sym, ...}}) and
%% the parameters. Self-adjoint words have real variables.
equations(Spaces, H, Ls, Ops, K) ->
    ?Q:with_spaces(Spaces, fun() ->
        Start = lists:usort([rep(W) || Op <- Ops, W := _ <- Op, W =/= []]),
        %% an average of more than K operators has no equation at order K (QuantumCumulants.jl refuses it too)
        case [W || W <- Start, length(elementary(W)) > K] of
            [] -> complete(Spaces, H, Ls, K, Start, #{}, []);
            Over -> error({pq_cumulant, {order_too_low, K, Over}})
        end
    end).

complete(_, _, _, _, [], _, Acc) -> lists:reverse(Acc);
complete(Spaces, H, Ls, K, [W | Rest], Done, Acc) when is_map_key(W, Done) -> complete(Spaces, H, Ls, K, Rest, Done, Acc);
complete(Spaces, H, Ls, K, [W | Rest], Done, Acc) ->
    D = ?Q:heisenberg(Spaces, H, Ls, #{W => 1}),
    {Rhs, Used} = average(D, K, Spaces),
    New = [U || U <- lists:usort(Used), not is_map_key(U, Done), U =/= W],
    complete(Spaces, H, Ls, K, Rest ++ New, Done#{W => true}, [{W, average_symbol(W), Rhs} | Acc]).

%% <E> at order K: {Expr, [representative words used]}
average(E, K, Spaces) ->
    {Terms, Used} = lists:foldl(fun({W, C}, {Ts, Us}) ->
                                        {M, U} = moment(W, K),
                                        {[{mul, [C, M]} | Ts], U ++ Us}
                                end, {[], []}, maps:to_list(E)),
    {can({add, Terms}, Spaces), Used}.

%% the moment of a word at order K, in the variables
moment([], _) -> {1, []};
moment(W, K) ->
    Xs = elementary(W),
    case length(Xs) =< K of
        true -> {var(W), [rep(W)]};
        false ->
            Idx = lists:seq(1, length(Xs)),
            Ps = [P || P <- partitions(Idx), lists:all(fun(B) -> length(B) =< K end, P)],
            Terms = [{mul, [cumulant(B, Xs) || B <- P]} || P <- Ps],
            Used = lists:usort([rep(word(Xs, C)) || P <- Ps, B <- P, S <- partitions(B), C <- S]),
            {{add, Terms}, Used}
    end.

%% kappa(B) = sum_sigma (|sigma| - 1)! (-1)^(|sigma| - 1) prod_C <C>
cumulant(B, Xs) ->
    {add, [{mul, [sign_fact(length(S)) | [var(word(Xs, C)) || C <- S]]} || S <- partitions(B)]}.

sign_fact(N) -> case (N - 1) rem 2 of 0 -> fact(N - 1); 1 -> -fact(N - 1) end.

%% the word of the elementary operators at positions C (in order)
word(Xs, C) ->
    Sub = [lists:nth(I, Xs) || I <- C],
    E = ?Q:mul([#{[X] => 1} || X <- Sub]),
    [W] = maps:keys(E),
    W.

%% a word as its elementary operators, in order: a boson's (a^dagger)^P a^Q as P a^dagger then Q a; an atom's |i><j|
elementary(W) ->
    lists:append([case F of
                      {P, Q} when is_integer(P) ->
                          case kind(N) of
                              boson -> lists:duplicate(P, {N, {1, 0}}) ++ lists:duplicate(Q, {N, {0, 1}});
                              _ -> [{N, F}]
                          end
                  end || {N, F} <- W]).

%% the order of a word (its number of elementary operators)
order(W, Spaces) -> ?Q:with_spaces(Spaces, fun() -> length(elementary(W)) end).

kind(N) ->
    case get({pq_qop, spaces}) of
        #{N := {atom, _}} -> atom;
        _ -> boson
    end.

%%% The variables: a word and its adjoint share one; the representative is the smaller of the two in term order

rep(W) -> min(W, dag(W)).

dag(W) -> [{N, {Y, X}} || {N, {X, Y}} <- W].

var(W) ->
    R = rep(W),
    S = {sym, average_symbol(R)},
    case R =:= W orelse R =:= dag(R) of
        true -> S;
        false -> {conj, S}
    end.

%% avg_<subsystem>_<i>_<j>__... (an atom the engine takes as a symbol)
average_symbol(W) ->
    list_to_atom("avg" ++ lists:append(["__" ++ atom_to_list(N) ++ "_" ++ integer_to_list(X) ++ "_" ++ integer_to_list(Y)
                                       || {N, {X, Y}} <- W])).

%% the word of a variable
word_of(Sym) ->
    "avg" ++ Rest = atom_to_list(Sym),
    [begin [N, X, Y] = string:split(F, "_", all), {list_to_atom(N), {list_to_integer(X), list_to_integer(Y)}} end
     || F <- string:split(string:trim(Rest, leading, "_"), "__", all)].

%% canonical, the averages of words that are not self-adjoint complex (the engine takes undeclared symbols as real),
%% the self-adjoint ones real
can(E, Spaces) ->
    Avgs = [S || S <- symbols(E), lists:prefix("avg", atom_to_list(S))],
    Doms = maps:merge(maps:get(domains, Spaces, #{}),
                      maps:from_list([{S, case word_of(S) of W -> case W =:= dag(W) of true -> real; false -> complex end end}
                                      || S <- Avgs])),
    case pq_symbolic:canonicalise(E, #{domains => Doms}) of
        {ok, C} -> C;
        {error, R} -> error({pq_cumulant, R})
    end.

symbols({sym, S}) -> [S];
symbols(T) when is_tuple(T) -> lists:append([symbols(X) || X <- tuple_to_list(T)]);
symbols(L) when is_list(L) -> lists:append([symbols(X) || X <- L]);
symbols(_) -> [].

%%% Set partitions of a list, each block in the list's order

partitions([]) -> [[]];
partitions([X | Rest]) ->
    lists:append([[[[X] | P] | [insert(X, P, I) || I <- lists:seq(1, length(P))]] || P <- partitions(Rest)]).

insert(X, P, I) ->
    {Before, [B | After]} = lists:split(I - 1, P),
    Before ++ [[X | B] | After].

fact(0) -> 1;
fact(N) -> N * fact(N - 1).

%%% Numerical solution: the variables split into real and imaginary parts, integrated by pq_jet

%% the variables at the times Ts (ascending, > 0), from the initial values Init (#{Symbol => number or {Re, Im}}; missing
%% ones 0), the parameters Params (#{symbol => number}): [{T, #{Symbol => {Re, Im}}}] (floats). Digits: pq_jet's.
solve(Eqs, Params, Init, Ts, Digits) ->
    Syms = [S || {_, S, _} <- Eqs],
    Real = maps:from_list([{S, W =:= dag(W)} || {W, S, _} <- Eqs]),
    Re = fun(S) -> list_to_atom(atom_to_list(S) ++ "_re") end,
    Im = fun(S) -> list_to_atom(atom_to_list(S) ++ "_im") end,
    Sub = maps:from_list([{S, case maps:get(S, Real) of
                                  true -> {sym, Re(S)};
                                  false -> {add, [{sym, Re(S)}, {mul, [{complex, 0, 1}, {sym, Im(S)}]}]}
                              end} || S <- Syms]),
    Doms = maps:from_list([{V, real} || S <- Syms, V <- [Re(S), Im(S)]] ++ [{P, real} || P <- maps:keys(Params)]),
    Can = fun(E) -> {ok, C} = pq_symbolic:canonicalise(E, #{domains => Doms}), C end,
    Pars = maps:map(fun(_, V) -> num(V) end, Params),
    Parts = lists:append([begin
                              E = Can(pq_symbolic:substitute(pq_symbolic:substitute(Rhs, Sub), Pars)),
                              ReE = Can({mul, [{rat, 1, 2}, {add, [E, {conj, E}]}]}),
                              case maps:get(S, Real) of
                                  true -> [{Re(S), ReE}];
                                  false -> [{Re(S), ReE},
                                            {Im(S), Can({mul, [{complex, 0, -1}, {rat, 1, 2}, {add, [E, {mul, [-1, {conj, E}]}]}]})}]
                              end
                          end || {_, S, Rhs} <- Eqs]),
    Vars = [V || {V, _} <- Parts],
    Rhss = [R || {_, R} <- Parts],
    Y0 = [begin
              Base = case atom_to_list(V) of L -> list_to_atom(lists:sublist(L, length(L) - 3)) end,
              {A, B} = case maps:get(Base, Init, 0) of {X, Y} -> {X, Y}; X -> {X, 0} end,
              case lists:suffix("_re", atom_to_list(V)) of true -> A; false -> B end
          end || V <- Vars],
    {Out, _, _} = lists:foldl(fun(T, {Acc, T0, Y}) ->
                                      {ok, Y2} = pq_jet:integrate(Rhss, t, Vars, T0, Y, T, #{digits => Digits}),
                                      Vals = maps:from_list(lists:zip(Vars, [pq_mp:to_float(X) || X <- Y2])),
                                      {[{T, maps:from_list([{S, {maps:get(Re(S), Vals), maps:get(Im(S), Vals, 0.0)}} || S <- Syms])} | Acc],
                                       T, Y2}
                              end, {[], 0, Y0}, Ts),
    lists:reverse(Out).

num({N, D}) -> {rat, N, D};
num(X) -> X.
