%% Dirac algebra in d dimensions (plan "theoretical physics instruments", F2): traces of products of gamma matrices,
%% natively (until now through FORM).
%%
%% A product is a list of {idx, Mu} (gamma^mu) and {vec, P} (p-slash = p_mu gamma^mu). In d dimensions, with Tr 1 = 4,
%%   Tr(g1 g2 ... gn) = sum_(k = 2..n) (-1)^k <g1, gk> Tr(g2 ... gk-1 gk+1 ... gn),
%% <idx mu, idx nu> = g^(mu nu), <idx mu, vec p> = p^mu, <vec p, vec q> = p.q: (n - 1)!! products of pairings for n
%% gamma matrices, n even. Repeated indices are removed first by d-dimensional identities (below), which keeps long
%% traces tractable; an index that appears once stays (a free index).
%%
%% The result is a sum of terms Coef * d^K * prod(p.q) * prod(free factors), returned as a map
%% {K, Dots, Free} => Coef (Dots sorted pairs {P, Q} with P =< Q; Free sorted factors {g, Mu, Nu} or {p, P, Mu}), or as
%% an expression of the engine in the symbols d and dot_P_Q (to_expr/1, for fully contracted traces).
-module(pq_dirac).
-export([trace/1, to_expr/1, dot_symbol/2]).

trace(Gammas) ->
    length(Gammas) rem 2 =:= 0 orelse error({pq_dirac, odd_number_of_gammas}),
    erase_memo(),
    R = tr(Gammas),
    erase_memo(),
    R.

%% Repeated indices first, each by the d-dimensional identity
%%   gamma^mu X S gamma_mu = 2 S X - X (gamma^mu S gamma_mu),   gamma^mu gamma_mu = d,
%% after rotating the trace so that the index's first occurrence leads; then the trace of what is left (only vectors
%% and free indices) by the pairing recursion, memoised on the sequence. The coefficients are polynomials in d.
tr(Gammas) ->
    Key = {?MODULE, tr, Gammas},
    case get(Key) of
        undefined -> V = tr1(Gammas), put(Key, V), V;
        V -> V
    end.

tr1(Gammas) ->
    case repeated(Gammas) of
        none -> pairs_trace(Gammas);
        Mu ->
            Rot = rotate_to(Gammas, {idx, Mu}),
            [{idx, Mu} | Rest] = Rot,
            {S, [{idx, Mu} | B]} = lists:splitwith(fun(X) -> X =/= {idx, Mu} end, Rest),
            lists:foldl(fun({Poly, L}, Acc) -> madd(Acc, scale(tr(L ++ B), Poly)) end, #{}, sandwich(S))
    end.

%% gamma^mu S gamma_mu as [{PolyInD, Product}], PolyInD a map K => integer (the coefficient of d^K)
sandwich([]) -> [{#{1 => 1}, []}];
sandwich([X | S]) -> [{#{0 => 2}, S ++ [X]} | [{maps:map(fun(_, C) -> -C end, P), [X | L]} || {P, L} <- sandwich(S)]].

repeated(Gammas) ->
    Counts = lists:foldl(fun({idx, M}, A) -> A#{M => maps:get(M, A, 0) + 1}; (_, A) -> A end, #{}, Gammas),
    case [M || M := 2 <- Counts] of [] -> none; Ms -> lists:min(Ms) end.

rotate_to(L, X) -> {Before, After} = lists:splitwith(fun(Y) -> Y =/= X end, L), After ++ Before.

%% the trace of a product without repeated indices: Tr(g1 ... gn) = sum_k (-1)^k <g1, gk> Tr(rest), memoised
pairs_trace([]) -> #{{0, [], []} => 4};
pairs_trace([X | Rest]) ->
    lists:foldl(fun({K, Y}, Acc) ->
                        Sub = tr(drop(Rest, K)),
                        madd(Acc, maps:fold(fun({Kd, Dots, Free}, C, M) ->
                                                    {D2, F2} = pair_factor(X, Y, Dots, Free),
                                                    add(M, {Kd, lists:sort(D2), lists:sort(F2)}, C * sign(K))
                                            end, #{}, Sub))
                end, #{}, lists:enumerate(Rest)).

%% <X, Y> joined to the factors already there (no repeated indices remain at this stage, so no contraction)
pair_factor({vec, P}, {vec, Q}, Dots, Free) -> {[{min(P, Q), max(P, Q)} | Dots], Free};
pair_factor({vec, P}, {idx, Mu}, Dots, Free) -> {Dots, [{p, P, Mu} | Free]};
pair_factor({idx, Mu}, {vec, P}, Dots, Free) -> {Dots, [{p, P, Mu} | Free]};
pair_factor({idx, Mu}, {idx, Nu}, Dots, Free) -> {Dots, [{g, min(Mu, Nu), max(Mu, Nu)} | Free]}.

%% sums of terms
madd(A, B) -> maps:fold(fun(K, C, M) -> add(M, K, C) end, A, B).

%% a trace times a polynomial in d
scale(T, Poly) ->
    maps:fold(fun({K, Dots, Free}, C, M) ->
                      maps:fold(fun(Kp, Cp, M2) -> add(M2, {K + Kp, Dots, Free}, C * Cp) end, M, Poly)
              end, #{}, T).

erase_memo() -> _ = [erase(K) || {{?MODULE, tr, _} = K, _} <- get()], ok.

%% the k-th of the rest is the (k+1)-th of the product: (-1)^(k+1)
sign(K) -> case K rem 2 of 1 -> 1; 0 -> -1 end.

drop(L, K) -> {A, [_ | B]} = lists:split(K - 1, L), A ++ B.

add(M, Key, C) ->
    case maps:get(Key, M, 0) + C of
        0 -> maps:remove(Key, M);
        V -> M#{Key => V}
    end.

%% a fully contracted trace as an expression of the engine in d and dot_P_Q
to_expr(T) ->
    Terms = [case Free of
                 [] -> {mul, [C, {pow, {sym, d}, K} | [{sym, dot_symbol(P, Q)} || {P, Q} <- Dots]]};
                 _ -> error({pq_dirac, {free_indices, Free}})
             end || {K, Dots, Free} := C <- T],
    {ok, E} = pq_symbolic:canonicalise({add, Terms}),
    E.

dot_symbol(P, Q) -> list_to_atom("dot_" ++ atom_to_list(P) ++ "_" ++ atom_to_list(Q)).
