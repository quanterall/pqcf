%% pq_symbolic on the tasks of the FORM comparison (gen.py): the same expressions, timed, terms counted.
-module(symbench).
-export([main/0]).

main() ->
    [application:set_env(pq, K, V) || {K, V} <- [{sym_timeout_ms, 900000}, {sym_max_heap_words, 4000000000},
                                               {max_poly_terms, 50000000}, {max_expand_power, 64},
                                               {max_expr_nodes, 500000000}, {max_symbolic_photons, 8}]],
    Big = #{sym_timeout_ms => 900000, sym_max_heap_words => 4000000000, max_poly_terms => 50000000,
            max_expand_power => 64, max_expr_nodes => 500000000},
    Sym = fun(N) -> {sym, list_to_atom(N)} end,
    Xs = [Sym("x" ++ integer_to_list(I)) || I <- lists:seq(1, 6)],
    [run(io_lib:format("t1_~b", [N]), fun() -> pq_symbolic:canonicalise({pow, {add, [1 | Xs]}, N}, Big) end)
     || N <- [8, 12, 16, 20]],
    {ok, P} = pq_symbolic:canonicalise({pow, {add, [1, Sym("x"), Sym("y"), Sym("z"), Sym("w")]}, 10}, Big),
    run("t2", fun() -> pq_symbolic:canonicalise({mul, [P, P]}, Big) end),
    [run(io_lib:format("t3_~b", [N]),
         fun() -> pq_symbolic:permanent({matrix, [[Sym(io_lib:format("a~b~b", [I, J])) || J <- lists:seq(1, N)]
                                                  || I <- lists:seq(1, N)]}) end) || N <- [5, 6, 7, 8]],
    [run(io_lib:format("t4_~b", [K]), fun() -> pq_symbolic:circuit_unitary(mesh(K), Big) end) || K <- [4, 6, 8]],
    halt().

mesh(K) ->
    Bs = lists:sublist(lists:append([lists:seq(L rem 2, K - 2, 2) || L <- lists:seq(0, K - 1)]), K * (K - 1) div 2),
    Is = [{beamsplitter, M, M + 1, {sym, list_to_atom("t" ++ integer_to_list(J))}, {sym, list_to_atom("p" ++ integer_to_list(J))}}
          || {J, M} <- lists:enumerate(0, Bs)],
    {ok, C} = pq_circuit:new(K),
    setelement(5, C, Is).

run(Name, F) ->
    {T, R} = timer:tc(F),
    case R of
        {ok, V} -> io:format("~s ~.2f s ~b terms~n", [Name, T / 1.0e6, terms(V)]);
        {error, E} -> io:format("~s ~.2f s ERROR ~p~n", [Name, T / 1.0e6, E])
    end.

terms({matrix, Rows}) -> lists:sum([terms(E) || Row <- Rows, E <- Row]);
terms({add, L}) -> length(L);
terms(0) -> 0;
terms(_) -> 1.
