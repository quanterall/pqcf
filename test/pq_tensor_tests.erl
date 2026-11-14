%% pq_tensor (decision record SY2, phase 4): canonical forms of Riemann monomials against SymPy 1.14's Butler-Portugal
%% implementation (test/oracle/gt_tensor.py, test/fixtures/gt_tensor.json), and the numbers of independent quadratic
%% and cubic Riemann invariants against the literature (S. A. Fulling, R. C. King, B. G. Wybourne, C. J. Cummins,
%% Class. Quantum Grav. 9 (1992) 1151, https://doi.org/10.1088/0264-9381/9/5/003: 3 and 8).
-module(pq_tensor_tests).
-include_lib("eunit/include/eunit.hrl").

fixture() ->
    Dir = filename:dirname(code:which(?MODULE)),
    Paths = [filename:join([Dir, "..", "test", "fixtures", "gt_tensor.json"]),
             filename:join([Dir, "..", "..", "..", "..", "test", "fixtures", "gt_tensor.json"]),
             "test/fixtures/gt_tensor.json"],
    [P | _] = [P || P <- Paths, filelib:is_regular(P)],
    {ok, B} = file:read_file(P),
    json:decode(B).

%% the same classes as SymPy: equal SymPy forms <=> equal pq_tensor forms, the same zeros, consistent relative signs
sympy_test() ->
    #{<<"cases">> := Cases} = fixture(),
    Rows = [begin
                Mono = [{riemann, [binary_to_atom(I) || I <- F]} || F <- Fs],
                Sy = binary_to_list(S),
                {Sign, Class} = case Sy of [$- | Rest] -> {-1, Rest}; _ -> {1, Sy} end,
                {Class, Sign, pq_tensor:canonical(1, Mono)}
            end || #{<<"factors">> := Fs, <<"sympy">> := S} <- Cases],
    % zeros
    [?assertEqual(Class =:= "0", Ours =:= 0, {Class, Ours}) || {Class, _, Ours} <- Rows],
    NonZero = [{Class, Sign, C, Key} || {Class, Sign, {C, Key}} <- Rows],
    % one pq_tensor form per SymPy class, and a different one for each class
    ByClass = maps:groups_from_list(fun({Class, _, _, _}) -> Class end, NonZero),
    Forms = maps:map(fun(_, L) -> lists:usort([Key || {_, _, _, Key} <- L]) end, ByClass),
    [?assertMatch([_], F) || F <- maps:values(Forms)],
    ?assertEqual(map_size(Forms), length(lists:usort([hd(F) || F <- maps:values(Forms)]))),
    % the sign relative to the class's form agrees: SymPy's sign times ours is constant within a class
    [?assertMatch([_], lists:usort([S * C || {_, S, C, _} <- L])) || L <- maps:values(ByClass)].

basics_test() ->
    % R_aacd = 0; R_abcd R_bacd = -R_abcd R_abcd
    ?assertEqual(0, pq_tensor:canonical(1, [{riemann, [a, a, c, d]}])),
    {C1, K1} = pq_tensor:canonical(1, [{riemann, [a, b, c, d]}, {riemann, [a, b, c, d]}]),
    ?assertEqual({-C1, K1}, pq_tensor:canonical(1, [{riemann, [a, b, c, d]}, {riemann, [b, a, c, d]}])),
    % free indices keep their names: R_abcd R_cdef with a, b, e, f free is R_efab R_ab.. reordered consistently
    {_, K2} = pq_tensor:canonical(1, [{riemann, [a, b, c, d]}, {riemann, [c, d, e, f]}]),
    {_, K3} = pq_tensor:canonical(1, [{riemann, [c, d, e, f]}, {riemann, [a, b, c, d]}]),
    ?assertEqual(K2, K3).

invariants_test_() ->
    {timeout, 300, fun() ->
        Q = pq_tensor:full_contractions(riemann, 2),
        C = pq_tensor:full_contractions(riemann, 3),
        % the classes under the slot symmetries alone, then the independent ones (the cyclic identity included)
        ?assertEqual(4, length(Q)),
        ?assertEqual(13, length(C)),
        ?assertEqual(3, pq_tensor:independent_invariants(Q, 4, 8)),
        ?assertEqual(8, pq_tensor:independent_invariants(C, 6, length(C) + 4))
    end}.
