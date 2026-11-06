%% Common subexpressions across straight-line programs (compile/3 with optimize => cse | horner_cse;
%% docs/decisions/H2-cse.md): the same values as the default program within rounding, never more operations than
%% the program the pass starts from, deterministic.
-module(pq_cse_tests).
-include_lib("eunit/include/eunit.hrl").

s(X) -> {sym, X}.

ops({slp, 1, _, Ops, _}) -> length([x || Op <- Ops, lists:member(element(1, Op), [mul, add, pow, exp, conj])]).

prog(Exprs, Syms, O) -> {ok, P} = pq_symbolic:compile(Exprs, Syms, #{optimize => O}), P.

close(A, B) -> abs(c(A) - c(B)) =< 1.0e-12 * max(1.0, abs(c(A))).
c({complex, R, I}) -> math:sqrt(R * R + I * I) * (if R < 0 -> -1; true -> 1 end) + I * 1.0e-3;
c(X) -> X.

%% every program of the four options has the default's values; cse never adds to the program it starts from
check(Exprs, Syms, B) ->
    [D, C, H, HC] = [prog(Exprs, Syms, O) || O <- [none, cse, horner, horner_cse]],
    ?assert(ops(C) =< ops(D)),
    ?assert(ops(HC) =< ops(H)),
    {ok, Vd} = pq_symbolic:eval_slp(D, B),
    [begin {ok, V} = pq_symbolic:eval_slp(P, B), ?assert(close_all(Vd, V)) end || P <- [C, H, HC]],
    {D, C, H, HC}.

close_all(A, B) when is_list(A) -> lists:all(fun({X, Y}) -> close(X, Y) end, lists:zip(A, B));
close_all(A, B) ->
    {ok, {cmat, Ra}} = pq_core:unpack(A),
    {ok, {cmat, Rb}} = pq_core:unpack(B),
    lists:all(fun({{complex, R1, I1}, {complex, R2, I2}}) -> abs(R1 - R2) + abs(I1 - I2) < 1.0e-12 end,
              lists:zip(lists:append(Ra), lists:append(Rb))).

%% the pair search itself: a + b + d and a + c + d share a + d, which the default's chains (a + b) + d and
%% (a + c) + d don't; products likewise
shared_pair_test() ->
    Es = [begin {ok, E} = pq_symbolic:parse(T), E end || T <- [<<"a + b + d">>, <<"a + c + d">>, <<"a*b*d">>, <<"a*c*d">>]],
    {D, C, _, _} = check(Es, [a, b, c, d], #{a => 0.5, b => 1.5, c => -2.0, d => 0.25}),
    ?assertEqual(8, ops(D)),
    ?assertEqual(6, ops(C)).

%% random polynomials in five symbols, several outputs sharing terms
random_polynomials_test() ->
    rand:seed(exsss, {3, 1, 4}),
    Xs = [s(a), s(b), s(c), s(d), s(e)],
    Term = fun() -> {mul, [rand:uniform(9) - 5 | [{pow, lists:nth(rand:uniform(5), Xs), rand:uniform(4)}
                                                 || _ <- lists:seq(1, rand:uniform(4))]]} end,
    B = #{a => 0.3, b => -1.2, c => 0.7, d => 1.9, e => -0.4},
    [begin
         Shared = [Term() || _ <- lists:seq(1, 3 + rand:uniform(10))],
         Es = [{add, Shared ++ [Term() || _ <- lists:seq(1, rand:uniform(20))]} || _ <- lists:seq(1, 1 + rand:uniform(3))],
         check(Es, [a, b, c, d, e], B)
     end || _ <- lists:seq(1, 40)].

%% a resultant (the benchmark of the FORM comparison): fewer operations than Horner alone
resultant_test() ->
    {ok, E} = pq_symbolic:parse(<<"a0^2*b2^3 - a0*a1*b1*b2^2 + a1^2*b0*b2^2 - 2*a0*a2*b0*b2^2 + a0*a2*b1^2*b2 - a1*a2*b0*b1*b2 + a2^2*b0^2*b2 + a0^2*b1*b2*b3 + a3*a1*b0*b3^2">>),
    Syms = pq_symbolic:free_symbols(E),
    B = maps:from_list([{X, 0.1 * K + 0.3} || {K, X} <- lists:enumerate(Syms)]),
    {D, C, H, HC} = check([E], Syms, B),
    ?assert(ops(C) < ops(D)),
    ?assert(ops(HC) < ops(H)),
    %% deterministic: the same program twice
    ?assertEqual(HC, prog([E], Syms, horner_cse)).

%% radicals, exponentials, conjugates, complex symbols, fractional powers
other_atoms_test() ->
    {ok, E} = pq_symbolic:parse(<<"sqrt(2)*x^2*y + exp(i*t)*x*y^2 + x*y*conj(z) + 3*x^2 + sqrt(x)*y + x*y*exp(i*t)">>,
                                #{domains => #{x => positive, y => real, t => real, z => complex}}),
    check([E], [t, x, y, z], #{t => 0.4, x => 1.3, y => -0.6, z => {complex, 0.2, 0.9}}).

%% outputs that are a symbol, a constant, zero, or the same expression twice
trivial_outputs_test() ->
    {ok, E} = pq_symbolic:parse(<<"x*y + x*z">>),
    {D, C, _, HC} = check([s(x), 3, 0, E, E], [x, y, z], #{x => 2.0, y => 3.0, z => 5.0}),
    {ok, V} = pq_symbolic:eval_slp(HC, #{x => 2.0, y => 3.0, z => 5.0}),
    ?assertEqual([2.0, 3.0, 0.0, 16.0, 16.0], V),
    ?assert(ops(C) =< ops(D)).

%% a matrix: the symbolic unitary of a mesh, every entry within rounding
matrix_test() ->
    {ok, C0} = pq_circuit:new(3),
    Is = [{beamsplitter, 0, 1, s(t1), s(p1)}, {beamsplitter, 1, 2, s(t2), s(p2)}, {beamsplitter, 0, 1, s(t3), s(p3)}],
    C = lists:foldl(fun(I, A) -> pq_circuit:add(A, I) end, C0, Is),
    {ok, U} = pq_symbolic:circuit_unitary(C),
    Syms = pq_symbolic:free_symbols(U),
    check([U], Syms, maps:from_list([{X, 0.2 * K} || {K, X} <- lists:enumerate(Syms)])).
