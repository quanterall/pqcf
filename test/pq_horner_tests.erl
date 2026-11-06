%% Horner schemes in straight-line programs (compile/3 with optimize => horner; docs/decisions/H1-horner.md): the
%% same values as the default program within rounding, never more operations, deterministic.
-module(pq_horner_tests).
-include_lib("eunit/include/eunit.hrl").

s(X) -> {sym, X}.

ops({slp, 1, _, Ops, _}) -> length([x || Op <- Ops, lists:member(element(1, Op), [mul, add, pow, exp, conj])]).

progs(Exprs, Syms) ->
    {ok, D} = pq_symbolic:compile(Exprs, Syms),
    {ok, H} = pq_symbolic:compile(Exprs, Syms, #{optimize => horner}),
    {D, H}.

close(A, B) -> abs(c(A) - c(B)) =< 1.0e-12 * max(1.0, abs(c(A))).
c({complex, R, I}) -> math:sqrt(R * R + I * I) * (if R < 0 -> -1; true -> 1 end) + I * 1.0e-3;
c(X) -> X.

%% random polynomials in five symbols: values within 1e-12, no more operations
random_polynomials_test() ->
    rand:seed(exsss, {7, 8, 9}),
    Xs = [s(a), s(b), s(c), s(d), s(e)],
    [begin
         Terms = [{mul, [rand:uniform(9) - 5 | [{pow, lists:nth(rand:uniform(5), Xs), rand:uniform(4)} || _ <- lists:seq(1, rand:uniform(4))]]}
                  || _ <- lists:seq(1, 5 + rand:uniform(40))],
         E = {add, Terms},
         {D, H} = progs([E], [a, b, c, d, e]),
         ?assert(ops(H) =< ops(D)),
         B = #{a => 0.3, b => -1.2, c => 0.7, d => 1.9, e => -0.4},
         {ok, [Vd]} = pq_symbolic:eval_slp(D, B),
         {ok, [Vh]} = pq_symbolic:eval_slp(H, B),
         ?assert(close(Vd, Vh))
     end || _ <- lists:seq(1, 60)].

%% a resultant (the benchmark of the FORM comparison): far fewer operations, the same value
resultant_test() ->
    {ok, E} = pq_symbolic:parse(<<"a0^2*b2^3 - a0*a1*b1*b2^2 + a1^2*b0*b2^2 - 2*a0*a2*b0*b2^2 + a0*a2*b1^2*b2 - a1*a2*b0*b1*b2 + a2^2*b0^2*b2 + a0^2*b1*b2*b3 + a3*a1*b0*b3^2">>),
    Syms = pq_symbolic:free_symbols(E),
    {D, H} = progs([E], Syms),
    ?assert(ops(H) < ops(D)),
    B = maps:from_list([{X, 0.1 * K + 0.3} || {K, X} <- lists:enumerate(Syms)]),
    {ok, [Vd]} = pq_symbolic:eval_slp(D, B),
    {ok, [Vh]} = pq_symbolic:eval_slp(H, B),
    ?assert(close(Vd, Vh)),
    %% deterministic: the same program twice
    ?assertEqual({ok, H}, pq_symbolic:compile([E], Syms, #{optimize => horner})).

%% radicals, exponentials, conjugates, complex symbols: the atoms Horner does not split go through unchanged
other_atoms_test() ->
    {ok, E} = pq_symbolic:parse(<<"sqrt(2)*x^2*y + exp(i*t)*x*y^2 + x*y*conj(z) + 3*x^2 + sqrt(x)*y + x*y*exp(i*t)">>,
                                #{domains => #{x => positive, y => real, t => real, z => complex}}),
    Syms = [t, x, y, z],
    {D, H} = progs([E], Syms),
    B = #{t => 0.4, x => 1.3, y => -0.6, z => {complex, 0.2, 0.9}},
    {ok, [Vd]} = pq_symbolic:eval_slp(D, B),
    {ok, [Vh]} = pq_symbolic:eval_slp(H, B),
    ?assert(close(Vd, Vh)),
    ?assert(ops(H) =< ops(D)).

%% a matrix: the symbolic unitary of a mesh, every entry within rounding
matrix_test() ->
    {ok, C0} = pq_circuit:new(3),
    Is = [{beamsplitter, 0, 1, s(t1), s(p1)}, {beamsplitter, 1, 2, s(t2), s(p2)}, {beamsplitter, 0, 1, s(t3), s(p3)}],
    C = lists:foldl(fun(I, A) -> pq_circuit:add(A, I) end, C0, Is),
    {ok, U} = pq_symbolic:circuit_unitary(C),
    Syms = pq_symbolic:free_symbols(U),
    {D, H} = progs([U], Syms),
    B = maps:from_list([{X, 0.2 * K} || {K, X} <- lists:enumerate(Syms)]),
    {ok, Md} = pq_symbolic:eval_slp(D, B),
    {ok, Mh} = pq_symbolic:eval_slp(H, B),
    {ok, {cmat, Rd}} = pq_core:unpack(Md),
    {ok, {cmat, Rh}} = pq_core:unpack(Mh),
    [?assert(abs(R1 - R2) + abs(I1 - I2) < 1.0e-12) || {{complex, R1, I1}, {complex, R2, I2}} <- lists:zip(lists:append(Rd), lists:append(Rh))],
    ?assert(ops(H) =< ops(D)).

%% the default is unchanged: compile/3 without the option is compile/2
default_unchanged_test() ->
    {ok, E} = pq_symbolic:parse(<<"x^2*y + 2*x*y + y^2">>),
    ?assertEqual(pq_symbolic:compile([E], [x, y]), pq_symbolic:compile([E], [x, y], #{})).
