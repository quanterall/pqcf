-module(pq_symbolic_tests).
-include_lib("proper/include/proper.hrl").
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

-define(S, pq_symbolic).

x() -> {sym, x}.
y() -> {sym, y}.
t() -> {sym, t}.

canon(E) -> {ok, C} = ?S:canonicalise(E), C.
same(A, B) -> canon(A) =:= canon(B).
val(E, B) -> {ok, V} = ?S:evaluate(E, B), V.

%%% GT-S01: identities (SY 9.3)

pythagoras_test() ->
    ?assertEqual(1, canon({add, [{pow, {sin, t()}, 2}, {pow, {cos, t()}, 2}]})).

angle_addition_test() ->
    A = {sym, a}, B = {sym, b},
    ?assert(same({cos, {add, [A, B]}}, {add, [{mul, [{cos, A}, {cos, B}]}, {mul, [-1, {sin, A}, {sin, B}]}]})),
    ?assert(same({sin, {mul, [2, t()]}}, {mul, [2, {sin, t()}, {cos, t()}]})).

special_values_test() ->
    Pi = {const, pi},
    ?assert(same({cos, {mul, [{rat, 1, 4}, Pi]}}, {mul, [{rat, 1, 2}, {pow, 2, {rat, 1, 2}}]})),
    ?assert(same({cos, {mul, [{rat, 1, 12}, Pi]}},
                 {mul, [{rat, 1, 4}, {add, [{pow, 6, {rat, 1, 2}}, {pow, 2, {rat, 1, 2}}]}]})),
    ?assertEqual(-1, canon({exp, {mul, [{complex, 0, 1}, Pi]}})),
    ?assertEqual({complex, 0, 1}, canon({exp, {mul, [{complex, 0, {rat, 1, 2}}, Pi]}})),
    %% exp(pi/4 + i pi/3): the real part stays a factor, the imaginary part is evaluated
    ?assert(same({exp, {mul, [{complex, {rat, 1, 4}, {rat, 1, 3}}, Pi]}},
                 {mul, [{exp, {mul, [{rat, 1, 4}, Pi]}}, {add, [{rat, 1, 2}, {mul, [{complex, 0, {rat, 1, 2}}, {pow, 3, {rat, 1, 2}}]}]}]})),
    %% every multiple of pi/12 against the floating-point value (SS4)
    [begin
         V = val({exp, {mul, [{complex, 0, {rat, K, 12}}, Pi]}}, #{}),
         {Re, Im} = case V of {complex, R, I} -> {R, I}; R -> {R, 0.0} end,
         ?assert(abs(Re - math:cos(K * math:pi() / 12)) < 1.0e-14),
         ?assert(abs(Im - math:sin(K * math:pi() / 12)) < 1.0e-14)
     end || K <- lists:seq(-30, 30)].

radicals_test() ->
    R = fun(N) -> {pow, N, {rat, 1, 2}} end,
    ?assertEqual(2, canon({mul, [R(2), R(2)]})),
    ?assert(same({mul, [R(6), R(2)]}, {mul, [2, R(3)]})),
    ?assert(same(R(12), {mul, [2, R(3)]})),
    ?assert(same(R({rat, 1, 2}), {mul, [{rat, 1, 2}, R(2)]})),
    ?assert(same({pow, R(2), -1}, {mul, [{rat, 1, 2}, R(2)]})),
    ?assert(same(R(-4), {complex, 0, 2})),
    ?assert(same({mul, [R(2), R(3), R(6)]}, 6)).

scope_fence_test() ->
    %% (x^2 - 1)/(x - 1) is not reduced to x + 1 (no polynomial gcd, SY 1.2, 10.3)
    Q = {mul, [{add, [{pow, x(), 2}, -1]}, {pow, {add, [x(), -1]}, -1}]},
    ?assertNot(same(Q, {add, [x(), 1]})),
    ?assert(abs(val(Q, #{x => 3}) - 4.0) < 1.0e-12).

primitive_bases_test() ->
    %% 1/(2 + 2x) and 1/(2(1 + x)) are the same form; so are 1/(-1 - x) and -1/(1 + x)
    ?assert(same({pow, {add, [2, {mul, [2, x()]}]}, -1}, {mul, [{rat, 1, 2}, {pow, {add, [1, x()]}, -1}]})),
    ?assert(same({pow, {add, [-1, {mul, [-1, x()]}]}, -1}, {mul, [-1, {pow, {add, [1, x()]}, -1}]})),
    %% sqrt(1 - eta) twice is 1 - eta
    Eta = {sym, eta},
    S = {pow, {add, [1, {mul, [-1, Eta]}]}, {rat, 1, 2}},
    ?assert(same({mul, [S, S]}, {add, [1, {mul, [-1, Eta]}]})).

domains_test() ->
    %% sqrt(x^2) is x only for a nonnegative x
    E = {pow, {pow, x(), 2}, {rat, 1, 2}},
    ?assertEqual({ok, {mul, [1, x()]}}, ?S:canonicalise(E, #{domains => #{x => nonnegative}})),
    ?assertNotEqual({mul, [1, x()]}, canon(E)),
    ?assert(abs(val(E, #{x => -3}) - 3.0) < 1.0e-12),
    %% conj: real symbols are fixed, complex ones are not
    ?assertEqual({mul, [1, x()]}, canon({conj, x()})),
    {ok, Z} = ?S:canonicalise({conj, {sym, z}}, #{domains => #{z => complex}}),
    ?assertEqual({mul, [1, {conj, {sym, z}}]}, Z),
    {ok, ZZ} = ?S:canonicalise({conj, {conj, {sym, z}}}, #{domains => #{z => complex}}),
    ?assertEqual({mul, [1, {sym, z}]}, ZZ).

conj_branch_cut_test() ->
    %% conj((-8)^(1/3)) = 1 - sqrt(3) i, not conj(-8)^(1/3) = 1 + sqrt(3) i
    E = {conj, {pow, {add, [x(), y()]}, {rat, 1, 3}}},
    {complex, Re, Im} = val(E, #{x => -5, y => -3}),
    ?assert(abs(Re - 1) < 1.0e-12 andalso abs(Im + math:sqrt(3)) < 1.0e-12),
    ?assert(same({conj, E}, {pow, {add, [x(), y()]}, {rat, 1, 3}})),
    ?assert(same({mul, [E, E, E]}, {add, [x(), y()]})).

%%% Worked examples (SY 10)

hong_ou_mandel_test() ->
    %% perm of B(theta, 0) is cos^2 - sin^2; p(1,1) is its squared modulus
    Perm = {add, [{pow, {cos, t()}, 2}, {mul, [-1, {pow, {sin, t()}, 2}]}]},
    P = {mul, [Perm, {conj, Perm}]},
    Expected = {add, [{rat, 1, 2}, {mul, [{rat, 1, 4}, {exp, {mul, [{complex, 0, 4}, t()]}}]},
                      {mul, [{rat, 1, 4}, {exp, {mul, [{complex, 0, -4}, t()]}}]}]},
    ?assert(same(P, Expected)),
    ?assert(same(P, {add, [{rat, 1, 2}, {mul, [{rat, 1, 2}, {cos, {mul, [4, t()]}}]}]})),
    ?assert(abs(val(P, #{t => math:pi() / 4})) < 1.0e-15),
    ?assertEqual(<<"1/2 + 1/4*exp(-4*i*t) + 1/4*exp(4*i*t)">>, ?S:pretty(P)).

mach_zehnder_test() ->
    Th = {sym, theta},
    P = {add, [{rat, 1, 2}, {mul, [{rat, -1, 2}, {cos, Th}]}]},
    ?assert(same(P, {add, [{rat, 1, 2}, {mul, [{rat, -1, 4}, {exp, {mul, [{complex, 0, 1}, Th]}}]},
                           {mul, [{rat, -1, 4}, {exp, {mul, [{complex, 0, -1}, Th]}}]}]})).

%%% Matrices

beamsplitter() ->
    {matrix, [[{cos, t()}, {mul, [-1, {exp, {mul, [{complex, 0, -1}, {sym, phi}]}}, {sin, t()}]}],
              [{mul, [{exp, {mul, [{complex, 0, 1}, {sym, phi}]}}, {sin, t()}]}, {cos, t()}]]}.

unitary_test() ->
    U = beamsplitter(),
    ?assertEqual({ok, {matrix, [[1, 0], [0, 1]]}}, ?S:matmul(U, {dagger, U})),
    ?assertEqual({ok, {matrix, [[1, 0], [0, 1]]}}, ?S:canonicalise({matmul, {dagger, U}, U})),
    %% one sign changed: no longer unitary
    {matrix, [[A, B], [C, D]]} = U,
    {ok, {matrix, [[_, Off], _]}} = ?S:matmul({matrix, [[A, B], [{mul, [-1, C]}, D]]},
                                             {dagger, {matrix, [[A, B], [{mul, [-1, C]}, D]]}}),
    ?assertNotEqual(0, Off).

matrix_ops_test() ->
    M = {matrix, [[1, x()], [0, 2]]},
    ?assertEqual({ok, {2, 2}}, ?S:shape(M)),
    ?assertEqual({ok, {4, 4}}, ?S:shape({kron, M, {identity, 2}})),
    ?assertEqual({ok, {matrix, [[1, 0], [{mul, [1, x()]}, 2]]}}, ?S:transpose(M)),
    ?assertEqual(?S:canonicalise(M), ?S:matmul({identity, 2}, M)),
    ?assertEqual({ok, {matrix, [[2, {mul, [2, x()]}], [0, 4]]}}, ?S:smul(2, M)),
    ?assertMatch({error, {bad_shape, _}}, ?S:shape({add, [M, 1]})),
    ?assertMatch({error, {bad_shape, _}}, ?S:shape({matmul, M, {matrix, [[1, 2, 3]]}})),
    ?assertMatch({error, {bad_shape, _}}, ?S:canonicalise({matmul, M, {matrix, [[1, 2, 3]]}})),
    {ok, Packed} = ?S:evaluate(M, #{x => 3}),
    ?assertEqual({ok, {cmat, [[{complex, 1.0, 0.0}, {complex, 3.0, 0.0}], [{complex, 0.0, 0.0}, {complex, 2.0, 0.0}]]}},
                 pq_core:unpack(Packed)).

%%% Evaluation, substitution, errors

evaluate_test() ->
    ?assertEqual(5.0, val({add, [x(), y()]}, #{x => 2, y => {rat, 6, 2}})),
    ?assertEqual({error, {unbound_symbols, [y]}}, ?S:evaluate({add, [x(), y()]}, #{x => 1})),
    ?assertEqual({error, {domain_violation, {x, positive}}},
                 ?S:evaluate(x(), #{x => -1}, #{domains => #{x => positive}})),
    ?assertMatch({error, {non_finite, _}}, ?S:evaluate({pow, x(), -1}, #{x => 0})),
    ?assertEqual({complex, 0.0, 1.0}, val({pow, x(), {rat, 1, 2}}, #{x => -1})),
    ?assertEqual({ok, {rat, 53, 44}}, ?S:evaluate({add, [{pow, x(), 2}, {pow, y(), -1}]}, #{x => {rat, 1, 2}, y => {rat, 22, 21}},
                                                #{exact => true})),
    ?assertEqual({error, {symbolic_unsupported, exact}}, ?S:evaluate({const, pi}, #{}, #{exact => true})).

subs_test() ->
    ?assertEqual({ok, {add, [1, {mul, [2, y()]}]}},
                 ?S:subs({add, [x(), y()]}, #{x => {add, [1, y()]}})),
    %% simultaneous: the replacements are not substituted again
    ?assertEqual(?S:canonicalise({add, [x(), y()]}), ?S:subs({add, [x(), y()]}, #{x => y(), y => x()})),
    ?assertEqual({error, {domain_violation, {x, positive}}},
                 ?S:subs(x(), #{x => y()}, #{domains => #{x => positive}})),
    ?assertMatch({ok, _}, ?S:subs(x(), #{x => {exp, y()}}, #{domains => #{x => positive}})),
    ?assertEqual({ok, 3}, ?S:canonicalise({subs, {add, [x(), 1]}, #{x => 2}})).

errors_test() ->
    ?assertEqual({error, {bad_expr, {foo}}}, ?S:canonicalise({add, [{foo}]})),
    ?assertEqual({error, {symbolic_limit, #{reason => poly_terms, terms => 400, max => 100}}},
                 ?S:canonicalise({pow, {add, [{sym, a}, {sym, b}, {sym, c}, {sym, d}, {sym, e}, {sym, f}, {sym, g},
                                               {sym, h}, {sym, k}, {sym, l}, {sym, m}, {sym, n}, {sym, o}, {sym, p},
                                               {sym, q}, {sym, r}, {sym, s}, {sym, u}, {sym, v}, {sym, w}]}, 2},
                                 #{max_poly_terms => 100})),
    %% a large power of a sum stays as Base(B)^N
    {ok, Big} = ?S:canonicalise({pow, {add, [x(), 1]}, 9}),
    ?assertMatch({mul, [1, {pow, {add, _}, 9}]}, Big).

%%% Printing and parsing

parse_test() ->
    {ok, E} = ?S:parse(<<"-x^2 + 2*sin(theta)/3 - sqrt(2) + e + 1.5e-3*i">>),
    ?assert(same(E, {add, [{mul, [-1, {pow, x(), 2}]}, {mul, [{rat, 2, 3}, {sin, {sym, theta}}]},
                           {mul, [-1, {pow, 2, {rat, 1, 2}}]}, {const, e}, {complex, 0, 1.5e-3}]})),
    ?assertEqual({ok, {pow, x(), {rat, -1, 3}}}, ?S:parse(<<"x^(-2/6)">>)),
    ?assertEqual({error, {sym_parse, 3}}, ?S:parse(<<"2 x">>)),
    ?assertEqual({error, {sym_parse, 7}}, ?S:parse(<<"(x + 1">>)),
    ?assertEqual({error, {sym_parse, 4}}, ?S:parse(<<"sin x">>)),
    ?assertEqual({error, {sym_parse, 1}}, ?S:parse(<<"$">>)),
    ?assertEqual({error, {sym_parse, 1}}, ?S:parse(list_to_binary(lists:duplicate(65, $a)))),
    ?assertMatch({error, {sym_parse, _}}, ?S:parse(list_to_binary(lists:duplicate(300, $() ++ "x" ++ lists:duplicate(300, $))))).

atom_safety_test() ->
    Fresh = <<"pqcf_never_an_atom_", (integer_to_binary(erlang:unique_integer([positive])))/binary>>,
    ?assertEqual({error, {sym_parse, 5}}, ?S:parse(<<"x + ", Fresh/binary>>, #{allowed_symbols => [x]})),
    ?assertError(badarg, binary_to_existing_atom(Fresh)),
    ?assertEqual({error, {sym_parse, 5}}, ?S:parse(<<"x + y">>, #{allowed_symbols => [x]})),
    ?assertMatch({error, {symbolic_limit, #{reason := max_symbols}}}, ?S:parse(<<"a + b + c">>, #{max_symbols => 2})).

latex_test() ->
    ?assertEqual(<<"\\frac{1}{2} \\cdot \\sqrt{2}">>, ?S:to_latex({cos, {mul, [{rat, 1, 4}, {const, pi}]}})).

%%% Parameter resolution

resolve_test() ->
    Defaults = #{a => {mul, [{rat, 1, 4}, {const, pi}]}, b => {cos, {sym, a}}, c => {mul, [2, {sym, b}]}},
    {ok, R} = ?S:resolve_params(Defaults, #{}),
    ?assert(abs(maps:get(c, R) - math:sqrt(2)) < 1.0e-15),
    {ok, R2} = ?S:resolve_params(Defaults, #{a => 0}),
    ?assertEqual(2.0, maps:get(c, R2)),
    ?assertEqual({error, {param_cycle, [p, q]}}, ?S:resolve_params(#{p => {sym, q}, q => {sym, p}}, #{})),
    ?assertEqual({error, {unbound_symbols, [z]}}, ?S:resolve_params(#{p => {sym, z}}, #{})).

%%% Bounded execution (SY 7.2, SS5)

bounded_test_() ->
    Big = {pow, {add, [{sym, A} || A <- [a, b, c, d, e, f, g, h, k, l, m, n, o, p, q, r]]}, 8},
    Wide = {add, [Big | lists:duplicate(100, 1)]},
    {timeout, 60,
     [?_assertEqual({error, {symbolic_limit, #{reason => timeout}}},
                    ?S:canonicalise(Wide, #{sym_timeout_ms => 50, max_poly_terms => 100000000})),
      ?_assertEqual({error, {symbolic_limit, #{reason => heap}}},
                    ?S:canonicalise(Wide, #{sym_max_heap_words => 200000, max_poly_terms => 100000000})),
      ?_assertEqual({error, {symbolic_limit, #{reason => max_expr_nodes, max => 50}}},
                    ?S:canonicalise(Wide, #{max_expr_nodes => 50})),
      ?_test(begin
                 Before = length(processes()),
                 _ = ?S:canonicalise(Wide, #{sym_timeout_ms => 20, max_poly_terms => 100000000}),
                 timer:sleep(50),
                 ?assert(length(processes()) =< Before)
             end)]}.

%%% Properties (SY 9.1)

%% expressions over x, y (real) and p (positive)
leaf() -> oneof([integer(-3, 3), ?LET({N, D}, {integer(-5, 5), integer(1, 4)}, {rat, N, D}),
                 {sym, x}, {sym, y}, {sym, p}, {const, pi}, {complex, 0, 1}]).

expr(0) -> leaf();
expr(N) ->
    S = expr(N - 1),
    frequency([{3, leaf()},
               {3, ?LET(L, vector(2, S), {add, L})},
               {3, ?LET(L, vector(2, S), {mul, L})},
               {2, ?LET({B, E}, {S, integer(-2, 2)}, {pow, B, E})},
               {1, {sin, S}}, {1, {cos, S}},
               {1, ?LET(A, S, {exp, {mul, [{rat, 1, 4}, A]}})},
               {1, {conj, S}},
               {1, ?LET(A, poly(2), {pow, A, {rat, 1, 2}})},
               {1, ?LET(A, poly(2), {pow, A, {rat, 1, 3}})}]).

%% real polynomials in x, y, p (a real argument for the roots)
poly(0) -> oneof([integer(-3, 3), {sym, x}, {sym, y}, {sym, p}]);
poly(N) -> ?LET({Op, A, B}, {oneof([add, mul]), poly(N - 1), poly(N - 1)}, {Op, [A, B]}).

binding() -> ?LET({X, Y, P}, {float(-2.0, 2.0), float(-2.0, 2.0), float(0.1, 2.0)}, #{x => X, y => Y, p => P}).

-define(DOM, #{domains => #{p => positive}, max_poly_terms => 3000, sym_timeout_ms => 2000}).

%% an evaluator of the raw tree that shares nothing with the canonical form
raw(N, _) when is_integer(N) -> {float(N), 0.0};
raw({rat, N, D}, _) -> {N / D, 0.0};
raw({complex, R, I}, B) -> {element(1, raw(R, B)), element(1, raw(I, B))};
raw({const, pi}, _) -> {math:pi(), 0.0};
raw({sym, X}, B) -> {maps:get(X, B), 0.0};
raw({add, L}, B) -> lists:foldl(fun(E, {R, I}) -> {R2, I2} = raw(E, B), {R + R2, I + I2} end, {0.0, 0.0}, L);
raw({mul, L}, B) -> lists:foldl(fun(E, A) -> cmul(A, raw(E, B)) end, {1.0, 0.0}, L);
raw({pow, A, N}, B) when is_integer(N) -> ipow(raw(A, B), N);
raw({pow, A, {rat, N, D}}, B) -> cpow(raw(A, B), N / D);
raw({exp, A}, B) -> cexp(raw(A, B));
raw({sin, A}, B) -> {R, I} = raw(A, B), {math:sin(R) * math:cosh(I), math:cos(R) * math:sinh(I)};
raw({cos, A}, B) -> {R, I} = raw(A, B), {math:cos(R) * math:cosh(I), -math:sin(R) * math:sinh(I)};
raw({conj, A}, B) -> {R, I} = raw(A, B), {R, -I}.

cmul({A, B}, {C, D}) -> {A * C - B * D, A * D + B * C}.
ipow(_, 0) -> {1.0, 0.0};
ipow(Z, N) when N < 0 -> {R, I} = ipow(Z, -N), M = R * R + I * I, {R / M, -I / M};
ipow(Z, N) -> cmul(Z, ipow(Z, N - 1)).
cpow({R, I}, _W) when R == 0, I == 0 -> {0.0, 0.0};
cpow({R, I}, W) -> L = {math:log(math:sqrt(R * R + I * I)) * W, math:atan2(I + 0.0, R) * W}, cexp(L).
cexp({R, I}) -> M = math:exp(R), {M * math:cos(I), M * math:sin(I)}.

as_pair({complex, R, I}) -> {R, I};
as_pair(R) when is_float(R) -> {R, 0.0}.

close({A, B}, {C, D}) ->
    Scale = max(1.0, max(abs(A) + abs(B), abs(C) + abs(D))),
    abs(A - C) + abs(B - D) =< 1.0e-9 * Scale.

%% P1: the canonical form has the value of the expression
prop_value() ->
    ?FORALL({E, B}, {expr(3), binding()},
            begin
                Raw = try raw(E, B) catch error:badarith -> skip end,
                case Raw of
                    skip -> true;
                    {R, I} when abs(R) > 1.0e6; abs(I) > 1.0e6 -> true;
                    _ ->
                        case ?S:evaluate(E, B, ?DOM) of
                            {ok, V} -> close(Raw, as_pair(V));
                            {error, {non_finite, _}} -> true;
                            {error, {symbolic_limit, _}} -> true
                        end
                end
            end).

%% P2: idempotence
prop_idempotent() ->
    ?FORALL(E, expr(3),
            case ?S:canonicalise(E, ?DOM) of
                {ok, C} -> {ok, C} =:= ?S:canonicalise(C, ?DOM);
                {error, _} -> true
            end).

%% P8: parse(pretty(E)) normalises to E
prop_round_trip() ->
    ?FORALL(E, expr(3),
            case ?S:canonicalise(E, ?DOM) of
                {ok, C} ->
                    Text = pq_sym_print:ascii(C),
                    {ok, P} = ?S:parse(Text),
                    ?WHENFAIL(io:format("~s~n", [Text]), {ok, C} =:= ?S:canonicalise(P, ?DOM));
                {error, _} -> true
            end).

%% P1 exact: rational bindings give the exact value of the algebraic fragment
afrag(0) -> oneof([integer(-3, 3), ?LET({N, D}, {integer(-5, 5), integer(1, 4)}, {rat, N, D}), {sym, x}, {sym, y}]);
afrag(N) ->
    S = afrag(N - 1),
    frequency([{2, afrag(0)}, {3, ?LET(L, vector(2, S), {add, L})}, {3, ?LET(L, vector(2, S), {mul, L})},
               {1, ?LET({B, E}, {S, integer(-2, 3)}, {pow, B, E})}]).

rraw(N, _) when is_integer(N) -> {N, 1};
rraw({rat, N, D}, _) -> pq_sym_exact:rat(N, D);
rraw({sym, X}, B) -> maps:get(X, B);
rraw({add, L}, B) -> lists:foldl(fun(E, A) -> pq_sym_exact:radd(A, rraw(E, B)) end, {0, 1}, L);
rraw({mul, L}, B) -> lists:foldl(fun(E, A) -> pq_sym_exact:rmul(A, rraw(E, B)) end, {1, 1}, L);
rraw({pow, A, N}, B) -> rpow(rraw(A, B), N).
rpow(_, 0) -> {1, 1};
rpow(R, N) when N < 0 -> pq_sym_exact:rinv(rpow(R, -N));
rpow(R, N) -> pq_sym_exact:rmul(R, rpow(R, N - 1)).

prop_exact() ->
    ?FORALL({E, X, Y}, {afrag(3), {integer(-7, 7), integer(1, 5)}, {integer(-7, 7), integer(1, 5)}},
            begin
                Xr = pq_sym_exact:rat(element(1, X), element(2, X)),
                Yr = pq_sym_exact:rat(element(1, Y), element(2, Y)),
                case try rraw(E, #{x => Xr, y => Yr}) catch error:badarith -> zero_division end of
                    {N, D} ->
                        Want = case D of 1 -> N; _ -> {rat, N, D} end,
                        {ok, Want} =:= ?S:evaluate(E, #{x => {rat, element(1, Xr), element(2, Xr)},
                                                        y => {rat, element(1, Yr), element(2, Yr)}}, #{exact => true});
                    zero_division -> true
                end
            end).

%% the property bites: a value off by 1e-6 is caught
prop_mutant() ->
    ?FORALL({E, B}, {expr(3), binding()},
            case try raw(E, B) catch error:badarith -> skip end of
                {R, I} when abs(R) < 1.0e6, abs(I) < 1.0e6 ->
                    case ?S:evaluate({add, [E, {rat, 1, 1000000}]}, B, ?DOM) of
                        {ok, V} -> close({R, I}, as_pair(V));
                        {error, _} -> true
                    end;
                _ -> true
            end).

%% the generator makes expressions worth testing
generator_test() ->
    Sizes = [?S:size(element(2, proper_gen:pick(expr(3), 20))) || _ <- lists:seq(1, 1000)],
    Mean = lists:sum(Sizes) / 1000,
    io:format(user, "~nexpr(3): mean ~.1f nodes, largest ~p~n", [Mean, lists:max(Sizes)]),
    ?assert(Mean > 8).

properties_test_() ->
    Q = fun(P, N) -> ?_assert(proper:quickcheck(P, [{numtests, N}, {to_file, user}])) end,
    {timeout, 600, [{timeout, 150, Q(prop_value(), 1000)}, {timeout, 150, Q(prop_idempotent(), 600)},
                    {timeout, 150, Q(prop_round_trip(), 600)}, {timeout, 150, Q(prop_exact(), 1000)},
                    {timeout, 150, ?_assertNot(proper:quickcheck(prop_mutant(), [{numtests, 1000}, quiet]))}]}.

%%% Differentiation (SY 5.3, P6)

d(E, X) -> {ok, D} = ?S:differentiate(E, X), D.

diff_examples_test() ->
    T = {sym, t},
    %% SY 10.1: d/dtheta of the Hong-Ou-Mandel probability is -2 sin(4 theta)
    Perm = {add, [{pow, {cos, T}, 2}, {mul, [-1, {pow, {sin, T}, 2}]}]},
    P = {mul, [Perm, {conj, Perm}]},
    ?assert(same(d(P, t), {mul, [-2, {sin, {mul, [4, T]}}]})),
    %% SY 10.2: d/dTheta (1/2 - cos(Theta)/2) = sin(Theta)/2
    ?assert(same(d({add, [{rat, 1, 2}, {mul, [{rat, -1, 2}, {cos, T}]}]}, t), {mul, [{rat, 1, 2}, {sin, T}]})),
    ?assert(same(d({pow, x(), 5}, x), {mul, [5, {pow, x(), 4}]})),
    ?assert(same(d({pow, {add, [1, x()]}, -1}, x), {mul, [-1, {pow, {add, [1, x()]}, -2}]})),
    ?assert(same(d({exp, {pow, x(), 2}}, x), {mul, [2, x(), {exp, {pow, x(), 2}}]})),
    ?assertEqual(0, d({add, [y(), {const, pi}, {pow, 2, {rat, 1, 2}}]}, x)),
    %% sqrt(1 - x) for a real x: -1/2 (1 - x)^(-1/2)
    ?assert(same(d({pow, {add, [1, {mul, [-1, x()]}]}, {rat, 1, 2}}, x),
                 {mul, [{rat, -1, 2}, {pow, {add, [1, {mul, [-1, x()]}]}, {rat, -1, 2}}]})),
    {ok, D2} = ?S:differentiate({sin, x()}, x, #{order => 2}),
    ?assert(same(D2, {mul, [-1, {sin, x()}]})),
    {ok, [Gx, Gy]} = ?S:gradient({mul, [x(), {pow, y(), 2}]}, [x, y]),
    ?assert(same(Gx, {pow, y(), 2}) andalso same(Gy, {mul, [2, x(), y()]})),
    {ok, J} = ?S:jacobian({matrix, [[{mul, [x(), y()]}], [{add, [x(), y()]}]]}, [x, y]),
    ?assertEqual(?S:canonicalise({matrix, [[y(), x()], [1, 1]]}), {ok, J}),
    {ok, DM} = ?S:differentiate({matrix, [[x(), {pow, x(), 2}]]}, x),
    ?assertEqual(?S:canonicalise({matrix, [[1, {mul, [2, x()]}]]}), {ok, DM}).

diff_errors_test() ->
    Z = #{domains => #{z => complex}},
    ?assertEqual({error, {differentiate_unsupported, {sym, z}}}, ?S:differentiate({sym, z}, z, Z)),
    %% a real variable with a complex symbol around: fine, conj(z) does not depend on x
    ?assertMatch({ok, _}, ?S:differentiate({mul, [x(), {conj, {sym, z}}]}, x, Z)),
    ?assertMatch({error, {bad_shape, _}}, ?S:jacobian({matrix, [[x(), y()]]}, [x])).

%% P6: product and chain rules as canonical identities; mixed second derivatives commute
diff_rules_test() ->
    F = {add, [{sin, {mul, [x(), y()]}}, {pow, x(), 3}]},
    G = {exp, {mul, [{rat, 1, 2}, y(), x()]}},
    ?assert(same(d({mul, [F, G]}, x), {add, [{mul, [d(F, x), G]}, {mul, [F, d(G, x)]}]})),
    ?assert(same(d({sin, F}, x), {mul, [{cos, F}, d(F, x)]})),
    ?assert(same(d(d(F, x), y), d(d(F, y), x))).

%% P6: first derivatives agree with central differences within 1e-6 (relative) on random expressions
prop_diff_numeric() ->
    ?FORALL({E, B}, {expr(3), binding()},
            case ?S:differentiate(E, x, ?DOM) of
                {ok, D} ->
                    H = 1.0e-6,
                    X0 = maps:get(x, B),
                    F = fun(X) -> ?S:evaluate(E, B#{x => X}, ?DOM) end,
                    case {?S:evaluate(D, B, ?DOM), F(X0 + H), F(X0 - H)} of
                        {{ok, DV}, {ok, Fp}, {ok, Fm}} ->
                            {R1, I1} = as_pair(Fp), {R2, I2} = as_pair(Fm),
                            Num = {(R1 - R2) / (2 * H), (I1 - I2) / (2 * H)},
                            {A, Bi} = as_pair(DV),
                            Big = abs(R1) + abs(I1) + abs(A) + abs(Bi) > 1.0e4,
                            %% near a branch cut or a pole a difference quotient means nothing
                            Big orelse abs(A - element(1, Num)) + abs(Bi - element(2, Num))
                                       =< 1.0e-6 * max(1.0, abs(A) + abs(Bi)) orelse near_cut(E, B);
                        _ -> true
                    end;
                {error, _} -> true
            end).

%% the expression has a fractional power whose argument is within 1e-3 of the negative real axis or of 0
near_cut(E, B) ->
    Args = [A || {pow, A, {rat, _, _}} <- subterms(E)],
    lists:any(fun(A) ->
                      case ?S:evaluate(A, B, ?DOM) of
                          {ok, V} -> {R, I} = as_pair(V), abs(I) < 1.0e-3 andalso R < 1.0e-3;
                          _ -> true
                      end
              end, Args).

subterms(T) when is_tuple(T) -> [T | lists:append([subterms(X) || X <- tuple_to_list(T)])];
subterms(L) when is_list(L) -> lists:append([subterms(X) || X <- L]);
subterms(_) -> [].

%% the property bites: a derivative off by 0.1 % (plus 1e-3) is caught
prop_diff_mutant() ->
    ?FORALL({E, B}, {expr(3), binding()},
            case ?S:differentiate(E, x, ?DOM) of
                {ok, D} ->
                    H = 1.0e-6,
                    X0 = maps:get(x, B),
                    F = fun(X) -> ?S:evaluate(E, B#{x => X}, ?DOM) end,
                    case {?S:evaluate({add, [{mul, [{rat, 1001, 1000}, D]}, {rat, 1, 1000}]}, B, ?DOM),
                          F(X0 + H), F(X0 - H)} of
                        {{ok, DV}, {ok, Fp}, {ok, Fm}} ->
                            {R1, I1} = as_pair(Fp), {R2, I2} = as_pair(Fm),
                            {A, Bi} = as_pair(DV),
                            abs(A - (R1 - R2) / (2 * H)) + abs(Bi - (I1 - I2) / (2 * H))
                                =< 1.0e-6 * max(1.0, abs(A) + abs(Bi));
                        _ -> true
                    end;
                {error, _} -> true
            end).

diff_property_test_() ->
    {timeout, 300, [?_assert(proper:quickcheck(prop_diff_numeric(), [{numtests, 800}, {to_file, user}])),
                    ?_assertNot(proper:quickcheck(prop_diff_mutant(), [{numtests, 800}, quiet]))]}.

%%% Straight-line programs (SY 5.5, FR-S10) and equivalence (SY 5.4, FR-S11)

prop_slp() ->
    ?FORALL({E, B}, {expr(3), binding()},
            case {?S:compile([E], [x, y, p]), ?S:evaluate(E, B, ?DOM)} of
                {{ok, P}, {ok, V}} ->
                    {ok, [W]} = ?S:eval_slp(P, B),
                    close(as_pair(V), as_pair(W));
                _ -> true
            end).

slp_test_() ->
    {timeout, 120,
     [?_assert(proper:quickcheck(prop_slp(), [{numtests, 600}, {to_file, user}])),
      ?_test(begin
                 M = {matrix, [[{cos, {sym, t}}, {mul, [-1, {sin, {sym, t}}]}], [{sin, {sym, t}}, {cos, {sym, t}}]]},
                 {ok, P} = ?S:compile([M], [t]),
                 {ok, [A, B]} = ?S:eval_batch(P, [#{t => 0.3}, #{t => 1.0}]),
                 {ok, {cmat, [[{complex, C, _}, _], _]}} = pq_core:unpack(A),
                 ?assert(abs(C - math:cos(0.3)) < 1.0e-15),
                 ?assertMatch({cmat, 2, 2, _}, B),
                 ?assertEqual({error, {unbound_symbols, [t]}}, ?S:eval_slp(P, #{}))
             end)]}.

equiv_test() ->
    X = {sym, x},
    ?assertEqual({ok, true}, ?S:equiv({pow, {add, [X, 1]}, 2}, {add, [{pow, X, 2}, {mul, [2, X]}, 1]}, #{})),
    ?assertEqual({ok, false}, ?S:equiv({pow, X, 2}, {pow, X, 3}, #{})),
    ?assertEqual({ok, true}, ?S:equiv({add, [{pow, {sin, X}, 2}, {pow, {cos, X}, 2}]}, 1, #{})),
    %% the scope fence (SY 10.3)
    ?assertEqual({ok, unknown}, ?S:equiv({mul, [{add, [{pow, X, 2}, -1]}, {pow, {add, [X, -1]}, -1}]}, {add, [X, 1]}, #{})),
    %% the counter-example to SY 3.5: equal, different forms, so not decidable here
    Cos5 = {cos, {mul, [{rat, 1, 5}, {const, pi}]}},
    Rad5 = {mul, [{rat, 1, 4}, {add, [1, {pow, 5, {rat, 1, 2}}]}]},
    ?assertEqual({ok, unknown}, ?S:equiv(Cos5, Rad5, #{})),
    ?assertEqual({ok, unknown}, ?S:equiv(Cos5, Rad5, #{numeric => 5})),
    %% sqrt(x^2) is |x|, not x: a negative point shows it
    ?assertMatch({ok, {false, #{x := V}}} when V < 0, ?S:equiv({pow, {pow, X, 2}, {rat, 1, 2}}, X, #{numeric => 20})).

%% SS3: the straight-line program against tree evaluation on a symbolic Clements mesh unitary
ss3_test_() ->
    {timeout, 300,
     [?_test(begin
                 [begin
                      Instrs = [{beamsplitter, J, J + 1, {sym, list_to_atom("t" ++ integer_to_list(I))},
                                 {sym, list_to_atom("p" ++ integer_to_list(I))}}
                                || {I, J} <- lists:enumerate([J0 || L <- lists:seq(0, K - 1), J0 <- lists:seq(L rem 2, K - 2, 2)])],
                      {ok, C0} = pq_circuit:new(K),
                      Syms = lists:usort(lists:append([[element(2, element(4, I)), element(2, element(5, I))] || I <- Instrs])),
                      {Tu, {ok, U}} = timer:tc(fun() -> pq_symbolic:circuit_unitary(C0#circuit{instrs = Instrs}) end),
                      {Tc, {ok, P}} = timer:tc(fun() -> ?S:compile([U], Syms) end),
                      B = maps:from_list([{S, 0.1 * N} || {N, S} <- lists:enumerate(Syms)]),
                      Reps = 20,
                      {Ts, _} = timer:tc(fun() -> [?S:eval_slp(P, B) || _ <- lists:seq(1, Reps)] end),
                      {Tt, _} = timer:tc(fun() -> [?S:evaluate(U, B) || _ <- lists:seq(1, Reps)] end),
                      {ok, A1} = ?S:eval_slp(P, B), {ok, A2} = ?S:evaluate(U, B),
                      {ok, {cmat, R1}} = pq_core:unpack(A1), {ok, {cmat, R2}} = pq_core:unpack(A2),
                      ?assert(lists:all(fun({{complex, Ar, Ai}, {complex, Br, Bi}}) -> abs(Ar - Br) + abs(Ai - Bi) < 1.0e-12 end,
                                        lists:zip(lists:append(R1), lists:append(R2)))),
                      {slp, 1, _, Ops, _} = P,
                      io:format(user, "~nSS3 k = ~p: unitary ~.2f s, compile ~.2f s, ~p ops; eval ~.3f ms (slp) vs ~.3f ms (tree): ~.1fx~n",
                                [K, Tu / 1.0e6, Tc / 1.0e6, length(Ops), Ts / Reps / 1000, Tt / Reps / 1000, Tt / Ts])
                  end || K <- [4, 6]]
             end)]}.
