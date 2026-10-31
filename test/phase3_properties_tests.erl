%% The property tests of FS 12.3 not covered elsewhere (the symbolic ones are in pq_symbolic_tests, the codec round
%% trip in pq_core_tests): composition associativity, lowering commutes with evaluation, Gaussian maps preserve
%% physicality, unitary Fock gates preserve the norm, loss preserves the trace, marginal probabilities sum to 1.
-module(phase3_properties_tests).
-include_lib("proper/include/proper.hrl").
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

-define(N, 3).

angle() -> ?LET(X, float(-3.2, 3.2), X).
mode() -> integer(0, ?N - 1).
pair() -> ?SUCHTHAT({A, B}, {mode(), mode()}, A =/= B).

passive() ->
    oneof([?LET({{A, B}, T, P}, {pair(), angle(), angle()}, {beamsplitter, A, B, T, P}),
           ?LET({M, P}, {mode(), angle()}, {phaseshift, M, P})]).

fock_gate() -> oneof([passive(), ?LET({M, K}, {mode(), angle()}, {kerr, M, K})]).

gaussian_gate() ->
    oneof([passive(),
           ?LET({M, R, P}, {mode(), float(-0.8, 0.8), angle()}, {squeezing, M, R, P}),
           ?LET({M, R, P}, {mode(), float(0.0, 1.0), angle()}, {displacement, M, R, P}),
           ?LET({M, E}, {mode(), float(0.0, 1.0)}, {loss, M, E})]).

fock_input() -> ?SUCHTHAT(In, vector(?N, integer(0, 2)), lists:sum(In) > 0 andalso lists:sum(In) =< 3).

circuit(In, Is) -> {ok, C} = pq_circuit:new(?N, In), C#circuit{instrs = Is}.

cvec(Packed) -> {ok, {cvec, L}} = pq_core:unpack(Packed), L.

norm2(L) -> lists:sum([Re * Re + Im * Im || {complex, Re, Im} <- L]).

%% compose(compose(A, B), C) = compose(A, compose(B, C))
prop_compose_associative() ->
    ?FORALL({A, B, C}, {list(passive()), list(passive()), list(passive())},
            begin
                [Ca, Cb, Cc] = [circuit(vacuum, X) || X <- [A, B, C]],
                {ok, AB} = pq_circuit:compose(Ca, Cb),
                {ok, BC} = pq_circuit:compose(Cb, Cc),
                {ok, L} = pq_circuit:compose(AB, Cc),
                {ok, R} = pq_circuit:compose(Ca, BC),
                L =:= R
            end).

%% the symbolic unitary evaluated at a binding equals the unitary of the circuit lowered at that binding
prop_lowering_commutes() ->
    %% at most 6 gates: longer ones reach max_expr_nodes (symbolic_limit, by design)
    ?FORALL({Is, Vals}, {?SUCHTHAT(L, non_empty(list(passive())), length(L) =< 6), vector(4, angle())},
            begin
                Syms = [a, b, c, d],
                %% every angle parameter becomes a + k (a symbol, so lowering has work to do)
                Sym = fun(K) -> {add, [{sym, lists:nth(K rem 4 + 1, Syms)}, K / 7]} end,
                {Si, _} = lists:mapfoldl(fun({beamsplitter, A, B, _, _}, K) -> {{beamsplitter, A, B, Sym(K), Sym(K + 1)}, K + 2};
                                            ({phaseshift, M, _}, K) -> {{phaseshift, M, Sym(K)}, K + 1}
                                         end, 0, Is),
                C = circuit(vacuum, Si),
                B = maps:from_list(lists:zip(Syms, Vals)),
                {ok, U} = pq_symbolic:circuit_unitary(C),
                {ok, Us} = pq_symbolic:evaluate(U, B),
                {ok, L} = pq_circuit:lower(C, B),
                Un = pq_sampling:passive_unitary(L),
                {ok, {cmat, X}} = pq_core:unpack(Us),
                {ok, {cmat, Y}} = pq_core:unpack(Un),
                lists:max([abs(Ra - Rb) + abs(Ia - Ib) || {{complex, Ra, Ia}, {complex, Rb, Ib}}
                                                             <- lists:zip(lists:append(X), lists:append(Y))]) < 1.0e-12
            end).

%% every Gaussian instruction maps a physical state to a physical state (sigma + i (hbar/2) Omega >= 0)
prop_gaussian_physical() ->
    ?FORALL(Is, list(gaussian_gate()),
            begin
                {ok, S, _} = pq_api:evolve(circuit(vacuum, Is), #{backend => gaussian}),
                pq_state:check_physical(S) =:= ok
            end).

%% beamsplitters, phase shifts and Kerr gates keep the norm of a Fock state (no truncation: they conserve photons)
prop_fock_norm() ->
    ?FORALL({In, Is}, {fock_input(), list(fock_gate())},
            begin
                {ok, #fock_state{amps = A}, _} = pq_api:evolve(circuit({fock, In}, Is), #{backend => fock}),
                abs(norm2(cvec(A)) - 1) < 1.0e-12
            end).

%% loss keeps the trace of a density matrix
prop_loss_trace() ->
    ?FORALL({In, Is}, {fock_input(), list(oneof([fock_gate(), ?LET({M, E}, {mode(), float(0.0, 1.0)}, {loss, M, E})]))},
            begin
                {ok, #fock_dm{rho = R}, _} = pq_api:evolve(circuit({fock, In}, Is), #{backend => fock_dm}),
                {ok, {cmat, Rows}} = pq_core:unpack(R),
                Tr = lists:sum([Re || {K, Row} <- lists:enumerate(1, Rows), {complex, Re, _} <- [lists:nth(K, Row)]]),
                abs(Tr - 1) < 1.0e-12
            end).

%% the probabilities of any set of measured modes sum to 1 (pure states, density matrices, permanents)
prop_marginals_sum_to_one() ->
    ?FORALL({In, Is, Modes, Backend}, {fock_input(), list(passive()), non_empty(list(mode())), oneof([fock, fock_dm, boson])},
            begin
                Ms = lists:usort(Modes),
                {ok, P, _} = pq_api:probabilities(circuit({fock, In}, Is ++ [{measure_pnr, Ms}]), #{backend => Backend}),
                abs(lists:sum(maps:values(P)) - 1) < 1.0e-12
            end).

properties_test_() ->
    Opts = [{numtests, 200}, {to_file, user}],
    [{timeout, 300, ?_assert(proper:quickcheck(P, Opts))} || P <- [prop_compose_associative(), prop_lowering_commutes(),
                                                    prop_gaussian_physical(), prop_fock_norm(), prop_loss_trace(),
                                                    prop_marginals_sum_to_one()]].
