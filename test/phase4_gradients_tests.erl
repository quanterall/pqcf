%% Phase 4, item 3 (FS 15 item 3): adjoint gradients and parameter-shift rules (docs/decisions/AD1-adjoint-gradients.md).
%% Exit: against symbolic gradients (passive circuits) and central differences (Kerr, displacement, squeezing),
%% one forward and one backward pass.
-module(phase4_gradients_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

c(N, In, Is) -> {ok, C} = pq_circuit:new(N, In), C#circuit{instrs = Is}.
s(X) -> {sym, X}.

grad(C, Obs, B, Opts) -> {ok, G} = pq_gradient:gradient(C, Obs, B, Opts), G.

close(A, B, Tol) ->
    [?assertEqual({K, within}, {K, case abs(maps:get(K, A) - maps:get(K, B)) of D when D =< Tol -> within; D -> D end})
     || K <- maps:keys(B)],
    ?assertEqual(lists:sort(maps:keys(A)), lists:sort(maps:keys(B))).

%% three modes, two photons, every passive parameter symbolic
passive() ->
    c(3, {fock, [1, 1, 0]}, [{beamsplitter, 0, 1, s(a), s(b)}, {phaseshift, 1, s(c)},
                             {beamsplitter, 1, 2, s(d), {rat, 2, 5}}, {beamsplitter, 0, 1, {mul, [2, s(a)]}, s(e)},
                             {measure_pnr, [0, 1, 2]}]).

passive_binding() -> #{a => 0.7, b => -0.3, c => 1.1, d => 0.45, e => 0.25}.

passive_against_symbolic_test() ->
    C = passive(),
    B = passive_binding(),
    Obs = [{probability, [0, 1], [1, 1]}, {probability, [2], [1]}, {mean_photon_number, [1, 2]},
           {weighted_sum, [{2, {probability, [0], [2]}}, {{sym, w}, {mean_photon_number, [0]}}]}],
    [begin
         Cw = C#circuit{params = #{w => 0.25}},
         Gs = grad(Cw, O, B, #{mode => symbolic}),
         close(grad(Cw, O, B, #{mode => adjoint}), Gs, 1.0e-12),
         close(grad(Cw, O, B, #{mode => parameter_shift}), Gs, 1.0e-10)
     end || O <- Obs].

%% a weight that depends on a differentiated symbol contributes dW/ds times its observable
weight_depends_on_symbol_test() ->
    C = c(2, {fock, [1, 0]}, [{beamsplitter, 0, 1, s(t), 0}]),
    Obs = {weighted_sum, [{{mul, [s(t), s(t)]}, {probability, [0], [1]}}]},
    T = 0.6,
    Want = 2 * T * math:pow(math:cos(T), 2) - T * T * math:sin(2 * T),
    [?assert(abs(maps:get(t, grad(C, Obs, #{t => T}, #{mode => M})) - Want) < 1.0e-12) || M <- [adjoint, parameter_shift]].

%% a default that uses the differentiated symbol: q = 3 s, the gate uses q
defaults_chain_test() ->
    C = (c(2, {fock, [1, 1]}, [{beamsplitter, 0, 1, s(q), s(s)}, {phaseshift, 0, s(s)}, {beamsplitter, 0, 1, 0.3, 0}]))
        #circuit{params = #{q => {mul, [3, s(s)]}}},
    Obs = {probability, [0], [2]},
    Gn = grad(C, Obs, #{s => 0.2}, #{}),
    close(grad(C, Obs, #{s => 0.2}, #{mode => adjoint}), Gn, 1.0e-8),
    close(grad(C, Obs, #{s => 0.2}, #{mode => parameter_shift}), grad(C, Obs, #{s => 0.2}, #{mode => adjoint}), 1.0e-12),
    %% bound explicitly, q no longer follows s
    Gb = grad(C, Obs, #{s => 0.2, q => 0.5}, #{symbols => [s]}),
    close(grad(C, Obs, #{s => 0.2, q => 0.5}, #{symbols => [s], mode => adjoint}), Gb, 1.0e-8).

%% the truncated gates: Kerr, displacement and squeezing against central differences, at a cutoff
nonlinear() ->
    c(2, vacuum, [{displacement, 0, s(r), s(p)}, {squeezing, 1, s(z), s(q)}, {beamsplitter, 0, 1, s(t), 0.2},
                  {kerr, 1, s(k)}, {displacement, 1, 0.2, {mul, [s(p), s(t)]}}, {phaseshift, 0, s(f)}]).

nonlinear_binding() -> #{r => 0.4, p => 0.3, z => 0.25, q => -0.6, t => 0.8, k => 0.3, f => 0.5}.

truncated_against_differences_test() ->
    C = nonlinear(),
    B = nonlinear_binding(),
    Opts = #{cutoff => 9},
    [begin
         Gn = grad(C, O, B, Opts#{step => 1.0e-5}),
         Ga = grad(C, O, B, Opts#{mode => adjoint}),
         close(Ga, Gn, 1.0e-8)
     end || O <- [{probability, [0, 1], [1, 0]}, {mean_photon_number, [0, 1]}]].

%% the parameter-shift rules of the periodic parameters agree with the adjoint, Kerr included (R = P^2)
shift_against_adjoint_test() ->
    C = c(2, {fock, [1, 1]}, [{beamsplitter, 0, 1, s(t), s(p)}, {kerr, 0, s(k)}, {phaseshift, 1, s(f)},
                              {beamsplitter, 0, 1, 0.7, s(p)}]),
    B = #{t => 0.4, p => 0.9, k => 0.35, f => -0.2},
    Obs = {probability, [0], [1]},
    close(grad(C, Obs, B, #{mode => parameter_shift}), grad(C, Obs, B, #{mode => adjoint}), 1.0e-10),
    %% the same rule passed as a strategy
    close(grad(C, Obs, B, #{strategy => pq_gradient_shift}), grad(C, Obs, B, #{mode => parameter_shift}), 0.0),
    close(grad(C, Obs, B, #{strategy => pq_gradient_adjoint}), grad(C, Obs, B, #{mode => adjoint}), 0.0).

%% the angles of displacement and squeezing under a cutoff are periodic too
shift_phases_test() ->
    C = c(2, vacuum, [{displacement, 0, 0.5, s(p)}, {squeezing, 1, 0.3, s(q)}, {beamsplitter, 0, 1, 0.6, 0},
                      {kerr, 0, 0.2}]),
    B = #{p => 0.4, q => 1.2},
    Opts = #{cutoff => 5},
    Obs = {mean_photon_number, [1]},
    close(grad(C, Obs, B, Opts#{mode => parameter_shift}), grad(C, Obs, B, Opts#{mode => adjoint}), 1.0e-10).

%% parameter shift with loss (the density-matrix backend): exact against central differences
shift_with_loss_test() ->
    C = c(2, {fock, [2, 1]}, [{beamsplitter, 0, 1, s(t), 0.3}, {loss, 0, 0.7}, {phaseshift, 0, s(f)},
                              {beamsplitter, 0, 1, 0.5, 0}]),
    B = #{t => 0.35, f => 0.8},
    Obs = {probability, [0], [1]},
    close(grad(C, Obs, B, #{mode => parameter_shift}), grad(C, Obs, B, #{step => 1.0e-5}), 1.0e-9).

%% squeezed vacuum: d<n>/dr = sinh 2r, within the truncated mass at cutoff 30
squeezed_vacuum_test() ->
    C = c(1, vacuum, [{squeezing, 0, s(r), 0.4}]),
    #{r := G} = grad(C, {mean_photon_number, [0]}, #{r => 0.4}, #{mode => adjoint, cutoff => 30}),
    ?assert(abs(G - math:sinh(0.8)) < 1.0e-10).

%% bit-identical on repetition (D1: the kernel's order is fixed)
deterministic_test() ->
    C = nonlinear(),
    G = grad(C, {mean_photon_number, [0]}, nonlinear_binding(), #{mode => adjoint, cutoff => 7}),
    ?assertEqual(G, grad(C, {mean_photon_number, [0]}, nonlinear_binding(), #{mode => adjoint, cutoff => 7})).

%% one forward and one backward pass: with 24 parameters, faster than the 48 evaluations of central differences
cost_test_() ->
    {timeout, 120, fun() ->
        Is = lists:append([[{beamsplitter, M, M + 1, s(list_to_atom("t" ++ integer_to_list(K))),
                             s(list_to_atom("p" ++ integer_to_list(K)))} || M <- [0, 1, 2]] || K <- lists:seq(1, 4)]),
        C = c(4, {fock, [1, 1, 1, 0]}, Is),
        B = maps:from_list([{list_to_atom(L ++ integer_to_list(K)), 0.1 * K + 0.05 * J}
                            || K <- lists:seq(1, 4), {J, L} <- [{1, "t"}, {2, "p"}]]),
        ?assertEqual(8, map_size(B)),
        Obs = {probability, [3], [1]},
        {Ta, Ga} = timer:tc(fun() -> grad(C, Obs, B, #{mode => adjoint}) end),
        {Tn, Gn} = timer:tc(fun() -> grad(C, Obs, B, #{step => 1.0e-5}) end),
        close(Ga, Gn, 1.0e-8),
        ?debugFmt("adjoint ~p us, central differences ~p us (24 gate parameters, 8 symbols)", [Ta, Tn]),
        ?assert(Ta < Tn)
    end}.

refusals_test() ->
    Obs = {mean_photon_number, [0]},
    Lossy = c(2, {fock, [1, 0]}, [{beamsplitter, 0, 1, s(t), 0}, {loss, 0, 0.5}]),
    %% loss runs on density matrices, within max_dm_dim
    ?assertMatch({error, {dimension_limit, 3, 2}}, pq_gradient:gradient(Lossy, Obs, #{t => 0.2}, #{mode => adjoint, max_dm_dim => 2})),
    ?assertMatch({error, {differentiate_unsupported, #{mode := adjoint, backend := gaussian}}},
                 pq_gradient:gradient(Lossy, Obs, #{t => 0.2}, #{mode => adjoint, backend => gaussian})),
    %% the magnitude of a displacement is not periodic; nor is anything on the Gaussian backend
    D = c(1, vacuum, [{displacement, 0, s(r), 0}, {kerr, 0, 0.1}]),
    ?assertMatch({error, {differentiate_unsupported, #{mode := parameter_shift, reason := <<"not periodic">>}}},
                 pq_gradient:gradient(D, Obs, #{r => 0.2}, #{mode => parameter_shift, cutoff => 4})),
    G = c(1, vacuum, [{displacement, 0, 0.3, s(p)}]),
    ?assertMatch({error, {differentiate_unsupported, #{mode := parameter_shift, backend := gaussian}}},
                 pq_gradient:gradient(G, Obs, #{p => 0.2}, #{mode => parameter_shift})),
    %% a symbolic unitary in the differentiated symbols
    U = {matrix, [[{cos, s(t)}, {neg, {sin, s(t)}}], [{sin, s(t)}, {cos, s(t)}]]},
    Cu = c(2, {fock, [1, 0]}, [{unitary, [0, 1], U}]),
    ?assertMatch({error, {differentiate_unsupported, {symbolic_unitary, 0}}},
                 pq_gradient:gradient(Cu, Obs, #{t => 0.2}, #{mode => adjoint})),
    %% feed-forward branches on outcomes
    Ff = c(2, {fock, [1, 0]}, [{beamsplitter, 0, 1, s(t), 0}, {measure_pnr, [0]},
                               {feed_forward, {phaseshift, 1, s(m)}, #{m => 0}}]),
    ?assertMatch({error, {differentiate_unsupported, {instructions, [feed_forward]}}},
                 pq_gradient:gradient(Ff, Obs, #{t => 0.2}, #{mode => adjoint})).

%% density matrices (the follow-up of item 3): loss differentiated, thermal loss passed through
lossy_against_differences_test() ->
    C = c(3, {fock, [2, 1, 0]}, [{beamsplitter, 0, 1, s(t), s(p)}, {loss, 0, s(e)}, {kerr, 1, s(k)},
                                 {beamsplitter, 1, 2, 0.6, s(p)}, {loss, 2, {mul, [s(e), s(e)]}}, {phaseshift, 2, s(f)},
                                 {thermal_loss, 1, 0.9, 0.2}, {beamsplitter, 0, 1, 0.5, 0.1}]),
    B = #{t => 0.4, p => 0.3, e => 0.8, k => 0.25, f => 0.6},
    Opts = #{cutoff => 5},
    [begin
         Gn = grad(C, O, B, Opts),
         Ga = grad(C, O, B, Opts#{mode => adjoint}),
         close(Ga, Gn, 1.0e-8)
     end || O <- [{probability, [0, 1], [1, 1]}, {mean_photon_number, [1, 2]},
                  {weighted_sum, [{s(e), {probability, [2], [0]}}, {1, {mean_photon_number, [0]}}]}]],
    %% the periodic parameters against the exact shift rules (density-matrix backend)
    Bp = maps:remove(e, B),
    Cp = C#circuit{params = #{e => 0.8}},
    close(grad(Cp, {probability, [0], [1]}, Bp, Opts#{mode => adjoint}),
          grad(Cp, {probability, [0], [1]}, Bp, Opts#{mode => parameter_shift}), 1.0e-10).

%% pq_stats:mean_photon_number/2 on a density matrix (it had no clause): one photon through loss 0.3 keeps 0.3
mean_photon_number_dm_test() ->
    {ok, S, _} = pq_api:evolve(c(1, {fock, [1]}, [{loss, 0, 0.3}]), #{}),
    {ok, M} = pq_stats:mean_photon_number(S, [0]),
    ?assert(abs(M - 0.3) < 1.0e-15).

%% the kernel directly: a malformed gate list is a bad_shape, not a crash
kernel_shape_test() ->
    {ok, #fock_state{amps = A}} = pq_state:fock_state([1, 0], 1),
    ?assertMatch({error, {bad_shape, _}}, pq_nif:fock_adjoint_gradient(A, [{phaseshift, 0.1, 5}], [], 2, 1)),
    {ok, {[V], [[G]]}} = pq_nif:fock_adjoint_gradient(A, [{phaseshift, 0.1, 0}], [{photons, 1.0, [0]}], 2, 1),
    ?assertEqual({1.0, 0.0}, {V, abs(G)}).
