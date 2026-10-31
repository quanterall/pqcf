%% FR-V03, Decision 6: numeric and symbolic gradients agree with each other and with the analytic derivatives.
-module(pq_gradient_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

c(N, In, Is) -> {ok, C} = pq_circuit:new(N, In), C#circuit{instrs = Is}.

mzi_test() ->
    C = c(2, {fock, [1, 0]}, pq_circuit:mzi(0, 1, {sym, th}, {rat, 3, 10})),   % exact: floats are contagious
    Obs = {probability, [0], [1]},
    Th = 0.8,
    Want = math:sin(Th) / 2,                         % d/dTheta (1 - cos Theta)/2
    {ok, #{th := Gn}} = pq_gradient:gradient(C, Obs, #{th => Th}, #{}),
    {ok, #{th := Gs}} = pq_gradient:gradient(C, Obs, #{th => Th}, #{mode => symbolic}),
    ?assert(abs(Gn - Want) < 1.0e-6),
    ?assert(abs(Gs - Want) < 1.0e-13),
    {ok, #{th := E}} = pq_gradient:symbolic_gradient(C, Obs, [th]),
    ?assertEqual(pq_symbolic:canonicalise({mul, [{rat, 1, 2}, {sin, {sym, th}}]}), {ok, E}).

hom_test() ->
    %% SY 10.1: p(1,1) = (cos^2 - sin^2)^2, derivative -2 sin(4 theta)
    C = c(2, {fock, [1, 1]}, [{beamsplitter, 0, 1, {sym, t}, 0}]),
    T = 0.3,
    {ok, #{t := Gs}} = pq_gradient:gradient(C, {probability, [0, 1], [1, 1]}, #{t => T}, #{mode => symbolic}),
    {ok, #{t := Gn}} = pq_gradient:gradient(C, {probability, [0, 1], [1, 1]}, #{t => T}, #{}),
    ?assert(abs(Gs + 2 * math:sin(4 * T)) < 1.0e-13 andalso abs(Gn - Gs) < 1.0e-6),
    %% mean photon number in mode 1: 2 sin^2... for |1,1>: <n_1> = 1 always, derivative 0
    {ok, #{t := Gm}} = pq_gradient:gradient(C, {mean_photon_number, [1]}, #{t => T}, #{mode => symbolic}),
    ?assert(abs(Gm) < 1.0e-14),
    %% a weighted sum of two observables
    W = {weighted_sum, [{2, {probability, [0, 1], [1, 1]}}, {{sym, c}, {probability, [0], [2]}}]},
    {ok, Gw} = pq_gradient:gradient(C#circuit{params = #{c => 3}}, W, #{t => T}, #{mode => symbolic}),
    {ok, Gwn} = pq_gradient:gradient(C#circuit{params = #{c => 3}}, W, #{t => T}, #{}),
    ?assert(abs(maps:get(t, Gw) - maps:get(t, Gwn)) < 1.0e-6).

gaussian_test() ->
    %% <n> = sinh^2 r for squeezed vacuum: d/dr = sinh(2r); numeric only (symbolic Gaussian evolution, FR-G04, is
    %% not built)
    C = c(1, vacuum, [{squeezing, 0, {sym, r}, 0}]),
    {ok, #{r := G}} = pq_gradient:gradient(C, {mean_photon_number, [0]}, #{r => 0.4}, #{}),
    ?assert(abs(G - math:sinh(0.8)) < 1.0e-6),
    ?assertMatch({error, {differentiate_unsupported, _}},
                 pq_gradient:gradient(C, {mean_photon_number, [0]}, #{r => 0.4}, #{mode => symbolic})).
