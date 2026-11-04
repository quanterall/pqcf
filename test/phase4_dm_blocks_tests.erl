%% Phase 4, item 6 (FS 15 item 6): density matrices stored by photon-number blocks inside the backend
%% (docs/decisions/B1-dm-blocks.md). Exit: every density-matrix result unchanged -- here bit for bit against the
%% dense storage -- and the memory at cutoff 20 on 4 modes.
-module(phase4_dm_blocks_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

c(N, In, Is) -> {ok, C} = pq_circuit:new(N, In), C#circuit{instrs = Is}.

blocks() -> #{dm_storage => blocks}.

setup_test_() ->
    {setup, fun() -> {ok, _} = application:ensure_all_started(pq) end, fun(_) -> ok end,
     [fun same_probabilities/0, fun same_states/0, fun same_samples/0, {timeout, 600, fun memory/0}, fun limits/0]}.

circuits() ->
    [%% number-conserving: Fock input, passive gates, Kerr, loss
     c(3, {fock, [2, 1, 1]}, [{beamsplitter, 0, 1, 0.7, 0.3}, {kerr, 1, 0.4}, {loss, 0, 0.6}, {beamsplitter, 1, 2, 1.1, -0.2},
                              {loss, 2, 0.8}, {phaseshift, 2, 0.3}, {measure_pnr, [0, 1, 2]}]),
     %% displacement, squeezing and thermal loss fill the off-diagonal blocks
     c(2, vacuum, [{displacement, 0, 0.4, 0.2}, {squeezing, 1, 0.3, 0.1}, {beamsplitter, 0, 1, 0.6, 0}, {loss, 1, 0.7},
                   {thermal_loss, 0, 0.8, 0.2}, {kerr, 0, 0.3}, {measure_pnr, [0, 1]}]),
     %% a pure state promoted at the first loss
     c(2, {fock, [1, 2]}, [{beamsplitter, 0, 1, 0.5, 0.1}, {unitary, [0, 1], mzi()}, {loss, 0, 0.5}, {measure_threshold, [0, 1]}])].

mzi() ->
    {ok, U} = pq_core:pack({cmat, [[{complex, 0.6, 0.0}, {complex, 0.0, 0.8}], [{complex, 0.0, 0.8}, {complex, 0.6, 0.0}]]}),
    U.

opts(#circuit{input = vacuum}) -> #{cutoff => 8, backend => fock_dm};
opts(_) -> #{backend => fock_dm}.

%% probabilities identical, float for float
same_probabilities() ->
    [begin
         O = opts(C),
         {ok, Pd, _} = pq_api:probabilities(C, O),
         {ok, Pb, _} = pq_api:probabilities(C, maps:merge(O, blocks())),
         ?assertEqual(Pd, Pb)
     end || C <- circuits()].

%% the public state (evolve) is the full matrix, equal to the dense storage's entry for entry; the only bit
%% differences are signed zeros (the dense kernel leaves -0.0 in entries that the blocks never store)
same_states() ->
    [begin
         O = opts(C),
         {ok, #fock_dm{rho = {cmat, _, _, _} = Rd} = Sd, Wd} = pq_api:evolve(C, O),
         {ok, #fock_dm{rho = {cmat, _, _, _} = Rb} = Sb, Wb} = pq_api:evolve(C, maps:merge(O, blocks())),
         ?assertEqual(Sd#fock_dm{rho = x}, Sb#fock_dm{rho = x}),
         {ok, {cmat, Ed}} = pq_core:unpack(Rd),
         {ok, {cmat, Eb}} = pq_core:unpack(Rb),
         Pairs = lists:zip(lists:append(Ed), lists:append(Eb)),
         ?assertEqual([], [P || {{complex, A, B}, {complex, X, Y}} = P <- Pairs, not (A == X andalso B == Y)]),
         ?assertEqual([], [P || {{complex, A, B}, {complex, X, Y}} = P <- Pairs, A =/= X orelse B =/= Y,
                                not (A == 0 andalso B == 0 andalso X == 0 andalso Y == 0)]),
         ?assertEqual(Wd, Wb)
     end || C <- circuits()].

%% mid-circuit measurements and terminal sampling draw the same samples
same_samples() ->
    C = c(3, {fock, [2, 1, 0]}, [{beamsplitter, 0, 1, 0.7, 0.3}, {loss, 0, 0.6}, {measure_pnr, [0]}, {beamsplitter, 1, 2, 0.9, 0},
                                 {loss, 2, 0.7}, {measure_pnr, [1, 2]}]),
    O = #{shots => 2000, seed => 13},
    {ok, #result{samples = Sd}} = pq_api:sample(C, O),
    {ok, #result{samples = Sb}} = pq_api:sample(C, maps:merge(O, blocks())),
    ?assertEqual(Sd, Sb),
    %% with promote_dm on the pure-state backend too
    {ok, #result{samples = Sp}} = pq_api:sample(C, maps:merge(O#{backend => fock, promote_dm => true}, blocks())),
    ?assertEqual(Sd, Sp).

%% cutoff 20 on 4 modes: the dense storage refuses (D = 10626 > max_dm_dim), the blocks hold the 11,480,634 entries
%% of the diagonal blocks (0.18 GB instead of 1.81 GB)
memory() ->
    C = c(4, {fock, [5, 5, 5, 5]}, [{beamsplitter, 0, 1, 0.7, 0.1}, {loss, 0, 0.8}, {beamsplitter, 2, 3, 0.4, 0.0},
                                    {loss, 3, 0.9}, {measure_pnr, [0]}]),
    ?assertMatch({error, {dimension_limit, 10626, 4096}}, pq_api:probabilities(C, #{})),
    {T, {ok, P, _}} = timer:tc(fun() -> pq_api:probabilities(C, blocks()) end),
    ?assert(abs(lists:sum(maps:values(P)) - 1) < 1.0e-12),
    %% the mean photon number of mode 0 against the closed form: 0.8 times the beamsplitter's mean
    Mean = lists:sum([V * N || [N] := V <- P]),
    Want = 0.8 * (5 * math:pow(math:cos(0.7), 2) + 5 * math:pow(math:sin(0.7), 2)),
    ?assert(abs(Mean - Want) < 1.0e-10),
    {ok, {S, _}, _} = pq_simulator:evolve_state(pq_fock_dm, C#circuit{instrs = lists:droplast(C#circuit.instrs)}, blocks()),
    #fock_dm{rho = {dm_blocks, B}} = S,
    Entries = lists:sum([byte_size(X) div 16 || {_, _, X} <- B]),
    ?debugFmt("cutoff 20, 4 modes: ~p stored entries (~.2f GB; dense ~.2f GB), ~p blocks, ~.1f s",
              [Entries, Entries * 16 / 1.0e9, 10626 * 10626 * 16 / 1.0e9, length(B), T / 1.0e6]),
    ?assert(Entries =< 11480634),
    ?assert(lists:all(fun({A, Cc, _}) -> A =:= Cc end, B)).

limits() ->
    C = hd(circuits()),
    ?assertMatch({error, {dimension_limit, #{max_dm_entries := 100}}},
                 pq_api:probabilities(C, maps:merge(opts(C), (blocks())#{max_dm_entries => 100}))),
    %% evolve presents the full matrix, so it stays within max_dm_dim
    Big = c(4, {fock, [5, 5, 5, 5]}, [{loss, 0, 0.8}]),
    ?assertMatch({error, {dimension_limit, 10626, 4096}}, pq_api:evolve(Big, blocks())).
