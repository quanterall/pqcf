%% Phase 4, item 5 (FS 15 item 5): GPU offload of the permanent (docs/decisions/G1-gpu-offload.md). Exit: GT-P01 on
%% the GPU path, results equal to the CPU kernel within the scale, reproducible bit for bit. Without a GPU (CI) the
%% GPU cases check the refusal instead: device => gpu is an error, never a silent CPU run.
-module(phase4_gpu_tests).
-include_lib("eunit/include/eunit.hrl").

-define(TOL, 1.0e-10).

fixture() ->
    {ok, Bin} = file:read_file(filename:join([code:lib_dir(pq), "..", "..", "..", "..", "test", "fixtures", "gt_p01.json"])),
    json:decode(Bin).

cmat(Rows) ->
    {ok, M} = pq_core:pack({cmat, [[{complex, Re, Im} || [Re, Im] <- Row] || Row <- Rows]}),
    M.

random(N, Seed) ->
    {ok, Us} = pq_nif:rng_u64(Seed, 0, 2 * N * N),
    F = [(U bsr 11) * math:pow(2, -53) - 0.5 || U <- Us],
    cmat(chunk([[Re, Im] || [Re, Im] <- pairs(F)], N)).

pairs([A, B | R]) -> [[A, B] | pairs(R)];
pairs([]) -> [].
chunk([], _) -> [];
chunk(L, N) -> {A, B} = lists:split(N, L), [A | chunk(B, N)].

gpu() ->
    case pq_nif:permanent(random(2, 1), #{device => gpu}) of
        {ok, _} -> true;
        {error, {unsupported_input, #{device := gpu, reason := Why}}} ->
            ?debugFmt("no GPU (~s): the GPU cases check the refusal", [Why]),
            false
    end.

gpu_test_() ->
    case gpu() of
        true -> [{"GT-P01 on the GPU", fun gt_p01/0}, {"against the CPU kernel", fun against_cpu/0},
                 {"auto and repeated rows", fun auto/0}, {timeout, 60, {"cancel", fun cancel/0}}];
        false -> [fun refused/0]
    end.

%% every permanent and probability of the golden file through the GPU path
gt_p01() ->
    #{<<"permanents">> := Ps, <<"probabilities">> := Pr} = fixture(),
    [begin
         {ok, {{Re, Im}, #{scale := S, device := gpu}}} = pq_nif:permanent(cmat(M), #{device => gpu}),
         ?assert(abs(Re - WRe) =< ?TOL * max(S, 1.0) andalso abs(Im - WIm) =< ?TOL * max(S, 1.0))
     end || #{<<"matrix">> := M, <<"perm">> := [WRe, WIm]} <- Ps],
    [begin
         {ok, {P, #{permanent_scale := S}}} = pq_nif:bs_probability(cmat(U), In, Out, #{device => gpu}),
         ?assert(abs(P - W) =< ?TOL * max(S * S, 1.0e-6))
     end || #{<<"unitary">> := U, <<"input">> := In, <<"output">> := Out, <<"probability">> := W} <- Pr].

%% n = 2..22: within 8 eps S of the CPU kernel, the scale too; bit-identical on repetition
against_cpu() ->
    [begin
         M = random(N, 1000 + N),
         {ok, {{Cr, Ci}, #{scale := S, device := cpu}}} = pq_nif:permanent(M, #{threads => 4}),
         {ok, {{Gr, Gi}, #{scale := Sg, device := gpu}} = G} = pq_nif:permanent(M, #{device => gpu}),
         Eps = math:pow(2, -52),
         ?assert(abs(Gr - Cr) =< 8 * Eps * S andalso abs(Gi - Ci) =< 8 * Eps * S),
         ?assert(abs(Sg - S) =< 8 * Eps * S),
         ?assertEqual({ok, G}, {ok, element(2, pq_nif:permanent(M, #{device => gpu}))})
     end || N <- [2, 3, 7, 12, 13, 16, 20, 22]].

%% auto: the CPU below n = 24, the GPU from it; repeated multiplicities stay on the CPU
auto() ->
    ?assertMatch({ok, {_, #{device := cpu}}}, pq_nif:permanent(random(10, 3), #{device => auto})),
    ?assertMatch({ok, {_, #{device := gpu}}}, pq_nif:permanent(random(24, 3), #{device => auto})),
    U = random(4, 5),
    ?assertMatch({ok, {_, #{device := cpu}}}, pq_nif:permanent_repeated(U, [2, 0, 1, 1], [1, 1, 1, 1], #{device => gpu})),
    {ok, {Vg, #{device := gpu}}} = pq_nif:permanent_repeated(U, [1, 0, 1, 1], [1, 1, 0, 1], #{device => gpu}),
    {ok, {Vc, #{device := cpu}}} = pq_nif:permanent_repeated(U, [1, 0, 1, 1], [1, 1, 0, 1], #{}),
    ?assert(abs(element(1, Vg) - element(1, Vc)) < 1.0e-15 andalso abs(element(2, Vg) - element(2, Vc)) < 1.0e-15),
    %% through the facade: boson probabilities with device => gpu equal the CPU's
    {ok, C0} = pq_circuit:new(4, {fock, [1, 1, 1, 0]}),
    C = lists:foldl(fun(I, A) -> pq_circuit:add(A, I) end, C0,
                    [{beamsplitter, 0, 1, 0.7, 0.2}, {beamsplitter, 1, 2, 0.5, -0.4}, {beamsplitter, 2, 3, 0.9, 0.1},
                     {measure_pnr, [0, 1, 2, 3]}]),
    {ok, Pc, _} = pq_api:probabilities(C, #{backend => boson}),
    {ok, Pg, _} = pq_api:probabilities(C, #{backend => boson, device => gpu}),
    [?assert(abs(V - maps:get(K, Pg)) < 1.0e-15) || K := V <- Pc].

%% a set cancel token stops the GPU path between launches
cancel() ->
    T = pq_nif:cancel_token_new(),
    ok = pq_nif:cancel_token_set(T),
    ?assertEqual({error, cancelled}, pq_nif:permanent(random(30, 9), #{device => gpu, cancel => T})).

refused() ->
    ?assertMatch({error, {unsupported_input, #{device := gpu}}}, pq_nif:permanent(random(5, 1), #{device => gpu})),
    %% auto falls back to the CPU
    ?assertMatch({ok, {_, #{device := cpu}}}, pq_nif:permanent(random(24, 1), #{device => auto})).
