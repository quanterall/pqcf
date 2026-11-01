%% NFR-06 (unitarity and symplecticity errors at most 1e-10 after 1,000 gates) and NFR-09 (ping latency on normal
%% schedulers below 5 ms at p99 under a saturating dirty-scheduler workload). The thresholds are the specification's;
%% NFR-09 is measured on the reference hardware as well (bench/), this test keeps it from regressing.
-module(nfr_tests).
-include_lib("eunit/include/eunit.hrl").
-include("../include/pq.hrl").

%% a deterministic stream of gates (no rand: the same 1,000 gates every run)
gates(N, Count, Kinds) ->
    {Gs, _} = lists:mapfoldl(fun(K, X) ->
                                     X1 = (X * 6364136223846793005 + 1442695040888963407) band (1 bsl 64 - 1),
                                     U = fun(Sh) -> ((X1 bsr Sh) band 16#FFFF) / 16#FFFF end,
                                     A = trunc(U(0) * (N - 1)),
                                     Kind = lists:nth(K rem length(Kinds) + 1, Kinds),
                                     G = case Kind of
                                             bs -> {beamsplitter, A, A + 1, 3.2 * U(16), 6.3 * U(32) - 3.15};
                                             ps -> {phaseshift, A, 6.3 * U(16)};
                                             kerr -> {kerr, A, 2 * U(16)};
                                             sq -> {squeezing, A, 0.1 * U(16) - 0.05, 6.3 * U(32)};
                                             disp -> {displacement, A, 0.2 * U(16), 6.3 * U(32)}
                                         end,
                                     {G, X1}
                             end, 12345, lists:seq(1, Count)),
    Gs.

c(N, In, Is) -> {ok, C} = pq_circuit:new(N, In), C#circuit{instrs = Is}.

unitarity_test() ->
    N = 6,
    {ok, L} = pq_circuit:lower(c(N, vacuum, gates(N, 1000, [bs, ps, bs])), #{}),
    {ok, {cmat, U}} = pq_core:unpack(pq_sampling:passive_unitary(L)),
    Err = lists:max([abs(dot(Ri, Rj) - case I =:= J of true -> 1; false -> 0 end)
                     || {I, Ri} <- lists:enumerate(U), {J, Rj} <- lists:enumerate(U)]),
    io:format(user, "~nNFR-06: |U U^dagger - I| = ~.3g after 1,000 gates~n", [Err]),
    ?assert(Err =< 1.0e-10).

%% |sum_k a_k conj(b_k)|
dot(A, B) ->
    {Re, Im} = lists:foldl(fun({{complex, Ar, Ai}, {complex, Br, Bi}}, {R, I}) -> {R + Ar * Br + Ai * Bi, I + Ai * Br - Ar * Bi} end,
                           {0.0, 0.0}, lists:zip(A, B)),
    math:sqrt(Re * Re + Im * Im).

fock_norm_test() ->
    {ok, #fock_state{amps = A}, _} = pq_api:evolve(c(3, {fock, [1, 1, 1]}, gates(3, 1000, [bs, ps, kerr])), #{}),
    {ok, {cvec, V}} = pq_core:unpack(A),
    Err = abs(lists:sum([X * X + Y * Y || {complex, X, Y} <- V]) - 1),
    io:format(user, "NFR-06: Fock norm error ~.3g after 1,000 gates~n", [Err]),
    ?assert(Err =< 1.0e-10).

%% a pure Gaussian state stays pure under symplectic maps: det(Sigma) = (hbar/2)^(2n) = 1 at hbar = 2
gaussian_symplectic_test() ->
    {ok, S, _} = pq_api:evolve(c(3, vacuum, gates(3, 1000, [bs, sq, ps, disp, bs])), #{}),
    {_, Cov} = pq_state:mean_cov(S),
    Err = abs(det(Cov) - 1),
    io:format(user, "NFR-06: |det(Sigma) - 1| = ~.3g after 1,000 gates~n", [Err]),
    ?assert(Err =< 1.0e-10),
    ?assertEqual(ok, pq_state:check_physical(S)).

%% Gaussian elimination with partial pivoting
det(M) -> det(M, 1.0).
det([], D) -> D;
det(Rows, D) ->
    %% the row with the largest first entry swapped to the top: one transposition, or none
    {_, K} = lists:max([{abs(hd(R)), I} || {I, R} <- lists:enumerate(Rows)]),
    {Before, [P | After]} = lists:split(K - 1, Rows),
    Rest = case Before of [] -> After; [H | T] -> T ++ [H | After] end,
    Swaps = case K of 1 -> 1; _ -> -1 end,
    [Pv | Pt] = P,
    Reduced = [[X - (F / Pv) * Y || {X, Y} <- lists:zip(Rt, Pt)] || [F | Rt] <- Rest],
    det(Reduced, D * Pv * Swaps).

%% NFR-09: every dirty CPU scheduler busy with long boson-sampling chunks; a process on the normal schedulers
%% answers pings within 5 ms at p99
responsiveness_test_() ->
    {setup, fun() -> {ok, _} = application:ensure_all_started(pq) end, fun(_) -> ok end,
     {timeout, 120, fun responsiveness/0}}.

responsiveness() ->
    Shared = os:getenv("PQ_SHARED_RUNNER") =:= "1",
    Dirty = erlang:system_info(dirty_cpu_schedulers_online),
    {ok, C0} = pq_circuit:new(12, {fock, lists:duplicate(8, 1) ++ lists:duplicate(4, 0)}),
    Big = C0#circuit{instrs = lists:append([[{beamsplitter, J, J + 1, 0.3 + 0.05 * J, 0.1 * L} || J <- lists:seq(L rem 2, 10, 2)]
                                            || L <- lists:seq(0, 11)]) ++ [{measure_pnr, lists:seq(0, 11)}]},
    Self = self(),
    spawn(fun() -> Self ! {done, pq_api:run(boson, Big, #{shots => 200 * Dirty * 20000, seed => 1, chunk => 20000,
                                                          max_inflight => Dirty, notify => Self})} end),
    Job = receive {pq_job, J} -> J after 10000 -> error(no_job) end,
    timer:sleep(200),
    {ok, #{state := running, inflight := Inflight}} = pq_parallel:status(Job),
    ?assertEqual(Dirty, Inflight),
    Echo = spawn(fun Loop() -> receive {ping, From, T} -> From ! {pong, T}, Loop(); stop -> ok end end),
    Tp0 = erlang:monotonic_time(millisecond),
    Rtts = [begin
                timer:sleep(1),
                T0 = erlang:monotonic_time(microsecond),
                Echo ! {ping, self(), T0},
                receive {pong, T0} -> erlang:monotonic_time(microsecond) - T0 end
            end || _ <- lists:seq(1, 1000)],
    Echo ! stop,
    %% still saturated at the end, or the measurement proves nothing; if not, say how the job ended
    T1 = erlang:monotonic_time(millisecond),
    Saturated = case pq_parallel:status(Job) of
                    {ok, #{state := running}} -> true;
                    _ when Shared -> io:format(user, "~nNFR-09 (shared runner): the job ended before the pings "
                                                     "(~p ms for 1,000 pings)~n", [T1 - Tp0]),
                                     false;
                    Other -> error({job_ended_before_the_pings, Other, receive {done, R} -> R after 1000 -> no_reply end,
                                    #{pings_ms => T1 - Tp0}})
                end,
    ok = pq_parallel:cancel(Job),
    receive {done, _} -> ok after 10000 -> ok end,
    Sorted = lists:sort(Rtts),
    P99 = lists:nth(990, Sorted),
    io:format(user, "~nNFR-09: ~p dirty schedulers busy; ping p50 ~p us, p99 ~p us, max ~p us~n",
              [Dirty, lists:nth(500, Sorted), P99, lists:last(Sorted)]),
    %% NFR-09 is a requirement on the reference hardware: on a shared CI runner (PQ_SHARED_RUNNER=1, a 4-vCPU
    %% virtual machine whose dirty schedulers take every vCPU) it is measured and reported, not enforced
    case Shared of
        true -> io:format(user, "NFR-09 (shared runner, not enforced): p99 ~p us, saturated to the end: ~p~n",
                          [P99, Saturated]);
        false -> ?assert(P99 < 5000)
    end.
