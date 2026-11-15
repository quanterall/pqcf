%% pq_dicke (quantum optics plan, Q5): N two-level atoms in the Dicke basis, against an exact solution of the
%% superradiant burst (mpmath, 40 digits: test/fixtures/superradiance_mp.json from comparison/dicke/superradiance_mp.py)
%% and against QuTiP's PIQS (test/fixtures/piqs_dicke.json from comparison/dicke/piqs_check.py).
-module(pq_dicke_tests).
-include_lib("eunit/include/eunit.hrl").

-define(D, pq_dicke).

fixture(Name) ->
    Dir = filename:dirname(code:which(?MODULE)),
    Paths = [filename:join([Dir, "..", "test", "fixtures", Name]),
             filename:join([Dir, "..", "..", "..", "..", "test", "fixtures", Name]), filename:join("test/fixtures", Name)],
    [P | _] = [P || P <- Paths, filelib:is_regular(P)],
    {ok, B} = file:read_file(P),
    json:decode(B).

%% "0.046" as {46, 1000}
dec(B) ->
    case binary:split(B, <<".">>) of
        [I, F] -> {binary_to_integer(<<I/binary, F/binary>>), round(math:pow(10, byte_size(F)))};
        [I] -> {binary_to_integer(I), 1}
    end.
rat(X) -> {round(X * 1000), 1000}.

%% the superradiant burst of 100 atoms (collective emission, all excited) equals the exact solution to 30 digits
superradiance_exact_test_() ->
    {timeout, 120, fun() ->
        #{<<"t">> := Ts, <<"jpjm">> := Want} = fixture("superradiance_mp.json"),
        P = 160,
        Tr = ?D:triples(100, diagonal),
        Sol = pq_lindblad:evolve(?D:liouvillian(100, #{collective_emission => 1}, Tr, P), ?D:state(Tr, {excited, 100}, P),
                                 [dec(T) || T <- Ts], P),
        [begin
             W = pq_mp:from_string(binary_to_list(Wb), P),
             D = pq_mp:sub(?D:expect(jp_jm, Tr, R, P), W, P),
             ?assert(element(2, D) =:= 0 orelse abs(pq_mp:to_float(D)) < 1.0e-30 * pq_mp:to_float(W))
         end || {{_, R}, Wb} <- lists:zip(Sol, Want)]
    end}.

%% the burst peaks near Gross and Haroche's t_D = ln N/(N gamma) (0.046 for N = 100): on PIQS's grid of 0.002, the
%% exact <J+ J-> is largest at 0.048, as in PIQS
superradiance_peak_test_() ->
    {timeout, 120, fun() ->
        #{<<"C">> := #{<<"t">> := Ts, <<"jpjm">> := Piqs}} = fixture("piqs_dicke.json"),
        P = 64,
        Tr = ?D:triples(100, diagonal),
        Sol = pq_lindblad:evolve(?D:liouvillian(100, #{collective_emission => 1}, Tr, P), ?D:state(Tr, {excited, 100}, P),
                                 [rat(T) || T <- Ts], P),
        Vals = [pq_mp:to_float(?D:expect(jp_jm, Tr, R, P)) || {_, R} <- Sol],
        Peak = fun(L) -> element(2, lists:max(lists:zip(L, Ts))) end,
        ?assertEqual(Peak(Piqs), Peak(Vals)),
        ?assert(abs(Peak(Vals) - math:log(100) / 100) < 0.004)
    end}.

%% against PIQS on its three problems: every process at N = 6; the driven N = 5 (half-integer ladders, J_x mixing
%% m - m'); the burst at N = 100. PIQS computes its rates in single precision (Cython float), so the agreement is
%% about 1e-8; the test allows 1e-6 relative
piqs_test_() ->
    {timeout, 300, fun() ->
        F = fixture("piqs_dicke.json"),
        P = 80,
        [begin
             #{<<"N">> := N, <<"rates">> := Rs, <<"t">> := Ts, <<"jz">> := Jz, <<"jpjm">> := Jj} = C = maps:get(K, F),
             Rates = maps:merge(maps:from_list([{binary_to_atom(Kk), rat(V)} || Kk := V <- Rs]),
                                maps:from_list([{A, rat(maps:get(Kb, C))} || {Kb, A} <- [{<<"omega">>, omega}, {<<"omega_x">>, omega_x}],
                                                                            is_map_key(Kb, C)])),
             Tr = ?D:triples(N, case is_map_key(omega_x, Rates) of true -> all; false -> diagonal end),
             Init = case K of <<"B">> -> {ground, N}; _ -> {excited, N} end,
             Sol = pq_lindblad:evolve(?D:liouvillian(N, Rates, Tr, P), ?D:state(Tr, Init, P), [rat(T) || T <- Ts], P),
             [begin
                  ?assert(abs(pq_mp:to_float(?D:expect(jz, Tr, R, P)) - Z) < 1.0e-6 * max(1.0, abs(Z))),
                  ?assert(abs(pq_mp:to_float(?D:expect(jp_jm, Tr, R, P)) - J) < 1.0e-6 * max(1.0, abs(J)))
              end || {{_, R}, Z, J} <- lists:zip3(Sol, Jz, Jj)]
         end || K <- [<<"A">>, <<"B">>, <<"C">>]]
    end}.

%% the trace is conserved by every process (to the working precision)
trace_test() ->
    P = 100,
    Rates = #{emission => {1, 2}, dephasing => {3, 10}, pumping => {1, 5}, collective_emission => 1,
              collective_dephasing => {2, 5}, collective_pumping => {1, 10}},
    Tr = ?D:triples(7, diagonal),
    [{_, R}] = pq_lindblad:evolve(?D:liouvillian(7, Rates, Tr, P), ?D:state(Tr, {excited, 7}, P), [{3, 2}], P),
    ?assert(abs(pq_mp:to_float(pq_mp:sub(?D:expect(trace, Tr, R, P), pq_mp:from(1, P), P))) < 1.0e-25).
