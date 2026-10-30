%% The hafnian family and the torontonian through the NIF (FS 8.4); their values are checked against thewalrus
%% in the Rust tests (gt_h01, gt_tor).
-module(hafnian_nif_tests).
-include_lib("eunit/include/eunit.hrl").

ones(N) -> {ok, M} = pq_core:pack({cmat, lists:duplicate(N, lists:duplicate(N, 1))}), M.

hafnian_test() ->
    {ok, {{V, Im}, #{scale := _, subsets := 8}}} = pq_nif:hafnian(ones(6), #{}),
    ?assert(Im == 0),
    ?assert(abs(V - 15) < 1.0e-12),                                  % 5!! perfect matchings
    {ok, {{Z1, Z2}, _}} = pq_nif:hafnian(ones(5), #{}),
    ?assert(Z1 == 0 andalso Z2 == 0),
    {ok, {{R, _}, _}} = pq_nif:hafnian_repeated(ones(2), [4, 2], #{}),
    ?assert(abs(R - 15) < 1.0e-12),
    {ok, D} = pq_core:pack({cvec, [1, 1, 1, 1]}),
    {ok, {{L, _}, _}} = pq_nif:loop_hafnian(ones(4), D, [1, 1, 1, 1], #{}),
    ?assert(abs(L - 10) < 1.0e-12).                                  % involutions of 4 elements

torontonian_test() ->
    %% one-mode squeezed vacuum, r = 0.5: P(click) = 1 - 1/cosh r (FS A.10)
    R = 0.5, Ch = math:cosh(R), Sh = math:sinh(R),
    Det = math:pow(Ch, 4) - math:pow(Sh * Ch, 2),
    A = Ch * Ch / Det, B = Sh * Ch / Det,
    {ok, O} = pq_core:pack({cmat, [[1 - A, -B], [-B, 1 - A]]}),
    {ok, {T, #{terms := 2}}} = pq_nif:torontonian(O, #{}),
    ?assert(abs(T / math:sqrt(Det) - (1 - 1 / Ch)) < 1.0e-14),
    {ok, Bad} = pq_core:pack({cmat, [[2, 0], [0, 2]]}),
    ?assertEqual({error, {non_physical_state, not_positive_definite}}, pq_nif:torontonian(Bad, #{})),
    ?assertMatch({error, {bad_shape, _}}, pq_nif:torontonian(ones(3), #{})).
