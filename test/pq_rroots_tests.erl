%% pq_rroots: real roots of rational polynomials in an interval, to any precision.
-module(pq_rroots_tests).
-include_lib("eunit/include/eunit.hrl").

r(Poly, Lo, Hi) -> [pq_mp:to_string(R, 25) || R <- pq_rroots:roots(Poly, Lo, Hi, 100)].

roots_test() ->
    ?assertEqual(["1.414213562373095048801689e0"], r([-2, 0, 1], 0, 2)),
    % (x - 1)(x - 2)(x - 3); in the open interval (1, 3) only 2
    ?assertEqual(["1.000000000000000000000000e0", "2.000000000000000000000000e0", "3.000000000000000000000000e0"],
                 r([6, -11, 6, -1], 0, 4)),
    ?assertEqual(["2.000000000000000000000000e0"], r([6, -11, 6, -1], 1, 3)),
    % x^3 - x, a root at a bisection point; a double root once
    ?assertEqual(["-1.000000000000000000000000e0", "0", "1.000000000000000000000000e0"], r([0, -1, 0, 1], -2, 2)),
    ?assertEqual(["1.000000000000000000000000e0"], r([1, -2, 1], 0, 2)),
    % x^4 + 7/4 x^3 + x^2 - 5/2 x + 1/3: its two real roots, as mpmath's polyroots gives them
    ?assertEqual(["1.438684097843289079045541e-1", "7.281827521581601138523503e-1"],
                 r([{1, 3}, {-5, 2}, 1, {7, 4}, 1], -3, 3)),
    ?assertEqual([], r([1, 0, 1], -5, 5)),
    ?assertEqual([], r([3], 0, 1)).
