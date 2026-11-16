%% pq_quad: tanh-sinh quadrature at arbitrary precision, on integrals with closed forms, with singular ends and a
%% breakpoint, and with a complex integrand.
-module(pq_quad_tests).
-include_lib("eunit/include/eunit.hrl").

-define(P, 200).

%% |x - y| =< 2^-(P - 8), x the real part of the integral
same({R, _}, Want) -> same(R, Want);
same(R, Want) ->
    D = pq_mp:sub(R, Want, ?P),
    pq_mp:cmp(abs_(D), pq_mp:ldexp({mp, 1, 0}, -(?P - 8))) =< 0.

abs_({mp, M, E}) -> {mp, abs(M), E}.

w() -> pq_quad:bits(?P).

closed_forms_test_() ->
    {timeout, 120, fun() ->
        W = w(),
        One = {mp, 1, 0},
        % int_0^1 1/(1 + x^2) = pi/4
        ?assert(same(pq_quad:integrate(fun(X) -> pq_mp:divide(One, pq_mp:add(One, pq_mp:mul(X, X, W), W), W) end, 0, 1, ?P),
                     pq_mp:ldexp(pq_mp:pi(?P), -2))),
        % int_0^1 log x = -1 (a logarithmic end)
        ?assert(same(pq_quad:integrate(fun(X) -> pq_mp:log(X, W) end, 0, 1, ?P), pq_mp:from(-1, ?P))),
        % int_0^1 x^-1/2 = 2 (an inverse square root at the end)
        ?assert(same(pq_quad:integrate(fun(X) -> pq_mp:divide(One, pq_mp:sqrt(X, W), W) end, 0, 1, ?P), pq_mp:from(2, ?P))),
        % int_0^1 log|x - 1/3| = (2/3) log(2/3) + (1/3) log(1/3) - 1, split at the singular point 1/3
        Third = pq_mp:from({1, 3}, W),
        Want = pq_mp:sub(pq_mp:add(pq_mp:mul(pq_mp:from({2, 3}, ?P), pq_mp:log(pq_mp:from({2, 3}, ?P), ?P), ?P),
                                   pq_mp:mul(pq_mp:from({1, 3}, ?P), pq_mp:log(pq_mp:from({1, 3}, ?P), ?P), ?P), ?P),
                         One, ?P),
        ?assert(same(pq_quad:integrate(fun(X) -> pq_mp:log(abs_(pq_mp:sub(X, Third, W)), W) end, 0, 1, [{1, 3}], ?P), Want))
    end}.

%% a complex integrand: int_0^pi e^(i x) = 2 i
complex_test_() ->
    {timeout, 120, fun() ->
        W = w(),
        {R, I} = pq_quad:integrate(fun(X) -> {pq_mp:cos(X, W), pq_mp:sin(X, W)} end, 0, pq_mp:pi(W), ?P),
        ?assert(same(R, pq_mp:from(0, ?P))),
        ?assert(same(I, pq_mp:from(2, ?P)))
    end}.
