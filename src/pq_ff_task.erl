%% The finite-field reconstruction's tasks (decision record R1) as a pq_task: the tuples pq_recon builds, run where
%% the executor sends them. Moved here unchanged from pq_rexec, so their results are what they were.
-module(pq_ff_task).
-behaviour(pq_task).
-export([run/2, describe/1]).

%% {column, ProgramKey, P, R, N, SupportKey, Suffix}: one column of Zippel's interpolation (pq_nif:ff_column/6):
%%   {ok, Coefficients} | bad | singular | bad_prime
%% {line_column, ProgramKey, P, Shift, R, N, SupportKey, Suffix, Ts, Dn, Dd}: one column of the coefficients of a
%%   rational function along lines (pq_nif:ff_line_column/10): {ok, [Coefficients per output]} | bad | singular
%% {linsys_all, BlackBox, P, Points}: a linear system's wanted component for every right-hand side at the points
%%   (pq_nif:ff_linsys_all/3): {ok, [Values | bad]}
%% {values, ProgramKey, P, Points}: the program's outputs at the points (pq_ff:eval/3)
-spec run(tuple(), #{home := node()}) -> term().
run({column, ProgKey, P, R, N, SupportKey, Suffix}, #{home := Home}) ->
    {ok, {ff, Ref, _, _}} = pq_rstore:prepared(ProgKey, Home),
    {ok, Support} = pq_rstore:fetch(SupportKey, Home),
    pq_nif:ff_column(Ref, P, R, N, Support, Suffix);
run({line_column, ProgKey, P, Shift, R, N, SupportKey, Suffix, Ts, Dn, Dd}, #{home := Home}) ->
    {ok, {ff, Ref, _, _}} = pq_rstore:prepared(ProgKey, Home),
    {ok, Support} = pq_rstore:fetch(SupportKey, Home),
    pq_nif:ff_line_column(Ref, P, Shift, R, N, Support, Suffix, Ts, Dn, Dd);
run({linsys_all, BB, P, Points}, #{home := Home}) ->
    {ok, {ff, Ref, _, _}} = pq_rstore:prepared(BB, Home),
    pq_nif:ff_linsys_all(Ref, P, lists:append(Points));
run({values, ProgKey, P, Points}, #{home := Home}) ->
    {ok, FF} = pq_rstore:prepared(ProgKey, Home),
    pq_ff:eval(FF, P, Points).

-spec describe(tuple()) -> map().
describe(Task) -> #{kind => element(1, Task), prime => element(3, Task)}.
