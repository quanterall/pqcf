%% A real task from FORM's literature, reproduced with pq_symbolic: the resultant benchmarks of Ruijl, Plaat, van den
%% Herik and Vermaseren, "Why Local Search Excels in Expression Simplification" (arXiv:1409.5223, table II; also
%% HEPGAME, arXiv:1405.6369): res(m, n) = res_x(sum_(i<=m) a_i x^i, sum_(i<=n) b_i x^i), expanded, and the number of
%% operations needed to evaluate it -- as expanded, after FORM's output optimisation (Format O1-O4) and after
%% pq_symbolic's straight-line programs (compile/2).
%%
%% The resultant is the determinant of the Sylvester matrix, expanded here row by row over the subsets of columns
%% used so far (each step pq_symbolic's exact arithmetic; the determinant itself is outside pq_symbolic's scope
%% fence, SY 1.2, so it lives in this study). Writes res_<m>_<n>.txt (pq_symbolic's ASCII, which FORM reads) and
%% prints one line per expression.
%%
%%   erl -pa _build/default/lib/*/ebin -pa <dir> -noshell -eval 'resultants:main("<dir>", [{7,4},{7,5},{7,6}])'
-module(resultants).
-export([main/2, resultant/2, ops_expanded/1, ops_slp/1, form_script/3]).

-define(LIMITS, #{sym_timeout_ms => 3600000, sym_max_heap_words => 4000000000, max_poly_terms => 50000000,
                  max_expand_power => 64, max_expr_nodes => 500000000}).

main(Dir, Cases) ->
    {ok, _} = application:ensure_all_started(pq),
    [begin
         {T, R} = timer:tc(?MODULE, resultant, [M, N]),
         Terms = terms(R),
         ok = file:write_file(filename:join(Dir, io_lib:format("res_~b_~b.txt", [M, N])), pq_symbolic:pretty(R)),
         Syms = pq_symbolic:free_symbols(R),
         {Tc, {ok, Prog}} = timer:tc(pq_symbolic, compile, [[R], Syms, ?LIMITS]),
         {Mul, Add, Pow} = ops_slp(Prog),
         {Th, {ok, Hp}} = timer:tc(pq_symbolic, compile, [[R], Syms, maps:put(optimize, horner, ?LIMITS)]),
         {HM, HA, HP} = ops_slp(Hp),
         %% the two programs at a random point: equal within rounding
         Bind = maps:from_list([{X, 0.5 + rand:uniform()} || X <- Syms]),
         {ok, [V1]} = pq_symbolic:eval_slp(Prog, Bind),
         {ok, [V2]} = pq_symbolic:eval_slp(Hp, Bind),
         io:format("res(~b,~b): ~b variables, ~b terms, built in ~.2f s; ops expanded ~b; SLP ~b (mul ~b, add ~b, pow ~b) in ~.2f s; "
                   "Horner ~b (mul ~b, add ~b, pow ~b) in ~.2f s; values ~w / ~w, relative difference ~.2e~n",
                   [M, N, length(Syms), Terms, T / 1.0e6, ops_expanded(R), Mul + Add + Pow, Mul, Add, Pow, Tc / 1.0e6,
                    HM + HA + HP, HM, HA, HP, Th / 1.0e6, V1, V2, abs(V1 - V2) / max(abs(V1), 1.0e-300)])
     end || {M, N} <- Cases],
    ok.

%% res_x(A, B) for A of degree M, B of degree N with symbolic coefficients a0..aM, b0..bN
resultant(M, N) ->
    A = [{sym, list_to_atom("a" ++ integer_to_list(I))} || I <- lists:seq(M, 0, -1)],     % a_M .. a_0
    B = [{sym, list_to_atom("b" ++ integer_to_list(I))} || I <- lists:seq(N, 0, -1)],
    D = M + N,
    %% Sylvester rows: N shifted copies of A, then M shifted copies of B; an entry is a coefficient or 0
    Row = fun(Coeffs, Shift) -> fun(C) -> J = C - Shift, case J >= 0 andalso J < length(Coeffs) of
                                                             true -> lists:nth(J + 1, Coeffs);
                                                             false -> 0
                                                         end end end,
    Rows = [Row(A, K) || K <- lists:seq(0, N - 1)] ++ [Row(B, K) || K <- lists:seq(0, M - 1)],
    %% det = sum over permutations sigma of sign(sigma) prod_r S[r][sigma(r)], built row by row: the state is the set of
    %% columns used, the sign of choosing column c the parity of the used columns above c
    Final = lists:foldl(
              fun(RowF, States) ->
                      Next = maps:fold(
                               fun(Used, Val, Acc) ->
                                       lists:foldl(
                                         fun(C, Acc2) ->
                                                 case (Used band (1 bsl C)) =:= 0 andalso RowF(C) =/= 0 of
                                                     false -> Acc2;
                                                     true ->
                                                         Above = popcount(Used bsr (C + 1)),
                                                         Sgn = case Above rem 2 of 0 -> 1; 1 -> -1 end,
                                                         Term = {mul, [Sgn, RowF(C), Val]},
                                                         maps:update_with(Used bor (1 bsl C), fun(L) -> [Term | L] end, [Term], Acc2)
                                                 end
                                         end, Acc, lists:seq(0, D - 1))
                               end, #{}, States),
                      maps:map(fun(_, Ts) -> {ok, P} = pq_symbolic:canonicalise({add, Ts}, ?LIMITS), P end, Next)
              end, #{0 => 1}, Rows),
    maps:get((1 bsl D) - 1, Final).

%% the same resultant in FORM, independently: the Sylvester determinant as the product of its rows' sums over an
%% function of the columns used, antisymmetrized after every row (sorting its arguments gives the permutation's sign,
%% a repeated column gives 0, equal terms merge), then D(1, ..., m + n) = 1. With Check, PQCF's expression (file Check) is subtracted: the output must be 0.
form_script(M, N, Check) ->
    A = fun(I) -> "a" ++ integer_to_list(I) end,
    B = fun(I) -> "b" ++ integer_to_list(I) end,
    D = M + N,
    Entry = fun(R, C) -> case R < N of
                             true -> J = C - R, case J >= 0 andalso J =< M of true -> A(M - J); false -> none end;
                             false -> J = C - (R - N), case J >= 0 andalso J =< N of true -> B(N - J); false -> none end
                         end end,
    Rows = [["id D(?x) = ", lists:join(" + ", [[E, "*D(?x,", integer_to_list(C + 1), ")"]
                                               || C <- lists:seq(0, D - 1), E <- [Entry(R, C)], E =/= none]), ";\nAntiSymmetrize D;\n.sort\n"]
            || R <- lists:seq(0, D - 1)],
    lists:flatten(
      ["Symbols ", lists:join(",", [A(I) || I <- lists:seq(0, M)] ++ [B(I) || I <- lists:seq(0, N)]), ";\n",
       "CFunction D;\nLocal F = D;\n.sort\n", Rows,
       "id D(", lists:join(",", [integer_to_list(C) || C <- lists:seq(1, D)]), ") = 1;\n.sort\n",
       case Check of
           none -> "#write \"terms %d\", termsin_(F)\n.end\n";
           File -> ["Local P =\n#include ", File, "\n;\nLocal Z = F - P;\n.sort\nDrop F,P;\n#write \"difference: %E\", Z\n.end\n"]
       end]).

popcount(0) -> 0;
popcount(X) -> (X band 1) + popcount(X bsr 1).

terms({add, L}) -> length(L);
terms(0) -> 0;
terms(_) -> 1.

%% the operation count of the expanded form, as the paper counts it: per term, its multiplications (a power x^e as
%% e - 1 multiplications, a coefficient other than +-1 as one), and terms - 1 additions
ops_expanded({add, L}) -> length(L) - 1 + lists:sum([term_muls(T) || T <- L]);
ops_expanded(T) -> term_muls(T).

term_muls({mul, Fs}) ->
    Factors = [F || F <- Fs, not is_number(F) orelse (F =/= 1 andalso F =/= -1)],
    max(0, lists:sum([degree(F) || F <- Factors]) - 1);
term_muls(F) -> max(0, degree(F) - 1).

degree({pow, {sym, _}, E}) when is_integer(E) -> E;
degree({sym, _}) -> 1;
degree(N) when is_number(N) -> 1.

%% the operations of a straight-line program by kind (inputs and constants are not operations)
ops_slp({slp, 1, _, Ops, _}) ->
    {length([x || {mul, _, _, _} <- Ops]), length([x || {add, _, _, _} <- Ops]), length([x || {pow, _, _, _} <- Ops])}.
