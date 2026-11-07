%% Solving a linear system with polynomial entries (decision record R1, phase 5): PQCF reconstructs one solution
%% component from eliminations modulo primes at sample points (pq_recon:linear/4); FORM 5 eliminates symbolically,
%% its coefficients rational functions under PolyRatFun, reduced by polynomial GCDs (FLINT) at every step. The
%% system is shaped like an integration-by-parts system: n x n, sparse, every entry linear in V variables, the
%% wanted unknown one "master". Two kinds of right-hand side: random small integers (the answer is then about as
%% large as the elimination's intermediates), or b = A u* for a compact u* (every component a polynomial of degree
%% 1), as in reductions whose reduced coefficients are far simpler than the way to them. Writes the FORM script
%% (with PQCF's answer, so that FORM checks U - PQCF = 0).
%%
%%   erl -pa _build/default/lib/*/ebin -pa <dir> -noshell -eval 'linsys:main("<dir>", [{8, 3}])'
-module(linsys).
-export([main/2, main/3, system/5]).

vars(V) -> [list_to_atom("x" ++ integer_to_list(I)) || I <- lists:seq(1, V)].

%% A and b: the diagonal always present, other entries with probability 1/Density; entries a + sum c_i x_i with small
%% integer coefficients; b small integers; seeded, so FORM and PQCF get the same system
system(N, V, Density, Seed, Kind) ->
    S0 = rand:seed_s(exsss, {Seed, N, V}),
    Xs = [{sym, X} || X <- vars(V)],
    Lin = fun(S) ->
                  {Cs, S1} = lists:mapfoldl(fun(_, Sx) -> {C, Sy} = rand:uniform_s(7, Sx), {C - 4, Sy} end, S, Xs),
                  {C0, S2} = rand:uniform_s(9, S1),
                  {{add, [C0 - 5 | [{mul, [C, X]} || {C, X} <- lists:zip(Cs, Xs), C =/= 0]]}, S2}
          end,
    {Rows, S3} = lists:mapfoldl(fun(I, S) ->
                                        lists:mapfoldl(fun(J, Sx) ->
                                                               {R, Sy} = rand:uniform_s(Density, Sx),
                                                               case I =:= J orelse R =:= 1 of
                                                                   true -> Lin(Sy);
                                                                   false -> {0, Sy}
                                                               end
                                                       end, S, lists:seq(1, N))
                                end, S0, lists:seq(1, N)),
    {B, _} = case Kind of
                 random -> lists:mapfoldl(fun(_, S) -> {C, S2} = rand:uniform_s(5, S), {C - 3, S2} end, S3, lists:seq(1, N));
                 compact ->
                     %% u*_j = c_j + x_(j mod V); b = A u*, expanded
                     Ustar = [{add, [J, lists:nth((J rem V) + 1, Xs)]} || J <- lists:seq(1, N)],
                     {[begin {ok, Bi} = pq_symbolic:canonicalise({add, [{mul, [Aij, Uj]} || {Aij, Uj} <- lists:zip(Row, Ustar)]}), Bi end
                       || Row <- Rows], S3}
             end,
    {Rows, B}.

main(Dir, Cases) -> main(Dir, Cases, random).

main(Dir, Cases, Kind) ->
    {ok, _} = application:ensure_all_started(pq),
    [begin
         {A, B} = system(N, V, 3, 1, Kind),
         Vars = vars(V),
         {T, {ok, [{0, {Num, Den}}]}} = timer:tc(pq_recon, linear, [A, B, Vars, #{components => [0]}]),
         P = fun(E) -> binary_to_list(pq_symbolic:pretty(E)) end,
         File = filename:join(Dir, io_lib:format("lin_~s_~b_~b.frm", [Kind, N, V])),
         ok = file:write_file(File, form_script(A, B, Vars, P(pq_sym_norm:to_ast(Num)), P(pq_sym_norm:to_ast(Den)))),
         io:format("linear system (~s) n=~b V=~b: u_1 reconstructed in ~.2f s: numerator ~b terms, denominator ~b terms~n",
                   [Kind, N, V, T / 1.0e6, maps:size(Num), maps:size(Den)])
     end || {N, V} <- Cases],
    ok.

%% FORM: E_i = sum_j rat(A_ij, 1) u(j) - rat(b_i, 1) u(0); for k = n .. 2: E_k divided by its coefficient of u(k),
%% then subtracted from every E_j (j < k) times their coefficient of u(k); at the end E_1 = c1 u(1) + c0 u(0), so
%% u(1) = -c0 / c1. Then U - rat(PQCF numerator, PQCF denominator).
form_script(A, B, Vars, PNum, PDen) ->
    P = fun(E) -> binary_to_list(pq_symbolic:pretty(E)) end,
    N = length(B),
    Eqs = [["Local E", integer_to_list(I), " = ",
            lists:join(" + ", [["rat(", P(Aij), ",1)*u(", integer_to_list(J), ")"] || {J, Aij} <- lists:enumerate(Row), Aij =/= 0]
                       ++ [["rat(", P({mul, [-1, Bi]}), ",1)*u(0)"]]), ";\n"]
           || {I, {Row, Bi}} <- lists:enumerate(lists:zip(A, B))],
    lists:flatten(
      ["#-\n#:MaxTermSize 40M\n#:WorkSpace 400M\n#:LargeSize 2G\n#:SmallSize 500M\n#:TermsInSmall 10M\n",
       "Symbols ", lists:join(",", [atom_to_list(X) || X <- Vars]), ";\nCFunction rat,u;\nPolyRatFun rat;\n",
       Eqs, ".sort\n",
       "#do k = ", integer_to_list(N), ", 2, -1\n",
       "  Bracket u;\n  .sort\n",
       "  #$p = E`k'[u(`k')];\n  #inside $p\n    id rat(?a) = rat(reverse_(?a));\n  #endinside\n",
       "  Local E`k' = E`k'*$p;\n  Bracket u;\n  .sort\n",
       "  #do j = 1, `k'-1\n    #$c`j' = E`j'[u(`k')];\n  #enddo\n",
       "  #do j = 1, `k'-1\n    Local E`j' = E`j' - $c`j'*E`k';\n  #enddo\n  .sort\n  Drop E`k';\n#enddo\n",
       "Bracket u;\n.sort\n#$c1 = E1[u(1)];\n#$c0 = E1[u(0)];\n#inside $c1\n  id rat(?a) = rat(reverse_(?a));\n#endinside\n",
       "Local U = -$c0*$c1;\n.sort\nDrop E1;\n",
       "Local Z = U - rat(\n", PNum, ",\n", PDen, ");\n.sort\nDrop U;\nPrint Z;\n.end\n"]).
