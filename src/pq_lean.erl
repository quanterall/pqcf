%% PQCF's results certified by Lean 4 (through lean-check, a separate application: Lean 4 + Mathlib behind a pool
%% of supervised workers). PQCF stays the fast, untrusted solver; for a result it writes a Lean statement a reader
%% can check by eye and a proof Lean checks by machine.
%%
%%   factor(E)            the factorisation of pq_symbolic:factor/1:  forall x y .. : Q,  E = c * F1^m1 * ..   (ring)
%%   gcd(A, B)            the gcd of pq_symbolic:gcd/2 with its cofactors:  A = G * A' /\ B = G * B'          (ring)
%%                        (G a common divisor with these cofactors; that it is the greatest is not certified)
%%   linear(Rows, X)      an exact solution of a linear system, Rows = [{[a_i1, .., a_in], b_i}]:
%%                        a_i1 x_1 + .. + a_in x_n = b_i for every row, in exact rationals                    (norm_num)
%%
%% Each returns {ok, #{source => LeanSource, verdict => Verdict}} with lean-check's verdict, or {ok, #{source => ..}}
%% when no checker is asked for (Opts: checker => local | {node, Node} | none, default local when lean_check is
%% loaded, else none), or {error, Why}. PQCF does not depend on lean-check at build time.
-module(pq_lean).
-export([factor/1, factor/2, gcd/2, gcd/3, linear/2, linear/3, to_lean/1, check/2]).

%%% Certificates

factor(E) -> factor(E, #{}).
factor(E, Opts) ->
    case pq_symbolic:factor(E) of
        {ok, {C, Fs}} ->
            Rhs = {mul, [C | [case M of 1 -> F; _ -> {pow, F, M} end || {F, M} <- Fs]]},
            Src = theorem(<<"pqcf_factor">>, vars([E | [F || {F, _} <- Fs]]), [eq(E, Rhs)], <<"ring">>),
            check(Src, Opts);
        Err -> Err
    end.

gcd(A, B) -> gcd(A, B, #{}).
gcd(A, B, Opts) ->
    case pq_symbolic:gcd(A, B) of
        {ok, G} ->
            case {cofactor(A, G), cofactor(B, G)} of
                {{ok, Ca}, {ok, Cb}} ->
                    Src = theorem(<<"pqcf_gcd">>, vars([A, B, G, Ca, Cb]),
                                  [eq(A, {mul, [G, Ca]}), eq(B, {mul, [G, Cb]})], <<"constructor <;> ring">>),
                    check(Src, Opts);
                _ -> {error, {pq_lean, cofactor}}
            end;
        Err -> Err
    end.

%% A / G as a polynomial
cofactor(_, 0) -> {error, zero_gcd};
cofactor(A, G) -> pq_symbolic:cancel({mul, [A, {pow, G, -1}]}).

linear(Rows, X) -> linear(Rows, X, #{}).
linear(Rows, X, Opts) ->
    Eqs = [eq({add, [{mul, [num(Aij), num(Xj)]} || {Aij, Xj} <- lists:zip(As, X)]}, num(B)) || {As, B} <- Rows],
    Src = theorem(<<"pqcf_linear">>, [], Eqs, <<"norm_num">>),
    check(Src, Opts).

num(N) when is_integer(N) -> N;
num({N, D}) -> {rat, N, D};
num({rat, _, _} = Q) -> Q.

%%% Lean source

%% theorem Name (x y : ℚ) : Eq1 ∧ Eq2 := by Tactic
theorem(Name, Vars, Eqs, Tactic) ->
    Binders = case Vars of
                  [] -> <<>>;
                  _ -> [" (", lists:join(" ", [ident(V) || V <- Vars]), " : ℚ)"]
              end,
    Stmt = lists:join(" ∧ ", Eqs),
    unicode:characters_to_binary(["theorem ", Name, Binders, " :\n    ", Stmt, " := by\n  ", Tactic, "\n"]).

eq(L, R) -> ["(", to_lean(L), " : ℚ) = ", to_lean(R)].

%% a polynomial expression with exact rational coefficients as Lean source over ℚ (character data: Unicode code points)
to_lean(N) when is_integer(N), N >= 0 -> integer_to_list(N);
to_lean(N) when is_integer(N) -> ["(", integer_to_list(N), ")"];
to_lean({rat, N, D}) -> ["((", integer_to_list(N), " : ℚ) / ", integer_to_list(D), ")"];
to_lean({sym, X}) -> ident(X);
to_lean({add, []}) -> "0";
to_lean({mul, []}) -> "1";
to_lean({add, L}) -> ["(", lists:join(" + ", [to_lean(T) || T <- L]), ")"];
to_lean({mul, L}) -> ["(", lists:join(" * ", [to_lean(T) || T <- L]), ")"];
to_lean({pow, X, K}) when is_integer(K), K >= 0 -> ["(", to_lean(X), " ^ ", integer_to_list(K), ")"];
to_lean(X) -> error({pq_lean, {not_a_polynomial, X}}).

%% a symbol as a Lean identifier: plain when it is one, else in «»
ident(X) when is_atom(X) -> ident(atom_to_binary(X));
ident(X) ->
    case re:run(X, <<"^[A-Za-z_][A-Za-z0-9_]*$">>) =/= nomatch andalso not lists:member(X, keywords()) of
        true -> X;
        false -> ["«", X, "»"]
    end.

keywords() -> [<<"theorem">>, <<"lemma">>, <<"def">>, <<"fun">>, <<"by">>, <<"have">>, <<"show">>, <<"from">>, <<"at">>,
               <<"if">>, <<"then">>, <<"else">>, <<"do">>, <<"let">>, <<"in">>, <<"end">>, <<"with">>, <<"match">>,
               <<"open">>, <<"namespace">>, <<"section">>, <<"variable">>, <<"universe">>, <<"where">>, <<"calc">>,
               <<"import">>, <<"instance">>, <<"structure">>, <<"class">>, <<"axiom">>, <<"example">>].

vars(Es) -> lists:usort(lists:append([syms(E) || E <- Es])).
syms({sym, X}) -> [X];
syms(T) when is_tuple(T) -> lists:append([syms(X) || X <- tuple_to_list(T)]);
syms(L) when is_list(L) -> lists:append([syms(X) || X <- L]);
syms(_) -> [].

%%% Checking

%% the source and, when a checker is available, lean-check's verdict
check(Src, Opts) ->
    Checker = maps:get(checker, Opts, default),
    Timeout = maps:get(timeout, Opts, 300000),
    case checker(Checker) of
        none -> {ok, #{source => Src}};
        % lean_check through a variable: PQCF does not depend on it, and xref must not look for it
        local -> M = lean_check, {ok, #{source => Src, verdict => M:check(Src, #{timeout => Timeout})}};
        {node, N} -> {ok, #{source => Src, verdict => erpc:call(N, lean_check, check, [Src, #{timeout => Timeout}],
                                                                Timeout + 10000)}}
    end.

checker(default) ->
    case code:ensure_loaded(lean_check) of
        {module, _} -> local;
        _ -> none
    end;
checker(C) -> C.
