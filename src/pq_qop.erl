%% Operator algebra for quantum optics (plan Q1, after decision record SY2): bosonic modes and N-level atoms,
%% products in a canonical form, commutators, adjoints, and the Heisenberg-Lindblad equations of motion.
%%
%% Spaces: #{Name => boson | {atom, N}} (and domains => #{Symbol => Domain} for the coefficients), the levels of an atom 1 (the ground state) .. N. An expression is a map
%% #{Word => Coef}: Coef an exact expression of the symbolic engine (canonical, never 0), Word a sorted list of
%% {Name, Factor} with one factor per subsystem that is not the identity:
%%   a boson:  {P, Q} for (a^dagger)^P a^Q (normal order), P + Q > 0;
%%   an atom:  {I, J} for |I><J|, never |1><1| (completeness: |1><1| = 1 - sum_(i>1) |i><i|).
%% Different subsystems commute, so the word is sorted by subsystem. This form is unique: equal operators have equal
%% expressions.
%%
%% Two-level atoms: sigma^- = |1><2|, sigma^+ = |2><1|, sigma_z = |2><2| - |1><1| = 2 |2><2| - 1.
-module(pq_qop).
-export([one/0, num/1, a/1, ad/1, proj/3, sm/1, sp/1, sz/1,
         add/1, add/2, sub/2, mul/2, mul/1, scale/2, pow/2, comm/2, dag/1, is_zero/1,
         heisenberg/4, with_spaces/2, to_ast/1, pretty/1]).

-export_type([expr/0]).
-type expr() :: #{list() => term()}.

%%% Construction

-spec one() -> expr().
one() -> #{[] => 1}.

%% a number or an expression of the engine, times the identity
-spec num(term()) -> expr().
num(C) -> mk([{[], C}]).

a(M) -> #{[{M, {0, 1}}] => 1}.
ad(M) -> #{[{M, {1, 0}}] => 1}.

%% |I><J| of atom At; |1><1| is expanded by completeness under the spaces (mul/2 and resolve/1)
proj(At, I, J) -> resolve(#{[{At, {I, J}}] => 1}).

sm(At) -> proj(At, 1, 2).
sp(At) -> proj(At, 2, 1).
sz(At) -> add(scale(2, proj(At, 2, 2)), num(-1)).

%%% Arithmetic

add(L) when is_list(L) -> lists:foldl(fun add/2, #{}, L).

add(A, B) ->
    maps:fold(fun(W, C, Acc) ->
                      case maps:find(W, Acc) of
                          {ok, D} -> put_coef(W, cadd(C, D), Acc);
                          error -> Acc#{W => C}
                      end
              end, A, B).

sub(A, B) -> add(A, scale(-1, B)).

scale(K, A) -> maps:filtermap(fun(_, C) -> nz(cmul(K, C)) end, A).

mul(L) when is_list(L) -> lists:foldl(fun(X, Acc) -> mul(Acc, X) end, one(), L).

%% the product A B in canonical form; atoms with their levels from the spaces (2 if not declared)
mul(A, B) ->
    maps:fold(fun(Wa, Ca, Acc) ->
                      maps:fold(fun(Wb, Cb, Acc2) ->
                                        lists:foldl(fun({W, K}, Acc3) -> add(Acc3, mk([{W, cmul(cmul(Ca, Cb), K)}])) end,
                                                    Acc2, word_mul(Wa, Wb))
                                end, Acc, B)
              end, #{}, A).

pow(_, 0) -> one();
pow(A, N) when N > 0 -> mul(A, pow(A, N - 1)).

%% [A, B] = A B - B A
comm(A, B) -> sub(mul(A, B), mul(B, A)).

%% the adjoint, coefficients conjugated. Factors of different subsystems commute, so each factor is taken alone:
%% (a^dagger)^P a^Q becomes (a^dagger)^Q a^P, still normal ordered, and |I><J| becomes |J><I|. Both are {X, Y} ->
%% {Y, X}, and the words stay canonical.
dag(A) -> maps:fold(fun(W, C, Acc) -> Acc#{[{N, {Y, X}} || {N, {X, Y}} <- W] => conj(C)} end, #{}, A).

is_zero(A) -> map_size(A) =:= 0.

%%% The product of two words (subsystem by subsystem)

%% [{Word, Coefficient}]: the canonical words of the product and their integer coefficients
word_mul(Wa, Wb) ->
    Names = lists:usort([N || {N, _} <- Wa ++ Wb]),
    lists:foldl(fun(N, Acc) ->
                        Fa = proplists:get_value(N, Wa, none),
                        Fb = proplists:get_value(N, Wb, none),
                        Parts = factor_mul(N, Fa, Fb),
                        [{W ++ Fw, K * Kf} || {W, K} <- Acc, {Fw, Kf} <- Parts]
                end, [{[], 1}], Names).

%% the product of two factors of subsystem N: [{FactorList, K}] (FactorList [] for the identity)
factor_mul(_, none, none) -> [{[], 1}];
factor_mul(N, F, none) -> single(N, F);
factor_mul(N, none, F) -> single(N, F);
factor_mul(N, Fa, Fb) ->
    case kind(N) of
        boson ->
            % (ad^p a^q)(ad^r a^s) = sum_k C(q,k) C(r,k) k! ad^(p+r-k) a^(q+s-k)
            {P, Q} = Fa, {R, S} = Fb,
            [{[{N, {P + R - K, Q + S - K}} || P + R - K + Q + S - K > 0], binom(Q, K) * binom(R, K) * fact(K)}
             || K <- lists:seq(0, min(Q, R))];
        {atom, L} ->
            {I, J} = Fa, {K, M} = Fb,
            case J =:= K of
                false -> [];
                true when I =:= 1, M =:= 1 ->
                    % |1><1| = 1 - sum_(i>1) |i><i|
                    [{[], 1} | [{[{N, {X, X}}], -1} || X <- lists:seq(2, L)]];
                true -> [{[{N, {I, M}}], 1}]
            end
    end.

%% a lone factor, |1><1| of an atom expanded
single(N, {1, 1} = F) ->
    case kind(N) of
        {atom, L} -> [{[], 1} | [{[{N, {X, X}}], -1} || X <- lists:seq(2, L)]];
        boson -> [{[{N, F}], 1}]
    end;
single(N, F) -> [{[{N, F}], 1}].

%%% Spaces: set for the calling process (with_spaces/2 in the tests and the API below), default: a name declared
%%% boson unless declared an atom

kind(N) ->
    case get({?MODULE, spaces}) of
        undefined -> default_kind(N);
        Sp -> maps:get(N, Sp, default_kind(N))
    end.

default_kind(_) -> boson.

%%% The Heisenberg-Lindblad equation of an operator O: dO/dt = i [H, O] + sum_k (L_k^dagger O L_k -
%%% (1/2) {L_k^dagger L_k, O}), the jump operators L_k including their rates (sqrt(gamma) a). Spaces declare the atoms.

-spec heisenberg(map(), fun(() -> expr()) | expr(), [expr()], expr()) -> expr().
heisenberg(Spaces, H, Ls, O) ->
    with_spaces(Spaces, fun() ->
                                Hh = resolve(H), Os = resolve(O), Lss = [resolve(L) || L <- Ls],
                                Unitary = scale(i(), comm(Hh, Os)),
                                Diss = add([begin
                                                Ld = dag(L),
                                                LdL = mul(Ld, L),
                                                sub(mul([Ld, Os, L]), scale({rat, 1, 2}, add(mul(LdL, Os), mul(Os, LdL))))
                                            end || L <- Lss]),
                                add(Unitary, Diss)
                        end).

%% Spaces may also give the domains of the coefficients' symbols (domains => #{kappa => positive, g => real}), which
%% the engine uses when it simplifies (conj(sqrt(kappa)) sqrt(kappa) = kappa for a positive kappa)
with_spaces(Spaces, F) ->
    Old = get({?MODULE, spaces}),
    OldD = get({?MODULE, domains}),
    put({?MODULE, spaces}, maps:remove(domains, Spaces)),
    put({?MODULE, domains}, maps:get(domains, Spaces, #{})),
    try F()
    after
        case Old of undefined -> erase({?MODULE, spaces}); _ -> put({?MODULE, spaces}, Old) end,
        case OldD of undefined -> erase({?MODULE, domains}); _ -> put({?MODULE, domains}, OldD) end
    end.

%% an expression built before the spaces were known, put in canonical form under them
resolve(A) -> mul(one(), A).

%%% Output

%% the expression as an engine-style sum: {mul, [Coef, {op, ...}...]}, the operators as {op, Name, {P, Q}} (boson)
%% or {op, Name, {I, J}} (atom)
to_ast(A) ->
    case lists:sort(maps:to_list(A)) of
        [] -> 0;
        Ts -> {add, [{mul, [C | [{op, N, F} || {N, F} <- W]]} || {W, C} <- Ts]}
    end.

pretty(A) ->
    case lists:sort(maps:to_list(A)) of
        [] -> "0";
        Ts ->
            lists:flatten(lists:join(" + ", [[coef_s(C), [[" ", op_s(N, F)] || {N, F} <- W]] || {W, C} <- Ts]))
    end.

coef_s(C) ->
    case catch pq_symbolic:pretty(C) of
        B when is_binary(B) -> "(" ++ binary_to_list(B) ++ ")";
        _ -> io_lib:format("~p", [C])
    end.

op_s(N, {P, Q}) ->
    case kind(N) of
        boson -> [[io_lib:format("~s^+", [N]) || _ <- lists:seq(1, P)], [io_lib:format("~s", [N]) || _ <- lists:seq(1, Q)]];
        _ -> io_lib:format("|~b><~b|_~s", [P, Q, N])
    end.

%%% Coefficients (expressions of the engine)

mk(Terms) -> lists:foldl(fun({W, C}, Acc) -> add(Acc, maps:filtermap(fun(_, X) -> nz(X) end, #{W => C})) end, #{}, Terms).

nz(C) -> case can(C) of 0 -> false; X -> {true, X} end.

put_coef(W, C, Acc) -> case can(C) of 0 -> maps:remove(W, Acc); X -> Acc#{W => X} end.

cadd(A, B) -> can({add, [A, B]}).
cmul(A, B) -> can({mul, [A, B]}).
conj(C) -> can({conj, C}).
i() -> {complex, 0, 1}.

can(X) when is_integer(X) -> X;
can(X) ->
    Domains = case get({?MODULE, domains}) of undefined -> #{}; D -> D end,
    case pq_symbolic:canonicalise(X, #{domains => Domains}) of
        {ok, C} -> C;
        {error, E} -> error({pq_qop, E})
    end.

binom(N, K) -> fact(N) div (fact(K) * fact(N - K)).
fact(0) -> 1;
fact(N) -> N * fact(N - 1).
