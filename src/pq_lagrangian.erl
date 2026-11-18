%% Feynman rules from a Lagrangian (plan "theoretical physics instruments", F4): the vertices by functional
%% derivatives in momentum space, in a small exact tensor algebra.
%%
%% A Lagrangian is a list of monomials {Coef, Factors, Fields}:
%%   Coef     an exact fraction {N, D}, an integer, or {complex, Re, Im} of those;
%%   Factors  commuting objects with indices: {c, Name} a coupling or mass, {g, Mu, Nu} the metric, {d, A, B} the
%%            Kronecker delta of an internal index, {f, A, B, C} totally antisymmetric structure constants,
%%            {gamma, Mu} a Dirac matrix and {t, A} a colour generator (both between the fermion legs; at most one
%%            fermion line per monomial);
%%   Fields   occurrences {Field, Indices, Derivs}: the field's own indices (e.g. [mu, a] for A^a_mu) and the Lorentz
%%            indices of the derivatives acting on it.
%% Repeated index names in a monomial are summed (Lorentz or internal alike).
%%
%% vertex(L, Legs, Kinds), Legs = [{Field, Indices}] numbered 1..n with incoming momenta p_1..p_n, Kinds a map from
%% each field to the kinds of its indices (lorentz or internal): i times the sum, over the monomials with exactly
%% these fields and over every assignment of legs to the occurrences of equal fields, of the coefficient and factors,
%% each field index I taken by a leg with index L joined to it by g(I, L) (lorentz) or d(I, L) (internal), each
%% derivative d_mu on a field taken by leg k giving -i p_(k,mu); the contractions then do the renaming (an index shared
%% by two fields becomes a metric between the two legs). The legs' index names must differ from the monomials'. Grassmann fields may occur once per monomial (no sign from reordering is
%% needed then). The result is a canonical list of terms {{Re, Im}, Factors}: metrics and deltas contracted, f's indices
%% in order with the sign of the permutation, dummy indices renamed, equal terms added. Momenta are {p, K, Mu}, their
%% contraction {dot, K1, K2}.
-module(pq_lagrangian).
-export([vertex/3, propagator/2, canonical/1, equal/2, scale/2, add/2]).

%%% Vertices

vertex(L, Legs, Kinds) ->
    LegFields = lists:sort([F || {F, _} <- Legs]),
    Numbered = lists:enumerate(Legs),
    Terms = lists:append([[term(Coef, Factors, Assign, Kinds) || Assign <- assignments(Occ, Numbered)]
                          || {Coef, Factors, Occ} <- L, lists:sort([F || {F, _, _} <- Occ]) =:= LegFields]),
    canonical(scale(Terms, {{0, 1}, {1, 1}})).

%% every assignment of the legs to the monomial's occurrences of the same field
assignments([], []) -> [[]];
assignments([{F, _, _} = O | Rest], Legs) ->
    [[{O, Leg} | A] || {_, {F2, _}} = Leg <- Legs, F2 =:= F, A <- assignments(Rest, Legs -- [Leg])].

term(Coef, Factors, Assign, Kinds) ->
    Links = lists:append([begin
                              length(Idx) =:= length(LegIdx) orelse error({pq_lagrangian, {index_count, F}}),
                              [{metric_kind(K), I, Li} || {K, I, Li} <- lists:zip3(maps:get(F, Kinds, []), Idx, LegIdx)]
                          end || {{F, Idx, _}, {_, {_, LegIdx}}} <- Assign]),
    Mom = lists:append([[{p, K, Mu} || Mu <- Der] || {{_, _, Der}, {K, _}} <- Assign]),
    NDer = length(Mom),
    C = lists:foldl(fun(_, Cx) -> cmul(Cx, {{0, 1}, {-1, 1}}) end, coef(Coef), lists:seq(1, NDer)),
    {C, Factors ++ Links ++ Mom}.

metric_kind(lorentz) -> g;
metric_kind(internal) -> d.

coef({complex, Re, Im}) -> {q(Re), q(Im)};
coef(X) -> {q(X), {0, 1}}.

rename({c, _} = X, _) -> X;
rename({g, A, B}, R) -> {g, R(A), R(B)};
rename({d, A, B}, R) -> {d, R(A), R(B)};
rename({f, A, B, C}, R) -> {f, R(A), R(B), R(C)};
rename({gamma, A}, R) -> {gamma, R(A)};
rename({t, A}, R) -> {t, R(A)};
rename({p, K, A}, R) -> {p, K, R(A)}.

%%% Propagators

%% The standard propagators, momentum p (= leg 1's) flowing through, as {Numerator terms, Denominator}: the
%% denominator a list of {Coef, Factors} terms in the same algebra ({dot, 1, 1} = p^2). They are -Gamma_2^-1 for
%% the quadratic Lagrangians below (checked numerically in the tests):
%%   scalar  1/2 (d phi)^2 - 1/2 m^2 phi^2                         i / (p^2 - m^2)
%%   dirac   psibar (i gamma.d - m) psi                            i (pslash + m) / (p^2 - m^2)
%%   vector  -1/4 F^2 - 1/(2 xi) (d.A)^2, massless, R_xi gauge    -i [g^(mu nu) - (1 - xi) p^mu p^nu / p^2] / p^2
%% The vector's numerator is returned over p^2 (its second term carries {inv, dot11} for 1/p^2).
propagator(scalar, M) ->
    {[{{{0, 1}, {1, 1}}, []}], [{{{1, 1}, {0, 1}}, [{dot, 1, 1}]}, {{{-1, 1}, {0, 1}}, [{c, M}, {c, M}]}]};
propagator(dirac, M) ->
    {[{{{0, 1}, {1, 1}}, [{gamma, {dummy, 1}}, {p, 1, {dummy, 1}}]}, {{{0, 1}, {1, 1}}, [{c, M}]}],
     [{{{1, 1}, {0, 1}}, [{dot, 1, 1}]}, {{{-1, 1}, {0, 1}}, [{c, M}, {c, M}]}]};
propagator({vector, Mu, Nu}, Xi) ->
    {[{{{0, 1}, {-1, 1}}, [{g, Mu, Nu}]},
      {{{0, 1}, {1, 1}}, [{c, {one_minus, Xi}}, {p, 1, Mu}, {p, 1, Nu}, {inv, dot11}]}],
     [{{{1, 1}, {0, 1}}, [{dot, 1, 1}]}]}.

%%% Canonical form

canonical(Terms) ->
    Collected = lists:foldl(fun({C, Fs}, M) ->
                                    {S, Key} = canon_term(Fs),
                                    case S of
                                        0 -> M;
                                        _ -> M#{Key => cadd(maps:get(Key, M, {{0, 1}, {0, 1}}), cmul(C, {{S, 1}, {0, 1}}))}
                                    end
                            end, #{}, Terms),
    lists:sort([{C, Fs} || Fs := C <- Collected, C =/= {{0, 1}, {0, 1}}]).

%% contract, order, rename dummies: {Sign, Factors} (Sign 0 if the term vanishes)
canon_term(Fs0) ->
    Fs1 = contract(Fs0),
    {S, Fs2} = order_f(Fs1),
    case S of
        0 -> {0, []};
        _ -> {S, rename_dummies(Fs2)}
    end.

%% g_(mu nu) X_nu -> X_mu and d_(a b) X_b -> X_a for an index that occurs elsewhere; p_(k mu) p_(l mu) -> dot
contract(Fs) ->
    case find_contraction(Fs) of
        none -> dots(Fs);
        {Fs2} -> contract(Fs2)
    end.

find_contraction(Fs) ->
    Metrics = [{X, A, B} || {Kind, A, B} = X <- Fs, Kind =:= g orelse Kind =:= d],
    find_c(Metrics, Fs).

find_c([], _) -> none;
find_c([{X, A, B} | Rest], Fs) ->
    Others = Fs -- [X],
    case {occurs(B, Others), occurs(A, Others)} of
        _ when A =:= B -> {[{trace, element(1, X)} | Others]};     % g^mu_mu = d (dimension), d_aa = N_A
        {true, _} -> {[rename(Y, fun(I) when I =:= B -> A; (I) -> I end) || Y <- Others]};
        {_, true} -> {[rename(Y, fun(I) when I =:= A -> B; (I) -> I end) || Y <- Others]};
        _ -> find_c(Rest, Fs)
    end.

occurs(I, Fs) -> lists:any(fun(F) -> lists:member(I, indices(F)) end, Fs).

indices({c, _}) -> [];
indices({trace, _}) -> [];
indices({dot, _, _}) -> [];
indices({g, A, B}) -> [A, B];
indices({d, A, B}) -> [A, B];
indices({f, A, B, C}) -> [A, B, C];
indices({gamma, A}) -> [A];
indices({t, A}) -> [A];
indices({p, _, A}) -> [A].

%% two momenta with the same index: their scalar product
dots(Fs) ->
    Ps = [X || {p, _, _} = X <- Fs],
    case [{X, Y} || {p, _, A} = X <- Ps, {p, _, B} = Y <- Ps, X < Y, A =:= B] of
        [] -> Fs;
        [{{p, K, _} = X, {p, L, _} = Y} | _] -> dots([{dot, min(K, L), max(K, L)} | (Fs -- [X]) -- [Y]])
    end.

%% f's indices in order (sign of the permutation), g and d symmetric
order_f(Fs) ->
    lists:foldl(fun({f, A, B, C}, {S, Acc}) ->
                        case length(lists:usort([A, B, C])) of
                            3 -> {Sx, [X, Y, Z]} = sort_sign([A, B, C]), {S * Sx, [{f, X, Y, Z} | Acc]};
                            _ -> {0, Acc}
                        end;
                   ({K, A, B}, {S, Acc}) when K =:= g; K =:= d -> {S, [{K, min(A, B), max(A, B)} | Acc]};
                   (X, {S, Acc}) -> {S, [X | Acc]}
                end, {1, []}, Fs).

sort_sign(L) -> sort_sign(L, 1).
sort_sign(L, S) ->
    case swap(L) of
        none -> {S, L};
        L2 -> sort_sign(L2, -S)
    end.
swap([A, B | T]) when A > B -> [B, A | T];
swap([A | T]) -> case swap(T) of none -> none; T2 -> [A | T2] end;
swap([]) -> none.

%% dummy indices (those that appear twice in a term) renamed {dummy, 1}, {dummy, 2}, .. to the smallest sorted form
%% over their orders; f's re-sorted after renaming
rename_dummies(Fs) ->
    All = lists:append([indices(F) || F <- Fs]),
    Dummies = lists:usort([I || I <- All, length([x || J <- All, J =:= I]) > 1]),
    case Dummies of
        [] -> lists:sort(Fs);
        _ ->
            Forms = [begin
                         Map = maps:from_list(lists:zip(Order, [{dummy, K} || K <- lists:seq(1, length(Order))])),
                         Rn = [rename(F, fun(I) -> maps:get(I, Map, I) end) || F <- Fs],
                         {S, Ord} = order_f(Rn),
                         {lists:sort(Ord), S}
                     end || Order <- perms(Dummies)],
            {Best, _} = lists:min(Forms),
            Best
    end.

perms([]) -> [[]];
perms(L) -> [[X | P] || X <- L, P <- perms(L -- [X])].

%% two lists of terms equal after canonicalisation
equal(A, B) -> canonical(A) =:= canonical(B).

%%% Coefficients: Gaussian rationals {Re, Im}

scale(Terms, C) -> [{cmul(C, Cx), Fs} || {Cx, Fs} <- Terms].
add(A, B) -> A ++ B.

q(N) when is_integer(N) -> {N, 1};
q({N, D}) -> qn(N, D).
qn(N, D) when D < 0 -> qn(-N, -D);
qn(0, _) -> {0, 1};
qn(N, D) -> G = gcd(abs(N), D), {N div G, D div G}.
gcd(A, 0) -> A;
gcd(A, B) -> gcd(B, A rem B).
qa({A, B}, {C, D}) -> qn(A * D + C * B, B * D).
qm({A, B}, {C, D}) -> qn(A * C, B * D).
qneg({A, B}) -> {-A, B}.
cadd({A, B}, {C, D}) -> {qa(A, C), qa(B, D)}.
cmul({A, B}, {C, D}) -> {qa(qm(A, C), qneg(qm(B, D))), qa(qm(A, D), qm(B, C))}.
