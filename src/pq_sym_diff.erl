%% Differentiation of the canonical form (SY 5.3, FS FR-S05): term by term, the product rule over each monomial's
%% factors, the result canonical.
%%
%%   x^e            e x^(e-1)                         (other symbols, pi, radicals: 0)
%%   Base(B)^e      e Base(B)^(e-1) dB/dx
%%   conj(B^e)      e conj(B^(e-1)) conj(dB/dx)       (x real)
%%   exp(c b)       c (db/dx) exp(c b)
%%
%% Only real-domain symbols (real, positive, nonnegative): a complex x, or z* in the expression with z the
%% variable, is {differentiate_unsupported, Node} (no Wirtinger calculus, FR-S05).
-module(pq_sym_diff).
-export([diff/3]).

-define(ONE, {[], []}).

diff(P, X, Ctx) ->
    case maps:get(X, maps:get(domains, Ctx, #{}), real) of
        complex -> throw({differentiate_unsupported, {sym, X}});
        _ -> ok
    end,
    maps:fold(fun(M, C, Acc) -> pq_sym_poly:add(Acc, pq_sym_poly:scale(mono(M, X, Ctx), C)) end,
              pq_sym_poly:zero(), P).

%% d/dx of one monomial (coefficient 1), as a polynomial
mono({Alg, Exp}, X, Ctx) ->
    AlgTerms = [begin
                    Rest = {lists:delete(F, Alg), Exp},
                    pq_sym_poly:mul(#{Rest => pq_sym_exact:c(1)}, atom(F, X, Ctx), Ctx)
                end || F <- Alg],
    ExpTerms = [begin
                    Rest = {Alg, lists:delete(F, Exp)},
                    pq_sym_poly:mul(#{Rest => pq_sym_exact:c(1)}, exp_factor(F, X, Ctx), Ctx)
                end || F <- Exp],
    pq_sym_poly:sum(AlgTerms ++ ExpTerms).

atom({{sym, X}, E}, X, _Ctx) ->
    Lower = pq_sym_exact:rsub(E, {1, 1}),
    Base = case Lower of
               {0, 1} -> pq_sym_poly:one();
               _ -> pq_sym_poly:atom({sym, X}, Lower)
           end,
    pq_sym_poly:scale(Base, pq_sym_exact:c(E));
atom({{conj, X}, _}, X, _Ctx) -> throw({differentiate_unsupported, {conj, {sym, X}}});
atom({{base, B}, E}, X, Ctx) ->
    DB = diff(B, X, Ctx),
    power_rule({base, B}, E, DB, Ctx);
atom({{cbase, B}, E}, X, Ctx) ->
    DB = pq_sym_poly:conj(diff(B, X, Ctx), Ctx),
    power_rule({cbase, B}, E, DB, Ctx);
atom(_, _X, _Ctx) -> pq_sym_poly:zero().

%% e Key^(e-1) dB, normalised through a product (so Base(B)^1 becomes B, and so on)
power_rule(_Key, _E, DB, _Ctx) when map_size(DB) =:= 0 -> pq_sym_poly:zero();
power_rule(Key, E, DB, Ctx) ->
    Lower = pq_sym_exact:rsub(E, {1, 1}),
    Pow = case Lower of
              {0, 1} -> pq_sym_poly:one();
              _ -> pq_sym_poly:atom(Key, Lower)
          end,
    pq_sym_poly:scale(pq_sym_poly:mul(Pow, DB, Ctx), pq_sym_exact:c(E)).

%% d/dx exp(c b) = c (db/dx) exp(c b)
exp_factor({B, C} = F, X, Ctx) ->
    DB = mono(B, X, Ctx),
    case map_size(DB) of
        0 -> pq_sym_poly:zero();
        _ -> pq_sym_poly:mul(#{{[], [F]} => pq_sym_exact:c(1)}, pq_sym_poly:scale(DB, C), Ctx)
    end.
