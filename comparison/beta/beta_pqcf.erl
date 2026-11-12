%% PQCF's side of the three-loop beta function (comparison/beta/README.md): the gluon and ghost self-energies'
%% integrands as FORM writes them before Forcer (se.frm -D DUMP: one line per graph, its vertices vx(...), its
%% scalar products sp(a, b, power) and a coefficient in ep and the colour factors), reduced by PQCF -- one integral
%% family per graph, one combination per colour monomial, all on one system (pq_ibp:reduce_sums/3) -- and composed
%% with the values of PQCF's masters as Forcer gives them, then expanded in MS-bar as Forcer's msbarexpand=3 does.
%%   reduce(Dir, Field, Loops): reduces Dir/dump_<Field>_<Loops>.txt; writes the masters for Forcer
%%     (masters_<Field>_<Loops>.frm) and the reductions (<Field>_<Loops>.red); {Seconds, Info}
%%   compose(Dir, Field, Loops, Tables): after Forcer has run on the masters, the MS-bar series of P_L
-module(beta_pqcf).
-export([reduce/3, compose/4, parse/1, split/1]).

-define(COLOUR, ['cA', 'cR', 'I2R', 'nf', 'NA']).

%% the dump: [{Vx, N, [{CoefficientText, Powers, Sps}]}] grouped by graph
parse(File) ->
    {ok, Bin} = file:read_file(File),
    Text = re:replace(Bin, "\\s+", "", [global, {return, list}]),
    Terms = [T || T <- re:split(Text, "(?=[+-]vx\\()", [{return, list}]), T =/= ""],
    Parsed = [parse_term(T) || T <- Terms],
    Groups = maps:groups_from_list(fun({Vx, _, _, _, _}) -> Vx end, fun({_, N, C, P, S}) -> {N, C, P, S} end, Parsed),
    [{Vx, N, [{C, P, S} || {_, C, P, S} <- Ts]} || Vx := Ts = [{N, _, _, _} | _] <- Groups].

parse_term([Sign | Rest]) ->
    %% Key*(Coefficient): the key is vx(...) and sp(...) factors; the coefficient the last bracket
    I = string:rstr(Rest, "*("),
    Key = string:slice(Rest, 0, I - 1),
    Coef = string:slice(Rest, I + 1, length(Rest) - I - 2),
    {match, Vxs} = re:run(Key, "vx\\([^)]*\\)", [global, {capture, all, list}]),
    Sps = case re:run(Key, "sp\\(([^,]+),([^,]+),(-?[0-9]+)\\)", [global, {capture, all_but_first, list}]) of
              {match, M} -> [{A, B, list_to_integer(E)} || [A, B, E] <- M];
              nomatch -> []
          end,
    Lines = lists:usort([list_to_integer(X) || {match, Ms} <- [re:run(Key, "p([0-9]+)", [global, {capture, all_but_first, list}])],
                                                [X] <- Ms]),
    N = lists:max(Lines),
    Powers = [lists:sum([-E || {A, B, E} <- Sps, A =:= "p" ++ integer_to_list(J), B =:= A]) || J <- lists:seq(1, N)],
    Num = lists:append([lists:duplicate(E, {name(A), name(B)}) || {A, B, E} <- Sps, A =/= B, E > 0]),
    [] = [x || {A, B, E} <- Sps, A =/= B, E < 0],
    {lists:flatten(lists:join("*", [V || [V] <- Vxs])), N, [Sign] ++ "(" ++ Coef ++ ")", Powers, Num}.

name("Q") -> q;
name(X) -> list_to_atom(X).

%% a coefficient (ep and the colour factors) as #{ColourMonomial => Ast in d}
split(Text) ->
    {ok, Ast} = pq_symbolic:parse(list_to_binary(re:replace("0" ++ Text, "\\(([-+])", "(0\\1", [global, {return, list}]))),
    {ok, C} = pq_symbolic:canonicalise(Ast),
    Terms = case C of {add, Ts} -> Ts; T -> [T] end,
    Eps = {mul, [{rat, 1, 2}, {add, [4, {mul, [-1, {sym, d}]}]}]},
    lists:foldl(fun(T, Acc) ->
                        Fs = case T of {mul, F} -> F; F -> [F] end,
                        {Col, Rest} = lists:partition(fun({sym, S}) -> lists:member(S, ?COLOUR);
                                                         ({pow, {sym, S}, _}) -> lists:member(S, ?COLOUR);
                                                         (_) -> false end, Fs),
                        Mono = lists:sort([case X of {sym, S} -> {S, 1}; {pow, {sym, S}, E} -> {S, E} end || X <- Col]),
                        R = pq_symbolic:substitute({mul, [1 | Rest]}, #{ep => Eps}),
                        R2 = case R of {ok, X2} -> X2; X2 -> X2 end,
                        maps:update_with(Mono, fun(A) -> {add, [A, R2]} end, R2, Acc)
                end, #{}, Terms).

reduce(Dir, Field, Loops) ->
    {ok, _} = application:ensure_all_started(pq),
    Graphs = parse(filename:join(Dir, io_lib:format("dump_~s_~b.txt", [Field, Loops]))),
    T0 = erlang:monotonic_time(millisecond),
    Results = [begin
                   P = forcer_cmp:problem(Vx, N, #{}),
                   ByMono = lists:foldl(fun({CT, Pw, Sps}, Acc) ->
                                                maps:fold(fun(Mono, A, Acc2) ->
                                                                  maps:update_with(Mono, fun(L) -> [{A, Pw, Sps} | L] end, [{A, Pw, Sps}], Acc2)
                                                          end, Acc, split(CT))
                                        end, #{}, Terms),
                   Monos = maps:keys(ByMono),
                   Combs = [forcer_cmp:targets(P, maps:get(M, ByMono)) || M <- Monos],
                   Tg = erlang:monotonic_time(millisecond),
                   {ok, Reds, Info} = pq_ibp:reduce_sums(maps:get(family, P), Combs, #{}),
                   io:format("  graph ~s: ~b lines, ~b terms, ~b colour monomials, ~b equations, seeds ~w, ~.2f s~n",
                             [Vx, N, length(Terms), length(Monos), maps:get(equations, Info, 0), maps:get(seeds, Info, none),
                              (erlang:monotonic_time(millisecond) - Tg) / 1000]),
                   {P, lists:zip(Monos, Reds), Info}
               end || {Vx, N, Terms} <- Graphs],
    Secs = (erlang:monotonic_time(millisecond) - T0) / 1000,
    %% every graph's masters, for Forcer
    Masters = [{G, M} || {G, {P, MR, _}} <- lists:enumerate(Results), M <- lists:usort([M || {_, R} <- MR, {M, _} <- R]), _ <- [P]],
    ok = file:write_file(filename:join(Dir, io_lib:format("masters_~s_~b.frm", [Field, Loops])),
                         ["#-\n#include- forcer.h\n",
                          [io:format("") || _ <- []],
                          [io_lib:format("L M~b = ~s;~n", [K, forcer_cmp:integral_text(element(1, lists:nth(G, Results)), M)])
                           || {K, {G, M}} <- lists:enumerate(Masters)],
                          "#call Forcer\nB ep;\nP +s;\n.end\n"]),
    Index = maps:from_list([{{G, M}, K} || {K, {G, M}} <- lists:enumerate(Masters)]),
    Red = [{Mono, maps:get({G, M}, Index), pq_symbolic:pretty({mul, [pq_sym_norm:to_ast(Nu), {pow, pq_sym_norm:to_ast(De), -1}]})}
           || {G, {_, MR, _}} <- lists:enumerate(Results), {Mono, R} <- MR, {M, {Nu, De}} <- R],
    ok = file:write_file(filename:join(Dir, io_lib:format("~s_~b.red", [Field, Loops])), io_lib:format("~p.~n", [Red])),
    Eqs = lists:sum([maps:get(equations, I, 0) || {_, _, I} <- Results]),
    io:format("~s ~b loops: ~b graphs, ~b masters, ~b equations, ~.2f s~n", [Field, Loops, length(Results), length(Masters), Eqs, Secs]),
    {Secs, #{graphs => length(Results), masters => length(Masters), equations => Eqs}}.

%% the MS-bar series of P_L (the gluon's still times D - 1), from the reductions and Forcer's master values
compose(Dir, Field, Loops, Tables) ->
    {ok, [Red]} = file:consult(filename:join(Dir, io_lib:format("~s_~b.red", [Field, Loops]))),
    MV = forcer_cmp:form_exprs(filename:join(Dir, io_lib:format("masters_~s_~b.out", [Field, Loops]))),
    N = 8,
    Val = fun(K) -> {_, F} = epsexp:series_of_form(maps:get("M" ++ integer_to_list(K), MV), Tables), F(N) end,
    Vals = maps:from_list([{K, Val(K)} || K <- lists:usort([K || {_, K, _} <- Red])]),
    G = lists:foldl(fun({Mono, K, CoefText}, Acc) ->
                            {ok, Ast} = pq_symbolic:parse(CoefText),
                            C = epsexp:eval(Ast, Tables, N),
                            M = epsexp:exact(#{0 => #{Mono => {1, 1}}}),
                            epsexp:sadd(Acc, epsexp:smul(M, epsexp:smul(C, maps:get(K, Vals))))
                    end, epsexp:exact(#{}), Red),
    %% Forcer's msbarexpand=3: times the MS-bar factor per loop, orders with loops + power > 3 dropped
    epsexp:drop(epsexp:smul(G, epsexp:spow(epsexp:msbar(12), Loops, 12)), 3 - Loops + 1).
