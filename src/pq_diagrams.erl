%% Feynman diagrams (plan "theoretical physics instruments", F4): every graph of a theory with given external legs and
%% number of loops, with its weight sign/S.
%%
%% A theory: #{vertices => [{Name, Fields}], propagators => [{F1, F2}], fermions => [F]}. Fields is the multiset of
%% field names on a vertex (e.g. {v4, [phi, phi, phi, phi]}, QED's {v, ['ELE', ele, pho]}); a propagator {F1, F2}
%% joins a leg F1 to a leg F2 ({phi, phi} for a real field, {ele, 'ELE'} for a charged one, as in FORM's and QGRAF's
%% models); fermions are the Grassmann fields (fermions and ghosts). External legs are nodes {ext, 1}, {ext, 2}, ..
%% carrying one field each (the field of the leg they join, as FORM writes them).
%%
%% Method (the default, orderly): for each multiset of vertices with sum_t n_t (k_t - 2) = 2 L - 2 + E, the multigraphs
%% are built by joining the first free leg to a partner, partners of one leg class in increasing order, and a vertex
%% not yet joined to anything only if it is the first such vertex of its type: each graph comes a few times instead of
%% once per labelling. Each is brought to a canonical form by individualisation and refinement (colour refinement of
%% the vertices by their type, legs and neighbours; a vertex of the first ambiguous class fixed in turn); the leaves of
%% that search that give the smallest edge list are its automorphisms. Then
%%   S = |Aut| prod (multiplicity of each multiple line)! 2^(lines from a vertex to itself with one field at both ends),
%%   sign = (-1)^(closed fermion loops)   (a vertex carries at most one fermion line).
%% method => wick counts every pairing of the legs instead (Wick's theorem, 1/S = pairings / prod_t (n_t! prod_f
%% m_f!^n_t)): slow, kept as an independent check for theories without fermions.
%% Filters, as in QGRAF and FORM: connected (always), one_pi, no_tadpoles (no part without external legs hanging on
%% one line), no_snails (no line from a vertex to itself).
-module(pq_diagrams).
-export([graphs/4, graphs/5, sum_weights/1, fermion_loops/2]).

-type field() :: atom().
-type vtype() :: {atom(), [field()]}.

%% graphs(Theory, Externals, Loops, Opts) -> [{Graph, {N, D}}], Graph = {Vertices, Edges}: Vertices a list of
%% {Index, TypeName}, Edges a sorted list of {{NodeA, FieldA}, {NodeB, FieldB}}
graphs(Theory, Ext, Loops, Opts) -> graphs(Theory, Ext, Loops, Opts, #{}).

-spec graphs(#{vertices := [vtype()], propagators := [{field(), field()}]}, [field()], non_neg_integer(), [atom()], map()) ->
          [{term(), {integer(), pos_integer()}}].
graphs(#{vertices := Vts, propagators := Props} = Th, Ext, Loops, Opts, More) ->
    E = length(Ext),
    Target = 2 * Loops - 2 + E,
    Counts = vertex_counts(Vts, Target),
    case maps:get(method, More, orderly) of
        wick ->
            maps:get(fermions, Th, []) =:= [] orelse error({pq_diagrams, wick_without_fermions_only}),
            lists:append([graphs_for(Vts, Props, Ext, Ns, Opts) || Ns <- Counts]);
        orderly ->
            Fer = sets:from_list(maps:get(fermions, Th, []), [{version, 2}]),
            lists:append([orderly(Vts, Props, Fer, Ext, Ns, Opts) || Ns <- Counts])
    end.

%% the vectors n_t >= 0 with sum n_t (k_t - 2) = Target (k_t > 2 for every type; k_t = 2 would give infinitely many)
vertex_counts(Vts, Target) ->
    Ks = [length(Fs) - 2 || {_, Fs} <- Vts],
    lists:all(fun(K) -> K > 0 end, Ks) orelse error({pq_diagrams, two_point_vertices_unsupported}),
    counts(Ks, Target).

counts([], 0) -> [[]];
counts([], _) -> [];
counts([K | Ks], T) when T >= 0 -> [[N | Rest] || N <- lists:seq(0, T div K), Rest <- counts(Ks, T - N * K)];
counts(_, _) -> [].

graphs_for(Vts, Props, Ext, Ns, Opts) ->
    % the internal vertices, numbered within the whole graph
    Vlist = lists:append([lists:duplicate(N, T) || {N, T} <- lists:zip(Ns, Vts)]),
    Verts = lists:enumerate(Vlist),
    Half = [{{ext, I}, F} || {I, F} <- lists:enumerate(Ext)] ++
           [{{v, I}, F} || {I, {_, Fs}} <- Verts, F <- Fs],
    Join = joinable(Props),
    Types = [{I, Name} || {I, {Name, _}} <- Verts],
    Pairings = pairings(Half, Join, #{}),
    Classes = maps:fold(fun(Edges, C, Acc) ->
                                {Key, _} = canonical(Types, Edges),
                                Acc#{Key => maps:get(Key, Acc, 0) + C}
                        end, #{}, Pairings),
    Den = lists:foldl(fun({N, {_, Fs}}, D) -> D * fact(N) * pow(mult_fact(Fs), N) end, 1, lists:zip(Ns, Vts)),
    [{{Types, Edges}, norm(C, Den)} || Edges := C <- Classes, keep({Types, Edges}, Opts)].

joinable(Props) ->
    lists:foldl(fun({A, B}, S) -> sets:add_element({A, B}, sets:add_element({B, A}, S)) end, sets:new([{version, 2}]), Props).

%% every complete pairing of the half-edges, counted by its edge multiset (sorted list of sorted pairs). Legs of one
%% vertex with the same field are interchangeable, so the pairing works on classes {Node, Field} with the number of
%% free legs in each: the first leg of the first class is paired with each joinable class, weighted by its count.
pairings(Half, Join, _) ->
    Free = lists:foldl(fun(H, M) -> M#{H => maps:get(H, M, 0) + 1} end, #{}, Half),
    pair_classes(lists:sort(maps:to_list(Free)), Join, [], 1, #{}).

pair_classes([], _, Edges, W, Acc) ->
    K = lists:sort(Edges),
    Acc#{K => maps:get(K, Acc, 0) + W};
pair_classes([{H, N} | Rest0], Join, Edges, W, Acc) ->
    {_, Fh} = H,
    Rest = case N - 1 of 0 -> Rest0; M -> [{H, M} | Rest0] end,
    lists:foldl(fun({I, {X, C}}, A) ->
                        {_, Fx} = X,
                        case sets:is_element({Fh, Fx}, Join) of
                            false -> A;
                            true ->
                                Rest2 = case C - 1 of
                                            0 -> delete_nth(Rest, I);
                                            C1 -> replace_nth(Rest, I, {X, C1})
                                        end,
                                pair_classes(Rest2, Join, [edge(H, X) | Edges], W * C, A)
                        end
                end, Acc, lists:enumerate(Rest)).

replace_nth(L, K, V) -> {A, [_ | B]} = lists:split(K - 1, L), A ++ [V | B].

edge(A, B) when A =< B -> {A, B};
edge(A, B) -> {B, A}.

%% the canonical form and the number of automorphisms: internal vertices coloured by colour refinement (type and legs,
%% then repeatedly the multiset of {own field, neighbour's field, neighbour's colour}; external legs keep their
%% numbers); while a colour class has more than one vertex, each of its vertices is individualised in turn and the
%% colours refined again. Each leaf orders the vertices; the smallest relabelled edge list is the form, and the leaves
%% that give it differ by automorphisms (colours are invariant, so equal graphs get equal forms).
canonical(Types, Edges) ->
    Col0 = maps:from_list([{{v, I}, {T, lists:sort(legs({v, I}, Edges))}} || {I, T} <- Types]),
    Leaves = search(refine(Col0, Edges, length(Types) + 1), Edges, length(Types) + 1),
    Forms = [relabel(Edges, maps:from_list([{I, K} || {K, {v, I}} <- lists:enumerate(Order)])) || Order <- Leaves],
    Min = lists:min(Forms),
    {Min, length([F || F <- Forms, F =:= Min])}.

search(Col, Edges, N) ->
    case [Ns || Ns <- group_by_colour(Col), length(Ns) > 1] of
        [] -> [lists:append(group_by_colour(Col))];
        [Cell | _] ->
            lists:append([search(refine(individualise(Col, V), Edges, N), Edges, N) || V <- Cell])
    end.

%% V gets a colour of its own, ahead of the rest of its class
individualise(Col, V) ->
    C = maps:get(V, Col),
    maps:map(fun(X, Cx) when Cx =:= C -> case X of V -> {Cx, 0}; _ -> {Cx, 1} end; (_, Cx) -> {Cx, 1} end, Col).

legs(N, Edges) ->
    [F || {{A, F}, _} <- Edges, A =:= N] ++ [F || {_, {B, F}} <- Edges, B =:= N].

refine(Col, _, 0) -> Col;
refine(Col, Edges, K) ->
    C = fun({ext, _} = X) -> X; (X) -> maps:get(X, Col) end,
    New = maps:map(fun(N, Old) ->
                           Nb = [{Fa, Fb, C(B)} || {{A, Fa}, {B, Fb}} <- Edges, A =:= N] ++
                                [{Fb, Fa, C(A)} || {{A, Fa}, {B, Fb}} <- Edges, B =:= N],
                           {Old, lists:sort(Nb)}
                   end, Col),
    case length(lists:usort(maps:values(New))) =:= length(lists:usort(maps:values(Col))) of
        true -> Col;
        false -> refine(New, Edges, K - 1)
    end.

%% the vertices grouped by colour, the groups in colour order
group_by_colour(Col) ->
    G = maps:fold(fun(N, C, M) -> M#{C => [N | maps:get(C, M, [])]} end, #{}, Col),
    [lists:sort(Ns) || {_, Ns} <- lists:sort(maps:to_list(G))].

%%% Orderly generation

orderly(Vts, Props, Fer, Ext, Ns, Opts) ->
    Vlist = lists:append([lists:duplicate(N, T) || {N, T} <- lists:zip(Ns, Vts)]),
    Verts = lists:enumerate(Vlist),
    Types = [{I, Name} || {I, {Name, _}} <- Verts],
    Nodes = [{ext, I} || {I, _} <- lists:enumerate(Ext)] ++ [{v, I} || {I, _} <- Verts],
    Legs0 = maps:from_list([{{ext, I}, count([F])} || {I, F} <- lists:enumerate(Ext)] ++
                           [{{v, I}, count(Fs)} || {I, {_, Fs}} <- Verts]),
    Full = maps:from_list([{{v, I}, {Name, count(Fs)}} || {I, {Name, Fs}} <- Verts]),
    Join = joinable(Props),
    Raw = build(Nodes, Legs0, Full, Join, none, [], []),
    Uniq = lists:foldl(fun(Edges, Acc) ->
                               G = {Types, Edges},
                               case keep(G, Opts) of
                                   false -> Acc;
                                   true ->
                                       {Key, Aut} = canonical(Types, Edges),
                                       case maps:is_key(Key, Acc) of
                                           true -> Acc;
                                           false -> Acc#{Key => Aut}
                                       end
                               end
                       end, #{}, Raw),
    [{{Types, Edges}, weight(Edges, Aut, Fer)} || Edges := Aut <- Uniq].

count(Fs) -> lists:foldl(fun(F, M) -> M#{F => maps:get(F, M, 0) + 1} end, #{}, Fs).

%% Last = {Node, Field, Partner}: the last partner chosen for that leg class (partners of one class go in order)
build(Nodes, Legs, Full, Join, Last, Edges, Acc) ->
    case first_free(Nodes, Legs) of
        none -> [Edges | Acc];
        {N, F} ->
            Min = case Last of {N, F, P} -> P; _ -> none end,
            Cands = [{M, G} || M <- lists:dropwhile(fun(X) -> X =/= N end, Nodes),
                               {G, C} <- lists:sort(maps:to_list(maps:get(M, Legs))), C > 0,
                               sets:is_element({F, G}, Join),
                               not (M =:= N andalso G =:= F andalso C < 2),
                               Min =:= none orelse {M, G} >= Min,
                               first_untouched(M, N, Nodes, Legs, Full)],
            lists:foldl(fun({M, G} = P, A) ->
                                L1 = take(take(Legs, N, F), M, G),
                                build(Nodes, L1, Full, Join, {N, F, P}, [edge({N, F}, {M, G}) | Edges], A)
                        end, Acc, Cands)
    end.

first_free([], _) -> none;
first_free([N | Rest], Legs) ->
    case lists:sort([F || {F, C} <- maps:to_list(maps:get(N, Legs)), C > 0]) of
        [] -> first_free(Rest, Legs);
        [F | _] -> {N, F}
    end.

take(Legs, N, F) -> M = maps:get(N, Legs), Legs#{N => M#{F => maps:get(F, M) - 1}}.

%% a vertex untouched so far may be joined only if no untouched vertex of its type comes before it
first_untouched({ext, _}, _, _, _, _) -> true;
first_untouched(M, M, _, _, _) -> true;
first_untouched(M, Cur, Nodes, Legs, Full) ->
    {Type, Legs0} = maps:get(M, Full),
    case maps:get(M, Legs) =:= Legs0 of
        false -> true;
        true ->
            % (the current node N does not count: its leg is being joined now)
            Before = [X || X <- lists:takewhile(fun(X) -> X =/= M end, Nodes), X =/= Cur],
            not lists:any(fun({v, _} = X) -> {T, L0} = maps:get(X, Full), T =:= Type andalso maps:get(X, Legs) =:= L0;
                             (_) -> false end, Before)
    end.

%% sign / S
weight(Edges, Aut, Fer) ->
    Mult = lists:foldl(fun(E, M) -> M#{E => maps:get(E, M, 0) + 1} end, #{}, Edges),
    S = maps:fold(fun({{A, F}, {B, G}}, K, P) ->
                          Self = case A =:= B andalso F =:= G of true -> pow(2, K); false -> 1 end,
                          P * fact(K) * Self
                  end, Aut, Mult),
    Sign = case fermion_loops(Edges, Fer) rem 2 of 0 -> 1; 1 -> -1 end,
    norm(Sign, S).

%% the closed fermion loops of a graph {Types, Edges} of a theory (exported for weights with N flavours)
fermion_loops({_, Edges}, Th) when is_map(Th) ->
    fermion_loops(Edges, sets:from_list(maps:get(fermions, Th, []), [{version, 2}]));

%% the closed fermion loops: components of the fermion lines (a vertex joins its two fermion legs) without an external
%% leg
fermion_loops(Edges, Fer) ->
    FE = [{A, B} || {{A, F}, {B, G}} <- Edges, sets:is_element(F, Fer), sets:is_element(G, Fer)],
    case FE of
        [] -> 0;
        _ ->
            Nodes = lists:usort(lists:append([[A, B] || {A, B} <- FE])),
            lists:foreach(fun(N) ->
                                  Deg = length([x || {A, _} <- FE, A =:= N]) + length([x || {_, B} <- FE, B =:= N]),
                                  case N of
                                      {v, _} when Deg =/= 2 -> error({pq_diagrams, {vertex_with_fermions, Deg}});
                                      _ -> ok
                                  end
                          end, Nodes),
            Comps = components(Nodes, [{{A, x}, {B, x}} || {A, B} <- FE]),
            length([C || C <- Comps, not lists:any(fun({ext, _}) -> true; (_) -> false end, sets:to_list(C))])
    end.

relabel(Edges, Map) ->
    R = fun({{v, I}, F}) -> {{v, maps:get(I, Map)}, F}; (X) -> X end,
    lists:sort([edge(R(A), R(B)) || {A, B} <- Edges]).

%%% Filters

keep(G, Opts) ->
    connected(G) andalso (not lists:member(one_pi, Opts) orelse one_pi(G))
        andalso (not lists:member(no_tadpoles, Opts) orelse no_tadpoles(G))
        andalso (not lists:member(no_snails, Opts) orelse no_snails(G)).

%% no snails (QGRAF's nosnail): no line from a vertex to itself
no_snails({_, Edges}) -> not lists:any(fun({{A, _}, {B, _}}) -> A =:= B end, Edges).

graph_nodes({Types, Edges}) ->
    lists:usort([{v, I} || {I, _} <- Types] ++ [N || {{N, _}, _} <- Edges] ++ [M || {_, {M, _}} <- Edges]).

connected({_, []}) -> true;
connected({_, Edges} = G) -> reach(Edges, [hd(graph_nodes(G))], sets:new([{version, 2}])) =:= length(graph_nodes(G)).

reach(_, [], Seen) -> sets:size(Seen);
reach(Edges, [N | Todo], Seen) ->
    case sets:is_element(N, Seen) of
        true -> reach(Edges, Todo, Seen);
        false ->
            Next = [B || {{A, _}, {B, _}} <- Edges, A =:= N] ++ [A || {{A, _}, {B, _}} <- Edges, B =:= N],
            reach(Edges, Next ++ Todo, sets:add_element(N, Seen))
    end.

%% one-particle irreducible: removing any one edge between internal vertices leaves the internal part connected
one_pi({Types, Edges}) ->
    Internal = [E || {{{v, _}, _}, {{v, _}, _}} = E <- Edges],
    Ivs = [{v, I} || {I, _} <- Types],
    Ivs =:= [] orelse
        lists:all(fun(K) ->
                          Rest = delete_nth(Internal, K),
                          reach(Rest, [hd(Ivs)], sets:new([{version, 2}])) =:= length(Ivs)
                  end, lists:seq(1, length(Internal))).

delete_nth(L, K) -> {A, [_ | B]} = lists:split(K - 1, L), A ++ B.

%% no tadpoles: no internal part connected to the rest by a single edge with no external leg on it, nor a self-loop
%% hanging from it... in the usual sense: no subgraph without external legs attached by one line
no_tadpoles({Types, Edges}) ->
    Internal = [E || {{{v, _}, _}, {{v, _}, _}} = E <- Edges],
    Ivs = [{v, I} || {I, _} <- Types],
    Exts = [E || E <- Edges, not lists:member(E, Internal)],
    lists:all(fun(K) ->
                      Rest = delete_nth(Internal, K),
                      Comps = components(Ivs, Rest),
                      length(Comps) =:= 1 orelse
                          lists:all(fun(C) -> has_external(C, Exts) end, Comps)
              end, lists:seq(1, length(Internal))).

components([], _) -> [];
components([N | _] = Ns, Edges) ->
    C = reach_set(Edges, [N], sets:new([{version, 2}])),
    [C | components([M || M <- Ns, not sets:is_element(M, C)], Edges)].

reach_set(_, [], Seen) -> Seen;
reach_set(Edges, [N | Todo], Seen) ->
    case sets:is_element(N, Seen) of
        true -> reach_set(Edges, Todo, Seen);
        false ->
            Next = [B || {{A, _}, {B, _}} <- Edges, A =:= N] ++ [A || {{A, _}, {B, _}} <- Edges, B =:= N],
            reach_set(Edges, Next ++ Todo, sets:add_element(N, Seen))
    end.

has_external(C, Exts) ->
    lists:any(fun({{A, _}, {B, _}}) -> sets:is_element(A, C) orelse sets:is_element(B, C) end, Exts).

%%% Weights

%% the sum of 1/S over a list of graphs
sum_weights(Gs) -> lists:foldl(fun({_, W}, Acc) -> qadd(Acc, W) end, {0, 1}, Gs).

mult_fact(Fs) ->
    M = lists:foldl(fun(F, A) -> A#{F => maps:get(F, A, 0) + 1} end, #{}, Fs),
    lists:foldl(fun(K, P) -> P * fact(K) end, 1, maps:values(M)).

fact(0) -> 1;
fact(N) -> N * fact(N - 1).

pow(_, 0) -> 1;
pow(X, K) -> X * pow(X, K - 1).

%% a fraction with its denominator positive
norm(N, D) -> G = gcd(abs(N), abs(D)), S = case D < 0 of true -> -1; false -> 1 end, {S * N div G, S * D div G}.
gcd(A, 0) -> A;
gcd(A, B) -> gcd(B, A rem B).
qadd({A, B}, {C, D}) -> norm(A * D + C * B, B * D).
