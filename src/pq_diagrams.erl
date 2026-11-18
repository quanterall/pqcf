%% Feynman diagrams (plan "theoretical physics instruments", F4, first part): every graph of a theory with given
%% external legs and number of loops, with its weight 1/S, as Wick's theorem produces them.
%%
%% A theory is a list of vertex types {Name, Fields} (Fields the multiset of field names on the vertex, e.g.
%% {v4, [phi, phi, phi, phi]}) and the propagators {F1, F2} (which field ends can be joined; {phi, phi} for a real
%% scalar). External legs are labelled nodes {ext, 1}, {ext, 2}, .. carrying a field each.
%%
%% Method: for each multiset of vertices with sum_t n_t (k_t - 2) = 2 L - 2 + E (L loops, E legs), every complete
%% pairing of the half-edges (Wick's contractions) is generated, the resulting graph brought to a canonical form
%% (internal vertices relabelled within each type to the smallest edge list), and the pairings of each graph counted.
%% With a vertex normalised as g / prod(m_f!) (m_f the multiplicity of field f on it, the usual g/4! phi^4),
%%   1/S = (number of pairings) / prod_t (n_t! prod_f m_f!^n_t),
%% exactly: no automorphism group is computed. Filters: connected (default), one_pi (no internal bridge), no_tadpoles
%% (no part without external legs hanging on one line) and no_snails (no line from a vertex to itself), as in QGRAF.
%% Reach: the pairings are enumerated with labelled vertices, so the work grows with the number of labelled graphs:
%% phi^4 vacuum graphs to six loops (5 vertices) in 0.3 s, phi^3 to four loops (6 vertices) in 4 s; eight cubic
%% vertices are out of reach. Orderly generation (as QGRAF does) is the next step.
-module(pq_diagrams).
-export([graphs/4, graphs/5, sum_weights/1]).

-type field() :: atom().
-type vtype() :: {atom(), [field()]}.

%% graphs(Theory, Externals, Loops, Opts) -> [{Graph, {N, D}}], Graph = {Vertices, Edges}: Vertices a list of
%% {Index, TypeName}, Edges a sorted list of {{NodeA, FieldA}, {NodeB, FieldB}}
graphs(Theory, Ext, Loops, Opts) -> graphs(Theory, Ext, Loops, Opts, #{}).

-spec graphs(#{vertices := [vtype()], propagators := [{field(), field()}]}, [field()], non_neg_integer(), [atom()], map()) ->
          [{term(), {integer(), pos_integer()}}].
graphs(#{vertices := Vts, propagators := Props}, Ext, Loops, Opts, _More) ->
    E = length(Ext),
    Target = 2 * Loops - 2 + E,
    Counts = vertex_counts(Vts, Target),
    lists:append([graphs_for(Vts, Props, Ext, Ns, Opts) || Ns <- Counts]).

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
                                Key = canonical(Types, Edges),
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

%% the canonical form: internal vertices coloured by colour refinement (their type and legs, then repeatedly the
%% multiset of {own field, neighbour's field, neighbour's colour} over their edges, external legs by number), ordered
%% by colour and renumbered 1..n; within a colour class every order is tried and the smallest edge list kept. Colours
%% are invariant terms, so equal graphs get equal forms, and the classes are small.
canonical(Types, Edges) ->
    Col0 = maps:from_list([{{v, I}, {T, lists:sort(legs({v, I}, Edges))}} || {I, T} <- Types]),
    Col = refine(Col0, Edges, length(Types) + 1),
    Classes = group_by_colour(Col),
    Orders = orders(Classes),
    lists:min([relabel(Edges, maps:from_list([{I, K} || {K, {v, I}} <- lists:enumerate(Order)])) || Order <- Orders]).

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

orders([]) -> [[]];
orders([Ns | Rest]) -> [P ++ R || P <- perms(Ns), R <- orders(Rest)].

perms([]) -> [[]];
perms(L) -> [[X | P] || X <- L, P <- perms(L -- [X])].

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

norm(N, D) -> G = gcd(N, D), {N div G, D div G}.
gcd(A, 0) -> A;
gcd(A, B) -> gcd(B, A rem B).
qadd({A, B}, {C, D}) -> norm(A * D + C * B, B * D).
