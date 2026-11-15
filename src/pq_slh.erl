%% Quantum input-output networks in the SLH formalism (quantum optics plan, Q6): Gough and James, IEEE TAC 54, 2530
%% (2009), arXiv:0707.0048; Combes, Kerckhoff, Sarovar, Adv. Phys. X 2, 784 (2017), arXiv:1611.00375.
%%
%% A component with n input and n output channels is G = (S, L, H): S an n x n scattering matrix, L an n-vector of
%% coupling operators, H a Hamiltonian; every entry a pq_qop expression (scalars as pq_qop:num/1), so networks reduce
%% exactly and symbolically. The compositions:
%%   concatenation  G1 + G2 = (diag(S1, S2), [L1; L2], H1 + H2)
%%   series         G2 <| G1 = (S2 S1, L2 + S2 L1, H1 + H2 + (1/2i)(L2^dagger S2 L1 - L1^dagger S2^dagger L2))
%%   feedback       [G]_(k -> l), output k fed into input l (S_kl a scalar with 1 - S_kl invertible):
%%       S' = S_[k][l] + S_[k]l (1 - S_kl)^-1 S_k[l],   L' = L_[k] + S_[k]l (1 - S_kl)^-1 L_k,
%%       H' = H + (1/2i) (sum_j L_j^dagger S_jl (1 - S_kl)^-1 L_k - h.c.)
%% ([k] leaves out row k, [l] column l). The reduced network's master equation has H and the L as jump operators.
%% pq_qop's spaces must be set by the caller (pq_qop:with_spaces/2) when atoms are involved. Channels are 1-based.
-module(pq_slh).
-export([new/3, identity/1, static/1, coherent/1, concat/2, series/2, feedback/3, master/1, size/1]).

-define(Q, pq_qop).

new(S, L, H) -> #{s => S, l => L, h => H}.

size(#{l := L}) -> length(L).

%% n channels passing through
identity(N) -> new([[?Q:num(case I =:= J of true -> 1; false -> 0 end) || J <- lists:seq(1, N)] || I <- lists:seq(1, N)],
                   lists:duplicate(N, ?Q:num(0)), ?Q:num(0)).

%% a static linear-optical element: S a matrix of scalars (engine expressions), no coupling, no Hamiltonian
static(Rows) -> new([[?Q:num(X) || X <- R] || R <- Rows], [?Q:num(0) || _ <- Rows], ?Q:num(0)).

%% a coherent input field of amplitude Alpha on one channel: (1, Alpha, 0)
coherent(Alpha) -> new([[?Q:num(1)]], [?Q:num(Alpha)], ?Q:num(0)).

concat(#{s := S1, l := L1, h := H1}, #{s := S2, l := L2, h := H2}) ->
    N1 = length(L1), N2 = length(L2),
    Z = ?Q:num(0),
    new([R ++ lists:duplicate(N2, Z) || R <- S1] ++ [lists:duplicate(N1, Z) ++ R || R <- S2], L1 ++ L2, ?Q:add(H1, H2)).

%% G2 after G1
series(#{s := S2, l := L2, h := H2}, #{s := S1, l := L1, h := H1}) ->
    S = matmul(S2, S1),
    L = vadd(L2, matvec(S2, L1)),
    Cross = ?Q:sub(dot(L2, matvec(S2, L1)), dot(L1, matvec(dagger(S2), L2))),
    new(S, L, ?Q:add([H1, H2, ?Q:scale({mul, [{rat, 1, 2}, {complex, 0, -1}]}, Cross)])).

%% output K fed into input Lc
feedback(#{s := S, l := L, h := H}, K, Lc) ->
    Skl = scalar(at(S, K, Lc)),
    Inv = {pow, {add, [1, {mul, [-1, Skl]}]}, -1},
    Rows = drop(S, K),
    ColL = [lists:nth(Lc, R) || R <- Rows],                         % S_[k] l
    RowK = drop(lists:nth(K, S), Lc),                               % S_k [l]
    S2 = [[?Q:add(X, ?Q:scale(Inv, ?Q:mul(C, Y))) || {X, Y} <- lists:zip(drop(R, Lc), RowK)] || {R, C} <- lists:zip(Rows, ColL)],
    Lk = lists:nth(K, L),
    L2 = [?Q:add(Lj, ?Q:scale(Inv, ?Q:mul(C, Lk))) || {Lj, C} <- lists:zip(drop(L, K), ColL)],
    T = ?Q:scale(Inv, ?Q:add([?Q:mul([?Q:dag(Lj), at(S, J, Lc), Lk]) || {J, Lj} <- lists:enumerate(L)])),
    new(S2, L2, ?Q:add(H, ?Q:scale({mul, [{rat, 1, 2}, {complex, 0, -1}]}, ?Q:sub(T, ?Q:dag(T))))).

%% the master equation of the network: {H, [L_k]} (the zero couplings left out)
master(#{l := L, h := H}) -> {H, [X || X <- L, not ?Q:is_zero(X)]}.

%%% matrices of operators

at(M, I, J) -> lists:nth(J, lists:nth(I, M)).

drop(List, I) -> {A, [_ | B]} = lists:split(I - 1, List), A ++ B.

scalar(E) ->
    case maps:to_list(E) of
        [] -> 0;
        [{[], C}] -> C;
        _ -> error({pq_slh, {not_a_scalar, E}})
    end.

matmul(A, B) ->
    Bt = transpose(B),
    [[?Q:add([?Q:mul(X, Y) || {X, Y} <- lists:zip(R, C)]) || C <- Bt] || R <- A].

matvec(A, V) -> [?Q:add([?Q:mul(X, Y) || {X, Y} <- lists:zip(R, V)]) || R <- A].

vadd(A, B) -> [?Q:add(X, Y) || {X, Y} <- lists:zip(A, B)].

%% sum_i A_i^dagger B_i
dot(A, B) -> ?Q:add([?Q:mul(?Q:dag(X), Y) || {X, Y} <- lists:zip(A, B)]).

dagger(M) -> [[?Q:dag(X) || X <- R] || R <- transpose(M)].

transpose([[] | _]) -> [];
transpose(M) -> [[hd(R) || R <- M] | transpose([tl(R) || R <- M])].
