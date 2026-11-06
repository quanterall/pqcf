%% Prime fields for reconstruction (decision record R1): primes below 2^63, compiled programs evaluated modulo them
%% (pq_nif:ff_prepare/3, ff_eval/3, native/pq_nif/src/ffslp.rs), and the number theory that brings residues back to
%% the rationals: inverses, Chinese remaindering, and rational reconstruction.
-module(pq_ff).
-export([prime/1, is_prime/1, prepare/1, eval/3, inv/2, crt/2, ratrecon/2]).

-define(BAD, 18446744073709551615).

%%% Primes

%% the K-th prime (K >= 1) below 2^63, from the top: the same list on every node (D1)
-spec prime(pos_integer()) -> pos_integer().
prime(K) when is_integer(K), K >= 1 ->
    Key = {?MODULE, primes},
    Known = try persistent_term:get(Key) catch error:badarg -> [] end,
    case length(Known) >= K of
        true -> lists:nth(K, Known);
        false ->
            Start = case Known of [] -> 1 bsl 63; _ -> lists:last(Known) end,
            More = next_primes(Start, K - length(Known)),
            persistent_term:put(Key, Known ++ More),
            lists:nth(K, Known ++ More)
    end.

next_primes(_, 0) -> [];
next_primes(Below, N) ->
    P = prev_prime(Below - 1),
    [P | next_primes(P, N - 1)].

prev_prime(N) ->
    case is_prime(N) of
        true -> N;
        false -> prev_prime(N - 1)
    end.

%% deterministic Miller-Rabin: these bases decide every number below 3.3 * 10^24
-spec is_prime(non_neg_integer()) -> boolean().
is_prime(N) when N < 2 -> false;
is_prime(N) ->
    Bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41],
    case [B || B <- Bases, N rem B =:= 0] of
        [B] when B =:= N -> true;
        [] ->
            {S, D} = split(N - 1, 0),
            lists:all(fun(A) -> witness(A, S, D, N) end, Bases);
        _ -> false
    end.

split(D, S) when D band 1 =:= 0 -> split(D bsr 1, S + 1);
split(D, S) -> {S, D}.

witness(A, S, D, N) ->
    X = powm(A, D, N),
    X =:= 1 orelse X =:= N - 1 orelse squares(X, S - 1, N).

squares(_, 0, _) -> false;
squares(X, K, N) ->
    Y = X * X rem N,
    Y =:= N - 1 orelse (Y =/= 1 andalso squares(Y, K - 1, N)).

powm(_, 0, M) -> 1 rem M;
powm(B, E, M) when E band 1 =:= 1 -> B * powm(B * B rem M, E bsr 1, M) rem M;
powm(B, E, M) -> powm(B * B rem M, E bsr 1, M).

%%% Programs over a prime field

%% A compiled program (pq_symbolic:compile/3 with exact_constants => true) prepared for prime fields: every
%% operation exact over the rationals -- sums, products, integer powers, rational constants. Refused: pi,
%% radicals and other float constants, complex constants, fractional powers, exponentials, conjugates.
-spec prepare(tuple()) -> {ok, {ff, reference(), [atom()], non_neg_integer()}} | {error, term()}.
prepare({slp, 1, Inputs, Ops, Outputs}) ->
    Index = maps:from_list([{X, I} || {I, X} <- lists:enumerate(0, Inputs)]),
    Outs = case Outputs of {matrix, _, _, Ss} -> Ss; Ss -> Ss end,
    try
        lists:foldl(fun(Op, N) -> element(2, Op) =:= N orelse throw({unsupported, slot_order}), N + 1 end, 0, Ops),
        FOps = [ff_op(Op, Index) || Op <- Ops],
        case pq_nif:ff_prepare(length(Inputs), FOps, Outs) of
            {ok, R} -> {ok, {ff, R, Inputs, length(Outs)}};
            Error -> Error
        end
    catch
        throw:{unsupported, _} = U -> {error, {symbolic_unsupported, U}}
    end.

ff_op({input, _, X}, Index) -> {input, maps:get(X, Index)};
ff_op({const, _, {q, {N, D}, {0, 1}}}, _) -> {const, sign(N), limbs(abs(N)), limbs(D)};
ff_op({const, _, Z}, _) -> throw({unsupported, {constant, Z}});
ff_op({mul, _, A, B}, _) -> {mul, A, B};
ff_op({add, _, A, B}, _) -> {add, A, B};
ff_op({sub, _, A, B}, _) -> {sub, A, B};
ff_op({pow, _, A, {E, 1}}, _) -> {pow, A, E};
ff_op(Op, _) -> throw({unsupported, element(1, Op)}).

sign(N) when N < 0 -> -1;
sign(_) -> 1.

%% a non-negative integer as 64-bit limbs, least significant first (any size: the constants of a result can be
%% far beyond 128 bits)
limbs(0) -> [0];
limbs(N) -> limbs1(N).
limbs1(0) -> [];
limbs1(N) -> [N band 16#FFFFFFFFFFFFFFFF | limbs1(N bsr 64)].

%% The outputs at every point (a list of coordinate lists, one residue per input, in the program's input order),
%% mod P: [Values | bad] per point, bad where the program divides by zero mod P.
-spec eval({ff, reference(), [atom()], non_neg_integer()}, pos_integer(), [[non_neg_integer()]]) ->
          {ok, [[non_neg_integer()] | bad]} | bad_prime | {error, term()}.
eval({ff, R, Inputs, NOut}, P, Points) ->
    NIn = length(Inputs),
    case lists:all(fun(X) -> length(X) =:= NIn end, Points) of
        false -> {error, {bad_shape, points}};
        true ->
            case pq_nif:ff_eval(R, P, lists:append(Points)) of
                {ok, Flat} -> {ok, rows(Flat, max(NOut, 1), NOut)};
                Other -> Other
            end
    end.

rows([], _, _) -> [];
rows(Flat, Step, NOut) ->
    {Row, Rest} = lists:split(Step, Flat),
    [case Row of [?BAD | _] -> bad; _ -> lists:sublist(Row, NOut) end | rows(Rest, Step, NOut)].

%%% Number theory

%% the inverse of A mod M (A and M coprime)
-spec inv(integer(), pos_integer()) -> non_neg_integer().
inv(A, M) ->
    {G, X, _} = egcd(mod(A, M), M),
    G =:= 1 orelse error({not_invertible, A, M}),
    mod(X, M).

egcd(0, B) -> {B, 0, 1};
egcd(A, B) ->
    {G, X, Y} = egcd(B rem A, A),
    {G, Y - (B div A) * X, X}.

mod(A, M) -> ((A rem M) + M) rem M.

%% the residue mod M1 M2 that is A1 mod M1 and A2 mod M2 (M1, M2 coprime)
-spec crt({integer(), pos_integer()}, {integer(), pos_integer()}) -> {non_neg_integer(), pos_integer()}.
crt({A1, M1}, {A2, M2}) ->
    K = mod((A2 - A1) * inv(M1, M2), M2),
    {mod(A1 + M1 * K, M1 * M2), M1 * M2}.

%% Wang's rational reconstruction: the fraction N / D with |N|, D =< sqrt(M / 2) that is A mod M, as {N, D} in
%% lowest terms with D > 0, or false if there is none
-spec ratrecon(integer(), pos_integer()) -> {integer(), pos_integer()} | false.
ratrecon(A, M) ->
    Bound = isqrt(M div 2),
    wang(M, 0, mod(A, M), 1, Bound, M).

wang(R0, T0, R1, T1, Bound, M) when R1 > Bound ->
    Q = R0 div R1,
    wang(R1, T1, R0 - Q * R1, T0 - Q * T1, Bound, M);
wang(_R0, _T0, R1, T1, Bound, _M) ->
    case T1 =/= 0 andalso abs(T1) =< Bound andalso pq_sym_exact:gcd(R1, abs(T1)) =:= 1 of
        true when T1 < 0 -> {-R1, -T1};
        true -> {R1, T1};
        false -> false
    end.

isqrt(0) -> 0;
isqrt(N) -> isqrt(N, N).
isqrt(N, X) ->
    Y = (X + N div X) div 2,
    case Y >= X of true -> X; false -> isqrt(N, Y) end.
