%% The JSON service (1.2; docs/decisions/C1-python-client.md): PQCF's facade behind JSON requests, for clients in
%% other languages (python/pqcf). One request is a JSON object
%%   {"id": any, "op": "version" | "validate" | "backend_for" | "run" | "sample" | "probabilities" | "evolve" |
%%           "gradient", "circuit": <Appendix B document>, "options": {...}, ...}
%% and one response {"id": the same, "ok": result} or {"id": ..., "error": {"category": ..., "detail": ...}}.
%%
%% Atom safety: option keys and enumerated values are matched against fixed tables, never made into atoms; symbol
%% names in bindings and observables are matched against the circuit's own symbols (atoms the JSON reader already
%% checked and created); predicates and circuits go through pq_json. Nothing a request contains can crash the
%% service: every failure is an error response.
%%
%% Transports: stdio/0, one request per line on standard input and one response per line on standard output (the
%% Python client starts `erl ... -s pq_service stdio` as a subprocess: no network); and start_http/1,2, POST /api
%% on OTP's inets httpd, bound to 127.0.0.1 unless told otherwise.
-module(pq_service).
-include("pq.hrl").
-include_lib("inets/include/httpd.hrl").
-export([handle/1, request/1, stdio/0, http/0, start_http/1, start_http/2, stop_http/1]).
%% inets httpd callback
-export([do/1]).

-define(SERVICE_VERSION, 1).
-define(MAX_REQUEST, 64 * 1024 * 1024).
-define(MAX_VALUES, 4000000).

%%% Requests

%% one encoded request to one encoded response
-spec handle(binary()) -> binary().
handle(Bin) when byte_size(Bin) > ?MAX_REQUEST ->
    encode(#{<<"id">> => null, <<"error">> => error_json({bad_term, request_too_large})});
handle(Bin) ->
    Resp = try json:decode(Bin) of
               #{} = Req -> request(Req);
               _ -> #{<<"id">> => null, <<"error">> => error_json({bad_term, not_an_object})}
           catch
               _:_ -> #{<<"id">> => null, <<"error">> => error_json({bad_term, invalid_json})}
           end,
    encode(Resp).

%% one decoded request to one decoded response
-spec request(map()) -> map().
request(Req) ->
    Id = maps:get(<<"id">>, Req, null),
    try op(maps:get(<<"op">>, Req, missing), Req) of
        Result -> #{<<"id">> => Id, <<"ok">> => Result}
    catch
        throw:{error, R} -> #{<<"id">> => Id, <<"error">> => error_json(R)};
        Class:Reason:Stack ->
            logger:warning("pq_service: ~p:~p ~p", [Class, Reason, Stack]),
            #{<<"id">> => Id, <<"error">> => error_json({pq_bug, #{class => Class}})}
    end.

op(<<"version">>, _Req) ->
    {Ma, Mi, Pa, Triple} = pq_nif:nif_version(),
    {ok, Vsn} = application:get_key(pq, vsn),
    #{<<"pq">> => list_to_binary(Vsn), <<"nif">> => [Ma, Mi, Pa], <<"triple">> => Triple,
      <<"otp">> => list_to_binary(erlang:system_info(otp_release)), <<"service">> => ?SERVICE_VERSION};
op(<<"validate">>, Req) ->
    C = circuit(Req),
    ok(pq_circuit:validate(C)),
    #{<<"valid">> => true, <<"free_symbols">> => [atom_to_binary(S) || S <- pq_circuit:free_symbols(C)]};
op(<<"backend_for">>, Req) ->
    C = circuit(Req),
    atom_to_binary(pq_api:backend_for(C, options(Req, C)));
op(Run, Req) when Run =:= <<"run">>; Run =:= <<"sample">> ->
    C = circuit(Req),
    Opts = options(Req, C),
    #result{kind = K, measured = M, samples = S, metadata = Meta} =
        ok(case Run of <<"run">> -> pq_api:run(C, Opts); <<"sample">> -> pq_api:sample(C, Opts) end),
    #{<<"kind">> => atom_to_binary(K), <<"measured">> => M, <<"samples">> => S, <<"metadata">> => term(Meta)};
op(<<"probabilities">>, Req) ->
    C = circuit(Req),
    {Map, Meta} = ok2(pq_api:probabilities(C, options(Req, C))),
    #{<<"probabilities">> => [[P, V] || {P, V} <- lists:sort(maps:to_list(Map))], <<"metadata">> => term(Meta)};
op(<<"evolve">>, Req) ->
    C = circuit(Req),
    {S, W} = ok2(pq_api:evolve(C, options(Req, C))),
    #{<<"state">> => state(S), <<"warnings">> => term(W)};
op(<<"gradient">>, Req) ->
    C = circuit(Req),
    Syms = symbol_table(C),
    Obs = observable(maps:get(<<"observable">>, Req, missing), Syms),
    Binding = binding(maps:get(<<"binding">>, Req, #{}), Syms),
    Opts = options(Req, C),
    G = ok(pq_gradient:gradient(C, Obs, Binding, Opts)),
    #{<<"gradient">> => maps:from_list([{atom_to_binary(K), V} || K := V <- G])};
op(Other, _Req) -> fail({bad_term, #{op => term(Other)}}).

%%% Circuits, options, observables

circuit(Req) ->
    case maps:find(<<"circuit">>, Req) of
        {ok, #{} = Doc} -> ok(pq_circuit:from_json(iolist_to_binary(json:encode(Doc))));
        {ok, Bin} when is_binary(Bin) -> ok(pq_circuit:from_json(Bin));
        _ -> fail({bad_term, circuit_required})
    end.

%% the circuit's symbols by name (atoms that exist: the JSON reader made them)
symbol_table(C) -> maps:from_list([{atom_to_binary(S), S} || S <- pq_circuit:all_symbols(C)]).

sym(Name, Syms) ->
    case maps:find(Name, Syms) of
        {ok, A} -> A;
        error -> fail({unbound_symbols, [Name]})
    end.

binding(#{} = B, Syms) -> maps:from_list([{sym(K, Syms), number(V, K)} || K := V <- B]);
binding(B, _) -> fail({bad_term, #{binding => term(B)}}).

number(V, _) when is_number(V) -> V;
number(V, Where) -> fail({bad_term, #{expected_number => Where, got => term(V)}}).

%% the observables of pq_gradient as JSON: {"probability": {"modes": [...], "pattern": [...]}},
%% {"mean_photon_number": {"modes": [...]}}, {"weighted_sum": [[weight, observable], ...]} with numeric weights or
%% expressions in the circuit's symbols
observable(#{<<"probability">> := #{<<"modes">> := Ms, <<"pattern">> := P}}, _) -> {probability, ints(Ms), ints(P)};
observable(#{<<"mean_photon_number">> := #{<<"modes">> := Ms}}, _) -> {mean_photon_number, ints(Ms)};
observable(#{<<"weighted_sum">> := Terms}, Syms) when is_list(Terms) ->
    {weighted_sum, [{weight(W, Syms), observable(O, Syms)} || [W, O] <- Terms]};
observable(O, _) -> fail({bad_term, #{observable => term(O)}}).

weight(W, _) when is_number(W) -> W;
weight(W, Syms) when is_binary(W) -> ok(pq_symbolic:parse(W, #{allowed_symbols => maps:values(Syms)}));
weight(W, _) -> fail({bad_term, #{weight => term(W)}}).

ints(L) when is_list(L) ->
    [case X of N when is_integer(N), N >= 0 -> N; _ -> fail({bad_term, #{expected_index => term(X)}}) end || X <- L];
ints(L) -> fail({bad_term, #{expected_list => term(L)}}).

%% run options: a fixed table of keys and value kinds
options(Req, C) ->
    case maps:get(<<"options">>, Req, #{}) of
        #{} = O -> maps:from_list([option(K, V, C) || K := V <- O]);
        O -> fail({bad_term, #{options => term(O)}})
    end.

-define(ENUM(K, Vals), option(<<K>>, V, _) -> {binary_to_atom(<<K>>), enum(<<K>>, V, Vals)}).
-define(INT(K), option(<<K>>, V, _) -> {binary_to_atom(<<K>>), int(<<K>>, V)}).
-define(NUM(K), option(<<K>>, V, _) -> {binary_to_atom(<<K>>), number(V, <<K>>)}).
-define(BOOL(K), option(<<K>>, V, _) -> {binary_to_atom(<<K>>), bool(<<K>>, V)}).

?INT("shots"); ?INT("seed"); ?INT("cutoff"); ?INT("max_photons"); ?INT("max_enumeration"); ?INT("max_inflight");
?INT("chunk"); ?INT("timeout"); ?INT("max_attempts"); ?INT("max_dm_dim"); ?INT("max_dm_entries"); ?INT("max_fock_dim");
?INT("burn_in"); ?INT("thin"); ?INT("chains");
?NUM("step"); ?NUM("truncation_tol");
?BOOL("strict"); ?BOOL("check_physical"); ?BOOL("collision_free"); ?BOOL("promote_dm"); ?BOOL("parallel");
?ENUM("backend", [<<"gaussian">>, <<"fock">>, <<"fock_dm">>, <<"boson">>]);
?ENUM("noise", [<<"none">>, <<"trajectory">>]);
?ENUM("shots_mode", [<<"attempts">>, <<"accepted">>]);
?ENUM("algorithm", [<<"clifford_clifford">>, <<"brute_force">>, <<"metropolis">>]);
?ENUM("device", [<<"cpu">>, <<"gpu">>, <<"auto">>]);
?ENUM("dm_storage", [<<"dense">>, <<"blocks">>]);
?ENUM("mode", [<<"numeric">>, <<"symbolic">>, <<"adjoint">>, <<"parameter_shift">>]);
option(<<"bindings">>, V, C) -> {bindings, binding(V, symbol_table(C))};
option(<<"symbols">>, V, C) when is_list(V) -> {symbols, [sym(S, symbol_table(C)) || S <- V]};
option(<<"visibility">>, V, _) -> {visibility, num_or_list(<<"visibility">>, V)};
option(<<"multi_photon">>, V, _) -> {multi_photon, num_or_list(<<"multi_photon">>, V)};
option(<<"dark_counts">>, V, _) when is_number(V) -> {dark_counts, V};
option(<<"dark_counts">>, #{} = V, _) ->
    {dark_counts, maps:from_list([{mode_key(K), number(X, <<"dark_counts">>)} || K := X <- V])};
option(<<"gram">>, Rows, _) when is_list(Rows) -> {gram, [[complex(X) || X <- list(R)] || R <- Rows]};
option(<<"postselect">>, P, _) -> {postselect, ok(pq_json:predicate_from_json(P))};
option(K, _, _) -> fail({bad_term, #{unknown_option => term(K)}}).

enum(K, V, Vals) ->
    case lists:member(V, Vals) of
        %% V is one of the fixed values above, so this makes no new atom from request data
        true -> binary_to_atom(V);
        false -> fail({bad_term, #{option => K, expected => Vals, got => term(V)}})
    end.

int(_, V) when is_integer(V), V >= 0 -> V;
int(K, V) -> fail({bad_term, #{option => K, expected => non_neg_integer, got => term(V)}}).

bool(_, V) when is_boolean(V) -> V;
bool(K, V) -> fail({bad_term, #{option => K, expected => boolean, got => term(V)}}).

num_or_list(_, V) when is_number(V) -> V;
num_or_list(K, V) when is_list(V) -> [number(X, K) || X <- V];
num_or_list(K, V) -> fail({bad_term, #{option => K, got => term(V)}}).

mode_key(K) ->
    try binary_to_integer(K) of
        N when N >= 0 -> N;
        _ -> fail({bad_term, #{mode => K}})
    catch
        error:_ -> fail({bad_term, #{mode => K}})
    end.

complex(X) when is_number(X) -> X;
complex([Re, Im]) when is_number(Re), is_number(Im) -> {complex, Re, Im};
complex(X) -> fail({bad_term, #{expected_complex => term(X)}}).

list(L) when is_list(L) -> L;
list(L) -> fail({bad_term, #{expected_list => term(L)}}).

%%% Results as JSON

%% a public state; large ones refused (MAX_VALUES numbers)
state(#gaussian_state{n_modes = N, hbar = H} = S) ->
    {Mu, Cov} = pq_state:mean_cov(S),
    #{<<"kind">> => <<"gaussian">>, <<"n_modes">> => N, <<"hbar">> => H, <<"mean">> => Mu, <<"cov">> => Cov};
state(#fock_state{n_modes = N, cutoff = C, amps = A}) ->
    {ok, {cvec, V}} = pq_core:unpack(A),
    length(V) =< ?MAX_VALUES orelse fail({dimension_limit, #{service_values => length(V), max => ?MAX_VALUES}}),
    #{<<"kind">> => <<"fock_state">>, <<"n_modes">> => N, <<"cutoff">> => C, <<"amps">> => [[R, I] || {complex, R, I} <- V]};
state(#fock_dm{n_modes = N, cutoff = C, rho = Rho}) ->
    {cmat, D, D, _} = Rho,
    D * D =< ?MAX_VALUES orelse fail({dimension_limit, #{service_values => D * D, max => ?MAX_VALUES}}),
    {ok, {cmat, Rows}} = pq_core:unpack(Rho),
    #{<<"kind">> => <<"fock_dm">>, <<"n_modes">> => N, <<"cutoff">> => C,
      <<"rho">> => [[[R, I] || {complex, R, I} <- Row] || Row <- Rows]}.

%% any Erlang term as JSON: maps as objects (keys as text), lists and tuples as arrays, atoms as strings, numbers
%% as numbers, anything else as its printed form
term(M) when is_map(M) -> maps:from_list([{key(K), term(V)} || K := V <- M]);
%% a list is always an array: a list of mode numbers such as [10] must not become a string
term(L) when is_list(L) -> [term(X) || X <- L];
term({complex, R, I}) when is_number(R), is_number(I) -> [R, I];
term(T) when is_tuple(T) -> [term(X) || X <- tuple_to_list(T)];
term(A) when is_boolean(A); A =:= null -> A;
term(A) when is_atom(A) -> atom_to_binary(A);
term(X) when is_number(X); is_binary(X) -> X;
term(X) -> iolist_to_binary(io_lib:format("~p", [X])).

key(K) when is_atom(K) -> atom_to_binary(K);
key(K) when is_binary(K) -> K;
key(K) -> iolist_to_binary(io_lib:format("~p", [K])).

error_json(R) ->
    Cat = case R of
              {C, _} when is_atom(C) -> C;
              {C, _, _} when is_atom(C) -> C;
              C when is_atom(C) -> C;
              _ -> unknown
          end,
    #{<<"category">> => atom_to_binary(Cat), <<"reason">> => term(R),
      <<"detail">> => iolist_to_binary(io_lib:format("~0p", [R]))}.

encode(Map) -> iolist_to_binary(json:encode(Map)).

ok(ok) -> ok;
ok({ok, V}) -> V;
ok({error, R}) -> fail(R).

ok2({ok, A, B}) -> {A, B};
ok2({error, R}) -> fail(R).

-spec fail(term()) -> no_return().
fail(R) -> throw({error, R}).

%%% Transport: standard input and output, one JSON document per line

-spec stdio() -> no_return().
stdio() ->
    %% standard output carries the protocol: log to standard error
    _ = logger:set_handler_config(default, config, #{type => standard_error}),
    {ok, _} = application:ensure_all_started(pq),
    ok = io:setopts(standard_io, [binary, {encoding, unicode}]),
    %% the client waits for this line before it sends anything
    io:put_chars(standard_io, [encode(#{<<"ready">> => ?SERVICE_VERSION}), $\n]),
    loop().

loop() ->
    case io:get_line(standard_io, "") of
        eof -> halt(0);
        {error, _} -> halt(1);
        Line ->
            Req = string:trim(Line, trailing, "\r\n"),
            case Req of
                <<>> -> ok;
                _ -> io:put_chars(standard_io, [handle(Req), $\n])
            end,
            loop()
    end.

%%% Transport: HTTP (POST /api with a JSON body), inets httpd

%% `erl ... -s pq_service http`: the HTTP transport as a program; PQ_HTTP_PORT (default 0, a free port) and
%% PQ_HTTP_BIND (default 127.0.0.1); prints {"listening": Port} on standard output, then serves until stopped
-spec http() -> no_return().
http() ->
    _ = logger:set_handler_config(default, config, #{type => standard_error}),
    Port = list_to_integer(os:getenv("PQ_HTTP_PORT", "0")),
    case start_http(Port, #{bind => os:getenv("PQ_HTTP_BIND", "127.0.0.1")}) of
        {ok, _, P} ->
            io:put_chars(standard_io, [encode(#{<<"listening">> => P}), $\n]),
            receive after infinity -> ok end;
        {error, R} ->
            io:put_chars(standard_error, io_lib:format("pq_service: ~p~n", [R])),
            halt(1)
    end.

%% Port 0 picks a free port; Opts: bind => "127.0.0.1" (the default; anything else is the caller's decision)
-spec start_http(non_neg_integer()) -> {ok, pid(), inet:port_number()} | {error, term()}.
start_http(Port) -> start_http(Port, #{}).

-spec start_http(non_neg_integer(), map()) -> {ok, pid(), inet:port_number()} | {error, term()}.
start_http(Port, Opts) ->
    {ok, _} = application:ensure_all_started(pq),
    {ok, _} = application:ensure_all_started(inets),
    {ok, Addr} = inet:parse_address(maps:get(bind, Opts, "127.0.0.1")),
    Root = filename:basedir(user_cache, "pqcf-service"),
    ok = filelib:ensure_path(Root),
    case inets:start(httpd, [{port, Port}, {bind_address, Addr}, {server_name, "pqcf"}, {server_root, Root},
                             {document_root, Root}, {modules, [?MODULE]}, {max_body_size, ?MAX_REQUEST},
                             {max_clients, 64}]) of
        {ok, Pid} ->
            [{port, P}] = httpd:info(Pid, [port]),
            {ok, Pid, P};
        E -> E
    end.

-spec stop_http(pid()) -> ok | {error, term()}.
stop_http(Pid) -> inets:stop(httpd, Pid).

%% inets callback: POST /api only
do(#mod{method = Method, request_uri = Uri, entity_body = Body}) ->
    {Status, Out} = case {Method, Uri} of
                        {"POST", "/api"} -> {200, handle(iolist_to_binary(Body))};
                        {_, "/api"} -> {405, encode(#{<<"error">> => error_json({bad_term, post_required})})};
                        _ -> {404, encode(#{<<"error">> => error_json({bad_term, not_found})})}
                    end,
    Head = [{code, Status}, {content_type, "application/json"}, {content_length, integer_to_list(byte_size(Out))}],
    {proceed, [{response, {response, Head, [Out]}}]}.
