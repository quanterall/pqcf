%% Noise (FS 5.8): loss channels inserted into a circuit. Loss is exact on Gaussian states and goes through
%% density matrices on Fock states (FR-N03, FR-F05).
-module(pq_noise).
-include("pq.hrl").
-export([insert_loss/3, detector_efficiency/2]).

%% FR-N02: after every instruction that is not a measurement, a loss on each mode it touches; Etas maps modes to
%% transmissivities (numeric or symbolic), modes without an entry get none. Opts is reserved.
-spec insert_loss(#circuit{}, #{non_neg_integer() => pq_symbolic:sym_expr()}, map()) -> {ok, #circuit{}}.
insert_loss(#circuit{instrs = Is} = C, Etas, _Opts) ->
    New = lists:append([[I | [{loss, M, maps:get(M, Etas)} || M <- touched(I), maps:is_key(M, Etas)]]
                        || I <- Is]),
    {ok, C#circuit{instrs = New}}.

%% FR-N04: a loss of transmissivity Eta (or Etas per mode) immediately before each measurement, on its modes
-spec detector_efficiency(#circuit{}, pq_symbolic:sym_expr() | #{non_neg_integer() => pq_symbolic:sym_expr()}) ->
          {ok, #circuit{}}.
detector_efficiency(#circuit{instrs = Is} = C, Eta) ->
    EtaOf = fun(M) when is_map(Eta) -> maps:get(M, Eta, none); (_) -> Eta end,
    New = lists:append([case measured(I) of
                            [] -> [I];
                            Ms -> [{loss, M, EtaOf(M)} || M <- Ms, EtaOf(M) =/= none] ++ [I]
                        end || I <- Is]),
    {ok, C#circuit{instrs = New}}.

touched({beamsplitter, A, B, _, _}) -> [A, B];
touched({unitary, Ms, _}) -> Ms;
touched({custom, _, Args}) -> maps:get(modes, Args, []);
touched({Tag, _}) when Tag =:= measure_pnr; Tag =:= measure_threshold -> [];
touched({measure_homodyne, _, _}) -> [];
touched({loss, _, _}) -> [];
touched(I) -> [element(2, I)].

measured({measure_pnr, Ms}) -> Ms;
measured({measure_threshold, Ms}) -> Ms;
measured({measure_homodyne, M, _}) -> [M];
measured(_) -> [].
