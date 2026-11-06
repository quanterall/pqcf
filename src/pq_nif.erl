%% The only module that calls Rust (FS 2.1). The functions below are replaced by the NIF library at load time.
-module(pq_nif).
-export([nif_version/0, check_version/0, cancel_token_new/0, cancel_token_set/1, rng_u64/3,
         permanent/2, permanent_repeated/4, bs_probability/4, sample_boson/4,
         hafnian/2, hafnian_repeated/3, loop_hafnian/4, torontonian/2,
         gaussian_apply_symplectic/5, gaussian_apply_passive/4, gaussian_apply_loss/5, gaussian_marginal/3,
         gaussian_physicality/2, sample_gaussian/6,
         fock_apply_beamsplitter/7, fock_apply_phaseshift/5, fock_apply_kerr/5, fock_apply_passive/5,
         fock_apply_displacement/6, fock_apply_squeezing/6, fock_pnr_probabilities/4,
         unitary_decompose/1, sample_categorical/3, gaussian_to_amat/3, gaussian_probability/5,
         fock_dm_apply_beamsplitter/7, fock_dm_apply_phaseshift/5, fock_dm_apply_kerr/5, fock_dm_apply_passive/5,
         fock_dm_apply_displacement/6, fock_dm_apply_squeezing/6, fock_dm_apply_loss/5,
         fock_dm_apply_thermal_loss/6, gaussian_apply_thermal_loss/6,
         gaussian_measure_homodyne/6, stream_new/2, fock_measure_pnr/5, fock_apply_loss_trajectory/6,
         fock_adjoint_gradient/5, boson_partial_probability/4,
         fock_dmb_from/3, fock_dmb_to_dense/3, fock_dmb_apply/4, fock_dmb_pnr_probabilities/4, fock_dmb_measure_pnr/5,
         fock_dm_adjoint_gradient/5, sym_poly_mul/6, sym_poly_pow/5,
         rpoly_new/4, rpoly_mul/5, rpoly_add/5, rpoly_pow/2, rpoly_scale/4, rpoly_info/1, rpoly_to_map/2, rpoly_to_ast/2,
         horner_prepare/2, horner_count/2, ff_prepare/3, ff_eval/3, ff_tvandermonde/3]).
-on_load(load/0).

%% The NIF {Major, Minor} this code was written against (FS 8.1).
-define(NIF_MAJOR, 0).
-define(NIF_MINOR, 2).

load() ->
    Dir = case code:priv_dir(pq) of
              {error, bad_name} -> filename:join(filename:dirname(filename:dirname(code:which(?MODULE))), "priv");
              D -> D
          end,
    erlang:load_nif(filename:join(Dir, "pq_nif"), 0).

-spec check_version() -> ok | {error, {nif_version_mismatch, map()}}.
check_version() ->
    case nif_version() of
        {?NIF_MAJOR, ?NIF_MINOR, _, _} -> ok;
        Found -> {error, {nif_version_mismatch, #{expected => {?NIF_MAJOR, ?NIF_MINOR}, found => Found}}}
    end.

-spec nif_version() -> {non_neg_integer(), non_neg_integer(), non_neg_integer(), binary()}.
nif_version() -> erlang:nif_error(nif_not_loaded).

-spec cancel_token_new() -> reference().
cancel_token_new() -> erlang:nif_error(nif_not_loaded).

-spec cancel_token_set(reference()) -> ok.
cancel_token_set(_Token) -> erlang:nif_error(nif_not_loaded).

%% The raw ChaCha8 stream (Seed, Stream) for the golden test GT-R01 (FS 8.3).
-spec rng_u64(non_neg_integer(), non_neg_integer(), non_neg_integer()) ->
          {ok, [non_neg_integer()]} | {error, {dimension_limit, #{count := non_neg_integer(), max := 65536}}}.
rng_u64(_Seed, _Stream, _Count) -> erlang:nif_error(nif_not_loaded).

%%% The permanent family (FS 8.4, AN 3; dirty CPU). Opts: cancel (a token), threads (default 1), device (Phase 4,
%%% FS 15 item 5: cpu, the default; gpu, the CUDA kernel of native/pq_nif/gpu/perm.cu, an unsupported_input error
%%% without a device; auto, the GPU from n = 24 when there is one). The GPU covers multiplicities 0 and 1; repeated
%%% rows or columns run on the CPU. Info carries the scale S of AN 2.4 (amendment A-1): the error bound is about
%%% n * eps * S; and the device that ran.

-type opts() :: #{cancel => reference(), threads => pos_integer(), device => cpu | gpu | auto}.
-type perm_info() :: #{scale := float(), positions := non_neg_integer(), blocks := non_neg_integer(),
                       device := cpu | gpu}.
-type kernel_error() :: {error, pq_core:reason()}.

-spec permanent(pq_core:packed_cmatrix(), opts()) -> {ok, {pq_core:nif_complex(), perm_info()}} | kernel_error().
permanent(_M, _Opts) -> erlang:nif_error(nif_not_loaded).

-spec permanent_repeated(pq_core:packed_cmatrix(), [non_neg_integer()], [non_neg_integer()], opts()) ->
          {ok, {pq_core:nif_complex(), perm_info()}} | kernel_error().
permanent_repeated(_M, _RowMult, _ColMult, _Opts) -> erlang:nif_error(nif_not_loaded).

%% FR-B02: probability of OutputOcc from InputOcc through the unitary U (U is not checked here).
-spec bs_probability(pq_core:packed_cmatrix(), [non_neg_integer()], [non_neg_integer()], opts()) ->
          {ok, {float(), #{permanent_scale := float(), positions := non_neg_integer(),
                           blocks := non_neg_integer(), device := cpu | gpu}}} | kernel_error().
bs_probability(_U, _InputOcc, _OutputOcc, _Opts) -> erlang:nif_error(nif_not_loaded).

%%% Boson sampling (FS 8.4, AN 5). Shot i uses stream (seed, stream_start + i), so samples do not depend on how
%%% shots are split between calls. Opts: algorithm (clifford_clifford, default | brute_force | metropolis), seed
%%% (default 0), stream_start (default 0), max_enumeration (brute_force, default 1,000,000), burn_in, thin, chains
%%% (metropolis: 100, 10, min(shots, 8); chain c on stream (seed, c), not shot-parallel), loss => {input | output,
%%% [Eta]} (FR-B05; input loss with clifford_clifford only), shared => [P] (Phase 4: per input mode the probability
%%% of the shared internal state, partial distinguishability; clifford_clifford only), cancel.

-spec sample_boson(pq_core:packed_cmatrix(), [non_neg_integer()], non_neg_integer(), map()) ->
          {ok, {pq_core:packed_imatrix(), #{algorithm := atom(), shots := non_neg_integer()}}} | kernel_error().
sample_boson(U, InputOcc, NShots, Opts) ->
    case sample_boson_raw(U, InputOcc, NShots, Opts) of
        {error, {probability_sum, Total}} -> error({pq_bug, {probability_sum, Total}});
        Other -> Other
    end.

sample_boson_raw(_U, _InputOcc, _NShots, _Opts) -> erlang:nif_error(nif_not_loaded).

%%% Hafnians and the torontonian (FS 8.4, AN 4; dirty CPU). Info carries the scale S (amendment A-1).
%%% The hafnians take device too (gpu/haf.cu; repeated indices on the GPU unless Kan's formula is cheaper; auto stays
%%% on the CPU, spike G1) and report it in info; the torontonian runs on the CPU.

-type haf_info() :: #{scale := float(), subsets := non_neg_integer(), device => cpu | gpu}.

-spec hafnian(pq_core:packed_cmatrix(), opts()) -> {ok, {pq_core:nif_complex(), haf_info()}} | kernel_error().
hafnian(_M, _Opts) -> erlang:nif_error(nif_not_loaded).

-spec hafnian_repeated(pq_core:packed_cmatrix(), [non_neg_integer()], opts()) ->
          {ok, {pq_core:nif_complex(), haf_info()}} | kernel_error().
hafnian_repeated(_M, _Mult, _Opts) -> erlang:nif_error(nif_not_loaded).

%% the hafnian of M with its diagonal replaced by Diag, indices repeated by Mult (FS A.9)
-spec loop_hafnian(pq_core:packed_cmatrix(), pq_core:packed_cvector(), [non_neg_integer()], opts()) ->
          {ok, {pq_core:nif_complex(), haf_info()}} | kernel_error().
loop_hafnian(_M, _Diag, _Mult, _Opts) -> erlang:nif_error(nif_not_loaded).

%% tor(O) for O = I - Q^-1 restricted to the clicked modes, laid out (modes, conjugate modes) (FS A.10)
-spec torontonian(pq_core:packed_cmatrix(), opts()) ->
          {ok, {float(), #{scale := float(), terms := non_neg_integer()}}} | kernel_error().
torontonian(_M, _Opts) -> erlang:nif_error(nif_not_loaded).

%%% Gaussian states (FS 8.4, A.2-A.6, A.12): Mean is a packed rvec of 2m, Cov a packed rmat 2m x 2m, xxpp.

-type gauss() :: {ok, {pq_core:packed_rvector(), pq_core:packed_rmatrix()}} | kernel_error().

%% mu' = S mu + D, Sigma' = S Sigma S^T on Modes; S must be symplectic within physicality_tol (not_symplectic)
-spec gaussian_apply_symplectic(pq_core:packed_rvector(), pq_core:packed_rmatrix(), pq_core:packed_rmatrix(),
                                pq_core:packed_rvector(), [non_neg_integer()]) -> gauss().
gaussian_apply_symplectic(Mean, Cov, S, D, Modes) ->
    gaussian_apply_symplectic_raw(Mean, Cov, S, D, Modes, float(pq_core:get(physicality_tol))).

gaussian_apply_symplectic_raw(_Mean, _Cov, _S, _D, _Modes, _Tol) -> erlang:nif_error(nif_not_loaded).

-spec gaussian_apply_passive(pq_core:packed_rvector(), pq_core:packed_rmatrix(), pq_core:packed_cmatrix(),
                             [non_neg_integer()]) -> gauss().
gaussian_apply_passive(_Mean, _Cov, _U, _Modes) -> erlang:nif_error(nif_not_loaded).

%% Phase 4: loss of transmissivity Eta to a thermal environment of mean occupation Nth
-spec gaussian_apply_thermal_loss(pq_core:packed_rvector(), pq_core:packed_rmatrix(), non_neg_integer(), float(), float(), float()) ->
          {ok, {pq_core:packed_rvector(), pq_core:packed_rmatrix()}} | kernel_error().
gaussian_apply_thermal_loss(_Mean, _Cov, _Mode, _Eta, _Nth, _Hbar) -> erlang:nif_error(nif_not_loaded).

-spec gaussian_apply_loss(pq_core:packed_rvector(), pq_core:packed_rmatrix(), non_neg_integer(), float(), float()) ->
          gauss().
gaussian_apply_loss(_Mean, _Cov, _Mode, _Eta, _Hbar) -> erlang:nif_error(nif_not_loaded).

-spec gaussian_marginal(pq_core:packed_rvector(), pq_core:packed_rmatrix(), [non_neg_integer()]) -> gauss().
gaussian_marginal(_Mean, _Cov, _Modes) -> erlang:nif_error(nif_not_loaded).

%% the smallest eigenvalue of Cov + i (hbar/2) Omega: physical when it is at least -physicality_tol
-spec gaussian_physicality(pq_core:packed_rmatrix(), float()) -> {ok, float()} | kernel_error().
gaussian_physicality(_Cov, _Hbar) -> erlang:nif_error(nif_not_loaded).

%% Detector {homodyne, Phi}: x_Phi on every mode, a packed rmat; pnr | threshold: the chain-rule samplers (AN 6.2,
%% 6.3), a packed imat, Info truncated_mass. Opts: seed, stream_start, cancel, max_photons, max_per_mode,
%% max_clicks.
-spec sample_gaussian(pq_core:packed_rvector(), pq_core:packed_rmatrix(), non_neg_integer(),
                      {homodyne, float()} | pnr | threshold, float(), opts() | map()) ->
          {ok, {pq_core:packed_rmatrix() | pq_core:packed_imatrix(), map()}} | kernel_error().
sample_gaussian(_Mean, _Cov, _NShots, _Detector, _Hbar, _Opts) -> erlang:nif_error(nif_not_loaded).

%%% Fock space, pure states (FS 8.4, AN 7): amplitudes as a packed cvec in the basis of FS 3.4; NModes and Cutoff
%%% explicit. Displacement and squeezing also return the norm deficit (FS 3.5).

-type amps() :: pq_core:packed_cvector().
-type fock() :: {ok, amps()} | kernel_error().

-spec fock_apply_beamsplitter(amps(), float(), float(), non_neg_integer(), non_neg_integer(), pos_integer(),
                              non_neg_integer()) -> fock().
fock_apply_beamsplitter(_A, _Theta, _Phi, _M1, _M2, _N, _C) -> erlang:nif_error(nif_not_loaded).

-spec fock_apply_phaseshift(amps(), float(), non_neg_integer(), pos_integer(), non_neg_integer()) -> fock().
fock_apply_phaseshift(_A, _Phi, _Mode, _N, _C) -> erlang:nif_error(nif_not_loaded).

-spec fock_apply_kerr(amps(), float(), non_neg_integer(), pos_integer(), non_neg_integer()) -> fock().
fock_apply_kerr(_A, _Kappa, _Mode, _N, _C) -> erlang:nif_error(nif_not_loaded).

-spec fock_apply_passive(amps(), pq_core:packed_cmatrix(), [non_neg_integer()], pos_integer(), non_neg_integer()) ->
          fock().
fock_apply_passive(_A, _U, _Modes, _N, _C) -> erlang:nif_error(nif_not_loaded).

-spec fock_apply_displacement(amps(), float(), float(), non_neg_integer(), pos_integer(), non_neg_integer()) ->
          {ok, {amps(), float()}} | kernel_error().
fock_apply_displacement(_A, _R, _Phi, _Mode, _N, _C) -> erlang:nif_error(nif_not_loaded).

-spec fock_apply_squeezing(amps(), float(), float(), non_neg_integer(), pos_integer(), non_neg_integer()) ->
          {ok, {amps(), float()}} | kernel_error().
fock_apply_squeezing(_A, _R, _Phi, _Mode, _N, _C) -> erlang:nif_error(nif_not_loaded).

-spec fock_pnr_probabilities(amps(), [non_neg_integer()], pos_integer(), non_neg_integer()) ->
          {ok, [{[non_neg_integer()], float()}]} | kernel_error().
fock_pnr_probabilities(_A, _Modes, _N, _C) -> erlang:nif_error(nif_not_loaded).

%% FR-C10: the Clements mesh, beamsplitters {M, M+1, Theta, Phi} in application order, then the phases
-spec unitary_decompose(pq_core:packed_cmatrix()) ->
          {ok, {[{non_neg_integer(), pos_integer(), float(), float()}], pq_core:packed_rvector()}} | kernel_error().
unitary_decompose(_U) -> erlang:nif_error(nif_not_loaded).

%% indices drawn from Probs, one uniform per shot from stream (seed, stream_start + i) (AN 5.1)
-spec sample_categorical(pq_core:packed_rvector(), non_neg_integer(), map()) -> {ok, [non_neg_integer()]} | kernel_error().
sample_categorical(_Probs, _NShots, _Opts) -> erlang:nif_error(nif_not_loaded).

%% FS 8.4, A.9: {A, Gamma, Norm} with P(S) = (loop) hafnian(A_S, Gamma_S) / (prod S! Norm)
-spec gaussian_to_amat(pq_core:packed_rvector(), pq_core:packed_rmatrix(), float()) ->
          {ok, {pq_core:packed_cmatrix(), pq_core:packed_cvector(), float()}} | kernel_error().
gaussian_to_amat(_Mean, _Cov, _Hbar) -> erlang:nif_error(nif_not_loaded).

%% the probability of one photon-number (pnr) or click (threshold) pattern of a Gaussian state
-spec gaussian_probability(pq_core:packed_rvector(), pq_core:packed_rmatrix(), float(), [non_neg_integer()],
                           pnr | threshold) -> {ok, float()} | kernel_error().
gaussian_probability(_Mean, _Cov, _Hbar, _Pattern, _Detector) -> erlang:nif_error(nif_not_loaded).

%%% Density matrices (FS 8.4, AN 7.5, A.11): rho' = G rho G^dagger for every gate, loss by its Kraus sum. Rho is a
%%% packed D x D cmat. Each returns {ok, Rho} except displacement and squeezing: {ok, {Rho, TraceDeficit}}.

-type rho() :: pq_core:packed_cmatrix().

dm(Rho, Op, N, C) ->
    case fock_dm_apply_raw(Rho, Op, N, C) of
        {ok, {R, _}} when element(1, Op) =/= displacement, element(1, Op) =/= squeezing, element(1, Op) =/= thermal_loss -> {ok, R};
        Other -> Other
    end.

-spec fock_dm_apply_beamsplitter(rho(), float(), float(), non_neg_integer(), non_neg_integer(), pos_integer(),
                                 non_neg_integer()) -> {ok, rho()} | kernel_error().
fock_dm_apply_beamsplitter(R, T, P, M1, M2, N, C) -> dm(R, {beamsplitter, T, P, M1, M2}, N, C).
-spec fock_dm_apply_phaseshift(rho(), float(), non_neg_integer(), pos_integer(), non_neg_integer()) -> {ok, rho()} | kernel_error().
fock_dm_apply_phaseshift(R, P, M, N, C) -> dm(R, {phaseshift, P, M}, N, C).
-spec fock_dm_apply_kerr(rho(), float(), non_neg_integer(), pos_integer(), non_neg_integer()) -> {ok, rho()} | kernel_error().
fock_dm_apply_kerr(R, K, M, N, C) -> dm(R, {kerr, K, M}, N, C).
-spec fock_dm_apply_passive(rho(), pq_core:packed_cmatrix(), [non_neg_integer()], pos_integer(), non_neg_integer()) ->
          {ok, rho()} | kernel_error().
fock_dm_apply_passive(R, U, Ms, N, C) -> dm(R, {passive, U, Ms}, N, C).
-spec fock_dm_apply_displacement(rho(), float(), float(), non_neg_integer(), pos_integer(), non_neg_integer()) ->
          {ok, {rho(), float()}} | kernel_error().
fock_dm_apply_displacement(R, Rr, P, M, N, C) -> dm(R, {displacement, Rr, P, M}, N, C).
-spec fock_dm_apply_squeezing(rho(), float(), float(), non_neg_integer(), pos_integer(), non_neg_integer()) ->
          {ok, {rho(), float()}} | kernel_error().
fock_dm_apply_squeezing(R, Rr, P, M, N, C) -> dm(R, {squeezing, Rr, P, M}, N, C).
%% FS A.11, intensity transmissivity Eta on Mode
-spec fock_dm_apply_loss(rho(), float(), non_neg_integer(), pos_integer(), non_neg_integer()) -> {ok, rho()} | kernel_error().
fock_dm_apply_loss(R, Eta, M, N, C) -> dm(R, {loss, Eta, M}, N, C).

%% Phase 4: loss to a thermal environment (the amplifier after pure loss, truncated at the cutoff)
-spec fock_dm_apply_thermal_loss(rho(), float(), float(), non_neg_integer(), pos_integer(), non_neg_integer()) ->
          {ok, {rho(), float()}} | kernel_error().
fock_dm_apply_thermal_loss(R, Eta, Nth, M, N, C) -> dm(R, {thermal_loss, Eta, Nth, M}, N, C).

fock_dm_apply_raw(_Rho, _Op, _N, _C) -> erlang:nif_error(nif_not_loaded).

%% FS A.7, FR-G08: measure x_Phi on Mode and condition the rest; Opts: fix_outcome, or seed and stream (the first
%% normal of the stream). {singular_matrix, Variance} when the quadrature has no variance.
-spec gaussian_measure_homodyne(pq_core:packed_rvector(), pq_core:packed_rmatrix(), non_neg_integer(), float(), float(),
                                map()) -> {ok, {float(), pq_core:packed_rvector(), pq_core:packed_rmatrix()}} | kernel_error().
gaussian_measure_homodyne(_Mean, _Cov, _Mode, _Phi, _Hbar, _Opts) -> erlang:nif_error(nif_not_loaded).

%% the stream (Seed, Stream) as a resource: kernels given it as stream_ref draw from it in order (a shot driven
%% step by step, FS 9.3)
-spec stream_new(non_neg_integer(), non_neg_integer()) -> reference().
stream_new(_Seed, _Stream) -> erlang:nif_error(nif_not_loaded).

%% FR-M07: {Prob, State projected onto Outcome on Modes and renormalised}; the measured modes stay in the state
-spec fock_measure_pnr(pq_core:packed_cvector() | pq_core:packed_cmatrix(), [non_neg_integer()], [non_neg_integer()],
                       pos_integer(), non_neg_integer()) -> {ok, {float(), term()}} | kernel_error().
fock_measure_pnr(_State, _Modes, _Outcome, _N, _C) -> erlang:nif_error(nif_not_loaded).

%% FR-F05 trajectory mode: a lost-photon count drawn (one uniform of stream_ref, or of stream (seed, stream)), the
%% pure state replaced by K_l psi renormalised; {Amps, Lost}
-spec fock_apply_loss_trajectory(pq_core:packed_cvector(), float(), non_neg_integer(), pos_integer(), non_neg_integer(),
                                 map()) -> {ok, {pq_core:packed_cvector(), non_neg_integer()}} | kernel_error().
fock_apply_loss_trajectory(_A, _Eta, _Mode, _N, _C, _Opts) -> erlang:nif_error(nif_not_loaded).

%% Adjoint gradients (Phase 4, FS 15 item 3; native/pq_nif/src/fockgrad.rs): the value of every observable term at
%% the output of Gates applied to the pure state Amps, and the gradient of the weighted sum of the terms with
%% respect to every gate parameter (one list per gate, in the instruction's parameter order; none for passive).
-type adjoint_gate() :: {beamsplitter, float(), float(), non_neg_integer(), non_neg_integer()}
                      | {phaseshift, float(), non_neg_integer()} | {kerr, float(), non_neg_integer()}
                      | {passive, pq_core:packed_cmatrix(), [non_neg_integer()]}
                      | {displacement, float(), float(), non_neg_integer()} | {squeezing, float(), float(), non_neg_integer()}.
-type adjoint_term() :: {probability, float(), [non_neg_integer()], [non_neg_integer()]}
                      | {photons, float(), [non_neg_integer()]}.
-spec fock_adjoint_gradient(amps(), [adjoint_gate()], [adjoint_term()], pos_integer(), non_neg_integer()) ->
          {ok, {[float()], [[float()]]}} | kernel_error().
fock_adjoint_gradient(_A, _Gates, _Terms, _N, _C) -> erlang:nif_error(nif_not_loaded).

%% Partial distinguishability (Phase 4, FS 15 item 4): the probability of Output through U for photons whose
%% internal states have the Gram matrix Gram (n x n, the photons in input-mode order); n! permanents, n <= 12.
-spec boson_partial_probability(pq_core:packed_cmatrix(), [non_neg_integer()], [non_neg_integer()], pq_core:packed_cmatrix()) ->
          {ok, float()} | kernel_error().
boson_partial_probability(_U, _In, _Out, _Gram) -> erlang:nif_error(nif_not_loaded).

%%% Photon-number blocks (Phase 4, FS 15 item 6; internal to pq_fock_dm, Decision 3: the NIFs of FS 8.4 and every
%%% public function present the full matrix). Blocks: [{N1, N2, Binary}], block (N1, N2) d_N1 x d_N2 complex
%%% row-major little-endian, zero blocks absent; every operation equals the dense kernel's entry for entry.
-type dm_blocks() :: [{non_neg_integer(), non_neg_integer(), binary()}].

%% blocks from amplitudes (cvec) or a full matrix (cmat)
-spec fock_dmb_from(pq_core:packed_cvector() | pq_core:packed_cmatrix(), pos_integer(), non_neg_integer()) -> {ok, dm_blocks()} | kernel_error().
fock_dmb_from(_State, _N, _C) -> erlang:nif_error(nif_not_loaded).

-spec fock_dmb_to_dense(dm_blocks(), pos_integer(), non_neg_integer()) -> {ok, pq_core:packed_cmatrix()} | kernel_error().
fock_dmb_to_dense(_Blocks, _N, _C) -> erlang:nif_error(nif_not_loaded).

%% the operations of fock_dm_apply_*, as their tuples; {Blocks, TraceDeficit, StoredEntries}
-spec fock_dmb_apply(dm_blocks(), tuple(), pos_integer(), non_neg_integer()) ->
          {ok, {dm_blocks(), float(), non_neg_integer()}} | kernel_error().
fock_dmb_apply(_Blocks, _Op, _N, _C) -> erlang:nif_error(nif_not_loaded).

-spec fock_dmb_pnr_probabilities(dm_blocks(), [non_neg_integer()], pos_integer(), non_neg_integer()) ->
          {ok, [{[non_neg_integer()], float()}]} | kernel_error().
fock_dmb_pnr_probabilities(_Blocks, _Modes, _N, _C) -> erlang:nif_error(nif_not_loaded).

-spec fock_dmb_measure_pnr(dm_blocks(), [non_neg_integer()], [non_neg_integer()], pos_integer(), non_neg_integer()) ->
          {ok, {float(), dm_blocks()}} | kernel_error().
fock_dmb_measure_pnr(_Blocks, _Modes, _Outcome, _N, _C) -> erlang:nif_error(nif_not_loaded).

%% The adjoint on density matrices (Phase 4 follow-up, fockgrad.rs adjoint_dm): from the pure input Amps through the
%% channels -- the gates of fock_adjoint_gradient/5, {loss, Eta, M} (one gradient, d/dEta) and
%% {thermal_loss, Eta, Nth, M} (none) -- the terms' values and the gradient per channel.
-spec fock_dm_adjoint_gradient(amps(), [adjoint_gate() | {loss, float(), non_neg_integer()} | {thermal_loss, float(), float(), non_neg_integer()}],
                               [adjoint_term()], pos_integer(), non_neg_integer()) -> {ok, {[float()], [[float()]]}} | kernel_error().
fock_dm_adjoint_gradient(_A, _Chans, _Terms, _N, _C) -> erlang:nif_error(nif_not_loaded).

%% The product of two plain polynomials of the symbolic engine (pq_sym_poly): each as its terms' exponents, one per
%% atom key of Vars (sorted), and its coefficients times L, integers, real and imaginary per term; the product as
%% pq_sym_poly's map, or overflow (a coefficient outside i128) or packing (the exponent ranges need more than 64
%% bits) for pq_sym_poly to multiply in Erlang.
-spec sym_poly_mul([term()], [integer()], [integer()], [integer()], [integer()], pos_integer()) ->
          {ok, map()} | overflow | packing.
sym_poly_mul(_Vars, _ExpsA, _CoefA, _ExpsB, _CoefB, _L) -> erlang:nif_error(nif_not_loaded).

%% A plain polynomial as for sym_poly_mul/6, its coefficients times L, to the power N; D = L^N, the denominator of
%% the result's coefficients; overflow or packing as there.
-spec sym_poly_pow([term()], [integer()], [integer()], non_neg_integer(), pos_integer()) -> {ok, map()} | overflow | packing.
sym_poly_pow(_Vars, _Exps, _Coef, _N, _D) -> erlang:nif_error(nif_not_loaded).

%% Polynomials kept in Rust (K2; pq_sym_poly holds them as {rpoly, Vars, Resource}): terms as for sym_poly_mul/6,
%% their coefficients over Denom; products, sums, powers and scalings of them (a variable map sends each operand's
%% variable j to the result's position), refused as overflow or packing; their size; back to pq_sym_poly's map;
%% pq_sym_norm:to_ast/1's result for variables that are all symbols.
-spec rpoly_new(non_neg_integer(), [integer()], [integer()], pos_integer()) -> reference().
rpoly_new(_NVars, _Exps, _Coef, _Denom) -> erlang:nif_error(nif_not_loaded).
-spec rpoly_mul(reference(), [non_neg_integer()], reference(), [non_neg_integer()], non_neg_integer()) -> {ok, reference()} | overflow | packing.
rpoly_mul(_A, _MapA, _B, _MapB, _N) -> erlang:nif_error(nif_not_loaded).
-spec rpoly_add(reference(), [non_neg_integer()], reference(), [non_neg_integer()], non_neg_integer()) -> {ok, reference()} | overflow | packing.
rpoly_add(_A, _MapA, _B, _MapB, _N) -> erlang:nif_error(nif_not_loaded).
-spec rpoly_pow(reference(), non_neg_integer()) -> {ok, reference()} | overflow | packing.
rpoly_pow(_A, _K) -> erlang:nif_error(nif_not_loaded).
-spec rpoly_scale(reference(), integer(), integer(), pos_integer()) -> {ok, reference()} | overflow.
rpoly_scale(_A, _Re, _Im, _Den) -> erlang:nif_error(nif_not_loaded).
-spec rpoly_info(reference()) -> {non_neg_integer(), non_neg_integer()}.
rpoly_info(_A) -> erlang:nif_error(nif_not_loaded).
-spec rpoly_to_map(reference(), [term()]) -> map().
rpoly_to_map(_A, _Vars) -> erlang:nif_error(nif_not_loaded).
-spec rpoly_to_ast(reference(), [term()]) -> term().
rpoly_to_ast(_A, _Vars) -> erlang:nif_error(nif_not_loaded).

%% Plain polynomials of symbols ({Exps, Coef, Den} each, over NVars variables) prepared once, and the operations of
%% pq_sym_slp's Horner programs of them for a variable order (indices, first to last), counted without building
%% them (H4).
-spec horner_prepare(non_neg_integer(), [{[integer()], [integer()], pos_integer()}]) -> {ok, reference()} | overflow | packing.
horner_prepare(_NVars, _Polys) -> erlang:nif_error(nif_not_loaded).
-spec horner_count(reference(), [non_neg_integer()]) -> {ok, non_neg_integer()} | overflow | packing.
horner_count(_Prepared, _Order) -> erlang:nif_error(nif_not_loaded).

%% A straight-line program over prime fields (R1; pq_ff): prepared from its operations in slot order ({input, I},
%% {const, Sign, NumLimbs, DenLimbs}, {mul | add | sub, A, B}, {pow, A, E}) and outputs, then evaluated mod a prime < 2^63 at a flat
%% list of points; a point dividing by zero has every output 18446744073709551615.
-spec ff_prepare(non_neg_integer(), [tuple()], [non_neg_integer()]) -> {ok, reference()} | {error, term()}.
ff_prepare(_NInputs, _Ops, _Outputs) -> erlang:nif_error(nif_not_loaded).
-spec ff_eval(reference(), pos_integer(), [non_neg_integer()]) -> {ok, [non_neg_integer()]} | bad_prime | {error, term()}.
ff_eval(_Program, _P, _Points) -> erlang:nif_error(nif_not_loaded).
-spec ff_tvandermonde(pos_integer(), [non_neg_integer()], [non_neg_integer()]) -> {ok, [non_neg_integer()]} | singular.
ff_tvandermonde(_P, _Nodes, _Values) -> erlang:nif_error(nif_not_loaded).
