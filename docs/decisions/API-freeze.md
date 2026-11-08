# API freeze: 1.0.0-rc.1 (1 November 2026)

**Decision: the public API of `pq` is frozen as of version 1.0.0-rc.1, by the project owner, 1 November 2026.** Until
2.0, changes to it are additive only:

- new functions, new options, new error details and new metadata keys are allowed;
- nothing listed here is removed or renamed, and no return shape changes;
- a function that should go is deprecated first (documented, a warning, kept until 2.0).

Phase 3 closed with this freeze. Every functional requirement of the specification is implemented and tested.
The two deviations accepted on 31 October are recorded in `API-review.md`. The NFR-02 benchmark instance was
fixed on 1 November, and its results are in `bench/HARDWARE.md`.

## Frozen

The functions below, taken from the compiled modules, with their documented arguments, options and return
shapes (`docs/guide.md`):

- `pq_api`: `backend_for/2`, `circuit/1`, `circuit/2`, `evolve/2`, `evolve/3`, `probabilities/2`, `run/2`, `run/3`, `run_batch/3`, `sample/2`, `sample/3`
- `pq_circuit`: `add/2`, `all_symbols/1`, `check_capabilities/2`, `check_capabilities/3`, `compose/2`, `compose/3`, `delete/2`, `effective_binding/2`, `free_symbols/1`, `from_binary/1`, `from_json/1`, `from_json/2`, `insert/3`, `lower/2`, `mesh/1`, `mzi/4`, `new/1`, `new/2`, `replace/3`, `to_binary/1`, `to_json/1`, `two_mode_squeezing/3`, `two_mode_squeezing/4`, `validate/1`
- `pq_symbolic`: `add/1`, `bs_probability/3`, `canonicalise/1`, `canonicalise/2`, `circuit_unitary/1`, `circuit_unitary/2`, `compile/2`, `compile/3`, `conj/1`, `cos/1`, `dagger/1`, `differentiate/2`, `differentiate/3`, `divide/2`, `domain_of/2`, `equiv/3`, `eval_batch/2`, `eval_slp/2`, `evaluate/2`, `evaluate/3`, `exp/1`, `free_symbols/1`, `gaussian_evolve/2`, `gradient/2`, `gradient/3`, `hash/1`, `i/0`, `is_unitary/1`, `is_unitary/2`, `jacobian/2`, `jacobian/3`, `kron/2`, `madd/1`, `matmul/2`, `matrix/1`, `mul/1`, `neg/1`, `num/1`, `parse/1`, `parse/2`, `permanent/1`, `pi/0`, `pow/2`, `pretty/1`, `pretty/2`, `resolve_params/2`, `resolve_params/3`, `shape/1`, `simplify/1`, `simplify/2`, `sin/1`, `size/1`, `smul/2`, `sqrt/1`, `sub/2`, `subs/2`, `subs/3`, `substitute/2`, `sym/1`, `to_latex/1`, `transpose/1`
- `pq_state`: `basis_size/2`, `check_physical/1`, `coherent/3`, `fock_dm/1`, `fock_state/2`, `fock_superposition/2`, `from_mean_cov/3`, `from_xpxp/1`, `from_xpxp/2`, `gaussian_to_fock_dm/2`, `index/3`, `mean_cov/1`, `pattern/3`, `reduce/2`, `squeezed_vacuum/3`, `thermal/3`, `to_xpxp/1`, `vacuum/2`
- `pq_sampling`: `passive_unitary/1`, `postselect/2`, `probabilities/3`, `sample/3`
- `pq_stats`: `fidelity/2`, `g2/2`, `histogram/2`, `mean_photon_number/2`, `tv_distance/2`
- `pq_noise`: `detector_efficiency/2`, `insert_loss/3`
- `pq_gradient`: `behaviour_info/1`, `gradient/4`, `supports/1`, `symbolic_gradient/3`, `value/3`
- `pq_convert`: `from_perceval/1`, `from_piquasso/1`, `perceval_unitary/2`, `to_perceval/2`, `to_piquasso/2`
- `pq_gaussian`: `lower/2`, `symbolic_evolve/2`, `to_symbolic/1` (its `pq_backend` callbacks are frozen as that behaviour)
- `pq_parallel`: `cancel/1`, `job_opts/1`, `register_node/0`, `run_job/2`, `status/1` (`chunks/3` and `to_rows/1` are internal)
- `pq_core`: `categories/0`, `config/0`, `defaults/0`, `get/1`, `get/2`, `is_reason/1`, `pack/1`, `unpack/1`, `validate/1`

Also frozen:

- **Records** (`include/pq.hrl`): `#result{}` (`kind`, `measured`, `samples`, `metadata`) and the states
  `#gaussian_state{}`, `#fock_state{}` and `#fock_dm{}`, read-only for callers. Circuits are built through
  `pq_circuit`.
- **Behaviours:** `pq_backend` (for backend authors) and `pq_gradient` (strategies).
- **Instructions and inputs:** the instruction tuples and input specifications of the guide (section 3).
- **Options:** the keys of `pq_core:defaults/0` and the run options of the guide. New keys may be added.
- **Errors:** the categories of `pq_core:categories/0`. New categories may be added; none is removed.
- **Determinism (D1):** shot k uses stream (seed, k), whatever the workers, chunks or nodes. The stream algorithm
  (ChaCha8, SplitMix64 key) and the stream orders of decision S7 change only in a major version.
- **Telemetry:** the events `[pq, job, start | stop | exception]` and `[pq, chunk, start | stop | retry]`, with
  their measurement and metadata keys.
- **Interchange:** the JSON circuit format (Appendix B, version 1) and the safe binary format of
  `pq_circuit:to_binary/1`.

## Added in 1.1.0 (5 November 2026)

These additions are frozen from 1.1.0 by the same rule. Phase 4 decision records: `L1`, `F1`, `AD1`, `N1`, `G1` and
`B1`.

- **Instructions:** `{feed_forward, Gate, #{Symbol => Mode}}`, `{conditional, Predicate, [Gate]}`,
  `{thermal_loss, Mode, Eta, Nth}`, and their JSON ops (Appendix B, version 1, extended).
- **Gradients:**
  - the modes `adjoint` and `parameter_shift` of `pq_gradient:gradient/4`;
  - the strategy modules `pq_gradient_adjoint` (`gradient/4`, `supports/1`) and `pq_gradient_shift`
    (`gradient/4`, `supports/1`, `rule/2`);
  - the type `pq_gradient:observable/0`.
- **Options:** `visibility`, `gram`, `multi_photon`, `dark_counts`, `device`, `dm_storage`, `max_dm_entries`.
- **Metadata and info keys:** `dark_count_truncated_mass`, and `device` in the permanent and hafnian info maps.
- **Results that used to be errors:** threshold detection of displaced Gaussian states, and
  `pq_stats:mean_photon_number/2` on a `#fock_dm{}`.
- **Not frozen (internal):**
  - the new helper exports of `pq_noise` (`shared_probabilities/2`, `visibility_gram/2`, `gram/2`, `dark_means/2`,
    `dark_rows/4`, `dark_distribution/4`, `probabilities/3`, `source_mixture/3`);
  - `pq_gradient_chain`;
  - `pq_sampling:gram_distribution/4`;
  - `pq_simulator:evolve_state/3`;
  - the NIF contract, now 0.2.

## Added in 1.2.0 (5 November 2026)

Frozen by the same rule (decision records `C1`, `B2`):

- **`pq_service`:** `handle/1`, `request/1`, `stdio/0`, `http/0`, `start_http/1`, `start_http/2`, `stop_http/1`.
- **Its JSON protocol, service version 1:**
  - the ops `version`, `validate`, `backend_for`, `run`, `sample`, `probabilities`, `evolve` and `gradient`;
  - the request fields `id`, `op`, `circuit`, `options`, `observable` and `binding`;
  - the response shapes `ok` and `error` (`category`, `reason`, `detail`);
  - the option names, as for the Erlang API.

  New ops, options and result keys may be added; none is removed.
- **The Python package `pqcf`:** `Circuit` and its methods, `connect`, `Client` (`version`, `validate`,
  `backend_for`, `run`, `sample`, `probabilities`, `probabilities_with_metadata`, `evolve`, `gradient`,
  `request`), `Result`, `PqcfError`, `probability`, `mean_photon_number`, `weighted_sum`,
  `pqcf.convert` (`from_perceval`, `to_perceval`, `to_piquasso`, `from_piquasso`) and `pqcf.linear.unitary`.
- **Not frozen (internal):** `fockplan.rs` (behind the NIF contract), and `pq_json:predicate_from_json/1`.

## Added in 1.3.0 (7 November 2026)

Frozen by the same rule (decision records `K1`, `K2`, `H1`–`H4`):

- **`pq_symbolic:compile/3` options:** `optimize => none | horner | cse | horner_cse`, and `order_search =>
  #{trials, seed, batch, nodes, strategy => local | mcts, exploration, direction => forward | backward}`.
- **The option `poly_kernel => boolean()`** of the symbolic calls (default `true`).
- **The program operation `{sub, Slot, A, B}`** in programs compiled with Horner (format `{slp, 1, ...}` unchanged;
  `eval_slp/2` evaluates it).
- **Not frozen (internal):** `pq_sym_cse`, `pq_sym_search`, the Rust kernels (`polymul.rs`, `rpoly.rs`,
  `hcount.rs`) and their NIFs, the `{rpoly, ...}` handles inside `pq_sym_poly`.

## Added in 1.4.0 (8 November 2026)

Frozen by the same rule (decision record `R1`):

- **`pq_ff`:** `prime/1`, `is_prime/1`, `prepare/1`, `eval/3`, `inv/2`, `crt/2`, `ratrecon/2`.
- **`pq_recon`:** `polynomial/3`, `rational/3`, `linear/4`. A's and B's rows are lists or sparse maps. The options
  are `executor => local | {nodes, Nodes}`, `workers`, `per_node`, `batch`, `seed`, `components`, `rhs`,
  `guess => boolean()`, and the expression limits of the symbolic calls.
- **`pq_simplify`:** `rational/4`, `polynomial/3`, `alphabet/1`, `adds/1`, `verify/4`, `format/1`, `divide/2`.
- **`pq_ibp_family`:** `new/1`, `example/1`, `identities/2`, `zero/2`, `seeds/3`, `canon/2`, `derivative/3`,
  `vars/1`, `size/1`, `kin/3`. A family's map is opaque beyond the keys `new/1` takes.
- **`pq_ibp`:** `reduce/3`, `system/4`, `seeds/3`, `auto_seeds/3`, `differential_equation/4`. The system's map is
  for inspection, and its keys may grow.
- **Not frozen (internal):** `pq_rexec`, `pq_rstore`, the finite-field NIFs (`ff_*`) and `ffslp.rs`, among them
  the learned elimination.

## Not frozen

- **Internal modules:** the `pq_sym_*` modules, `pq_job`, `pq_worker`, the supervisors, `pq_node_registry`,
  `pq_spec_cache`, `pq_json`, `pq_telemetry`, the backend modules' callbacks outside `pq_backend`, and `pq_nif`
  (the Rust contract has its own version, `nif_version/0`).
- **Performance figures and limits' defaults:** these may improve; their keys stay.
