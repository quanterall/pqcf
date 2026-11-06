# Changelog

PQCF follows semantic versioning from the API freeze (`docs/decisions/API-freeze.md`). Until 2.0, every change to
the public API is additive.

## Unreleased (1.3.0)

### Added

- **Horner schemes for compiled expressions** (`H1`):
  - `pq_symbolic:compile/3` with `#{optimize => horner}`: every polynomial is compiled as a multivariate Horner
    scheme, in a fixed occurrence order, with identical sub-polynomials shared;
  - on the resultant benchmarks of FORM's expression-simplification papers: 2.6–3.0× fewer operations than the
    default program, 4–6 % fewer than FORM 5's `Format O1`, and evaluation 3.4–3.6× faster;
  - the default is unchanged, bit for bit.

## 1.2.0 (5 November 2026)

Beyond the specification's roadmap: PQCF from Python, and faster density matrices. Everything is additive and
results are bit for bit the same as 1.1.0's. Upgrade every node of a cluster together, as for 1.1.0.

### Added

- **A JSON service and a Python client** (`C1`):
  - `pq_service`: version, validate, backend_for, run, sample, probabilities, evolve and gradient as JSON
    requests;
  - two transports: stdio (`-s pq_service stdio`), and HTTP (`-s pq_service http`, `start_http/1,2`, 127.0.0.1 by
    default);
  - atom-safe: option keys, enum values and symbol names are matched against fixed tables or the circuit's own
    symbols.
- **The `pqcf` package** (`python/`, standard library only):
  - `Circuit`, a builder of the JSON circuit format;
  - `connect()`, the client;
  - observables;
  - `pqcf.convert`, to and from Perceval and Piquasso;
  - `pqcf.linear`, a passive circuit's unitary.
- `pq_json:predicate_from_json/1` (internal).

### Changed

- **Faster density-matrix gates** (`B2`): each gate is prepared once (index groups, sector matrices, the Clements
  mesh, the truncated matrices) instead of for every column and row. On photon-number blocks, number-conserving
  gates touch only the stored sectors.
  - Dense ρ: 5–23× per gate.
  - Blocks: 25× per beamsplitter.
  - Cutoff 20 on 4 modes, end to end: 12 s → 1.2 s.

  Results are identical bit for bit.

## 1.1.0 (5 November 2026)

The first release after 1.0.0-rc.1: Phase 4 of the specification's roadmap (FS 15). Everything is additive. Code
written for 1.0.0-rc.1 runs unchanged, and with the defaults its results are bit for bit the same, including the
sampling streams. Each item has a decision record in `docs/decisions/`.

### Added

- **Threshold detection of displaced Gaussian states** (the loop torontonian, `L1`). Threshold probabilities and
  the threshold sampler now accept non-zero means instead of returning `nonzero_mean_threshold`.
- **Feed-forward** (`F1`):
  - the instructions `{feed_forward, Gate, #{Symbol => Mode}}` and `{conditional, Predicate, [Gate]}`;
  - gate parameters and conditions may depend on earlier measurement outcomes in the same shot;
  - the JSON ops `feed_forward` and `conditional`.
- **Gradients** (`AD1`):
  - `mode => adjoint`: one forward and one backward pass over the Fock state, on density matrices when there is
    loss, with the loss transmissivity differentiated too;
  - `mode => parameter_shift`: exact shift rules with frequencies bounded by the cutoff, also with loss;
  - the strategy modules `pq_gradient_adjoint` and `pq_gradient_shift`, and the type `pq_gradient:observable/0`.
- **Noise models** (`N1`):
  - `visibility` and `gram`: partially distinguishable photons in boson sampling;
  - `multi_photon`: multi-photon source terms;
  - `dark_counts`: detector dark counts, also in exact probabilities;
  - the instruction `{thermal_loss, Mode, Eta, Nth}`, exact on the Gaussian backend and truncated on density
    matrices, with the JSON op `thermal_loss`.
- **GPU offload** (`G1`): `device => cpu | gpu | auto` for permanents, boson probabilities and hafnians, on NVIDIA
  GPUs through the CUDA driver.
  - The GPU computes every term bit-identically to the CPU kernel, and the results are reproducible bit for bit.
  - `auto` uses the GPU for permanents from n = 24 and keeps hafnians on the CPU (spike G1 on the GB10).
  - Building needs no CUDA.
- **Density matrices by photon-number blocks** (`B1`): `dm_storage => blocks` and `max_dm_entries`. Results are the
  same, with far less memory for number-conserving circuits: cutoff 20 on 4 modes takes 0.18 GB instead of 1.81 GB.
- **New metadata and info keys:** `dark_count_truncated_mass`, and `device` in the permanent and hafnian kernels'
  info maps.

### Fixed

- `pq_stats:mean_photon_number/2` had no clause for density matrices, so the numeric gradient of a mean photon
  number through loss crashed with `function_clause`.

### Changed

- The NIF contract is 0.2: new kernels only, and `pq_nif` remains internal.
- **Upgrade every node of a cluster together.** The node handshake (unchanged) shares work only between nodes with
  the same `pq` version and NIF minor version. A 1.1.0 node and a 1.0.0-rc.1 node do not split jobs: each runs its
  own jobs locally, with the same results.
- `pq_api:backend_for/2` chooses `boson` when `visibility`, `gram` or `multi_photon` is given, and `fock_dm` for
  circuits with `thermal_loss`. These are the only backends that model them.

## 1.0.0-rc.1 (1 November 2026)

- The API frozen (`docs/decisions/API-freeze.md`).
- Every functional requirement of the specification implemented and tested.
- Benchmarks on a provisional reference machine (`bench/HARDWARE.md`).
