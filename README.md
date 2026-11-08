# PQCF — Photonic Quantum Computing Framework

A simulator for photonic quantum computing: Gaussian and Fock-space backends, boson sampling and Gaussian boson
sampling, a symbolic engine for parameterised circuits, and parallel, reproducible sampling across nodes. The
top layer is an Erlang/OTP application (`pq`); the numerical kernels are a Rust library loaded as a NIF
(`pq_nif`, through [Rustler](https://github.com/rusterlium/rustler)).

PQCF stands for *Photonic Quantum Computing Framework*. It is not a post-quantum cryptography library and
implements no cryptographic algorithms.

**Status: 1.4.0 (`CHANGELOG.md`).**

- The API is frozen (`docs/decisions/API-freeze.md`): changes to the public API are additive only until 2.0.
- Every functional requirement of the specification is implemented and tested.
- Phase 4, the specification's roadmap, is complete:
  - threshold detection of displaced states (the loop torontonian);
  - feed-forward;
  - adjoint gradients (through loss too) and parameter-shift rules;
  - noise models: partial distinguishability, multi-photon sources, dark counts, thermal loss;
  - GPU offload of permanents and hafnians;
  - photon-number-block storage for density matrices.
- Beyond the roadmap:
  - 1.2: a JSON service and a Python client (`python/`); density-matrix gates 5–25 times faster;
  - 1.3: the symbolic engine measured against FORM 5;
  - 1.4: reducing large symbolic results by reconstruction over finite fields, with Feynman-integral reduction
    measured against FIRE (below).
- Benchmarks are from a provisional reference machine (`bench/HARDWARE.md`).

## What is done

### Core and kernels
- **Infrastructure:** the NIF with loading, version check and cancel tokens; the ChaCha8 random stream (GT-R01);
  core types, configuration, the error taxonomy and the packed binary codec (`pq_core`).
- **Permanents:** `permanent`, `permanent_repeated`, `bs_probability`, deterministic by fixed blocks, checked
  against thewalrus (GT-P01).
- **Hafnians:** the hafnian, loop hafnian and torontonian, within 1e-15 of thewalrus up to M = 40 and 20 clicks
  (S3, S4); Kan's repeated-index hafnian; the loop torontonian for threshold detection of displaced states.
- **Fock matrices:** beamsplitter, displacement and squeezing, checked against 80-digit closed forms (S5, S8).

### Circuits and states
- **Circuits** (`pq_circuit`): construction and editing; validation of modes, parameters, numeric domains,
  unitarity, inputs, measured-mode reuse and custom instructions; composition with mode maps and renames; symbols;
  safe binary serialisation; the gate library; the Clements mesh (`mesh/1`); JSON interchange (Appendix B).
- **States** (`pq_state`): the Fock basis index (FS 3.4), Gaussian preparation, reduced states, physicality;
  Fock state constructors; Gaussian to Fock conversion.
- **Converters** to and from Perceval and Piquasso (`pq_convert`), conventions asserted against the libraries.

### Backends
- **Gaussian** (`pq_gaussian`): evolution, terminal and mid-circuit homodyne sampling.
- **Fock, pure states** (`pq_fock`): beamsplitter sectors, phase, Kerr, passive unitaries, displacement and
  squeezing with the norm deficit; PNR, threshold and mixed sampling.
- **Fock, density matrices** (`pq_fock_dm`): the Kraus sum of loss, photon statistics after loss equal to the
  Gaussian backend's (GT-N01); storage by photon-number blocks (`dm_storage => blocks`, B1); gate plans, 5–25×
  faster (B2).
- **Measurements:** mid-circuit homodyne and PNR with post-selection (D1); feed-forward.
- **Noise** (`pq_noise`, N1): loss (at inputs or outputs, or by trajectories), partial distinguishability (HOM
  visibility or a Gram matrix), multi-photon sources, dark counts, thermal loss.

### Sampling
- **Boson sampling:** exact `sample_boson` (Clifford–Clifford in O(n 2^n), and brute force), chi-square tested (S6),
  the stream order pinned by an independent reference (S7); lossy and collision-free variants; a Metropolis
  sampler.
- **Gaussian boson sampling:** photon-counting and threshold probabilities, chain-rule samplers with the truncated
  mass, equal to thewalrus within 1e-14; prefix caching.
- **Statistics:** mean photon number, g2, total variation distance, histograms, fidelity (GT-X01).
- **GPU:** Glynn's permanent and the hafnian on an NVIDIA GPU (`device => gpu | auto`), bit-identical to the CPU
  kernels (G1).

### Jobs and distribution
- **The job server** (`pq_parallel`, `pq_job`, `pq_worker`): chunks, at most `max_inflight` in flight, retries,
  cancellation, timeouts, partial results; identical samples for any number of workers and any chunk size.
- **Several nodes** (`pq_node_registry`, `pq_spec_cache`): discovery and handshake, chunks on other nodes
  retried elsewhere when a node is lost, the memory pre-check, kernel isolation on peer nodes.
- **The facade** `pq_api`: backend choice, run, sample, evolve, exact probabilities, `run_batch`; telemetry events.

### Symbolic engine (`pq_symbolic`)
- **Canonical form:** exact Gaussian rationals, radicals and exponentials, trigonometric functions as exponentials;
  substitution with domains; differentiation (gradients, Jacobians, higher orders); exact or floating evaluation;
  matrices; ASCII and LaTeX printing; a parser with atom safety; bounded execution.
- **Physics:** the symbolic unitary of passive circuits with a decisive unitarity check; symbolic permanents;
  symbolic Gaussian evolution (FR-G04).
- **Compiled programs:** straight-line programs and `equiv`.
- **Speed:** compared with FORM every week (`comparison/form`); reductions compared with FIRE on FIRE's own
  examples, answers checked equal (`bench/ibp.sh`).

### Reducing results (`pq_recon`, `pq_ibp`, `pq_simplify`)
- **Reconstruction over finite fields:** polynomials, rational functions and linear systems reconstructed exactly
  from evaluations modulo primes (`pq_recon`), with denominator guessing and a learned, replayed elimination; the
  samples spread over local workers or Erlang nodes (`pq_rexec`), the same answer on any of them.
- **Feynman-integral reduction:** a family of integrals as data (`pq_ibp_family`), its IBP identities, zero sectors
  and symmetries; `pq_ibp:reduce/3` reduces targets to master integrals with automatically chosen seeds;
  `pq_ibp:differential_equation/4` derives the masters' differential equations. On FIRE's own double box, 3–7×
  faster than FIRE on one core, with the same answers (`comparison/fire`, decision record R1).
- **Simplification by verified rewriting** (`pq_simplify`): factors from an alphabet, coefficients collected by
  monomials, every result checked against the input.

### Gradients
- Central differences; symbolic gradients for passive circuits and Gaussian inputs; adjoint gradients on the Fock
  backend, through loss too, and exact parameter-shift rules (AD1).

### Python and services
- A JSON service over stdio and HTTP (`pq_service`) and the Python client `pqcf` (`python/`, C1).

### Quality
- The verification plan: million-shot chi-square tests of every sampler, property tests, xref, Dialyzer,
  AddressSanitizer runs.
- The user guide `docs/guide.md`, whose examples run as tests; the API review and freeze.

## One computational foundation

PQCF is one framework, not several projects under one name. Its photonic simulator and its symbolic work, Feynman
integrals included, rest on the same four layers. Each layer has its own checks.

| Layer | What it is | Used by the photonic simulator for | Used by the symbolic work for | How it is validated |
| --- | --- | --- | --- | --- |
| **Symbolic engine** (`pq_symbolic`) | Exact canonical expressions, matrices, derivatives | Parameterised circuits: the exact unitary of an interferometer whose angles are symbols, symbolic permanents and Gaussian states, exact gradients | Large expressions: the FORM benchmarks; the coefficients of integration-by-parts identities | Term-for-term agreement with FORM 5, run weekly (`comparison/form`); property tests |
| **Compiled programs** (`pq_symbolic:compile/3`) | An expression compiled once into a straight-line program, evaluated many times | Evaluating a circuit's symbolic results (amplitudes, probabilities, gradients) under many parameter values, for the user (guide, section 11) | The black boxes that reconstruction samples | Fewer operations than FORM's best level on its resultant benchmarks; values against the expression tree |
| **Reconstruction over finite fields** (`pq_ff`, `pq_recon`, `pq_simplify`) | Exact results rebuilt from values modulo primes, never building the large intermediate | Exact closed forms of photonic quantities, such as output amplitudes as polynomials of the circuit's parameters (guide, section 17) | Reduced results; linear systems; Feynman integrals (`pq_ibp`) | Checks at fresh points in fresh primes; the double box equal to FIRE's, master by master; textbook reductions in the tests |
| **Execution on the BEAM** (the job layer for sampling; `pq_rexec` and `pq_rstore` for reconstruction; Rust NIFs) | Two modules on one design: work spread over processes, cores and Erlang nodes, the same result for any number of them, a failed node's work retried elsewhere; numerical kernels in Rust on dirty schedulers | Boson-sampling and Gaussian-sampling shots (`pq_job`, `pq_parallel`) | Reconstruction samples, the order search | The same answer on any number of workers and nodes, and node failures, in the tests; NFR benchmarks (`bench/HARDWARE.md`, provisional) |

Why integral reduction lives here: reducing large exact results is the step after computing them, for photonic
expressions as for any other. Integral reduction is the field's hardest public benchmark for that machinery. It
has established tools to measure against (FORM, FIRE), so it tests the shared layers far harder than the photonic
expressions do today. The photonic simulator itself is validated separately, against The Walrus, Strawberry Fields
and Piquasso (below).

## 1.4.0: reducing large results

Reducing results by reconstruction over finite fields (`CHANGELOG.md`, decision record R1):

- **Reconstruction:** polynomials, rational functions and linear systems from evaluations modulo primes
  (`pq_recon`), the samples spread over workers or nodes (`pq_rexec`, `pq_rstore`), the same answer on any of them.
- **Feynman integrals** (`pq_ibp`, `pq_ibp_family`):
  - IBP identities, scaleless sectors, sector symmetries;
  - automatically chosen seeds; only the equations the targets need (mark and sweep);
  - denominator guessing; a learned, replayed elimination;
  - the masters' differential equations.
- **Against FIRE 6.5.2** on its own double box, one core: 1.21 s against 7.99, and 3.73 s against 11.99 for a
  harder target, the same answers (`comparison/fire`).
- **Simplification** by verified rewriting (`pq_simplify`).
- **A negative result,** recorded: a physics-informed neural network evaluating the box's masters from their
  differential equation reached 1e-4 after 52 s of training, against a classical solver's 1e-13 in 19 ms
  (`comparison/pinn`).

## 1.3.0: the symbolic engine against FORM

The symbolic engine against FORM 5 (`CHANGELOG.md`, `docs/decisions/` K1, K2, H1–H4):

- **Expansion:**
  - polynomial products and powers in Rust (K1);
  - polynomials kept in Rust between operations (K2);
  - dense products about 4× faster than FORM; sparse products, powers and the symbolic mesh 1.5–4.3× FORM's time
    (from 1.3–32× on 2 Nov).
- **Compiled programs:**
  - Horner schemes with primitive nodes and subtraction (H1, H4);
  - common subexpressions across them (H2);
  - a parallel, deterministic search over variable orders, local or Monte Carlo tree search (H3, H4);
  - on the resultant benchmarks: below FORM's `O2` at the default order, and below FORM's `O4` on all three with
    the tree search.

## Next

- The multi-node NFR-02 benchmark on separate machines, and the reference hardware for the benchmarks
  (`bench/HARDWARE.md` is provisional: a shared ARM machine, no x86-64 reference yet).
- Feynman-integral reduction spread across many machines, and on FIRE's larger examples (nbox2w; families with
  linear propagators).
- Regular comparisons with PennyLane and Strawberry Fields (Xanadu) and Perceval (Quandela).

## Documentation

The symbolic engine is compared with FORM every week: [`comparison/form`](comparison/form) (latest results in
`comparison/form/RESULTS.md`).

From Python: the `pqcf` package in [`python/`](python/README.md), a client of PQCF's JSON service (guide section 16).

The user guide is [`docs/guide.md`](docs/guide.md); every example in it is a test (`test/guide_examples_tests.erl`).
Design decisions and spike results are in `docs/decisions/`.

## Building and testing

Everything builds in a container with OTP 27 and the pinned Rust toolchain (`native/pq_nif/rust-toolchain.toml`):

```sh
docker build -t pqcf-dev .
docker run --rm pqcf-dev            # EUnit (incl. property and 1,000,000-shot statistical tests), Dialyzer, xref, cargo test, clippy
docker run --rm --gpus all pqcf-dev # the same with an NVIDIA GPU: the GPU kernels' tests run instead of their refusals
```

Sanitizer runs (AddressSanitizer over the Rust tests, nightly compiler in its own container): `sh tools/sanitize.sh`.

Without the container: OTP 27 or later, rebar3 and Rust 1.90, then `rebar3 eunit`. Supported targets are
`x86_64-unknown-linux-gnu` and `aarch64-unknown-linux-gnu`; no kernel uses a system BLAS or LAPACK.

## How results are checked

Every kernel is tested against an independent implementation: the golden vectors in `test/fixtures/` are
generated by the scripts in `test/oracle/`, which share no code with the kernels (for example
`test/oracle/chacha8_ref.py`, ChaCha8 written from its definition and checked against the RFC 7539 test
vector). Numerical oracles for the physics are thewalrus, Strawberry Fields and Piquasso, with their versions
recorded in `test/oracle/VERSIONS`; the oracle container is `test/oracle/Dockerfile`.

The code is written with an AI coding assistant (Claude); the checks above are how its results are verified.

## License

Apache License 2.0 — see `LICENSE` and `NOTICE`. Copyright 2026 Quanterall Ltd.
