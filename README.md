# PQCF — Photonic Quantum Computing Framework

A simulator for photonic quantum computing: Gaussian and Fock-space backends, boson sampling and Gaussian boson
sampling, a symbolic engine for parameterised circuits, and parallel, reproducible sampling across nodes. The
top layer is an Erlang/OTP application (`pq`); the numerical kernels are a Rust library loaded as a NIF
(`pq_nif`, through [Rustler](https://github.com/rusterlium/rustler)).

**Status: Phases 0-2 done except for the benchmarks that need the reference machines; Phase 3 (hardening) in progress.** Nothing here is ready for use yet.

| Done | Next |
| --- | --- |
| Project skeleton; NIF loading and version check; cancel tokens; the ChaCha8 random stream `rng_u64` (GT-R01); core types, configuration, error taxonomy and the packed binary codec (`pq_core`); the permanent family — `permanent`, `permanent_repeated`, `bs_probability` — with fixed-block determinism (GT-P01 against thewalrus); the core of the symbolic engine `pq_symbolic` — exact canonical form (Gaussian rationals, radicals, exponentials; trigonometric functions as exponentials), substitution with domains, differentiation (gradients, Jacobians, higher orders), evaluation (complex doubles, or exact), matrices, ASCII and LaTeX printing, a parser with atom safety, parameter resolution, bounded execution; circuits `pq_circuit` — construction and editing, validation (modes, parameters, numeric domains and unitarity, inputs, measured-mode reuse, custom instructions), effective bindings, composition with mode maps and renames, symbols, safe binary serialisation, the gate library; the Fock-space matrices of the beamsplitter, displacement and squeezing, checked against 80-digit closed forms (spikes S5, S8); exact boson sampling — `sample_boson` with `clifford_clifford` and `brute_force`, chi-square tested (S6), with the stream order pinned by an independent reference (S7); the hafnian, loop hafnian and torontonian kernels, equal to thewalrus within 1e-15 of their scale up to M = 40 and 20 clicks (S3, S4). Phase 1: states (`pq_state`: the Fock basis index of FS 3.4, Gaussian preparation, reduced states, physicality), lowering and capability checks, the Gaussian backend with terminal homodyne sampling (`pq_gaussian`, `pq_simulator`); the Clements mesh (`pq_circuit:mesh/1`); the pure-state Fock backend (`pq_fock`: beamsplitter sectors, phase, Kerr, passive unitaries, displacement and squeezing with the norm deficit, PNR, threshold and mixed sampling); `pq_stats:mean_photon_number`; boson sampling of passive circuits through the job server (`pq_sampling`, `pq_parallel`, `pq_job`, `pq_worker`: chunks, at most max_inflight in flight, retries, cancellation, timeout; identical samples for any number of workers and chunk size). Phase 2: Gaussian boson sampling (`gaussian_to_amat`, photon-counting and threshold probabilities, the chain-rule samplers with the truncated mass), equal to thewalrus within 1e-14; density matrices and loss (`pq_fock_dm`: the Kraus sum of A.11, `promote_dm`), photon statistics after loss equal to the Gaussian backend's (GT-N01); mid-circuit homodyne and PNR measurements (GT-G02) with post-selection (`docs/decisions/D1-midcircuit-streams.md`); noise insertion (`pq_noise`); statistics — g2, total variation distance, histograms, fidelity (GT-X01); lossy boson sampling (loss at the inputs or the outputs) and the Metropolis sampler; the symbolic unitary of passive circuits with a decisive unitarity check; Kan's repeated-index hafnian; distribution across nodes (`pq_node_registry`: pg scope pq, the handshake; `pq_spec_cache`; chunks on other nodes through erpc, retried elsewhere when a node is lost; the job memory pre-check; kernel isolation on peer nodes); the facade `pq_api` (backend choice, run, sample, evolve, exact probabilities, run_batch). Phase 3: JSON interchange of circuits (Appendix B); gradients (`pq_gradient`: central differences, and symbolic for passive circuits through the symbolic unitary and permanent); converters to and from Perceval and Piquasso (`pq_convert`, conventions asserted against the libraries in `test/oracle/converter_conventions.py`); the trajectory loss mode (pure states, the density matrix's distribution) | The multi-node NFR-02 benchmark on separate machines; Phase 3: performance, documentation, API freeze; the rest of the gate library; later: the rest of the symbolic engine; the benchmark harness and reference hardware; S1 on the reference machine and S7 on x86_64; the decision records are in `docs/decisions/` |

## Building and testing

Everything builds in a container with OTP 27 and the pinned Rust toolchain (`native/pq_nif/rust-toolchain.toml`):

```sh
docker build -t pqcf-dev .
docker run --rm pqcf-dev            # EUnit tests and clippy
```

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
