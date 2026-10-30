# S2: the dense linear-algebra crate

**Decision: `nalgebra` 0.35.** Status: closed, 30 October 2026. Reversible by a later spike if profiling shows the
dense operations dominate (Phase 3).

## Question (Algorithms and Numerics Specification, section 10)

Which pure-Rust crate, `nalgebra` or `faer`, for the Hermitian eigendecomposition, complex Cholesky and complex LU
(determinant, inverse) at sizes 2 to 200? Exit test: agrees with a numpy oracle to 1e-12; builds on both target
triples with no system library; bitwise repeatable single-threaded; speed within 3x of the faster; `nalgebra` on a
tie.

## Method

`spikes/s2_linalg`: deterministic test matrices for each size (a random Hermitian `H`, a Hermitian positive definite
`A = B B^dagger + n I`, a general complex `G`), the four operations in both crates, each run twice and compared bit
for bit, timed (median of 7 to 201 runs) on one core with faer set to sequential (`Par::Seq`). `compare.py` checks
every result against numpy 2.4.6 (LAPACK) in the oracle container. Errors are relative to the largest entry of the
reference; the inverse error is also shown divided by the condition number `kappa` of `G`.

    cargo build --release && ./target/release/s2_linalg gen > matrices.json && ./target/release/s2_linalg run > results.json
    python compare.py          # in the oracle container

Machine: aarch64 (a 20-core Arm machine shared with other work), one pinned core; nalgebra 0.35.0, faer 0.24.4,
Rust 1.90.0. Not yet run on x86_64 or on the reference machine.

## Results

| n | crate | eig | Cholesky | log det | inverse | inverse / kappa | kappa | repeat | eig, us | Cholesky, us | LU, us | inverse, us |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | nalgebra | 7.7e-16 | 3.2e-18 | 2.2e-16 | 8.4e-17 | 6.4e-17 | 1 | yes | 0.2 | 0.1 | 0.1 | 0.1 |
| 2 | faer | 2.6e-16 | 1.5e-16 | 0.0e+00 | 1.7e-16 | 1.3e-16 | 1 | yes | 0.4 | 0.1 | 0.2 | 0.3 |
| 4 | nalgebra | 2.4e-16 | 1.3e-17 | 2.1e-16 | 3.1e-16 | 6.0e-17 | 5 | yes | 0.7 | 0.1 | 0.1 | 0.3 |
| 4 | faer | 8.3e-16 | 2.0e-16 | 0.0e+00 | 2.8e-16 | 5.4e-17 | 5 | yes | 1.1 | 0.2 | 0.2 | 0.6 |
| 8 | nalgebra | 1.1e-15 | 1.4e-16 | 4.4e-16 | 1.8e-15 | 3.7e-17 | 48 | yes | 2.6 | 0.3 | 0.3 | 1.4 |
| 8 | faer | 9.0e-16 | 1.4e-16 | 4.4e-16 | 1.6e-15 | 3.3e-17 | 48 | yes | 3.3 | 0.3 | 0.6 | 1.5 |
| 16 | nalgebra | 2.2e-15 | 2.0e-16 | 1.1e-16 | 1.5e-15 | 7.0e-17 | 21 | yes | 13.2 | 0.9 | 1.6 | 8.8 |
| 16 | faer | 2.0e-15 | 2.0e-16 | 3.3e-16 | 1.1e-15 | 5.1e-17 | 21 | yes | 12.9 | 0.7 | 1.4 | 4.0 |
| 32 | nalgebra | 2.7e-15 | 2.9e-16 | 1.5e-16 | 2.5e-15 | 6.3e-17 | 40 | yes | 81.4 | 4.8 | 10.8 | 48.0 |
| 32 | faer | 1.8e-15 | 2.9e-16 | 4.5e-16 | 3.6e-15 | 9.1e-17 | 40 | yes | 58.0 | 2.9 | 5.1 | 15.0 |
| 64 | nalgebra | 2.9e-15 | 4.1e-16 | 1.6e-16 | 1.2e-14 | 9.5e-17 | 124 | yes | 558.3 | 31.8 | 82.2 | 311.2 |
| 64 | faer | 3.9e-15 | 6.1e-16 | 3.1e-16 | 9.5e-15 | 7.7e-17 | 124 | yes | 293.5 | 17.3 | 25.0 | 74.4 |
| 100 | nalgebra | 4.9e-15 | 4.9e-16 | 3.2e-16 | 4.0e-14 | 6.7e-17 | 598 | yes | 1951.4 | 113.9 | 309.6 | 1099.0 |
| 100 | faer | 1.0e-14 | 3.3e-16 | 4.8e-16 | 2.7e-14 | 4.6e-17 | 598 | yes | 880.5 | 45.2 | 75.9 | 227.7 |
| 128 | nalgebra | 5.5e-15 | 7.2e-16 | 2.1e-16 | 3.1e-14 | 9.9e-17 | 308 | yes | 4023.9 | 235.8 | 655.0 | 2246.5 |
| 128 | faer | 9.7e-15 | 7.2e-16 | 6.4e-16 | 2.4e-14 | 7.7e-17 | 308 | yes | 1671.9 | 85.6 | 144.1 | 428.2 |
| 200 | nalgebra | 9.3e-15 | 5.8e-16 | 3.4e-16 | 3.2e-14 | 5.3e-17 | 605 | yes | 14661.1 | 882.6 | 2462.2 | 8163.0 |
| 200 | faer | 8.0e-15 | 4.6e-16 | 1.1e-16 | 2.3e-14 | 3.9e-17 | 605 | yes | 4601.3 | 260.3 | 477.5 | 1440.5 |

## Against the exit test

| Criterion | nalgebra | faer |
| --- | --- | --- |
| Agrees with numpy to 1e-12 | yes: at most 4e-14 (inverse), 1e-14 elsewhere | yes: at most 3e-14 |
| Builds for x86_64 and aarch64 with no system library | yes (cross-built for x86_64 without a C toolchain) | yes with `default-features = false, features = ["std"]`; the default `npy` feature pulls in a C build step (`psm` through `pest`) |
| Bitwise repeatable, single thread | yes | yes, on one machine |
| Within 3x of the faster | no at n >= 64: LU 3.3-5.2x, inverse 4.2-5.7x, Cholesky 3.4x and eigendecomposition 3.2x at n = 200 | is the faster |

## Why nalgebra anyway

The exit test does not name a criterion that decides against faer, but the specification's numerical policy does:

- **D1 across machines of one target triple.** faer's SIMD layer (`pulp`) chooses its instruction set at run time
  (`is_x86_feature_detected!`, `is_aarch64_feature_detected!`). Two x86_64 machines with and without AVX-512 run
  different code paths, so the same call can differ in the last bits. The node handshake (FS 9.5) compares the
  target triple only, so such a cluster would silently break D1. faer offers no feature to pin the path.
- **No fused multiply-add (AN 2.2).** The x86-v3 path that faer dispatches to includes FMA; the policy forbids
  `mul_add` so that the two target triples differ only through the math library (D2).
- nalgebra is plain Rust: the operation order is fixed at compile time and LLVM does not contract `a * b + c`
  into FMA without fast-math flags.

The cost is speed on the dense O(n^3) operations at n >= 64. These run once per call (Gaussian maps, the per-call
setup of the samplers, conditioning), not inside the exponential kernels: the torontonian builds its own Cholesky
factors by bordering (AN 4.2), and the permanent and hafnian do not use the crate. At n = 200 an nalgebra
inverse takes about 8 ms on this machine.

## Follow-ups

- Run on x86_64 and on the reference machine; record cross-triple differences (D2).
- Revisit in Phase 3 if a profile shows these operations dominate a workload; faer would need its dispatch pinned.
