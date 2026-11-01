# S7: stream order and fixtures

**Decision: the draws and stream orders of FS 8.5 as implemented in `native/pq_nif/src/rng.rs` and
`boson.rs`, fixed by `test/fixtures/gt_s7.json`.** Closed on aarch64, 30 October 2026, and on x86_64, 31 October
2026: the first GitHub Actions run (commit 390a50a, runners ubuntu-24.04 and ubuntu-24.04-arm) passed the whole
suite, `gt_s7` included, on both targets.

## What is fixed

- The stream (Seed, Stream): ChaCha8, the key from SplitMix64 of the seed, 64-bit block counter and 64-bit
  stream (GT-R01).
- Uniform: `(x >> 11) * 2^-53` from one u64.
- Bounded integer in `0..s`: Lemire's multiply-and-reject on 64 bits (one u64 per try).
- Normal pair: Box-Muller, u1 rejected while 0, then u2; `r = sqrt(-2 ln u1)`, cosine first.
- Shuffle: Fisher-Yates from the last index down, `j = below(i + 1)`.
- `clifford_clifford`, per shot: the photons' input columns listed in mode order and shuffled, then one uniform
  per photon for the output draws. `brute_force`: one uniform per shot over the patterns in lexicographic order
  (mode 0 most significant). Shot i of a call uses stream (seed, stream_start + i).

## Exit test

`test/oracle/gt_s7.py` reimplements all of this in Python on the independent ChaCha8 of `chacha8_ref.py`
(including both samplers, with permanents by their definition). `native/pq_nif/tests/gt_s7.rs` reproduces it:
uniforms, bounded integers (bounds 1 to 2^64 - 1), shuffles and 40 samples per sampler for three inputs
(collision-free and bunched) bit for bit; normals within 1e-15 (they go through the math library: D2). The
NIF-level test `test/sample_boson_tests.erl` checks that splitting shots between calls changes no sample.

Both targets reproduce the fixtures bit for bit (CI, every push).
