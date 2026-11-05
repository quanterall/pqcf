# B2: Fock gate plans

**Decision: a Fock-space gate applied to many vectors is prepared once as a plan (`native/pq_nif/src/fockplan.rs`)
and then applied. On photon-number blocks, a number-conserving gate is applied only to the sectors a column
actually has.** The project owner chose this as the second item beyond the roadmap (5 Nov 2026). It is internal:
no interface changes, and every result is the same bit for bit.

## Why

A density-matrix step ρ → GρG† applies G to all D columns and then all D rows, so 2D vector applications. The
per-vector functions of `fockstate.rs` redo their setup on every call:

- the basis indices of every group, O(D m) per call;
- the beamsplitter's sector matrices, one eigendecomposition per photon number;
- for a passive unitary, the whole Clements decomposition;
- for displacement and squeezing, the truncated matrix (the Laguerre or Bargmann recurrences).

That setup dominated. A plan holds it:

- the index groups and their small matrices;
- the diagonal factors of a phase shift or Kerr gate;
- the mesh's beamsplitters and phases.

`apply` then performs exactly the per-vector function's arithmetic, in the same order.

**Sectors.** The basis of FS 3.4 orders states by total photon number. A beamsplitter, phase shift, Kerr gate or
passive unitary maps each sector into itself, so a column whose only stored blocks are in sectors S needs only the
groups of S (`apply_sectors`). A block-diagonal column then costs d_n instead of D. Displacement and squeezing mix
sectors, so they are applied in full.

## Where plans are used

- the dense density-matrix kernel, `fock_dm_apply_*`;
- the photon-number blocks, `fockdmb::Blocks::unitary_plan`;
- the density-matrix adjoint: forward, inverse and adjoint passes.

The pure-state backend applies each gate once per state, so it gains nothing and is unchanged.

## Measured (`fockplan::tests::speed`, the X925 cores of the provisional reference)

| Case | Per-vector | Plan | Speed-up |
| --- | --- | --- | --- |
| dense, 3 modes, cutoff 10 (D = 286): beamsplitter | 0.017 s | 0.001 s | 22× |
| the same: passive 3 × 3 | 0.051 s | 0.002 s | 23× |
| the same: displacement | 0.004 s | 0.001 s | 5× |
| the same: Kerr | 0.001 s | < 0.001 s | 7× |
| blocks, 4 modes, cutoff 16 (D = 4845), diagonal blocks: beamsplitter | 1.33 s | 0.053 s | 25× |

End to end (`phase4_dm_blocks_tests`, cutoff 20 on 4 modes: two beamsplitters and two losses through
`pq_api:probabilities`): 10–12 s → 1.2 s.

## Tests

- **`fockplan::tests`:** every plan (beamsplitter, phase shift, Kerr, passive, displacement, squeezing) bit-identical
  to its per-vector function on a random vector; sector-restricted application equal to full application on a
  vector supported in those sectors.
- **`fockdmb::tests`:** the block storage with plans equal to the dense kernel entry for entry, through
  beamsplitters, Kerr and a displacement.
- **The benchmark:** asserts equality of the results it times.
- **Unchanged and passing:** every density-matrix test of the suite, the density-matrix adjoint's
  finite-difference check, and the block tests.
