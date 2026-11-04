# B1: density matrices stored by photon-number blocks

**Decision: with `dm_storage => blocks` the density-matrix backend keeps ρ inside as a grid of photon-number
blocks ρ[n][n′] and leaves out the blocks that are zero.** Every public result is unchanged:

- probabilities and samples are identical float for float;
- `evolve` returns the full matrix, as Decision 3 requires;
- the NIFs of FS 8.4 still take full matrices.

The default stays `dense` (Decision 3: "v1 default is the full matrix"). Phase 4, item 6 (FS 15 item 6). Released
after 1.0.0-rc.1, purely additive under the freeze: two options (`dm_storage`, `max_dm_entries`) and internal
kernels.

## Why it saves memory

The basis of FS 3.4 orders states by total photon number, so sector n is a contiguous index range and a block is a
plain sub-matrix. Number-conserving gates (beamsplitters, phase shifts, Kerr, passive unitaries) act inside each
sector, and loss maps block (n, n′) to (n − l, n′ − l). Both keep n − n′, so a Fock input under such a circuit keeps
only its diagonal blocks: Σ d_n² entries instead of D². Displacement, squeezing and thermal loss fill other blocks,
and then the saving shrinks towards none.

| Modes, cutoff | D | Dense | Diagonal blocks |
| --- | --- | --- | --- |
| 4, 16 | 4,845 | 0.38 GB | 0.046 GB |
| 4, 20 | 10,626 | 1.81 GB (beyond `max_dm_dim` = 4096) | 0.18 GB |
| 3, 30 | 5,456 | 0.48 GB | 0.027 GB |

## Identical results

Each block operation (`native/pq_nif/src/fockdmb.rs`) is the dense kernel's operation in the same order:

- **Gates:** G on every column, then on every conjugated row, visiting only the sectors with a stored block. The
  others are zero, and a gate maps a zero vector to exactly zero.
- **Loss and thermal loss:** a source entry is visited in the dense order (row i, then column j, sector by sector)
  and sends the same weights to the same targets. Accumulations into shared targets happen in the same order.

So every stored entry is the dense entry bit for bit. The only bit difference in the full matrix is the sign of
zeros: the dense kernel leaves −0.0 in entries that the blocks never store.

**Limits:**

- `init` checks D against `max_fock_dim`, as the amplitudes do, instead of `max_dm_dim`.
- After every operation the stored entries are checked against `max_dm_entries` (default 4096², the dense limit's
  memory), and exceeding it is `dimension_limit`.
- `evolve` presents the full matrix, so it stays within `max_dm_dim`.
- Probabilities, sampling and mid-circuit measurements read the blocks directly
  (`pq_simulator:evolve_state/3` internally), so they never build the full matrix.

**The state between instructions** is an Erlang term `{dm_blocks, [{N1, N2, Binary}]}` in the `rho` field. A plain
term rather than a NIF resource, because the state before the first measurement travels in the job specification to
workers and other nodes. The block NIFs (`fock_dmb_*`) are additional internal entry points. The FS 8.4 contract and
its full-matrix arguments are untouched, which is how this record reads Decision 3's "every NIF MUST present the
full matrix": the contract's NIFs, not internal ones.

**Time:** the same as dense when the blocks are full. Fewer columns and rows are visited when they are not. A gate
still acts on full-length vectors (O(D) each), so the gain is memory first. Cutoff 20 on 4 modes, two beamsplitters
and two losses, took 10–12 s here.

## Tests (`test/phase4_dm_blocks_tests.erl`, Rust `fockdmb::tests`)

- **Rust:** every operation against the dense kernel, bit-identical:
  - two beamsplitters, Kerr, loss, thermal loss, a displacement, a measurement;
  - the Fock input keeping only diagonal blocks (fewer than a third of the entries).
- **Erlang, identical results:**
  - exact probabilities identical for three circuits: number-conserving with loss; displacement, squeezing and
    thermal loss; a pure state promoted at its first loss;
  - `evolve`'s full matrices equal entry for entry (signed zeros apart), with the same warnings;
  - samples with mid-circuit measurements identical, on `fock_dm` and on `fock` with `promote_dm`.
- **Cutoff 20 on 4 modes** (20 photons, two beamsplitters, two losses): the dense storage refuses
  (`dimension_limit`, 10626 > 4096). The blocks hold 11,480,634 entries in 21 diagonal blocks (0.18 GB against
  1.81 GB); the probabilities sum to 1, and the mean photon number of mode 0 matches the closed form within 1e-10.
- **Limits:** `max_dm_entries` refused; `evolve` beyond `max_dm_dim` refused.
- **The whole EUnit suite with `blocks` as the default** (a scratch copy): 312 of 314 pass. The two others are the
  tests that assert the dense limit itself.
