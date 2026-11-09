# Changelog

PQCF follows semantic versioning from the API freeze (`docs/decisions/API-freeze.md`). Until 2.0, every change to
the public API is additive.

## Unreleased

### Added

- Feynman-integral families with **massive propagators** and **any number of loops** (`pq_ibp_family`):
  - a propagator may be `{Sign, Momentum, MassSquared}`, D = Sign (q^2 - m^2), the mass squared an expression in
    the family's variables; the masses enter the scalar products' constants and the Symanzik polynomial F (zero
    sectors and symmetries see them);
  - U and F by determinant and adjugate for any number of loops (the two-loop special case is gone);
  - derivatives in a variable of the masses (the propagator's index raised), and vacuum families (no external
    momenta);
  - tests: the massive tadpole's I(2) = (d - 2)/(2 m^2) I(1) and its differential equation in m^2, a bubble with
    one massive line (zero sectors), and the three-loop massless banana (four lines, five numerators), whose
    differential equation in s is (3d/2 - 4)/s as dimensional analysis requires.
- **Sector symmetries with numerators** (`pq_ibp_family:expand/2`, `mapped/2`): for a sector equivalent to its
  representative, a change of the loop momenta k -> A k + B p (A integer, det A = +-1) carrying its propagators onto
  the representative's, found from the permutation (Pak) and an exact linear solve over the sign choices. Every
  integral of the sector, numerators included, then expands into integrals of the representative and its
  subsectors. Such sectors get no identities of their own, and a target in one is defined by its expansion. The
  kite's identities 2,274 -> 1,579; the double box 7,296 -> 6,807 and 18,366 -> 16,967, answers unchanged (equal to
  FIRE's). The maps may also move the external momenta (p1 <-> p2 and the like), when the move preserves every
  scalar product of the kinematics: the double box has 31 sector maps, its example target needs 3,018 equations
  (0.27 s on one core, FIRE 8.0 s). A sector whose dotted integrals expand into ones its representative's
  identities do not reduce is handled by the `sector_maps => false` option (every sector's identities, scalar
  symmetries only), which the automatic seeds try after `true` at each budget: I(1,1,1,1,1,1,1,-2,-1) takes
  1.28 s on one core this way (18,366 equations; FIRE 12.0 s).
- **Comparison with FORM's Forcer on FORM's own four-loop checks** (`comparison/forcer`; Forcer:
  [arXiv:1704.06650](https://arxiv.org/abs/1704.06650)): PQCF reduces the integrals of `forcer.frm`
  (`Forcer_example`, `Forcer_1`, `Forcer_2`: the four-loop non-planar topologies no1-no3 with dots and numerators,
  and a two-loop kite with insertions), and every answer is exact: equal to Forcer's results as epsilon expansions,
  and equal in MS-bar to the numbers FORM's checks assert (zeta_3 to zeta_7). Forcer is 17-25x faster on the
  four-loop checks (`Forcer_1`: 10.3 s against 175 s, from 918 s before the changes below); PQCF 13x faster on the
  three-loop example.
- **Reduction overheads:** the seed search's completeness check is one solve on a random combination of the
  targets (four loops: 37.7 s -> 5.2 s); `reduce_sum/3` reconstructs only the masters found non-zero;
  `pq_recon:linear/4` stores its program once per call and runs univariate columns concurrently. The seed search
  takes the global budgets first while they stay within 3x the smallest, then the sector-by-sector search: the
  double box's harder target is back to 1.4 s (3.0 s with the sector search first).
- **Seeds chosen sector by sector** (`pq_ibp:auto_seeds/3`): each master that is not a seed raises only its own
  sector's budget (`pq_ibp_family:sector_seeds/4`); the global budgets stay as the fallback, `sector_seeds =>
  false` restores them. Four loops: 374,416 equations -> 97,949 for one case.
- **`pq_ibp:reduce_sum/3`**: a combination of targets reduced as one result (one extra unknown and equation, one
  reconstruction); **`pq_recon:linear/4`** leaves out the zero right-hand sides (19,841 -> 21 on a four-loop
  system).
- **`pq_ibp_family`: insertions and numerators** -- `eps => #{Line => M}` (a line's power index + M epsilon, as a
  one-loop insertion leaves it) and `numerator/3` (an integral times scalar products, as integrals of the family).
- **Four-loop symmetry classes 75 s -> 0.7 s**: sector maps found by backtracking over lines with equal invariants
  instead of every permutation (11! for a four-loop top sector).
- **A reduction benchmark suite** (`bench/ibp.sh`, `bench/ibp_suite.erl`, `bench/fire_run.sh`): PQCF and FIRE
  6.5.2 on FIRE's own examples (the one-loop box, the two-loop vertex v2, the double box) and harder double-box
  targets, on the same core back to back. Every answer is checked equal to FIRE's (`fire_check:compare/5`,
  `tables/1`), and each run is compared with the previous one. PQCF-only cases: the kite, the massive tadpole and
  the three-loop banana. The first run: PQCF at 0.03–0.78 of FIRE's time on one core, every answer equal.
- **Smaller results** (decision record R1, "Smaller results"):
  - `pq_ibp:choose_masters/3`: a greedy search over master bases (each master swapped for a neighbour in its
    sector, kept when the simplified coefficients shrink), every candidate an exact reduction on the default
    reduction's seeds; `reduce/3` takes `masters => [Integral]` (preferred masters) and `pq_ibp:size_of/2` measures
    a reduction. The double box: 18 -> 12, 56 -> 32 and 102 -> 69 additions (12.7 s, 53 s and 295 s of search on
    16 cores); the expanded coefficients 25 -> 16, 79 -> 48 and 148 -> 103.
  - `pq_simplify`: multivariate partial fractions (Leinartas's decomposition for affine denominator factors, then
    the numerators reduced modulo the denominators; nested: in d, then in the kinematics) as a fourth form, chosen
    when smaller or forced by `form => partial`; `apart => [[Var]]` the variable sets tried.
  - `pq_simplify:discover/2` and `gcd_poly/2`: factors outside the alphabet, repeated in a coefficient or shared
    by several, found by exact multivariate gcds (primitive pseudo-remainder sequences) and a coprime refinement,
    to pass as `extra`.
- **Large systems: Laporta's elimination equation by equation** (`pq_nif:ff_pivots/4`): the equations from the
  simplest to the most complex, each reduced only while its leading integral already has a pivot, its other entries
  left as they are, the redundant ones dropped at once. The pivot columns (the masters) are the same; mark and
  sweep is unchanged. The column-by-column elimination it replaces carried every row along: on a double-box system
  of 194,739 equations (I(1,1,1,1,1,1,2,-2,-2): seeds with a dot and four numerators), 8.8 billion merged entries
  against a final 654,082, 31 s; now 3.2 s, and the target needs 18,543 equations instead of 67,663. That reduction:
  10.7 s on one core (FIRE 6.5.2: 12.8 s on one thread), 9.6 s on 20 (FIRE: 6.9 s on 20 threads), peak memory
  1.9-3.7 GB; FIRE's example targets 0.31 s and 0.79 s on one core (FIRE 8.0 s and 12.2 s), answers equal to FIRE's.
- The systems' coefficients are evaluated mod p in the processes that generate them (one copy of the equations
  instead of two): the large double-box system 6.45 -> 5.28 s on 20 cores, its whole reduction about 7.1 s (FIRE: 6.9
  s on 20 threads). A batched parallel learning elimination (equations reduced in parallel against earlier pivots,
  then inserted in order: exact, the same result as sequential) was measured and dropped: no faster, since
  neighbouring equations reduce mostly against each other's pivots.
- Fixed: `pq_ibp_family:canon/2` mapped integrals even in a family with `symmetries => false`; the automatic seeds
  could then never be complete for such a family. Found by the sector-map test, which compares mapped reductions
  with unmapped ones, master by master.

### Documentation

- After an external review of the repository (8 Nov 2026):
  - the user guide's introduction said PQCF was not ready and its API not frozen (stale since 1 Nov); it now
    states the freeze and the release;
  - a new guide section, 17, on reducing results: a photonic quantity (the permanent of a symbolic interferometer
    block) reconstructed over finite fields, a Feynman-integral reduction, a differential equation; every example
    is a test;
  - the README says that PQCF is not post-quantum cryptography;
  - "One computational foundation": the four layers (the symbolic engine, compiled programs, reconstruction over
    finite fields, execution on the BEAM) that the photonic simulator and the symbolic work share, what each serves
    on either side, and how each is validated;
  - "Next" brought up to date.

### Changed

- Programs are evaluated over finite fields in Montgomery arithmetic (`Mont` in `ffslp.rs`: three 64-bit products per
  multiplication, no 128-bit division, which is a library call on aarch64): a double-box sample 2.86 -> 2.14 ms,
  the same values.
- `pq_ibp:system/4` builds its systems from the identities' coefficients as exact polynomials
  (`pq_ibp_family:identities_poly/2`, the zero test and symmetry name of each integral computed once per call),
  evaluated mod p directly, only the needed rows turned into ASTs; Shoup multiplication in the search-based
  eliminations: the harder double-box target's seeds and system 1.12 -> 0.63 s, the whole reduction on one worker
  2.04 -> 1.63 s, the same answers (equal to FIRE's).
- Reductions use the cores: the systems' identities and their values are generated in parallel chunks (the
  family's caches computed once and handed to the chunks, `pq_ibp_family:warm/2`, `load/1`); the search-based
  eliminations update their column index only where a row changes (fill-in, cancellations), which made the
  learning elimination about 4x faster; a system's solved points are cells that concurrent requests wait on, and
  its learned elimination is shared by its black boxes; the guessing slices run as parallel tasks. Two double-box
  targets with four numerators each: 9.4 s on one core before, now 6.3 s on one and 4.4 s on 20 (the harder single
  target 1.14 s on 20); answers unchanged (equal to FIRE's). Reconstructing the masters' coefficients concurrently
  was measured slower than one after another (each then repeats the others' setup, and points are shared anyway),
  so they stay sequential.

## 1.4.0 (8 November 2026)

Reducing large symbolic results, additive under the API freeze (decision record `R1`). Results are reconstructed
exactly from evaluations modulo primes, spread over local workers or Erlang nodes with the same answer on any of
them. Feynman integrals reduce to master integrals through `pq_ibp`; on FIRE's own two-loop double box, PQCF is
3-7x faster than FIRE 6.5.2 on one core (1.21 s against 7.99, 3.73 s against 11.99), with answers checked equal
master by master. Results are simplified by verified rewriting (`pq_simplify`), and master-integral differential
equations are derived exactly.

### Added

- **Feynman-integral reduction as part of the library** (`pq_ibp`, `pq_ibp_family`; R1):
  - `pq_ibp_family:new/1`: a family as data (loop and external momenta, propagators with numerators, external
    scalar products as a map, variables, top sector), validated; `example/1` gives the benchmarks' box, kite and
    FIRE's double box; identities, scaleless sectors (Lee-Pak), sector symmetries (Pak), seeds and derivatives with
    respect to the invariants, cached per process by the family's hash, random points from local states;
  - `pq_ibp:reduce/3`: several targets at once, seeds automatic (`seeds => auto`, the default) or given
    (`{dots, D, S}`, Kira's `{R, S}`, a list), the masters' order (`order => numerators | dots`), pq_recon's options
    (executors, denominator guessing); a zero target maps to `[]`, a master to itself; `system/4`, `seeds/3` and
    `auto_seeds/3` for the system itself; `differential_equation/4`: the masters' equation in an invariant;
  - the double box's two targets in one call: 2.7 s on one core; the benchmark drivers (`comparison/form/ibp.erl`,
    `comparison/pinn/box_de.erl`, `comparison/fire/fire_check.erl`) now use the library, with the same results;

- **Reducing results by finite-field reconstruction, phase 1** (`R1`; the project owner's decision, 7 Nov 2026):
  - `pq_symbolic:compile/3` with `exact_constants => true`: constants known exactly stay exact in the program
    (float evaluation unchanged);
  - `pq_ff`: compiled programs evaluated modulo primes below 2^63 in Rust (`ffslp.rs`), a batch of points at a
    time, a division by zero marking a point as bad; deterministic primes (Miller-Rabin); inverses, Chinese
    remaindering and Wang's rational reconstruction;
  - `pq_recon:polynomial/3`: a polynomial reconstructed from its program's values (Zippel over the first prime,
    the support reused in the next, Chinese remaindering and rational reconstruction, a final check); the resultants
    res(7,4) and res(7,5) exactly, in 12.6 s and 248 s on one core;
  - the store and the executors (phase 3): `pq_rstore` (content-addressed blobs, fetched once by key, programs
    kept prepared) and `pq_rexec` (`local`, the default, or `{nodes, Nodes}` over erpc, a failed node's tasks run
    elsewhere); a reconstruction's columns run in parallel with the same result for any executor (res(7,4) in
    2.5 s, res(7,5) in 42 s on 16 local workers);
  - `pq_recon:rational/3` (phase 4): a rational function N / D reconstructed from its program, along lines with a
    shift, the coefficients by one multi-output Zippel, the shift taken out by a Taylor recursion; checked against
    FORM 5's `PolyRatFun` on sums of fractions (`comparison/form/fractions.erl`): the same results exactly, FORM
    300-12,000x faster on that task (an explicit sum, where adding fractions is the right tool);
  - `pq_recon:linear/4` (phase 5): the components of the solution of a linear system with polynomial entries, each
    reconstructed from eliminations mod p at the sample points (`pq_nif:ff_linsys/3`), never solved symbolically;
    against FORM 5's symbolic elimination under `PolyRatFun` (`comparison/form/linsys.erl`), the same answers, FORM
    12-18x faster when the answer is about as large as the elimination, PQCF 25 to 35,000x faster when it is
    compact (as in integration-by-parts reductions; n = 24: FORM 697 s, PQCF 0.02 s);
  - several right-hand sides at once (one per master integral): sparse elimination per point in Laporta order, all
    right-hand sides tested for 0 together; on the integration-by-parts identities of the one-loop massless box
    (`comparison/form/ibp.erl`, 86 to 497 equations), the same reductions as FORM 5 and, core for core, 1.9-4.5x
    faster than its symbolic elimination from about 260 equations up; on a two-loop family, the massless propagator "kite" (222 to 2,102
    equations, the kite reduced to the bubble squared and the sunset), 3-4.8x faster; `linear/4` passes its
    options' expression limits to the compiler, and takes A's rows sparse (maps Column => entry), only the non-zero
    entries compiled; B's rows may be sparse too (`rhs => M`), the right-hand sides carried as sparse columns of the
    elimination (`pq_nif:ff_linsys/7` with both patterns, a blob in the store), so thousands of masters cost their
    non-zero entries: the kite 11-38x faster than FORM on one core; `pq_nif:ff_pivots/3` (the pivots of a sparse
    matrix mod p, shortest row first); `comparison/form/ibp_family.erl`, IBP identities for any family from its
    propagators (scalar products solved, zero sectors by the Lee-Pak criterion, Kira-style seeds; it found and fixed
    a gap in the kite's zero rule); against FIRE 6.5.2 (FLINT) on its double box, FIRE reduces in 16 s what costs
    PQCF 56 s per sample point -- the gap is the Laporta front end (no symmetries, no equation selection), recorded
    in R1 with the steps that would close it; then closed: mark and sweep (`pq_nif:ff_pivots/4` records which
    rows each row was combined with, and only the equations the target needs are solved: 0.19 s per point instead
    of 56 s), sector symmetries by Pak's criterion (`ibp_family:canon/2`: the double box reduces to its 8 masters), a
    cache of solved points shared by a system's right-hand sides (prefix-stable sample points), and dot-free seeds
    (`{dots, D}`); the double box checked equal to FIRE's at a random point mod p, and reduced in 3.5 s on one core
    against FIRE's 7.7-8.4 s; seeds chosen automatically (`{auto, Family[, Target]}`: the smallest seeds whose
    reduction is complete, every non-zero master itself a seed): the double box in 3.95 s (FIRE 7.7-8.4 s), and
    I(1,1,1,1,1,1,1,-2,-1) in 32.2 s (FIRE 11.9-12.2 s), both checked equal to FIRE's;
  - `pq_simplify`: results simplified by verified rewriting (exact division by an alphabet of candidate factors,
    collecting by monomials, the smallest of the expanded, factored and collected forms by the additions written,
    checked at random points; `format/1` prints the structure); the double box's coefficients 25 -> 18 additions
    and 79 -> 58, in the textbook form (factors in d, monomials in s, t); the remaining numerators are irreducible
    (checked with FORM), so a network proposing factors has nothing to add on these results;
  - master-integral differential equations derived by PQCF (`ibp_family:derivative/3`, external-momentum
    operators reduced by IBP; `comparison/pinn/box_de.erl`): the one-loop box's textbook equation, checked to 1e-8
    against direct Feynman-parameter integration; a physics-informed neural network on it
    (`comparison/pinn/pinn_box.py`, after TechUPGRADE's PINN-PONN, in its JAX image): 1e-4 accuracy after 52 s of
    training against an ODE solver's 1e-13 in 19 ms, and 0.3-5.7% as a surrogate over unseen dimensions; recorded in
    R1 with what was positive and what was negative; `ibp:system/4` with `order => numerators | dots`;
  - denominator guessing in `pq_recon:linear/4` (default on, `guess => false` to turn off): univariate denominators
    on one line per variable matched against the factor alphabet, F * D reconstructed as a polynomial through a
    scaled black box (`pq_nif:ff_linsys/8`), falling back to the rational reconstruction when the guess fails its
    check; bit-identical results; the double box I(1,1,1,1,1,1,1,-2,-1) 32 -> 9 s reconstruction, whole runs on one
    core 3.16 s and 11.87 s against FIRE's 8.10 s and 12.21 s;
  - a learned elimination for linear systems (`Plan` in `ffslp.rs`: recorded at the first point, replayed on a
    flat array at the others, the wanted component through the transposed triangular system), `invmod` by the
    extended Euclidean algorithm instead of Fermat (every inversion in PQCF), Shoup multiplication in the replay:
    a double-box sample point 31 -> 2.9 ms, bit-identical results; whole runs on one core 1.21 s and 3.73 s against
    FIRE's 7.99 s and 11.99 s (6.6x and 3.2x); three executors
    behind one task interface (this node, Erlang nodes, FVLLMP federation peers), from a single node to 100+.

## 1.3.0 (7 November 2026)

The symbolic engine against FORM 5. Everything is additive; the default programs and canonical forms are the same as
1.2.0's, term for term. Expansion: dense products about 4 times faster than FORM, sparse products, powers and the
symbolic mesh 1.5–4.3 times FORM's time (from 1.3–32 times on 2 November). Compiled programs: on FORM's resultant
benchmarks, fewer operations than FORM 5's `O4` on all three (decision records K1, K2, H1–H4).

### Added

- **Horner schemes for compiled expressions** (`H1`):
  - `pq_symbolic:compile/3` with `#{optimize => horner}`: every polynomial is compiled as a multivariate Horner
    scheme, in a fixed occurrence order, with identical sub-polynomials shared;
  - on the resultant benchmarks of FORM's expression-simplification papers: 2.6–3.0× fewer operations than the
    default program, 4–6 % fewer than FORM 5's `Format O1`, and evaluation 3.4–3.6× faster;
  - the default is unchanged, bit for bit.
- **Common subexpressions across compiled programs** (`H2`):
  - `optimize => horner_cse` (Horner, then the pass) and `optimize => cse` (the default program, then the pass);
  - sums and products are treated as n-ary nodes, and the operand pair that most nodes share becomes one
    operation, repeatedly (FORM's `O2` step);
  - on the resultants: 4,285 / 18,166 / 61,404 operations, against Horner's 4,678 / 19,490 / 67,928 and FORM 5's
    `O2` 3,969 / 16,398 / 55,685.
- **A search over Horner variable orders** (`H3`):
  - `order_search => #{trials, seed, batch, nodes}` with `optimize => horner | horner_cse`;
  - batches of candidate orders are scored in parallel, one process each, on this node or spread over nodes, with
    the same program for any number of cores or nodes;
  - on the resultants, 160 trials give 4,273 / 16,408 / 54,643 operations, below FORM 5's `O2` on the two larger
    ones; 16 cores are 6.3 times faster than one.

- **Toward FORM's `O3`** (`H4`):
  - Horner nodes are computed as their primitive parts (content and sign taken out), so negated and scaled copies of a
    sub-polynomial share one slot;
  - a subtraction operation `{sub, Slot, A, B}` in Horner programs;
  - `order_search => #{strategy => mcts}`: a parallel, deterministic Monte Carlo tree search over order prefixes,
    with `exploration` and `direction`;
  - the search scores an order by counting its Horner operations in Rust (`pq_nif:horner_prepare/2`,
    `horner_count/2`), equal to the built program's count;
  - on the resultants: below FORM 5's `O2` at the occurrence order (3,771 / 15,875 / 52,423); with the tree search
    2,741 / 10,243 / 34,799 (5,000 / 5,000 / 10,000 trials), below FORM's `O4` (2,936 / 10,722 / 35,419) on all
    three.

### Changed

- **Polynomials that stay in Rust** (`K2`):
  - inside `canonicalise`, large plain products and powers stay in Rust between operations, and the result's
    expression tree is built there in print order;
  - the canonical form is the same, term for term; `poly_kernel => false` turns it off with the kernel;
  - against FORM 5: powers 4.3×, sparse products 2.4× (from 25–32× on 2 Nov), dense products about 4× faster.
  - exponentials exp(k i x) of symbols count as kernel variables, and `circuit_unitary` keeps its entries in Rust:
    the 8-mode symbolic mesh in 0.6 s (1.5× FORM 5, from 5×).
- **Products of polynomials in Rust** (`K1`):
  - `pq_sym_poly` multiplies plain polynomials in a NIF: symbols, conjugates and π with integer exponents, and
    exact coefficients;
  - exponents are packed into 64 bits and coefficients are exact i128, with a fallback to the Erlang product on
    overflow, on exponents that don't pack, or when no NIF is loaded;
  - powers of plain polynomials too (`pq_nif:sym_poly_pow/5`), with the multinomial enumeration as the fallback;
  - the result is the same polynomial, term for term (D1); `poly_kernel => false` turns the kernel off;
  - the canonical print order sorts once instead of twice, and by binary keys when every monomial is a product of
    symbols with integer exponents: 2.4 times faster on 591k terms, with the same output.

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
