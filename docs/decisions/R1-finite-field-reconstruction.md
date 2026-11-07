# R1: reducing large results by finite-field reconstruction, distributed

**Status: accepted, in progress (the project owner's decision, 7 Nov 2026: "go with B and the minimal store").
Phases 1–5 done.**

## Why

In large symbolic calculations, physicists often spend more time reducing a result than computing it:
- **the coefficients are rational functions:** of kinematic invariants, masses and ε;
- **reducing them needs heavy algebra:** combining many fractions into one reduced fraction needs multivariate GCDs;
- **intermediate expressions swell:** they grow far beyond the size of the final answer.

FORM normalises as it goes (`PolyRatFun`). The modern alternative (FiniteFlow, FireFly, Kira) never builds the
intermediate expressions:
1. Evaluate the result numerically, modulo large primes, at many points.
2. Reconstruct the reduced rational function from the values.

The samples are independent of each other, so the work spreads naturally over many machines. That is where PQCF's
Erlang distribution and deterministic jobs fit, and where a FORM-based workflow can be outperformed.

**Scale.** The work will run on many (100+) DGX Spark nodes of the FVLLMP federation. PQCF doesn't depend on
FVLLMP: it treats the federation as a pool of Erlang nodes it is given, and the network and node provisioning
belong to the federation.

## Design

### 1. The black box: compiled programs over a finite field
- **The function:** a result is given as a program (`pq_symbolic:compile/3`). Every operation is exact over ℚ: sums,
  products, integer powers (negative ones are divisions), and rational constants.
- **Exact constants:** compiling with `exact_constants => true` keeps each constant as a rational instead of a float.
- **Rust evaluation:** a program is prepared once, then evaluated modulo a 63-bit prime p at a batch of points.
- **Bad points:** a point where the program divides by zero mod p is reported as bad, and the reconstruction skips
  it.
- **Refused:** π, radicals, exponentials, conjugates and complex constants have no meaning in 𝔽_p.

### 2. Reconstruction
- **Polynomials (phase 1):** Zippel's sparse interpolation adds one variable at a time.
  - **First prime:** each stage interpolates the coefficients of the support found so far, from transposed
    Vandermonde systems, then interpolates each coefficient densely in the new variable. Degrees are found by early
    termination: interpolation stops when extra points agree.
  - **Later primes:** these reuse the support, so they need only one Vandermonde solve per prime.
  - **Back to ℚ:** Chinese remaindering, then Wang's rational reconstruction. More primes are added until one extra
    prime confirms the result.
- **Rational functions (phase 2):** the denominator is normalised, a univariate Thiele interpolation runs along a
  scaled direction, and Zippel then reconstructs the coefficient polynomials (the FireFly method, Klappert and
  Lange 2019).

### 3. Where the work runs: three executors, one task interface

The FVLLMP federation (fvllmp repo, `docs/FEDERATION-DESIGN.md`, read 7 Nov 2026) is not a cluster:
- **No central point:** independent DGX Sparks, with no coordinator or registry.
- **Pairwise trust:** owners pair by exchanging pinned self-signed certificates and member tokens.
- **One port:** nodes reach each other only through each node's HTTPS gate (`:8402`, mutual TLS).
- **Local work first:** a node lends only what its owner sets aside.
- **Accounting:** borrowed and lent work is counted on a two-sided ledger of credits.

Erlang distribution (a shared cookie, epmd and its ports, implicit full trust) therefore fits only inside one
owner's machines. So a task is plain data, and three executors run it:

1. **Local (the default; one node, no cluster):** processes on this machine, as many as there are cores. Everything
   works on a single node, with nothing to configure.
2. **Erlang nodes:** one owner's trusted machines, coordinated over `erpc` (described below).
3. **Federation peers (later):** the same tasks as requests to a paired node's PQCF service over HTTPS, the
   requesting node acting as coordinator. The tasks suit a federation:
   - **small:** a program hash, a prime and a range of point indices;
   - **cheap results:** a list of numbers;
   - **checkable:** every sample is a pure function of its inputs, so the requester spot-checks a peer by recomputing
     a few points, and the reconstruction rejects inconsistent data.

   **Open, outside PQCF:** reaching PQCF through the gate needs a route on FVLLMP's side, and its credits need a
   unit for PQCF work. That is the FVLLMP team's to decide and build; PQCF doesn't change FVLLMP.

The result is the same whichever executor runs the tasks, and however many nodes they use (D1).

### 4. Distribution and aggregation over Erlang nodes (built for 100+ nodes)
- **Topology:** a coordinator, a star of `erpc` calls to the worker nodes. Workers need no connections to each
  other (`-connect_all false` is enough), so there's no full mesh of 100² connections.
- **Nodes:** an explicit list (`nodes => [...]`), or the registry's members on a small cluster.
- **Tasks:** a task is (program hash, prime index, a range of sample indices, the point generator). The points
  are derived from (seed, prime index, sample index), so a worker generates them itself, and the messages stay small.
- **Determinism (D1):** the values are exact, and the points don't depend on which node computes them. So the
  result is the same for any number of nodes, any assignment of tasks, and any failures.
- **Faults:** a failed or timed-out task is retried on another node.
- **Aggregation:** the samples are assembled by index. Each prime's reconstruction is independent, so primes can be
  reconstructed on different nodes in parallel, and the coordinator combines them by Chinese remaindering.

### 5. The minimal store
- **Content-addressed:** blobs (programs, sample blocks, reconstructed results) are keyed by their SHA-256.
- **Placement:** in memory (ETS), or also on disk (`store_dir`).
- **Programs:** a worker fetches a program from the coordinator once, by hash, and caches it.
- **Resuming:** sample blocks are stored as they arrive, so a coordinator that restarts resumes without
  re-evaluating them.

## Phases

1. Finite-field evaluation of programs in Rust; prime and modular utilities; Chinese remaindering; rational
   reconstruction.
2. Zippel reconstruction of polynomials over ℚ, validated against the canonical expansion (the resultants).
3. The store and the executors (local, Erlang nodes), validated on peer nodes: the same result for any node count,
   and with a node lost.
4. Rational functions (Thiele + Zippel), and a benchmark: reducing a sum of many fractions, against FORM's
   `PolyRatFun`.
5. The federation executor, once FVLLMP provides the gate route.

## Phase 2, measured (7 Nov 2026)

`pq_recon:polynomial/3` reconstructs a polynomial from a compiled program (`exact_constants => true`):
- **First prime:** Zippel's interpolation, one variable per stage. Each stage solves transposed Vandermonde systems
  in Rust (`pq_nif:ff_tvandermonde/3`, O(n²)), then runs a Newton interpolation per coefficient, stopping when two
  further points add nothing.
- **Further primes:** one Vandermonde system each, on the known support.
- **Back to ℚ:** Chinese remaindering plus Wang's reconstruction, until two prime counts agree. The result is then
  checked at a fresh point in a fresh prime.
- **Points r¹ … rⁿ,** never r⁰ = (1, …, 1): there a fraction such as 1/(x − y) is undefined, which the first version
  hit (fixed by dividing the solution by the nodes).
- **Big constants:** passed to Rust as 64-bit limbs, of any size.

| From its `horner_cse` program (one core) | terms | variables | time | equal to the expansion |
| --- | --- | --- | --- | --- |
| res(7,4) | 2,562 | 13 | 12.6 s | yes |
| res(7,5) | 11,380 | 14 | 248 s | yes |

These are slow next to expanding them directly (0.65 s and 4.2 s): the resultants are polynomials, so expansion is
the right tool for them, and they serve here as a correctness check at scale. Reconstruction pays where expansion
swells, in rational functions, which is phase 4. The time goes into the Erlang side: the points, the Newton steps,
and big-integer arithmetic. Moving those into Rust, and the evaluations onto the executors (phase 3), is next.

Tests (`test/pq_recon_tests.erl`):
- 25 random sparse polynomials with rational coefficients, from default and `horner_cse` programs;
- a fraction that cancels to a polynomial, (x² − y²)/(x − y);
- coefficients of 2²⁰⁰, which need several primes;
- the zero polynomial, a constant, a power;
- a resultant;
- determinism, with another seed giving the same polynomial.

A Newton step that forgets the interpolant's value is caught (5 of 6 tests fail).

## Phase 3, measured (7 Nov 2026): the store and the executors

- **The task is a Zippel column.** For one value of the new variable, one task evaluates the program at |S| points
  and solves the O(|S|²) Vandermonde system, all in Rust (`pq_nif:ff_column/6`, with the points generated natively).
  It returns |S| coefficients. The columns of a stage are independent, so they are what the executors spread.
- **The coordinator works ahead, consumes in order.** `pq_recon` computes a batch of columns ahead (as many as the
  executor has slots) and takes them in order for Newton's method and the termination check. So the result
  doesn't depend on the batch, the executor or the number of nodes.
- **`pq_rexec`** runs tasks with one interface:
  - `local`, the default: a pool with one slot per dirty CPU scheduler;
  - `{nodes, Nodes}`: `erpc` from this node, as many slots per node as it has dirty CPU schedulers, interleaved
    across the nodes. A node that can't be reached, or fails, is dropped and its tasks are requeued; with every node
    gone, the work comes back here. A task's own exception is raised with its remote stack trace, never retried.
  - `slots => [Node]` gives the pool outright.
- **`pq_rstore`:** blobs keyed by SHA-256 of their deterministic form, in ETS and also on disk with `store_dir`.
  A worker fetches a program or support from the coordinator once, checks that it hashes to its key, and keeps the
  program prepared. On a node without the `pq` application, the store starts on first use. Two tasks arriving
  together on a fresh node raced on that start (a table read before `init` had created it); the start now waits for
  `init`.

| Reconstruction (local executor; the machine shared) | 1 worker | 16 workers | speed-up |
| --- | --- | --- | --- |
| res(7,4), 2,562 terms | 12.4 s | 2.48 s | 5.0× |
| res(7,5), 11,380 terms | 246 s | 42.4 s | 5.8× |

The same polynomials in every case. The speed-up is below 16 because the coordinator's share, Newton's steps in
Erlang, doesn't spread, and because the machine was shared.

Tests (`test/pq_rexec_tests.erl`, on three peer nodes):
- the same polynomial on this node with 1, 3 or 7 workers and other batches, on two peers, and on this node plus a
  peer;
- with an unreachable node in the list;
- with a peer stopped during the run;
- with every node gone (the work comes back here);
- the requeue made deterministic: unreachable nodes hold slots in the pool, and their tasks run on the others;
- the store: a peer fetches a blob by key, keeps it, and gets `not_found` for an unknown key;
- a failing task raises.

Recording a failed node's task as a result instead of requeueing it is caught.

**Not measured yet:** separate machines. Peer nodes on one host prove the protocol, not the network. The DGX Sparks
are the next measurement.

## Phase 4, measured (7 Nov 2026): rational functions, and the first benchmark against FORM

`pq_recon:rational/3` returns {N, D}: coprime, with D's coefficient of its last monomial (in Erlang's order of
exponent lists) equal to 1.
- **Lines:** along x = s + t (1, z₂, …, zₙ), for a random shift s (integers, the same in every prime), the function
  is n(t)/d(t) with d(0) = 1. Every coefficient of tʳ is homogeneous of degree r in z.
- **Degrees:** found by Thiele's continued fraction along one line, which grows until two further points are
  predicted.
- **The coefficients:** one multi-output Zippel in z₂ … zₙ reconstructs them all from the same points. The engine is
  shared with `polynomial/3`, with Newton interpolants per output and a union support. Each column is a Rust task
  (`pq_nif:ff_line_column/10`): for every point, the function along t, a small linear system for n and d, then one
  Vandermonde solve per output on the same nodes.
- **The shift** is taken out degree by degree from the top: N_k = n_k − Σ_{j>k} (s·∇)^(j−k) N_j / (j−k)!, which
  keeps the structure sparse instead of expanding N(x − s).

Tests (`test/pq_rational_tests.erl`), each result cross-multiplied exactly against the input:
- a fraction;
- 1/(xy), where the denominator vanishes at the origin and the shift is needed;
- 1/(x + y) + 1/(x − y) = 2x/(x² − y²);
- (x² − y²)/(x − y) = (x + y)/1;
- one variable;
- 12 random N/D in three variables with rational coefficients;
- determinism across seeds.

Rust tests cover the line solve.

**The benchmark: a sum of fractions** (`comparison/form/fractions.erl`). K terms, each a polynomial over a product of
up to three linear factors from a small set, in V variables, the shape of a physics coefficient. FORM 5 sums them as
`rat(n, d)` under `PolyRatFun`, adding with polynomial GCDs (FLINT). PQCF compiles the sum and reconstructs N/D.
Both results agree exactly, checked by cross-multiplication:

| | result | FORM 5 `PolyRatFun` | PQCF reconstruction (20 cores) |
| --- | --- | --- | --- |
| K = 20, V = 3 | numerator 787 terms, denominator 536 | 0.01 s | 5.5 s |
| K = 100, V = 3 | 1,136 / 810 | 0.03 s | 10.0 s |
| K = 100, V = 4 | 10,370 / 7,037 | 0.07 s (with a larger `MaxTermSize`) | 830 s |

**FORM is 300–500 times faster here, and that is the right result for this task.** The input is the sum written out,
with a modest common denominator. FORM sees every factor and adds fractions directly. Reconstruction samples a
black box whose degree is that of the whole result (here 13 to 20 along a line), so it does far more work than
needed. At V = 4, FORM first stopped with "output exceeds MaxTermSize" (under `PolyRatFun` the reduced sum is one term
rat(N, D)) and finished in 0.07 s with larger buffers. There FORM itself checked PQCF's result: the sum minus
rat(N, D) is 0. The smaller cases were checked by cross-multiplication in PQCF.

**Where reconstruction is the right tool:** when the reduced result is far smaller than the way to it. The textbook
case is solving a linear system whose entries are polynomials (integration-by-parts reductions, the bulk of
multi-loop calculations). Symbolic elimination swells, while an elimination modulo a prime at each point is cheap
and parallel, and only the solution is reconstructed. That is the method of FiniteFlow, FireFly and Kira. A
black-box task that solves such a system mod p per point is the next step, with that as the benchmark against
FORM.

## Phase 5, measured (7 Nov 2026): linear systems, never solved symbolically

`pq_recon:linear(A, B, Vars, Opts)` returns the reduced rational function {N, D} of each wanted component of the
solution of A u = b, whose entries are expressions.
- **The black box:** all entries compile into one program with n² + n outputs. At each sample point, the black box
  (`ffslp.rs` `BlackBox::Lin`, `pq_nif:ff_linsys/3`) evaluates them modulo p, solves by elimination, and returns
  u[want], or a bad point where A is singular there.
- **The same machinery:** the rational reconstruction, executors and store run unchanged on top. The column and
  line functions now take a `BlackBox` (a program's first output, or a linear system), and tasks name theirs.
  `pq_rstore:prepared/2` caches `{linsys, Key, N, Want}` once per node.

Tests (`test/pq_linear_tests.erl`) check A u = b exactly at random rational points, since the canonical form can't
combine fractions with different denominators:
- a 2 × 2 system (u₀ = y/(xy − 1));
- a sparse 5 × 5 system in three variables, every component;
- only the components asked for, the same as from a full solve;
- a singular and a malformed system refused.

Returning the wrong component is caught (3 of 4 tests fail).

**The benchmark** (`comparison/form/linsys.erl`). An IBP-shaped system, n × n: sparse, a full diagonal, every entry
linear in V variables. PQCF reconstructs u₁. FORM 5 eliminates symbolically under `PolyRatFun`, pivoting unknown k
on row k, k = n … 2, with the inverted pivot (`reverse_`). FORM itself checks every pair of results: U − rat(N, D)
is 0 in every case.

| right-hand side | n | V | answer (numerator / denominator terms) | FORM 5 (wall) | PQCF (20 cores) |
| --- | --- | --- | --- | --- | --- |
| random | 14 | 3 | 559 / 679 | 0.29 s | 3.9 s |
| random | 18 | 3 | 1,140 / 1,330 | 0.86 s | 15.3 s |
| random | 14 | 4 | 2,380 / 3,060 | 4.2 s | 50 s |
| b = A u*, u* compact | 14 | 3 | 2 / 1 | 0.52 s | 0.02 s |
| b = A u*, u* compact | 18 | 3 | 2 / 1 | 0.90 s | 0.01 s |
| b = A u*, u* compact | 14 | 4 | 2 / 1 | 4.2 s | 0.01 s |
| b = A u*, u* compact | 24 | 4 | 2 / 1 | 697 s | 0.02 s |

**What it shows:**
- **FORM's time follows the elimination,** whatever the answer: 4.2 s at n = 14, V = 4 in both families. Its
  intermediate fractions swell with the size of the system.
- **PQCF's time follows the answer.** With a random right-hand side, the answer is about as large as the
  intermediates (degree about n, nearly dense), and FORM is 12–18 times faster. When the answer is compact, PQCF is
  25 to 35,000 times faster: from n = 18 to 24 (V = 4), FORM went from 4.2 to 697 s, PQCF stayed at 0.02 s.
- **Real reductions sit between the two, closer to the compact family:** integration-by-parts reductions are known
  for reduced coefficients far simpler than their elimination. That is why reconstruction became the method there
  (FiniteFlow, FireFly, Kira).
- **Scale:** every sample is independent, so the executors spread them over cores or nodes (phase 3), and the work
  grows with the answer, not with the elimination.

**Next:**
- a real IBP system from the literature, or a public reduction benchmark, instead of synthetic ones;
- several components at once (one multi-output reconstruction);
- the federation executor (phase 6, once FVLLMP provides the gate route).

## Phase 5 on a real reduction: the one-loop massless box (7 Nov 2026)

`comparison/form/ibp.erl` generates the integration-by-parts identities of the one-loop massless box.
- **The family:** D₁ = k², D₂ = (k+p₁)², D₃ = (k+p₁+p₂)², D₄ = (k−p₄)², with massless external momenta and the
  invariants s, t and d.
- **The identities:** ∫ ∂/∂k·v / ∏ Dⱼ^aⱼ = 0 for v = k, p₁, p₂, p₄, Laporta-style over seed ranges, with the
  scaleless sectors dropped.
- **The square system:** one elimination modulo a prime, in Laporta order, picks the independent equations and the
  masters. A second elimination pivots the square system so that the target is eliminated last.
- **The target:** I(2,1,1,1) reduces to the box and the t-bubble with the textbook box coefficient (5 − d)/s, and
  the other masters drop out.
- **The masters' coefficients are rational functions of (d, s, t) only.** The answer is linear in the masters, so
  every master is a right-hand side: `pq_recon:linear/4` takes A u = B with several right-hand sides.

What made it fast:
- **One sparse elimination per point** (`solve_sparse` in `ffslp.rs`) tests every right-hand side for 0 at once, at
  three points. Only the non-zero coefficients are reconstructed, each with its own right-hand side alone.
- **The order of the elimination: Laporta's, the most complex column first.** The first sparse version eliminated
  in the opposite order: 26 ms per point at n = 480 from fill-in, 3.5 ms in the right order. FORM eliminates in
  the same order.

FORM 5 solves the same square systems symbolically under `PolyRatFun`, the masters as symbols, and checks
U − Σ rat(Nⱼ mⱼ, Dⱼ) = 0 (it is, in every case):

| seeds (max power, max numerator) | system | masters | FORM 5, 1 core | PQCF, 1 core | PQCF, 20 cores |
| --- | --- | --- | --- | --- | --- |
| 4, 0 | 86 × 86 | 39 | 0.06 s | — | 0.15 s |
| 5, 1 | 260 × 260 | 79 | 0.54 s | 0.28 s | 0.51 s |
| 6, 1 | 304 × 304 | 87 | 1.45 s | — | 0.69 s |
| 6, 2 | 480 × 480 | 103 | 2.4 s (4.4 s in another run) | 0.80 s | 1.29 s |
| 7, 2 | 493 × 493 | 103 | 3.8 s | 0.84 s | 1.37 s |
| 8, 2 | 497 × 497 | 103 | 2.5 s | — | 0.93 s |

- **Core for core, PQCF reduces these systems 1.9–4.5 times faster than FORM from about 260 equations up.** The
  smallest system goes to FORM.
- **On 20 cores PQCF is slower than on one here.** The answer is small, so each Zippel stage needs a few columns,
  and computing 20 ahead mostly wastes them. **Fixed by an adaptive batch:** the first round of a stage computes
  Terminate + 2 columns, each further round twice as many, up to the slots. On the 480 × 480 system the column
  tasks fell from 89 to 33 for the same result. While another job held 12 of the 20 cores, 20 cores still didn't
  beat one on a problem this small (1.22 s against 0.82 s): the parallel gain needs larger answers, or an idle
  machine to show.
- **These systems are small next to real reductions** (10⁴–10⁶ equations). The method's advantage grows with them:
  FORM's symbolic elimination swells, while PQCF's cost per point is a sparse elimination of numbers.
- **What wasn't compared:** FORM's threaded TFORM, and the dedicated reduction programs (Kira, FIRE) that physicists
  use for IBP. Those are the real competition at scale, and the next comparison.

## A two-loop family: the massless propagator "kite" (7 Nov 2026)

`comparison/form/ibp.erl` now takes the family with the seeds: `{kite, MaxPower, MaxNumerator}`.
- **The family:** D₁ = k₁², D₂ = k₂², D₃ = (k₁−p)², D₄ = (k₂−p)², D₅ = (k₁−k₂)², with p² = s. These five form a
  complete basis of the scalar products, so the family has no irreducible numerators. The invariants are d and s.
- **The identities:** ∂/∂kᵢ·v for i = 1, 2 and v = k₁, k₂, p, six per seed.
- **The zero sectors:** decided by the Lee–Pak criterion (`comparison/form/ibp_family.erl`, below). The integral
  is scaleless if it is invariant under some rescaling of the Feynman parameters, which shows on the support of
  U + F. A first, graph-based rule missed one case: a loop momentum left with no propagators, as in I(−1,1,−1,1,−1).
  That integral is 0 but stayed an unknown of the system. The answers were the same either way, and the table below
  is for the corrected systems.
- **The target:** I(1,1,1,1,1) reduces to the product of two one-loop bubbles, I(1,1,1,1,0), and the sunset. The
  sunset appears twice, as I(1,0,0,1,1) and I(0,1,1,0,1): the generator doesn't apply the family's symmetries, so
  two equal integrals stay separate masters.

`pq_recon:linear/4` passes its options' expression limits (`max_expr_nodes`, ...) to the compiler. A is passed
dense, and at 2,102 equations it has 4.4 million entries, far past the default limit of 10⁶ nodes.

FORM 5 checks every answer, as for the box: U − Σ rat(Nⱼ mⱼ, Dⱼ) = 0.

| seeds (max power, max numerator) | system | masters | FORM 5, 1 core | PQCF, 1 core | ratio |
| --- | --- | --- | --- | --- | --- |
| 5, 0 | 198 × 198 | 47 | 0.74 s | 0.07 s | 11× |
| 6, 0 | 365 × 365 | 67 | 2.3 s | 0.14 s | 16× |
| 6, 1 | 712 × 712 | 73 | 9.8 s | 0.50 s | 20× |
| 7, 1 | 916 × 916 | 88 | 16.3 s | 0.68 s | 24× |
| 7, 2 | 1,408 × 1,408 | 94 | 39.8 s | 1.06 s | 38× |

- **Core for core, PQCF is 11–38× faster than FORM on every kite system,** and the gap grows with the size. This is
  with the sparse input below: only the non-zero entries of A and B are compiled, and the right-hand sides are
  sparse columns of the elimination. The box systems above now take 0.22 s (260 equations), 0.38 s (480) and 0.40 s
  (493) on one core, 2.5–9.5× faster than FORM. The first kite runs had A compiled dense and the zero-rule gap:
  3–4.8× (2,102 equations: FORM 60.8 s, PQCF 12.65 s).
- **The masters column counts unreduced boundary integrals too.** These are the integrals at the edge of the seed
  range. Their coefficients in the target's reduction are 0, which one elimination per point finds at once.
- **Sparse input (done):** `pq_recon:linear/4` takes A's rows as lists or as maps Column => entry. The program
  outputs only the non-zero entries, a dense row's literal zeros included. The linear-system black box (`Lin` in
  `ffslp.rs`, `pq_nif:ff_linsys/6`) gets the sparsity pattern, each row's columns. The pattern is a blob in the
  store, so a remote node fetches it by key like the program. A system with 10⁴–10⁵ equations and a handful of
  entries per row now compiles to about 10⁵ outputs, not 10⁸–10¹⁰.
- **Sparse right-hand sides (done):** B's rows may be maps too (Rhs ⇒ entry, with `rhs => M`). The elimination
  (`solve_aug` in `ffslp.rs`) carries the right-hand sides as extra columns it never pivots on, so each costs only
  its non-zero entries. A reduction with thousands of boundary "masters" no longer needs n × M values per point.

## Against FIRE: the two-loop double box (7 Nov 2026)

**Kira couldn't be used.** Kira 3.1 builds from source on this machine (aarch64, with FireFly), but it needs Fermat
even with FireFly, and Fermat's license forbids use "as part of any business or for-profit endeavor" without its
author's written permission. Fermat is also x86-only. FIRE 6.5.2 has a fully open path: its FLINT simplifier
replaces Fermat.

Two changes were needed to build FIRE here:
- skip the PARI submodule, whose git server hangs;
- build with `-fsigned-char`. FIRE's `t_index` is a plain `char`, unsigned on ARM, and without the flag it aborts
  with "Bad inverse!".

**The problem:** FIRE's own `examples/doublebox`, the massless planar double box. It has 7 propagators and 2
numerators in d, s and t, and the target is I(1,1,1,1,1,1,1,−1,−1). `comparison/form/ibp_family.erl` builds the
same family for PQCF from its definition:
- the scalar products are solved from the propagators;
- the zero sectors come from the Lee–Pak criterion on the support of U + F, built mod p;
- the seeds follow Kira's (r, s) form.

The module reproduces `ibp.erl`'s hand-written box and kite identities exactly (372 and 1,026 identities, compared
as sets). It also found the kite's zero-rule gap.

| | FIRE 6.5.2 (FLINT), 1 thread | PQCF, as it stands |
| --- | --- | --- |
| equations | 70,476 generated, 58,349 used | 198,160 (r = 7, s = 2), 108,031 integrals |
| masters | 12 (with the family's symmetries) | 10,223 columns without a pivot, most of them boundary integrals |
| time | **16.3 s**, the whole reduction | **56 s for one elimination mod p** (`pq_nif:ff_pivots/3`, in Rust) |

- **FIRE wins by orders of magnitude here.** PQCF's reconstruction repeats that elimination at every sample
  point, and the answer, rational in d, s, t of degree up to about 5, needs hundreds of points. That is hours on one
  core. (This is an estimate from the measured per-point cost; the full run wasn't made.)
- **The gap is in setting up the system, not in the finite fields.** Generating PQCF's system takes 4 s; the
  elimination is the cost. FIRE doesn't solve one global Laporta system:
  - it works sector by sector, with the family's symmetries (12 masters, not 10,223 unknowns left over);
  - it uses only the equations the target needs.
- **What would close it, in order of effect:**
  1. the symmetries of the sectors (LiteRed-style), which shrink the system and the masters;
  2. "mark and sweep", as in FiniteFlow: one numeric elimination finds the equations the target needs, and every
     sample point solves only those;
  3. a learned elimination, recorded once and replayed per point without the search for pivots;
  4. the points spread over cores and nodes (`pq_rexec`), which helps only once a point costs milliseconds, not a
     minute.
- **Where PQCF stands:** on the systems both programs solve as given, PQCF is 11–38× faster than FORM's symbolic
  elimination. Against a dedicated IBP program on a real two-loop problem, its Laporta front end is the
  bottleneck. FIRE's sector-wise reduction with symmetries is the bar to meet.
