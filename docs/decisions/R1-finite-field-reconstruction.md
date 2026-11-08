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
- **Where PQCF stood:** on the systems both programs solve as given, PQCF was 11–38× faster than FORM's symbolic
  elimination. Against a dedicated IBP program on a real two-loop problem, its Laporta front end was the
  bottleneck. The next section closes the gap.

## Closing the gap: mark and sweep, symmetries, a shared cache (7 Nov 2026)

The steps listed above, in the order they paid off:

1. **Mark and sweep** (`pq_nif:ff_pivots/4`). The learning elimination records, for every row, the pivot rows it
   was combined with. The target needs its pivot row, every row that row was combined with, and the pivot row of
   every unknown left in a needed row, recursively. The other rows never touch the target. The needed rows, in the
   pivot order of their columns, form the square system. Its unknowns stay in Laporta order, so eliminating column
   k on row k repeats the learned elimination, fill-in included (`ibp:system/3`).
   - On the double box at Kira's r = 7, s = 2, 54,549 of 198,160 equations are needed.
   - A sample point costs **0.19 s instead of 56 s**.
2. **Sector symmetries** (`ibp_family:canon/2`, Pak's criterion). Two sectors are equivalent if U + F of one is U +
   F of the other with the Feynman parameters permuted, compared mod p at a random point of the invariants.
   - An integral without numerators then equals its images in the equivalent sectors and under the sector's own
     automorphisms, and is named by the smallest image.
   - The double box's 43 non-zero sectors fall into 18 classes; its top sector has 4 automorphisms.
   - The kite's two sunsets become one master, and the double box's target reduces to **8 masters**, the known
     number.
3. **A cache shared by the right-hand sides** (`LinCache` in `ffslp.rs`). One elimination yields every right-hand
   side, but each non-zero coefficient was reconstructed separately, so each one eliminated at its own sample points.
   - The points of a line along t are now a prefix-stable stream (`distinct/3`), so reconstructions of different
     degrees sample the same points.
   - The black boxes of one system share a cache of solved points: each point is eliminated once.
   - This changes no result. The double box's reduction fell from 135 s to 38 s.
   - Caveat for benchmarks: a second run of the identical system in the same VM is served from the cache. The timings
     below are fresh VMs.
4. **Seeds without dots** (`{dots, D}`: every sector's powers at most D above 1). Kira's r = 7 lets a 4-propagator
   sector carry 3 dots, and a reduction of a dot-free target rarely needs them. The double box at D = 0, s = 2 has
   7,300 equations instead of 198,160, and still reduces completely.

A pre-filter that dropped equations whose integral appears nowhere else removed only 3% of them, so it was left out.

**Checked against FIRE** (`comparison/fire/fire_check.erl`). The two programs pick different masters. FIRE's
example applies fewer symmetries, so it keeps 12 masters, which are 8 up to the symmetries. It also takes a dotted
top integral where PQCF takes one with a numerator. FIRE reduces PQCF's numerator master into its own basis. At a
random point mod p, PQCF's answer then equals FIRE's on every master, both for D = 0 and for D = 1.

**One core each, back to back on the same core, the machine otherwise loaded:**

| | equations generated / needed | masters | time |
| --- | --- | --- | --- |
| FIRE 6.5.2 (FLINT), `examples/doublebox` | 70,476 / 58,349 | 12 (8 up to symmetries) | 7.7–8.4 s |
| PQCF, seeds D = 0, s = 2 | 7,300 / 4,044 | 8 | **3.50 s** (0.3 s system, 3.15 s reconstruction) |
| PQCF, seeds D = 1, s = 2 | 39,790 / 15,379 | 8 | 24.3 s |
| PQCF, Kira's r = 7, s = 2 | 198,160 / 52,997 | 8 | about 100 s (learning 60 s) |

(FIRE's first runs took 16 s on a busier core.)

- **On this double box PQCF is now about 2.3× faster than FIRE on one core**, with the same answer.
- **The seeds matter as much as anything.** FIRE chooses its equations itself. PQCF needs a seed choice that fits the
  target; with Kira's standard seeds it is about 12× slower than FIRE. Choosing seeds automatically (growing them
  until the target reduces) is the next step for the front end.
- **The systems the box and kite benchmarks solve shrank with mark and sweep and the symmetries.** The kite at 7, 2
  keeps 261 of 2,310 equations. FORM, solving the same small systems, is now as fast as PQCF on the smallest:

| family, seeds | needed | FORM 5, 1 core | PQCF, 1 core (reconstruction) |
| --- | --- | --- | --- |
| box 5, 1 | 39 | 0.03 s | 0.05 s |
| box 6, 2 | 121 | 0.12 s | 0.11 s |
| kite 6, 1 | 112 | 0.09 s | 0.05 s |
| kite 7, 1 | 127 | 0.11 s | 0.06 s |
| kite 7, 2 | 261 | 0.43 s | 0.11 s |

  The large ratios over FORM in the sections above came from solving the full systems. Once only the needed
  equations are solved, symbolic elimination is fine at this size; the reconstruction wins on the large systems.

## Choosing the seeds automatically (7 Nov 2026)

`ibp:main/2` takes `{auto, Family}` or `{auto, Family, Target}`. `ibp:auto_seeds/2` tries seeds `{dots, D}` with
numerators up to S. It starts from the target's own dots and numerator powers and tries the cheaper candidates
first (D + S, then D).

A candidate is accepted when the reduction is **complete**. After one numeric solve at a random point mod p, every
master with a non-zero coefficient must itself be a seed, its own identities generated. A master at the edge of the
seed range, whose identities were never written, means the seeds are too small. No learning is involved: on these
families the first or second candidate is the answer, and the test costs one numeric solve (0.4–1.1 s on the
double box).

| family, target | seeds chosen | equations / needed | PQCF, 1 core | FIRE, same core |
| --- | --- | --- | --- | --- |
| box I(2,1,1,1) | D = 1, S = 0 | 108 / 5 | 0.08 s | — |
| kite I(1,1,1,1,1) | D = 0, S = 1 (after S = 0) | 102 / 15 | 0.05 s | — |
| double box I(1,1,1,1,1,1,1,−1,−1) | D = 0, S = 2 | 7,300 / 4,044 | 3.95 s | 7.7–8.4 s |
| double box I(1,1,1,1,1,1,1,−2,−1) | D = 0, S = 3 | 18,370 / 9,895 | 32.2 s | 11.9–12.2 s |

Both double-box results are checked equal to FIRE's at a random point mod p (`fire_check.erl`), all 8 masters.

- **The criterion picks the smallest working seeds**, the ones that were chosen by hand before, and the answers stay
  exact and verified.
- **On the harder target FIRE is 2.7× faster.** PQCF's time is now the reconstruction (30 of 32 s): the coefficients
  are larger (up to 15 / 10 terms), so they need more sample points. Two ways to cut it: fewer points, by guessing
  the denominators' factors (the factors of IBP coefficients come from a small set of polynomials in d, s, t), and a
  cheaper point, by replaying the learned elimination instead of searching for pivots again.
- **A learned model for the seeds isn't needed yet.** The simple growth finds the smallest seeds at the cost of one
  or two numeric solves; a model would have nothing to gain until families with many more candidates.

## Simplifying the results: verified rewriting (7 Nov 2026)

The question was whether neural networks could simplify symbolic results. A network's output is never exact, but
PQCF can check any rewrite exactly and cheaply: evaluate the original and the rewritten form at random points. So
the design is *propose, then verify*. The proposer can be a fixed rule, a heuristic or, later, a trained model; a
proposal can make a result smaller, never wrong.

`pq_simplify` is that environment, with no learning in it yet:
- **Moves.**
  - Exact division by candidate factors over the rationals. The default alphabet is the variables, a·x + b
    (a = 1..3, |b| ≤ 12) and x ± y, x + 2y, 2x + y; `extra` adds proposals.
  - Collecting by monomials in some variables (default: all but d), each coefficient factored and cancelled
    against the denominator's factors in the other variables.
- **Size: the additions written** (`adds/1`; a sum of k terms has k − 1). A product of linear factors counts one per
  factor, its expansion one per term. `rational/4` returns the smallest of the expanded, factored and collected
  forms, the more structured among equals: never larger than the input.
- **Verification.** Division is exact by construction, and `verify/4` checks the result against the input at random
  points anyway.
- **Printing:** `format/1` keeps the structure. `pq_symbolic:pretty/1` canonicalises, which would expand the
  products again.

On the double box (`ibp:main/2` writes `<name>.simplified`, every coefficient verified):

| target | additions, expanded | simplified |
| --- | --- | --- |
| kite I(1,1,1,1,1) | 6 | 5 |
| double box I(1,1,1,1,1,1,1,−1,−1) | 25 | 18 |
| double box I(1,1,1,1,1,1,1,−2,−1) | 79 | 58 |

For example, the coefficient of I(1,1,0,0,0,1,1,0,0) goes from (−255 + 323/2 d − 51/2 d²) / (16s − 8ds + d²s) to
−(17/2)(3d − 10)(d − 3) / (s (d − 4)²).

**Is there anything left for a network to find?** For factors, no. FORM's complete multivariate factorisation (FLINT)
of the remaining unfactored numerators finds them irreducible over the rationals: only content (2, −1/6, 1/4) comes
out. The alphabet already found every factor of these results, so a model proposing factors would have nothing to
add here. Where a network could still help:
- **The choice of masters.** The basis decides how large the coefficients are. Choices that separate d from the
  kinematics (d-factorising bases) are known to shrink coefficients a lot, and the space of choices grows quickly
  with the number of sectors: a search a model can steer, every candidate verified exactly by reducing again.
- **Larger results** (more loops, more scales), where factors outside a fixed alphabet appear.

## The masters' differential equations, and a physics-informed neural network on them (7 Nov 2026)

**Differential equations from PQCF** (`ibp_family:derivative/3`, `comparison/pinn/box_de.erl`).
- **The derivative.** The derivative of an integral with respect to an invariant is written with external-momentum
  operators p_i·∂/∂p_j acting on the propagators, like an IBP identity. A combination Σ c_ij p_i·∂/∂p_j that moves
  every external scalar product exactly as the invariant does (the others fixed, the legs on shell) is that
  derivative. The c_ij come from a small exact linear solve (`pq_recon:linear/4`).
- **The equation.** Each master's derivative is reduced to the masters (automatic seeds), and every entry of the
  matrix is reconstructed as one reduced rational function and simplified (`pq_simplify`).
- **The one-loop massless box gives the textbook equation**, dM/dt = A M for M = (box, s-bubble, t-bubble):
  - dBox/dt = [(d − 6)s − 2t] / (2t(s + t)) · Box − 2(d − 3) / (s t (s + t)) · Bub_s + 2(d − 3) / (t² (s + t)) · Bub_t;
  - dBub_t/dt = (d − 4)/(2t) · Bub_t, and dBub_s/dt = 0.
- **Checked independently.** The box, computed directly over its Feynman parameters (tensor tanh-sinh quadrature,
  about 1e-8), and carried by an ODE solver along PQCF's equation from t = −1, agree to 1e-8 at t = −0.1, −0.3, −3
  and −10 (d = 5.3, s = −1). The bubbles agree with their closed form to 1e-13.

**The PINN** (`comparison/pinn/pinn_box.py`). It learns M(t) from the equation and the value at t = −1 alone, with no
solution data.
- **Setup.** The network and training follow TechUPGRADE's PINN-PONN reproduction (TU-TR-2026-002; reused with its
  owners' permission): a tanh layer, three residual blocks of two tanh layers, AdamW inside a look-ahead (k = 5,
  α = 0.5), then L-BFGS. The boundary value is built in exactly. It ran in TechUPGRADE's JAX + CUDA image
  (phnn-trainer, read-only) on the GB10, with memory capped at 5% of the GPU because the machine was short.
- **Results:**

| task | PINN | classical |
| --- | --- | --- |
| M(t) on t ∈ [−10, −0.1], d = 5.3 | 52 s of training; max relative error 0.9e-4 (box), 1.0e-4 (s-bubble), 3.7e-4 (t-bubble) | ODE (DOP853): 19 ms at a tolerance of 1e-13 |
| one network over (t, d), d ∈ [4.6, 5.8], boundary values given at 5 values of d | 143 s of training; at unseen d (4.75, 5.05, 5.3, 5.65) max relative error 0.3–5.7% | an ODE solve per d (about 20 ms) after that d's boundary value (1–3 s of integration) |

(The first single-d run, with a dense Jacobian for the residual, took 398 s for the same 1e-4. Forward-mode with a
tangent of ones gives the same derivative in 52 s.)

**Negative:**
- **Accuracy and cost.** For one evaluation the PINN is about 2,700× slower and 10⁹× less accurate than an ODE
  solver on the same equation.
- **Exactness.** A PINN can't give the exact or verified results PQCF is built on. Its errors are unknown without a
  reference solution, and computing that reference is the classical work it was meant to replace.
- **The surrogate over d** reaches only percent-level accuracy away from the dimensions whose boundary values it
  saw. It still needs those boundary values, the expensive part: 1–3 s of integration each, against 20 ms for the
  equation.
- **It doesn't simplify anything.** A PINN evaluates integrals numerically; it doesn't make symbolic results smaller.

**Positive:**
- **No solution data.** The PINN learned the solution from the equation and one boundary value alone, to 1e-4,
  with no tuning beyond the borrowed recipe. The method works as advertised at this scale.
- **A differentiable model.** The result is smooth in t (and in d for the surrogate), differentiable, and cheap to
  evaluate in bulk. That could be useful where a downstream fit or optimisation needs derivatives of the integrals.
- **The setup was cheap.** Reusing TechUPGRADE's stack meant no new infrastructure, and it fit a tight memory
  budget.
- **The real gain of this step is the differential equation itself.** PQCF now derives master-integral equations
  exactly from IBP, verified against direct integration. That's the input every modern method needs: canonical
  bases, series solutions, numerical evaluation with classical solvers.

**Conclusion:** for evaluating master integrals, use PQCF's equation with a classical ODE or series solver. A PINN
is worth reconsidering only as a differentiable surrogate over many parameters, and only once its error is
acceptable, which here it wasn't. For simplifying results it is the wrong tool. The steerable search (the choice of
masters) stays the place for learning.

**The choice of masters, first look.** Ordering dots before numerators (`order => dots`) keeps the same 8 masters on
the double box. The dot-free seeds the automatic choice picks contain no dotted integrals to swap in. It still
shrinks the needed system (3,717 equations instead of 4,044). A different basis needs seeds with dots, chosen
together with the order: the next step for the basis search.

## Guessing the denominators: level with FIRE on the harder target (8 Nov 2026)

On I(1,1,1,1,1,1,1,−2,−1) the reconstruction took 30 of PQCF's 32 s: about 1,000 sample points at 0.031 s each,
almost all of it elimination (the program costs 0.002 s). Two measurements decided the fix:
- **Dropping the zero right-hand sides doesn't help.** Keeping only the 8 non-zero of 1,424 right-hand sides saves
  5% (0.0297 s against 0.0313 s): the cost is the elimination itself.
- **So, fewer points.** A rational function in three variables is reconstructed along lines, numerator and
  denominator together. A polynomial the size of these numerators needs far fewer points.

**Denominator guessing** (`pq_recon`, `guess => true` by default), as the reduction programs do it:
1. For every variable, a random line along it, 48 points. Every non-zero right-hand side's values come from the
   same eliminations, through the linear system's cache.
2. On each line, each right-hand side's univariate denominator by Thiele's continued fraction. A candidate factor's
   power is how often its root on that line divides the denominator; the candidates are `pq_simplify:alphabet/1`:
   the variables, small linear forms, x ± y.
3. The guess stands only if the factors found account for the whole univariate denominator on every line.
   Otherwise the column falls back to the rational reconstruction.
4. **F·D as a polynomial.** The linear-system black box multiplies its value by the guess (a `Scale` of linear
   factors in `ffslp.rs`; `pq_nif:ff_linsys/8`; the scale travels in the store's pattern blob, so a remote node
   applies it too). Zippel reconstructs F·D as a polynomial, from cached, unscaled solutions.
5. **N / D reduced and checked.** Each factor is divided out of the numerator as often as it goes (exact division,
   `pq_simplify:divide/2`), and D is normalised as before. A wrong guess fails the polynomial's final check at a
   fresh point in a fresh prime, and the column falls back to the rational reconstruction: the result stays exact.

**The results are bit-identical with and without guessing** (the same coefficients, the same hash). Reconstruction
alone, one core: 32.0 → **9.0 s** on I(1,1,1,1,1,1,1,−2,−1), and 3.3 → 2.4 s on I(1,1,1,1,1,1,1,−1,−1).

**Whole runs** (automatic seeds, system, reconstruction, simplification), one core, back to back with FIRE 6.5.2
(FLINT) on the same core:

| double box target | PQCF | FIRE |
| --- | --- | --- |
| I(1,1,1,1,1,1,1,−1,−1) | **3.16 s** | 8.10 s |
| I(1,1,1,1,1,1,1,−2,−1) | **11.87 s** | 12.21 s |

- **PQCF is now level with FIRE on the harder target and 2.6× faster on the other,** core for core, with answers
  checked equal to FIRE's.
- **What remains is the elimination per point** (0.031 s on the harder target). A learned elimination, recorded
  once and replayed per point without the search for pivots and the set bookkeeping, is the next lever. After
  that comes spreading the points over cores and nodes, which neither program has been given here.
