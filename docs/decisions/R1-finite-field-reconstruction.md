# R1: reducing large results by finite-field reconstruction, distributed

**Status: accepted, in progress (the project owner's decision, 7 Nov 2026: "go with B and the minimal store").
Phases 1–4 done.**

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
