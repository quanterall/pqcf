# R1: reducing large results by finite-field reconstruction, distributed

**Status: accepted, in progress (the project owner's decision, 7 Nov 2026: "go with B and the minimal store").**

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
