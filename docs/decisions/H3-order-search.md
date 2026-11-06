# H3: searching Horner variable orders, in parallel and deterministically

**Decision: `pq_symbolic:compile/3` with `optimize => horner | horner_cse` takes
`order_search => #{trials => N, seed => S, batch => B, nodes => Ns}`. It searches N variable orders, one order for
all the polynomials of the program, and keeps the order with the fewest operations. The candidates of a batch are
scored in parallel, one process each, spread over the given nodes. The program is the same whatever the number of
cores or nodes (D1).** This is step 3 of the plan agreed with the project owner on 6 Nov 2026 (after K1 and H2).
It is the step FORM takes in `O3` and `O4` (Monte Carlo tree search, local search), using Erlang processes and
distribution instead of one thread. It is additive; without `order_search`, or with `trials => 0`, nothing changes.

## The search (`pq_sym_search.erl`)

- **Start:** the occurrence order over all the outputs (most terms first, ties broken by canonical order).
- **Rounds:** each round makes B candidates from the current order (16 by default). Candidate I of round R gets one
  to three random moves (swap two atoms, or move one), drawn from a generator seeded with {S, R, I} alone.
- **Scoring:** a candidate's score is the operation count of the program compiled with that order, by Horner (and
  the pass of H2 for `horner_cse`). Every polynomial uses the order, followed by its own occurrence order.
- **Acceptance:** the candidate with the lowest {score, index} replaces the current order if it is no worse (equal
  scores move along plateaus). The best order seen is kept. It is used only if it beats each polynomial's own
  order, so the search never adds operations.
- **Parallel and deterministic:** a batch is scored by one linked process per candidate, the candidates dealt over
  `nodes` in turn (default: this node only). The results are collected in candidate order. Candidates depend only
  on (S, R, I) and scores are deterministic, so neither the cores nor the nodes nor the completion order change
  anything. The batch size does change the search, so it is an option, not derived from the cores.
- **Timeouts:** the workers are linked to the compiling process, so `sym_timeout_ms` stops them too.

**Annealing, tried and dropped:** accepting worse candidates with a falling temperature (0.2 % and 1 % of the
cost) gave the same result on res(7,4) as plain acceptance. There, 30 random orders score 4,538–5,811 against the
occurrence order's 4,285: the order hardly matters.

## Measured (`horner_cse`, batch 16, seed 1; the resultants of `SY-FORM-comparison.md`)

| | res(7,4) | res(7,5) | res(7,6) |
| --- | --- | --- | --- |
| `horner_cse` | 4,285 | 18,166 | 61,404 |
| **+ search, 160 trials** | **4,273** | **16,408** | **54,643** |
| + search, 800 trials | 4,273 | 16,327 | — |
| FORM 5 `O2` | 3,969 | 16,398 | 55,685 |
| FORM 5 `O3` | 3,458 | 11,006 | 35,945 |
| FORM 5 `O4` | 2,936 | 10,722 | 35,419 |
| time, 160 trials, 16 cores | 2.7 s | 29 s | 132 s (20 cores) |
| FORM `O4` (wall) | 2.6 s | 9.8 s | 57 s |

- **Below FORM's `O2`** on res(7,5) and res(7,6), 8 % above it on res(7,4). Still 1.5 times `O3` and `O4`.
- **Parallel:** res(7,5) with 160 trials took 184 s on one core and 29 s on 16, 6.3 times faster, with the same
  program.
- **What's left is in the scheme, not the order.** On res(7,4), no order gets near `O3`. FORM's remaining lead is in
  its Horner and greedy steps (its `O2` beats our `horner_cse` at the same, occurrence, order), and in its tree
  search, which builds orders prefix by prefix instead of moving around complete ones.

## Tests (`test/pq_order_search_tests.erl`)

- **The search on a known optimum** (inversions against the sorted order): it finds orders with fewer than a
  quarter of the start's inversions, the same order twice, and the start back for no trials or a single atom.
- **`compile/3`:** never more operations than without the search, the same values, the program without the search
  for `trials => 0`, and the option ignored without Horner.
- **Nodes:** the same program scored on this node, on a peer node, and on both.
- **Mutations:** never accepting a candidate is caught. Candidates that ignore the round number are not, rightly:
  they still differ once the current order moves.
