# H2: common subexpressions across straight-line programs

**Decision: `pq_symbolic:compile/3` takes `optimize => horner_cse` (the Horner scheme of H1, then a pass of common
subexpressions) and `optimize => cse` (the default program, then the same pass). The pass is `pq_sym_cse`. Values
agree with the default program's within rounding; within one option the program is deterministic (D1).** This is
step 2 of the plan agreed with the project owner on 6 Nov 2026 (after K1). It is FORM's `O2` step: `O1` is Horner,
and `O2` adds greedy common-subexpression elimination. It is additive.

## The pass (`pq_sym_cse.erl`)

1. **n-ary nodes.** A sum used once, and only by a sum, merges into it; likewise products. A chain
   ((a + b) + c) + d becomes one node {a, b, c, d}. The merged node leaves the map. Keeping it made a chain of n
   sums into n nodes of 1..n operands: quadratic memory, which on res(7,6) exhausted the heap.
2. **The pair search.**
   - Index every pair of operands of every node: {a, b}, or {a, a} for an operand that occurs twice.
   - Take the pair that the most nodes share, with ties broken by term order.
   - Make it one operation t, and substitute t into every node that has the pair.
   - Repeat until no pair is shared by two nodes. Each step saves (users − 1) operations, so the program never
     grows.
   - The index changes by the pairs a node lost or gained, and the priority queue (a `gb_sets` of
     {−count, pair}) by the same pairs. Re-indexing whole nodes, or a lazy queue, cost 15 s and the heap on res(7,5)
     and res(7,6); the diff costs 0.5 and 3 s.
   - Nodes with more than 64 distinct operands, such as the default program's top-level sum, take no part in the
     search (their pairs are quadratic), but they do gain from pairs found elsewhere.
3. **Back to a program.** Inputs come first in their order. Then every output is emitted depth first, an n-ary node
   as a chain of binary operations over its operands in slot order, with identical operations shared. Dead nodes
   drop out.

Sums and products of complex numbers are commutative and associative, so the result is the same function. The
operation order changes, and values agree within rounding (≤ 3e-13 relative on the resultants).

## Measured (`comparison/form/resultants.erl` expressions, core 8)

| | res(7,4) | res(7,5) | res(7,6) |
| --- | --- | --- | --- |
| default program | 12,375 | 54,288 | 203,775 |
| `cse` | 8,220 | 32,975 | 117,086 |
| Horner (H1) | 4,678 | 19,490 | 67,928 |
| **`horner_cse`** | **4,285** | **18,166** | **61,404** |
| FORM 5 `O1` | 4,968 | 20,210 | 71,262 |
| FORM 5 `O2` | 3,969 | 16,398 | 55,685 |
| FORM 5 `O4` | 2,936 | 10,722 | 35,419 |
| compile, Horner / `horner_cse` | 0.06 / 0.06 s | 0.46 / 0.54 s | 2.6 / 3.0 s |
| compile, `cse` | 0.07 s | 0.72 s | 5.6 s |
| 200 evaluations, Horner / `horner_cse` | 0.24 / 0.20 s | 0.68 / 0.63 s | — |

- `horner_cse` is 7–10 % below Horner, 8–10 % above FORM's `O2`, and 1.5–1.7 times `O4`.
- `cse` alone takes the default program 1.5–1.7 times lower. Use it where Horner doesn't apply.
- A Python estimate with the same greedy rule, made before this work, gave 4,285 and 18,164: the same pass.

**Open:** the rest of the way to `O2` is FORM's wider set of greedy steps. The way to `O3` and `O4` is the search
over Horner orders (step 3 of the plan: parallel and deterministic, with Erlang processes and nodes).

## Tests (`test/pq_cse_tests.erl`)

- **The pair search on a small case:** a + b + d and a + c + d share a + d, and the products likewise. 8 operations
  become 6.
- **40 random sets of polynomials** with shared terms: for every option, values within 1e-12 of the default's;
  `cse` never adds to the default, and `horner_cse` never to Horner.
- **A resultant:** fewer operations than Horner, and the same program twice.
- **Radicals, exponentials, conjugates, fractional powers.**
- **Trivial outputs:** a symbol, a constant, zero, and the same expression twice.
- **A matrix:** the symbolic unitary of a mesh, every entry.
- **Mutations:** substituting for only one operand of the pair, and merging operands that are used more than once,
  are each caught (3 and 5 failing tests).
