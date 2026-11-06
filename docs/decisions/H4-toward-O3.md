# H4: toward FORM's `O3`: primitive Horner nodes, subtraction, and tree search over orders

**Decision: three changes to the optimised programs of `compile/3` (all within the 1.3 additions, none to the default
program):**

1. **Primitive Horner nodes.** Every node of a Horner scheme is computed as its primitive part: its content, the sign
   included, is taken out (`pq_sym_poly:content/2`). So −(a − b), 2(a − b) and a − b share one slot, and the parent
   puts the sign into an addition or subtraction, multiplying by a content only where it isn't ±1.
2. **Subtraction.** The program gets the operation `{sub, Slot, A, B}` (A − B). Horner leaf sums add their positive
   terms and subtract the negative ones by magnitude, so a coefficient −1 no longer costs a multiplication. Only
   Horner programs contain it; `pq_sym_slp:eval/2` and `pq_sym_cse` handle it.
3. **Tree search.** `order_search` takes `strategy => mcts`: a Monte Carlo tree search over order prefixes, with
   `exploration` (the UCT constant, default 1.0; 0.25 did best on the resultants) and `direction => forward |
   backward` (the tree chooses the first or the last atoms of the order). It is parallel and deterministic like the
   local search (`strategy => local`, the default).

This is the project owner's "O3 research" (7 Nov 2026), after H3 found the order worth at most about 10 % with the
old scheme while FORM's `O3` gained 33 % over `O2`.

## What the research found

FORM's sources (`sources/optimize.cc`, read to understand the method, not copied: FORM is GPL) and paper
([Code Optimization in FORM, arXiv:1310.7007](https://arxiv.org/abs/1310.7007)) show:

- **The Horner trees take out the GCD of every coefficient group.** With our identical-polynomial sharing only, a
  sub-polynomial and its negative were computed twice. Taking the content out moved every program below `O2` at the
  occurrence order. It also made the order matter: 64 random orders of res(7,5) now span 15,966–21,342 operations,
  where before they were all worse than the occurrence order.
- **`O3`'s tree search scores orders by Horner plus subtree sharing only**, and runs its greedy step once at the end.
  So the order's value comes from the Horner tree itself, which our local search (swaps from the occurrence order)
  explored poorly. A greedy prefix construction did no better.
- **FORM's greedy step** (x^n, x·y, c·x, x + c, x ± y patterns over all intermediate expressions, the best 5 % per
  round, repeated) gains 20 % over Horner at the same order. Our pair search gains 11–13 %. That, and longer
  searches, are what remain.

## Measured (`comparison/form/resultants.erl`; operations: mul + add + sub + pow; FORM counts a subtraction as an addition)

| | res(7,4) | res(7,5) | res(7,6) |
| --- | --- | --- | --- |
| H2 `horner_cse` (old scheme) | 4,285 | 18,166 | 61,404 |
| H3, + local search, 160 trials | 4,273 | 16,408 | 54,643 |
| **`horner_cse`, primitive nodes + subtraction** | **3,823** | **15,974** | **52,700** |
| **+ tree search, 1,000 trials, exploration 0.25** | **3,140–3,277** | **13,617–13,678** | — |
| FORM 5 `O1` | 4,968 | 20,210 | 71,262 |
| FORM 5 `O2` | 3,969 | 16,398 | 55,685 |
| FORM 5 `O3` | 3,458 | 11,006 | 35,945 |
| FORM 5 `O4` | 2,936 | 10,722 | 35,419 |

- **Below FORM's `O2` at the occurrence order** on all three, and below `O3` on res(7,4) with the tree search
  (forward: 3,140 and 3,277 in two runs of the code as it developed).
- **res(7,5)** with the tree search: 13,617 forward, 13,678 backward, 24 % above `O3`.
- **Time:** a trial on res(7,5) compiles the whole program (about 0.7 s), so 1,000 trials take 200–220 s on 20 cores.
  FORM's `O3` scores by counting, without building the program.
- **Values** agree with the default program within 6.3e-15 (res(7,4)) to 5.5e-12 (res(7,5)): subtractions round
  differently from additions of negated products.

## Tests

- **`pq_horner_tests`:** a content-sharing case by hand (x(a − b) + y(b − a) + 2z(a − b): x − y + 2z once, two
  subtractions, 6 operations). The 60 random polynomials, a resultant, the other atoms and a matrix pass against the
  default program, with the operation counts now including `sub`.
- **`pq_order_search_tests`:** the tree search, forward and backward, on a known optimum (better than a quarter of the
  start's inversions, deterministic); `compile/3` with it never worse than without; and the same program on this
  node and on a peer.
- **`pq_cse_tests`** pass with `sub` in the programs.

## Next

FORM's greedy step at the same order (the 7–9 points that our pair search lacks), and a cheaper score for the
search: counting a Horner tree's operations without building its program.
