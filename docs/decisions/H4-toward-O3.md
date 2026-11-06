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

## Signed sums in the subexpression pass, and FORM's own order (7 Nov 2026)

`pq_sym_cse` now treats additions and subtractions as one signed n-ary sum. A pair records its relative sign, so
a + b and −a − b match, as do a − b and b − a, and an occurrence keeps its sign. `horner_cse` at the occurrence order
went to 3,771 / 15,875 / 52,423 (about −1 %).

FORM's `stats=on` breakdown shows where the totals stand. On res(7,4), FORM's `O2` uses 21 powers, 2,130
multiplications and 1,793 additions; our `horner_cse` uses 28 powers, 2,118 multiplications and 1,625 additions or
subtractions. At the occurrence order our scheme and pass now match or beat FORM's Horner and greedy step.

**The decisive test.** With `printscheme=on`, FORM 5 reports the order its `O3` chose for res(7,5):
`b5, a7, b3, a6, b4, a5, a4, b2, a3, a2, b1, a1, b0, a0`. With that order, our `horner` gives 12,374 operations and
`horner_cse` **10,565**, below FORM's own `O3` result (11,006). So what remains is the search, not the scheme. Our
tree search scores a candidate by building its whole program in Erlang (about 0.7 s on res(7,5)), so 1,000 trials
take 200–330 s and reach 13,334–13,678. FORM's 1,000 iterations take 16 s, because they only count.

## The score in Rust (7 Nov 2026)

`pq_nif:horner_count/3` (`native/pq_nif/src/hcount.rs`) counts the operations of the Horner program for an order
without building it. It follows `pq_sym_slp`'s Horner step for step: primitive nodes and their content's sign,
add / sub folding, leaf sums, and every atom power, monomial, product and sum keyed by structure. It applies to
plain polynomials of symbols, and a test requires it to equal the built program's count on 180 random polynomials
and orders. When it applies, the search scores with it (Horner without the pass, as FORM scores without its greedy
step). The final program is then built with the best order, and is kept only if it beats the occurrence order.

| `horner_cse` + tree search (forward, exploration 0.25, batch 16; 20 cores) | res(7,4) | res(7,5) | res(7,6) |
| --- | --- | --- | --- |
| 1,000 trials, scored by building (before) | 3,140–3,277 (29 s) | 13,617 (196 s) | — |
| 1,000 trials, scored in Rust | — | 11,815 (27 s) | — |
| **5,000 trials, scored in Rust** | **2,741** (22 s) | **10,243** (97 s) | — |
| 2,000 trials, scored in Rust | — | — | 40,842 (176 s) |
| FORM 5 `O3` | 3,458 | 11,006 | 35,945 |
| FORM 5 `O4` | 2,936 | 10,722 | 35,419 |
| FORM `O4` time (wall) | 2.6 s | 9.8 s | 57 s |

- **Below FORM's `O4` on res(7,4) and res(7,5)**, by 7 % and 4 %.
- **res(7,6) not yet:** 40,842 after 2,000 trials, 15 % above `O4`. It has the most variables and needs a larger
  budget.
- **Backward** did worse on res(7,5) at 5,000 trials (13,065), so forward stays the default.
- **FORM is still faster** for a given quality: its `O4` takes 2.6–57 s, ours 22–176 s on 20 cores.

## Tests

- **`pq_horner_tests`:** a content-sharing case by hand (x(a − b) + y(b − a) + 2z(a − b): x − y + 2z once, two
  subtractions, 6 operations). The 60 random polynomials, a resultant, the other atoms and a matrix pass against the
  default program, with the operation counts now including `sub`.
- **`pq_order_search_tests`:** the Rust score equal to the built program's operation count (60 random sets of
  polynomials, 3 orders each, several outputs; rationals, Gaussian coefficients, negative terms), and no fast score
  for polynomials with π. The tree search, forward and backward, on a known optimum (better than a quarter of the
  start's inversions, deterministic); `compile/3` with it never worse than without; and the same program on this
  node and on a peer.
- **`pq_cse_tests`** pass with `sub` in the programs.

## Next

res(7,6) at FORM's level: a larger budget, or a faster count (the Rust score still copies the polynomial's terms
into every trial).
