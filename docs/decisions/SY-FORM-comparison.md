# pq_symbolic and FORM 5: a comparison (2 November 2026)

## What `pq_symbolic` is for

The symbolic specification (SY 1.1) gives the engine six callers:

- circuit parameters: default expressions and lowering;
- the transfer matrix of a passive circuit, with a decisive unitarity check;
- symbolic Gaussian means and covariances (FR-G04);
- exact gradients of probabilities and mean photon numbers;
- the JSON parser;
- one expression evaluated under many bindings, through straight-line programs.

Its scope fence (SY 1.2) excludes:

- factorisation and gcd;
- solving, integration and series;
- symbolic determinants and inverses;
- symbolic Fock evolution.

Its expressions are those of a photonic circuit: thousands to a few hundred thousand terms, built inside the
Erlang VM, bounded (term counts, a timeout, a heap limit) because they may come from untrusted input. It has exact
radicals and π-exponentials, symbol domains and branch cuts, and compiles to fast numeric programs.

FORM (Vermaseren et al.; version 5, arXiv:2601.19982, January 2026) is a different kind of tool. It is a
stand-alone system for very large expressions in high-energy physics: terms streamed through sorts that spill to
disk, a parallel version (TFORM), and in version 5 a Feynman-diagram generator, arbitrary-precision floats and
FLINT for polynomial gcd and factorisation.

## The benchmark

`bench/form/` contains:

- `Dockerfile`: FORM from its repository (git fb12f2e, 13 Oct 2026) with FLINT 3.2.2;
- `gen_form.py`: the FORM scripts;
- `symbench.erl`: the same expressions through `pq_symbolic`, with its limits raised.

Both run single-threaded on the provisional reference (Cortex-X925).

| Task | Terms (both) | FORM 5 | pq_symbolic | ratio |
| --- | --- | --- | --- | --- |
| T1 (1 + x1 + … + x6)^8 | 3,003 | < 0.01 s | 0.09 s | |
| T1 the same to the power 12 | 18,564 | 0.01 s | 3.4 s | ~300× |
| T1 the same to the power 16 | 74,613 | 0.03 s | 61 s | ~2,000× |
| T1 the same to the power 20 | 230,230 | 0.12 s | refused (an intermediate of 64 M terms) | |
| T2 ((1 + x + y + z + w)^10)² | 10,626 | 0.14 s | 3.4 s | ~24× |
| T3 permanent of a 6 × 6 matrix of symbols | 720 | < 0.01 s | 12.6 s | |
| T3 the same, 7 × 7 | 5,040 | 0.03 s | 429 s | ~14,000× |
| T3 the same, 8 × 8 | 40,320 | 0.16 s | refused (heap limit) | |
| T4 Clements mesh unitary, 4 modes, all angles symbolic | 400 | < 0.01 s | 0.05 s | ~1× |
| T4 the same, 6 modes (30 symbols) | 11,616 | 0.03 s | 0.05 s | ~2× |
| T4 the same, 8 modes (56 symbols) | 276,800 | 0.40 s | 2.2 s | ~5× |

## What it shows

- **The answers agree:** the term counts are identical wherever both finish.
- **On the engine's own workload, the symbolic transfer matrix (T4), it is within 1–5× of FORM.** That covers 276,800
  terms for 8 modes, the `max_symbolic_modes` default.
- **On bulk polynomial algebra FORM is 25 to 2,000 times faster.** Part of that is FORM's C engine and its sorting.
  Part is our algorithms:
  - powers are built by repeated multiplication rather than the multinomial expansion;
  - the symbolic permanent uses Glynn's formula, which expands 2^(n−1) products of n linear forms and relies on
    massive cancellation. FORM, as written here, multiplies the rows one by one and drops every term that uses a
    column twice, so its intermediates stay small.

  Both are fixable inside the scope fence.
- **FORM 5's FLINT interface did not matter here:** FLINT speeds up gcd and factorisation (rational functions),
  which `pq_symbolic` does not do.
- **Beyond about 10⁶ terms, or for gcds and rational functions, FORM is the tool.** `pq_symbolic` refuses such
  sizes by design.

## Options (both taken, 2 November 2026)

1. Keep `pq_symbolic` for its role and fix the two algorithms: the multinomial power, and the pruned
   row-by-row permanent. Both are internal changes with no change to the API.
2. For large-scale work (symbolic permanents of large n, symbolic Fock evolution, rational functions), add an
   optional FORM backend. It would export the expression as a FORM script, run FORM as a separate process, and
   read the result back. FORM is GPL-3 and PQCF Apache-2.0, so FORM would be called as an external program, never
   linked or bundled.

## Taken

The project owner chose both options, and a permanent comparison:

**1. The two algorithms.**

- **Products of plain monomials** (symbols, `conj` and π, without radicals, opaque bases or exponentials) only
  merge their exponents. Before, about 31 times as much work went into the general normaliser.
- **Powers of sums of plain monomials** are enumerated multinomially, about one step per result term.
- **The permanent** is built row by row over the subsets of columns already used. Terms that would use a column
  twice are never made, which is the same pruning as the FORM script.
- **Tests:** each new method is exactly equal to the old one, in canonical form, on random and mixed inputs. Glynn's
  formula is kept as the test oracle.

Before and after, single-threaded, against FORM 5:

| Task | FORM 5 | before | after |
| --- | --- | --- | --- |
| (1 + x1 + … + x6)^16 | 0.03 s | 61 s | 0.64 s |
| (1 + x1 + … + x6)^20 | 0.12 s | refused | 2.4 s |
| dense product | 0.14 s | 3.4 s | 1.3 s |
| permanent 7 × 7 | 0.02–0.03 s | 429 s | 0.02 s |
| permanent 8 × 8 | 0.16 s | refused | 0.26 s |

**2. The bridge.** `pq_form` covers expansion and the symbolic permanent through FORM as an external process, on
the polynomial fragment. Its results equal `pq_symbolic`'s exactly. It is an addition under the freeze.

**The standing comparison.** `comparison/form/` holds the reference tasks, the latest results and their history.
It runs weekly on GitHub Actions until the end of the project.


## The literature's benchmarks (5 November 2026)

The project owner asked for a comparison on the problems that recent papers use to test FORM. Three sources:

- **FORM 5's own paper** (arXiv:2601.19982, section 3.3, table 1):
  - Its tests are `minceex` (three-loop moments), `forcer` (four-loop integration by parts) and `mbox1l` (a
    one-loop massive box), 24 TFORM workers on a Ryzen 9 7900X, with and without FLINT.
  - They run through FORM's own Mincer and Forcer packages, thousands of lines of FORM code. Their arithmetic is in
    rational functions: polynomial GCDs, which are outside `pq_symbolic`'s scope fence (SY 1.2).
  - They cannot be run here, and a reimplementation of the packages would compare our reimplementation, not the
    engine.
- **polybench** (T. Ueda, one of FORM's authors; github.com/tueda/polybench): multivariate GCD and factorisation in
  FORM, FLINT, Mathematica, Singular, Symbolica and others. Again outside the scope fence: PQCF has no
  polynomial GCD.
- **The classic multiplication benchmarks** used for FORM and other systems in the computer-algebra literature are
  within scope, and are now standing tasks of the comparison:
  - **T5,** Fateman's dense f (f + 1), with f = (1 + x + y + z + t)^n;
  - **T6,** Monagan and Pearce's sparse f g, with f = (1 + x + y + 2z² + 3t³ + 5u⁵)^n and g its mirror.

  The sizes are reduced from the classic n = 20 and n = 12, so that our engine's run stays within minutes.

**Results** (run of 5 Nov 2026, core 8 of the provisional reference, single-threaded, every check `ok`):

| Task | Terms | FORM 5 | native | ratio | bridge |
| --- | --- | --- | --- | --- | --- |
| T5 Fateman n = 10 | 10,626 | 0.14 s | 1.34 s | 9.7 | 0.24 s |
| T5 Fateman n = 15 | 46,376 | 2.17 s | 33.1 s | 15 | 3.1 s |
| T6 sparse n = 6 | 114,000 | 0.051 s | 1.28 s | 25 | 2.6 s |
| T6 sparse n = 8 | 591,235 | 0.38 s | 12.0 s | 32 | 20 s |

The bulk products sit where T2 already showed: 10 to 32 times FORM's time natively. That is the cost of a sparse
product in a map, without FORM's sorted streams. The bridge brings the dense ones within 1.1 to 1.7 times FORM
alone.

**Found and fixed on the way:**

- **The bridge timed out on T6** (600 s, where FORM alone takes 0.05 s): `pq_form:expand` wrote a product of
  powers as one expression, `Local F = (...)^n*(...)^n`, and FORM expands such a product inside every term before
  any sort. It now writes each factor as its own expression, sorts, then multiplies, as FORM practice does.
  - T6 n = 6: timeout → 2.6 s;
  - T2: 0.60 → 0.31 s;
  - T5 n = 15: 10.6 → 3.1 s.
- **Reading FORM's output back** spent most of its time on monomial powers such as u¹⁷, computed by
  square-and-multiply. `pq_sym_poly:pow_int` now raises a plain monomial directly (exponents times n, coefficient
  to the n), which halved the read-back at n = 4.

Results are unchanged: the bridge's 15 tests and the symbolic test suite pass.

## A real task from FORM's literature: the resultant benchmarks (5 November 2026)

The project owner asked for a real task from a paper that uses FORM, with a published description, done with PQCF
and compared with FORM. Most recent work that uses FORM is outside `pq_symbolic`'s scope: Dirac traces, integral
reductions, rational functions. The expression-simplification papers of Ruijl, Plaat, van den Herik and Vermaseren
(HEPGAME, arXiv:1405.6369; arXiv:1409.5223, table II) are inside it.

Their benchmarks include res(m, n) = res_x(Σ_(i≤m) a_i x^i, Σ_(i≤n) b_i x^i), the resultant of two polynomials
with symbolic coefficients. The task is to expand it, then find a cheap way to evaluate it. The papers' other
expressions (HEP processes) are not published; the resultants are fully specified.

**Method** (`comparison/form/resultants.erl`):

- **PQCF:** the Sylvester determinant is expanded row by row over the subsets of columns used, with
  `pq_symbolic`'s exact arithmetic at every step. The determinant itself is outside the scope fence (SY 1.2), so it
  lives in the study.
- **FORM:** the same determinant, built independently as a product of row sums over a column-recording function,
  antisymmetrised after every row.
- **The check:** FORM subtracts PQCF's expression and must get exactly 0. For res(7,4), sympy's resultant also
  agrees.
- **Evaluation cost:** FORM 5's output optimisation (`Format O1`–`O4`, `#Optimize`), against PQCF's straight-line
  programs (`compile/3`), both counted the way FORM counts: powers, multiplications and additions. Everything ran
  single-threaded on core 8 of the provisional reference.

**Results:**

| | res(7,4) | res(7,5) | res(7,6) |
| --- | --- | --- | --- |
| variables | 13 | 14 | 15 |
| terms (PQCF = FORM; = sympy for res(7,4)) | 2,562 | 11,380 | 43,166 |
| FORM − PQCF | 0 | 0 | 0 |
| build: PQCF | 0.68 s | 4.6 s | 26.5 s |
| build: FORM 5 | 0.16 s | 1.17 s | 4.8 s |
| operations as expanded, by FORM's counter | 29,163 | 142,711 | 587,880 |
| the paper's "original" | 29,163 | 142,711 | 587,880 |
| PQCF straight-line program | 12,375 | 54,288 | 203,775 |
| FORM 5 `O1` | 4,968 | 20,210 | 71,262 |
| FORM 5 `O2` | 3,969 | 16,398 | 55,685 |
| FORM 5 `O3` | 3,458 | 11,006 | 35,945 |
| FORM 5 `O4` | 2,936 | 10,722 | 35,419 |
| the paper's best (MCTS/SLS, 2014) | ≈3,820 | ≈13,770 | ≈45,500 |
| time: PQCF `compile/3` | 0.02 s | 0.17 s | 0.92 s |
| time: FORM `O1` / `O4` (wall, with start-up) | 0.6 / 2.6 s | 0.8 / 9.8 s | 1.4 / 57 s |

**What this shows:**

- **The same expressions.** FORM's own counter, applied to PQCF's result, gives the paper's "original" operation
  counts exactly. The paper's term counts (2,561, 11,379, 43,165) are one less than the actual number of terms:
  they equal FORM's count of additions.
- **Building the expanded resultant:** PQCF is 4 to 6 times FORM's time, as on the other bulk-product tasks, and
  about 800 times faster than sympy on res(7,4) (0.68 s against 548 s).
- **Evaluation cost:**
  - PQCF's straight-line programs only share atom powers, monomials and identical polynomials. They do no Horner
    schemes and no common-subexpression search. That reduces the expanded form's operations 2.4 to 2.9 times, where
    FORM 5's `O4` reduces them 9.9 to 16.6 times, so PQCF needs 4.2 to 5.8 times `O4`'s operations.
  - FORM 5's `O4` is 20–25 % better than the 2014 paper's best results, which it grew out of.

  This is the largest gap the comparison has found that matters to PQCF. Compiled expressions are what PQCF
  evaluates under many bindings (SY 8.1).
- **Next step, a candidate:** multivariate Horner schemes with a greedy variable order, and common subexpressions
  across them, in `pq_sym_slp`. That is behind an option, because the operation order and so the results'
  rounding would change (D1 holds within one choice).
