# H1: Horner schemes in straight-line programs

**Decision: `pq_symbolic:compile/3` takes `optimize => horner`, which compiles every polynomial as a multivariate
Horner scheme. The variables are taken in an occurrence order fixed for the polynomial, and identical
sub-polynomials are shared. The default, `optimize => none`, is the program of `compile/2`, unchanged operation for
operation.** The project owner asked for this on 6 Nov 2026, after the resultant study of `SY-FORM-comparison.md`
found PQCF's compiled programs needing 4.2–5.8 times the operations of FORM 5's `Format O4`. It is additive.

## The scheme (`pq_sym_slp.erl`)

For a polynomial P with more than one term:

1. **Choose the variable:** the first atom of the order that at least two terms of P contain with a positive
   integer exponent. The order is computed once for the whole polynomial, by the number of terms an atom occurs in,
   most first, with ties broken by canonical order.
2. **Split:** P = R + x^k Q, where Q collects the terms containing x, divided by x^k, and k is the smallest such
   exponent.
3. **Recurse** on R and Q.

A single term becomes its monomial as before, and a polynomial no atom splits becomes the plain sum.

Every node is a slot keyed by its structure, so identical sub-polynomials, and identical products and sums of
slots, are computed once. Radicals, exponentials, conjugates, opaque bases and rational exponents are not split;
they stay inside monomials and go through the default path.

**Why a fixed order.** Choosing the variable afresh at every node (greedily) gave 4,968 / 20,334 / 72,334
operations on the three resultants. The fixed order gives 4,678 / 19,490 / 67,928, 4–6 % fewer. This is the order
FORM's `O1` uses; its `O3` and `O4` search for better orders (Monte Carlo tree search, local search), which is
left open.

## Determinism

The program depends only on the expression: canonical forms, a fixed order, and ties broken canonically. It is
reproducible bit for bit (D1) within the option. Against the default program the operation order differs, so the
values agree within rounding. That is why it is an option and not the default: existing results keep their bits.

## Measured (`comparison/form/resultants.erl`, core 8 of the provisional reference)

| | res(7,4) | res(7,5) | res(7,6) |
| --- | --- | --- | --- |
| default program | 12,375 | 54,288 | 203,775 |
| **Horner** | **4,678** | **19,490** | **67,928** |
| FORM 5 `O1` | 4,968 | 20,210 | 71,262 |
| FORM 5 `O2` | 3,969 | 16,398 | 55,685 |
| FORM 5 `O4` | 2,936 | 10,722 | 35,419 |
| compile time, default / Horner | 0.02 / 0.06 s | 0.16 / 0.44 s | 0.87 / 2.5 s |
| default against Horner at a random point (relative difference) | 0 | 1.6e-13 | 1.2e-12 |
| 200 evaluations, default / Horner | 0.80 / 0.24 s | 3.6 / 0.99 s | — |

- Horner is 2.6–3.0 times smaller than the default program, below FORM's `O1`, and 1.6–1.9 times FORM's `O4`.
- Evaluation is 3.4–3.6 times faster, so the operation count carries over to the time.
- **Open:** common-subexpression elimination across different Horner branches (FORM's `O2`), and a search over
  variable orders (`O3`, `O4`). Both have measured targets in the table.

## Tests (`test/pq_horner_tests.erl`)

- **Random polynomials** (60, in five symbols): values within 1e-12 of the default program's, and never more
  operations.
- **A resultant:** fewer operations, the same value, the same program twice.
- **The other atoms:** radicals, an exponential, a conjugate of a complex symbol and a fractional power, values
  within 1e-12.
- **A symbolic 3-mode mesh unitary:** every entry within 1e-12.
- **The default unchanged:** `compile/3` without the option is `compile/2`.
- **Mutation check:** a Horner split that lowers the exponent by one too many, caught by two tests.
