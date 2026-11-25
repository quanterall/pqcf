# PQCF results certified by Lean 4

## Goals

PQCF's solvers are fast and not formally verified: polynomial factorisation, gcds, exact linear algebra,
reconstruction, integral reduction. Each is tested against independent programs, but a test covers the cases it
runs. For a single result, a stronger guarantee is possible: PQCF writes the result as a **Lean 4 statement** that a
reader can check by eye, with a **proof** that Lean's kernel checks by machine. PQCF stays the solver; Lean is the
judge of each result it is asked about.

- **Soundness of individual results**: the factors multiply back to the polynomial; the gcd divides both
  polynomials with the given cofactors; the solution satisfies every equation of the system. Exactly, over ℚ.
- **A small, named trusted base**: Lean's kernel, Mathlib and the pinned toolchain, and Lean's three standard axioms;
  each verdict lists the axioms its theorems use.
- **No build dependency**: PQCF writes the certificates itself (`pq_lean`, tested without Lean); the checking is
  done by [lean-check](#lean-check), a separate application ([quanterall/erlang-lean4](https://github.com/quanterall/erlang-lean4)), in the same node or on another.

Not certified (yet): that a result is complete or optimal where only its soundness is checked: that the gcd is
the *greatest* common divisor, that factors are irreducible, that a linear system has no other solution.

## lean-check

lean-check ([quanterall/erlang-lean4](https://github.com/quanterall/erlang-lean4), Apache 2.0) runs a supervised pool of Lean 4 workers with
Mathlib and returns a verdict per check: `valid` or `invalid` (with Lean's errors), `timeout` or `unavailable`, and
the axioms each theorem depends on. A submission is screened for constructs that would run code, and workers run
in a container without network access.

## Using it

`pq_lean` writes a certificate and, when a checker is available, returns lean-check's verdict:

| function | certifies | Lean proof |
| --- | --- | --- |
| `pq_lean:factor(E, Opts)` | `pq_symbolic:factor/1`'s result: E equals the product of the factors | `ring` |
| `pq_lean:gcd(A, B, Opts)` | `pq_symbolic:gcd/2`'s result G: A = G A' and B = G B' with the cofactors | `ring` |
| `pq_lean:linear(Rows, X, Opts)` | X solves the system, Rows = `[{[a_i1, .., a_in], b_i}]` | `norm_num` |

`Opts`: `checker => local | {node, Node} | none` (default: local when lean-check is loaded in the node, else none),
`timeout => Ms`. The result is `{ok, #{source => LeanSource, verdict => Verdict}}` (no verdict with `none`).

## Example (real output)

`examples/lean/run.sh` starts a lean-check node and a PQCF node on an internal Docker network (no outside access)
and certifies three results, then a wrong one. This is its output from 25 Nov 2026:

```erlang
{ok, E} = pq_symbolic:parse(<<"x^4 - y^4">>),
pq_lean:factor(E, #{checker => {node, 'lean_check@lc'}}).
```
```lean
theorem pqcf_factor (x y : ℚ) :
    (((x ^ 4) + ((-1) * (y ^ 4))) : ℚ) = (1 * ((1 * x) + ((-1) * y)) * ((1 * x) + (1 * y)) * ((1 * (x ^ 2)) + (1 * (y ^ 2)))) := by
  ring
```
```
-> valid (450 ms), axioms [{<<"pqcf_factor">>, [<<"propext">>, <<"Classical.choice">>, <<"Quot.sound">>]}]
```

The gcd of x^3 - x and x^2 + 2x + 1, with its cofactors:

```lean
theorem pqcf_gcd (x : ℚ) :
    (((x ^ 3) + ((-1) * x)) : ℚ) = ((1 + (1 * x)) * (((-1) * x) + (1 * (x ^ 2)))) ∧ (((x ^ 2) + (2 * x) + 1) : ℚ) = ((1 + (1 * x)) * (1 + (1 * x))) := by
  constructor <;> ring
```
```
-> valid (380 ms)
```

The solution x = 9/14, y = -13/28 of 3x + 2y = 1, x - 4y = 5/2:

```lean
theorem pqcf_linear :
    (((3 * ((9 : ℚ) / 14)) + (2 * ((-13 : ℚ) / 28))) : ℚ) = 1 ∧ (((1 * ((9 : ℚ) / 14)) + ((-4) * ((-13 : ℚ) / 28))) : ℚ) = ((5 : ℚ) / 2) := by
  norm_num
```
```
-> valid (604 ms)
```

A wrong solution (y = 13/28) is rejected: Lean reduces the statement to `False`.

```
-> invalid, errors
  unsolved goals
⊢ False
```

## Running the example

```sh
docker build -t pqcf-dev .
git clone https://github.com/quanterall/erlang-lean4 ../erlang-lean4 && (cd ../erlang-lean4 && docker build -t lean-check docker)
LEAN_CHECK_SRC=../erlang-lean4 sh examples/lean/run.sh     # about 20 s once both images exist
```

The lean-check image is large (11.8 GB: Lean 4.36.0, the Lean REPL and Mathlib's prebuilt cache) and takes 10 to
30 minutes to build the first time.

## Reading a verdict

`valid` means Lean accepted the proof of the statement shown, from the axioms listed. It does not mean the
statement is the claim one had in mind: read the statement. The certificates above are written by PQCF from its own
inputs and results, so the statement is exactly "this input equals that output" in each case.

## Next

The gcd's greatestness (a Bezout certificate u A + v B = G), the identities PQCF already compares with other
programs (Wigner symbols, colour factors, beta-function coefficients), exact steady states of master equations,
and interval certificates for high-precision numerical results.
