# K2: polynomials that stay in Rust between operations

**Decision: inside `pq_symbolic:canonicalise/1,2`, a polynomial product or power that the kernel (K1) computes stays
in Rust as `{rpoly, Vars, Resource}`. Sums, negations, scalings, products and powers of such polynomials run in Rust.
`pq_sym_norm:to_ast/1` builds the final expression tree in print order straight from Rust. Nothing outside
`pq_sym_poly` ever sees a handle, and the canonical form is identical, term for term, to the maps-only path's (D1).**
The project owner asked for this on 7 Nov 2026 ("go with the Rust-resident polynomials"), after K1 found that
building, hashing, ordering and copying large results as Erlang terms was the floor of the old design.

## Design

- **Where handles exist.** `canonicalise` sets the context flag `rpoly` (only with `poly_kernel`, the default). Under
  it, `pq_sym_poly:mul/3` makes a handle when both factors are plain with at least two terms and at least 64 term
  products, or when either factor already is a handle. `pow_int/3` does so for a plain sum within
  `max_expand_power`. Every other caller of `pq_sym_norm:norm/2` (compile, evaluate, gradients, equivalence,
  hashing, the permanent, …) runs without the flag and gets maps exactly as before.
- **Handle-aware functions.** `add`, `sum`, `neg`, `scale` (by an exact coefficient), `mul`, `pow_int`, `terms`,
  `is_const` and `single` work on handles. A plain map operand is converted (`pq_nif:rpoly_new`), and a non-plain
  one turns the handle back into a map. Every other function of `pq_sym_poly` (`exp`, `conj`, `pow_rat`,
  `base_power`, `content`, `sorted`, `const_value`) converts first (`pq_sym_poly:m/1`, `pq_nif:rpoly_to_map`).
  So the rest of the engine, and anything that pattern-matches polynomials as maps, never meets one.
- **In Rust** (`native/pq_nif/src/rpoly.rs`):
  - A term is its exponents plus a Gaussian-integer numerator over one common denominator, kept in lowest terms.
  - The variables stay on the Erlang side, sorted. An operation on two handles gets a map of each operand's
    variables into their union.
  - Products and powers use K1's packed kernel; sums use a packed hash merge. Each result is exact, or refused
    (overflow, packing), and then the maps do the work.
- **Output.** `rpoly_to_ast` sorts the terms by `print_key`'s order: the constant first, then the degree, then the
  names present as a list, then the exponents. It builds the terms `{mul, [Coeff | Factors]}`, `{pow, {sym, X}, E}`,
  `{rat, N, D}` and `{complex, Re, Im}` directly. It applies when every variable is a symbol; with π or `conj` the
  handle becomes a map and takes the Erlang path.
- **Limits.** A handle's memory is outside the process heap, so each new one is checked against
  `sym_max_heap_words` (terms × (variables + 4) words) and `max_poly_terms`. Resources are freed by the garbage
  collector with their last reference, at the latest when the bounded process ends.

## Measured (`comparison/form/compare.erl`, core 8, single scheduler; FORM's times as on a quiet machine)

| Task | FORM 5 | K1 | **K2** | K2 / FORM |
| --- | --- | --- | --- | --- |
| (1 + x1 + … + x6)^20, 230,230 terms | 0.111 s | 1.94 s | **0.48 s** | 4.3× |
| (1 + x1 + … + x6)^16 | 0.031 s | 0.50 s | 0.12 s | 4.0× |
| sparse n = 8, 591,235 terms | 0.382 s | 2.83 s | **0.92 s** | 2.4× |
| sparse n = 6 | 0.051 s | 0.40 s | 0.14 s | 2.7× |
| Fateman n = 15 | 2.14 s | 0.51 s | **0.46 s** | 0.22× |
| Fateman n = 10 | 0.137 s | 0.052 s | 0.035 s | 0.26× |
| dense product | 0.136 s | 0.072 s | 0.037 s | 0.27× |
| permanent 8 × 8 | 0.162 s | 0.31 s | 0.35 s | 2.2× (maps; no handles) |
| mesh k = 8 | 0.46 s | 2.2 s | **0.67 s** | 1.5× (with exponentials, below) |

Since the 2 Nov comparison: dense products went from 10–25× slower than FORM to about 4× faster, and sparse
products and powers from 25–32× and 10–25× slower to 2.4–4.3×. What remains on them is building the expression tree
of the result and copying it out of the process that bounds the heap. The permanent and the mesh don't use handles:
their products are monomial × polynomial (the mesh, with its exponentials, now uses handles: below).

## Exponentials of symbols, and circuit unitaries (7 Nov 2026)

An exponential factor exp(k i x) of a symbol x, with k an integer, only ever merges with another one by adding the k
(`normalize/3` keeps such factors as they are, and drops a zero one). So the kernel takes it as a variable
`{expi, X}` with exponent k. The Rust side builds the monomial's exponential list `[{{[{{sym, X}, {1, 1}}], []},
{q, {0, 1}, {K, 1}}}, …]` and, in `to_ast`, the factor `{exp, {mul, [{complex, 0, K}, {sym, X}]}}` after the atoms.
The print order compares the list of exponential factors (symbol, then k) after the atoms, as `print_key` does
through the monomial.

`pq_symbolic:circuit_unitary/1,2` now also sets the flag: its entries go straight to `value_to_ast`. A mesh's entries
are sums of products of exp(±i θ) and exp(±i φ), so every product of a gate entry with a large entry stays in Rust.

| mesh, `circuit_unitary` | FORM 5 | maps only | **kernel + handles** |
| --- | --- | --- | --- |
| k = 6 | 0.023 s | 0.061 s | **0.027–0.030 s** |
| k = 8, 276,800 terms | 0.46 s | 1.6–2.8 s | **0.58–0.67 s** (1.3–1.5× FORM) |

The same unitary, term for term, for meshes of 3–5 modes with phase shifts, for powers and products of sums of
exponentials (`exponentials_test_`), and for the random expressions, which include exp(i x). The k order of the
exponential list reversed is caught (9 failing tests).

## Tests (`test/pq_rpoly_tests.erl`)

- **200 random expressions**, `canonicalise` with handles against `poly_kernel => false`, `=:=`. They mix sums,
  products, powers of sums, rationals, Gaussian coefficients, π, `i`, conjugates, exponentials, square roots and
  negative powers of sums.
- **The benchmark shapes**, including a product that cancels to zero, and F × F⁻¹.
- **Matrices:** a matrix product and a Kronecker product of polynomial entries.
- **Every handle operation** against its map operation (20 random pairs): product, sums both ways, cancellation to
  zero, exact and float scaling, powers of a map and of a handle, conjugation, `to_ast`, `terms`, `is_const`.
- **The memory limit:** a handle beyond `sym_max_heap_words` gives `{symbolic_limit, #{reason := heap}}`.
- **Mutations caught:** a sum that forgets to rescale the second operand's denominator (1 failing test), and the
  print order with the name comparison reversed (4).
