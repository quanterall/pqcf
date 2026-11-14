# K1: products of polynomials in Rust

**Decision: `pq_sym_poly:mul/3` multiplies plain polynomials in Rust (`pq_nif:sym_poly_mul/6`,
`native/pq_nif/src/polymul.rs`). The result is the same polynomial, term for term, as the product in Erlang (D1).
It is internal: no public function changes, and `poly_kernel => false` in the options turns it off.** Decided on
6 Nov 2026, as the first step of the plan to close the gap to FORM. The step comes from the analysis of how close PQCF can
come to FORM: expansion was 10–30× slower than FORM 5, and a 40-line Rust prototype multiplied Fateman's n = 15
benchmark in 0.13 s against FORM's 2.2 s.

## What goes to Rust

A product goes to the kernel when all of the following hold:

- both factors are *plain*: every monomial is a product of symbols, `conj` atoms and π with integer exponents (no
  radicals, opaque bases, exponentials or fractional exponents), and every coefficient is exact (`{q, _, _}`, not
  float);
- both factors have at least two terms, and the product has at least 64 term products. A monomial times a
  polynomial stays in Erlang, where it is faster: the conversion costs more than the product (1 × 4,096 terms:
  0.8×; 2 × 4,096: 1.7×).

## How

1. **Denominators.** Erlang scales each factor by the least common multiple L of its coefficients' denominators,
   so the coefficients become Gaussian integers. Each must fit in i128, and so must L_P · L_Q.
2. **Exponents.** For every variable, the range of exponents in the product is known from the two factors. Each
   variable gets just the bits its range needs, offset from its minimum, so negative exponents work. If all ranges
   fit in 64 bits, a monomial is one u64 and a product of monomials is one addition, with no carries.
3. **Product.** A hash map (a deterministic multiplicative hash) accumulates the coefficients with checked i128
   arithmetic. Real products skip the imaginary parts.
4. **Back to Erlang.** The NIF builds `pq_sym_poly`'s map directly: atom keys from the variable list, exponents
   `{E, 1}`, coefficients over L_P · L_Q in lowest terms with a positive denominator, exactly as
   `pq_sym_exact:rat/2` would. Building the map in Erlang from flat lists cost 1.2 s for 591k terms; built in the
   NIF it costs about 0.4 s.

The kernel refuses (`overflow` or `packing`) instead of approximating, and the product is then computed in
Erlang. A missing NIF falls back the same way. It runs on a dirty CPU scheduler. The process's heap limit still
applies to the result (`sym_max_heap_words`).

**Also changed: the canonical print order.** `pq_sym_norm:to_ast/1` sorted the terms twice: once in canonical key
order, then by the print key. The print key ends with the monomial, which is unique, so the first sort changed
nothing. It is gone, and the output is identical. On 591k terms, `to_ast` went from 7.1 to 2.4 s.

## Measured (comparison/form/compare.erl, core 8, single scheduler, 6 Nov 2026)

| Task | FORM 5 | before | after | after / FORM |
| --- | --- | --- | --- | --- |
| T2 dense product | 0.137 s | 1.42 s | **0.072 s** | 0.53× |
| T5 Fateman n = 10 | 0.136 s | 1.34 s | **0.052 s** | 0.38× |
| T5 Fateman n = 15 | 2.13 s | 33.1 s | **0.51 s** | 0.24× |
| T6 sparse n = 6 | 0.051 s | 1.28 s | 0.40 s | 7.8× |
| T6 sparse n = 8 | 0.388 s | 12.0 s | 2.83 s | 7.3× |
| T1 (1 + x1 + … + x6)^20 | 0.111 s | 2.59 s | 1.94 s | 17× |
| T3 permanent 8 × 8 | 0.162 s | 0.333 s | 0.305 s | 1.9× |
| T4 mesh k = 8 | 0.46 s | 2.2 s | 2.2–2.7 s | ~5× (no plain products) |
| resultant res(7,6), built | 4.8 s | 26.5 s | 23.9 s | 5× |

- **Dense products:** PQCF is now faster than FORM, 1.9–4.2×.
- **Sparse products:** 4.2× faster than before, still 7–8× FORM. The product itself takes 0.5 s of the 2.8 s; the
  rest is the canonical print order (`to_ast`'s sort of 591k terms) and converting to the expression tree.
- **Unchanged:** powers (multinomial enumeration in Erlang), the permanent and the resultant (sums of monomial ×
  polynomial products), and the mesh (exponentials, not plain).

The machine is shared, and times vary by up to 2× between runs on the small tasks. The rows above are from one run
of the whole comparison.

**The print order by binary keys** (7 Nov 2026). When every monomial is a product of symbols with integer
exponents, `to_ast` sorts binary keys whose byte order is `print_key`'s order:
`<<NotOne, Degree, Name, 1, ..., 0, Exponents>>`, with names as UTF-8 text and integers offset by 2^62 big-endian.
Otherwise it takes the general path, so the output is identical either way. In one process, alternating, on the
591k-term sparse product: 1.43 / 1.48 s against 3.35 / 3.65 s for the general path, 2.4 times faster. The machine
had a load average near 38 then, so the absolute times are high. `test/pq_print_order_tests.erl` checks that the
order equals the general path's, through a π term that forces it (80 random polynomials over names with prefixes,
digits, capitals, underscores and a non-ASCII one, negative exponents, a constant). Exponents without the offset are
caught.

After it, sparse n = 8 spends about 0.5 s multiplying, 1.3 s ordering, and 0.3–0.4 s copying the result out of
the process that bounds the heap (by design).

**Powers in Rust** (7 Nov 2026, the project owner's choice). `pq_sym_poly:pow_int` sends P^N for a plain P of
at least two terms to `pq_nif:sym_poly_pow/5`: repeated products in the packed layout of P^N's exponent ranges, exact
in i128, with the coefficients over L^N. Overflow, exponents that don't pack, L^N beyond i128, or no NIF: the
multinomial enumeration in Erlang, the same polynomial. On (1 + x1 + … + x6)^20 (230,230 terms) the expansion went
from 0.49 to 0.40 s, and `canonicalise` from about 1.7 to 1.2–1.4 s (load average near 30). The arithmetic takes a
few hundredths of a second. The rest is building the Erlang map of 230k terms, and the print order and expression
tree (`to_ast`, about 1 s).

**The floor of this design.** The polynomials are Erlang maps and results are Erlang expression trees, so every
large result pays for building, hashing, sorting and copying its terms as Erlang terms. That is now most of the time
left on the bulk tasks: FORM's 0.11 s for the power against about 1.2 s here. Going further means polynomials that
stay in Rust (resource handles) between operations, with only the final expression tree built in Erlang. That is a
redesign of `pq_sym_poly`, not a kernel. Done in K2 (7 Nov 2026): the power 0.48 s, sparse n = 8 0.92 s.

## Open (measured targets)

- **The print-order sort** of large results: now 1.3–1.5 s at 591k terms with binary keys. Sorting in Rust would
  have to return the terms, which costs about as much as it saves.
- **Sums** of many monomial × polynomial products, the resultant's and permanent's inner loop: an n-ary sum in the
  kernel, or polynomials kept in Rust as handles between operations.
- **Powers** of sums: the multinomial enumeration in Rust.

## Tests

- **Rust** (`polymul.rs`): 200 random products against a schoolbook product over exponent vectors, covering
  negative exponents, Gaussian coefficients and up to five variables; overflow, packing and cancellation.
- **EUnit** (`test/pq_poly_kernel_tests.erl`):
  - 120 random products of plain polynomials (rational and Gaussian coefficients, π, `conj`, negative exponents),
    kernel against Erlang, `=:=`;
  - the NIF called directly, with coefficients over a denominator;
  - the refusals: coefficients beyond i128 (scaled or in the product), exponents of 2^40, a radical, a float, a
    fractional exponent;
  - small and empty products;
  - end to end, Fateman's f(f + 1) and a Laurent product with Gaussian rationals, with and without the kernel.
- **Mutations:** a wrong denominator on the imaginary part, and a wrong sign in the complex product, are each
  caught (4 of 5 tests fail).
