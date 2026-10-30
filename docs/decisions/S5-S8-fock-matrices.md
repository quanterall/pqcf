# S5 and S8: displacement, squeezing and beamsplitter matrices in Fock space

**Decisions** (closed 30 October 2026):

- **S8 (beamsplitter sectors):** AN 7.2 as written. Each sector matrix agrees with a 80-digit matrix exponential
  to 3.5e-15 (N = 0..20, three angle pairs), is unitary to 1e-12 up to N = 100, and the eigenvalue snap never
  exceeded its 1e-10 tolerance. The single-photon case reproduces the beamsplitter matrix of FS 3.3.
- **S5 (displacement and squeezing):** the column recurrences of AN 7.4 are **not used**: their error grows
  exponentially with the cutoff (amendment A-7 below). Instead:
  - displacement: the Laguerre form, by the three-term recurrence along each diagonal of the lower triangle and
    the symmetry `D[n][m] = (-1)^(m-n) conj(D[m][n])`; at most 5.5e-14 on the whole grid;
  - squeezing: the recurrence of the Bargmann generating function (Miatto and Quesada 2020; every coefficient
    of modulus at most 1) where it is accurate: `C <= 20`, or `r <= 0.1`, or `r >= 2.5`, or `C <= 40` with
    `r <= 0.25` or `r >= 1.5`. Elsewhere the exponential of the generator in a padded space, by a tridiagonal
    eigendecomposition per parity, doubling the padding until the block stops changing (by less than
    `max(2e-13, 1e-15 * size)`; at most 4,096 states per parity). At most 3.9e-13 on the whole grid.

## Method

`native/pq_nif/examples/s5_s8.rs` prints the matrices of every method; `test/oracle/s5_s8_compare.py`
compares them with closed forms in mpmath 1.3.0 at 80 digits: displacement by the generalised Laguerre
polynomials (as finite sums), squeezing by its normal-ordered form
`S = exp(-t a^dag^2 / 2) sech(r)^(n + 1/2) exp(conj(t) a^2 / 2)`, `t = e^(i phi) tanh r`, the beamsplitter by
`mpmath.expm` of the sector generator. The oracle was checked against itself at 40, 80 and 200 digits on the
worst entry. Grid: cutoffs C in (10, 20, 40, 80); r in (0.05 .. 3); |alpha| in (0.1 .. 10); phase 0.7. Error:
the largest absolute entry error of the (C+1) x (C+1) matrix (every entry has modulus at most 1).
`examples/s5_conv.rs` traces the doubling of the padding; `examples/s5_pad.rs` the error against the padding.

## Squeezing

| C | r | AN 7.4 recurrence | normal-ordered sum | Bargmann recurrence | chosen (Bargmann or padded eigen) |
|---|---|---|---|---|---|
| 10 | 0.05 | 6.7e-16 | 7.8e-16 | 2.2e-16 | 2.2e-16 |
| 10 | 0.1 | 8.9e-16 | 1.7e-15 | 5.6e-16 | 5.6e-16 |
| 10 | 0.25 | 1.0e-15 | 1.5e-15 | 3.7e-16 | 3.7e-16 |
| 10 | 0.5 | 3.4e-14 | 2.5e-15 | 5.8e-16 | 5.8e-16 |
| 10 | 0.75 | 7.3e-13 | 2.1e-15 | 3.9e-16 | 3.9e-16 |
| 10 | 1 | 5.8e-11 | 1.2e-15 | 3.7e-16 | 3.7e-16 |
| 10 | 1.5 | 2.3e-08 | 5.6e-16 | 1.9e-16 | 1.9e-16 |
| 10 | 2 | 3.8e-06 | 7.8e-16 | 1.4e-16 | 1.4e-16 |
| 10 | 2.5 | 5.3e-04 | 3.2e-16 | 1.1e-16 | 1.1e-16 |
| 10 | 3 | 7.9e-02 | 2.2e-16 | 6.8e-17 | 6.8e-17 |
| 20 | 0.05 | 1.3e-15 | 3.8e-15 | 3.3e-16 | 3.3e-16 |
| 20 | 0.1 | 1.2e-15 | 3.9e-15 | 8.9e-16 | 8.9e-16 |
| 20 | 0.25 | 1.9e-14 | 8.3e-15 | 1.1e-15 | 1.1e-15 |
| 20 | 0.5 | 3.4e-11 | 5.3e-14 | 5.9e-15 | 5.9e-15 |
| 20 | 0.75 | 6.8e-08 | 8.0e-14 | 3.9e-15 | 3.9e-15 |
| 20 | 1 | 5.0e-05 | 1.3e-13 | 4.4e-15 | 4.4e-15 |
| 20 | 1.5 | 1.7e+01 | 2.1e-14 | 5.3e-16 | 5.3e-16 |
| 20 | 2 | 2.9e+05 | 3.9e-15 | 2.0e-16 | 2.0e-16 |
| 20 | 2.5 | 1.1e+10 | 2.3e-15 | 1.7e-16 | 1.7e-16 |
| 20 | 3 | 3.6e+14 | 4.8e-16 | 1.1e-16 | 1.1e-16 |
| 40 | 0.05 | 1.7e-15 | 1.2e-14 | 8.7e-16 | 8.7e-16 |
| 40 | 0.1 | 8.6e-15 | 3.8e-14 | 1.7e-15 | 1.7e-15 |
| 40 | 0.25 | 8.6e-12 | 1.6e-12 | 5.6e-14 | 5.6e-14 |
| 40 | 0.5 | 4.1e-05 | 5.9e-11 | 1.6e-12 | 2.2e-14 |
| 40 | 0.75 | 1.1e+02 | 2.6e-10 | 6.1e-12 | 4.9e-14 |
| 40 | 1 | 1.1e+08 | 2.5e-10 | 2.7e-12 | 7.1e-14 |
| 40 | 1.5 | 7.6e+17 | 1.7e-11 | 1.2e-13 | 1.2e-13 |
| 40 | 2 | 9.7e+26 | 7.1e-13 | 9.8e-15 | 9.8e-15 |
| 40 | 2.5 | 1.3e+36 | 4.6e-14 | 3.0e-16 | 3.0e-16 |
| 40 | 3 | 6.5e+44 | 5.2e-15 | 1.4e-16 | 1.4e-16 |
| 80 | 0.05 | 5.1e-15 | 1.8e-13 | 3.3e-15 | 3.3e-15 |
| 80 | 0.1 | 3.7e-13 | 3.1e-12 | 4.2e-14 | 4.2e-14 |
| 80 | 0.25 | 1.2e-06 | 9.8e-09 | 4.9e-11 | 1.2e-14 |
| 80 | 0.5 | 9.6e+06 | 2.0e-05 | 7.9e-08 | 5.5e-14 |
| 80 | 0.75 | 1.0e+20 | 3.6e-04 | 8.6e-07 | 1.6e-13 |
| 80 | 1 | 2.6e+31 | 1.8e-04 | 8.8e-07 | 3.2e-13 |
| 80 | 1.5 | 8.1e+51 | 2.1e-06 | 2.0e-08 | 2.1e-13 |
| 80 | 2 | 5.7e+70 | 2.0e-09 | 1.0e-11 | 3.9e-13 |
| 80 | 2.5 | 2.6e+89 | 9.7e-12 | 5.7e-14 | 5.7e-14 |
| 80 | 3 | 2.0e+106 | 2.6e-13 | 7.2e-16 | 7.2e-16 |

## Displacement

| C | |alpha| | AN 7.4 recurrence | Laguerre diagonals (chosen) |
|---|---|---|---|
| 10 | 0.1 | 2.2e-16 | 6.7e-16 |
| 10 | 0.5 | 6.1e-16 | 1.6e-15 |
| 10 | 1 | 2.6e-15 | 2.8e-16 |
| 10 | 2 | 1.1e-13 | 2.8e-16 |
| 10 | 3 | 1.0e-12 | 1.9e-16 |
| 10 | 4 | 8.9e-13 | 2.4e-16 |
| 10 | 5 | 3.6e-13 | 3.2e-16 |
| 10 | 6 | 5.3e-14 | 2.0e-16 |
| 10 | 8 | 7.8e-18 | 6.0e-19 |
| 10 | 10 | 9.6e-24 | 1.5e-24 |
| 20 | 0.1 | 3.3e-16 | 5.2e-15 |
| 20 | 0.5 | 2.0e-15 | 2.2e-15 |
| 20 | 1 | 2.9e-14 | 1.3e-15 |
| 20 | 2 | 1.2e-11 | 7.6e-16 |
| 20 | 3 | 5.5e-10 | 5.2e-16 |
| 20 | 4 | 5.6e-08 | 7.1e-16 |
| 20 | 5 | 1.4e-07 | 8.6e-16 |
| 20 | 6 | 6.0e-08 | 5.6e-16 |
| 20 | 8 | 4.5e-10 | 5.2e-16 |
| 20 | 10 | 4.1e-14 | 6.7e-18 |
| 40 | 0.1 | 3.7e-16 | 5.8e-15 |
| 40 | 0.5 | 6.9e-15 | 3.1e-15 |
| 40 | 1 | 1.4e-12 | 1.6e-15 |
| 40 | 2 | 4.8e-08 | 1.0e-15 |
| 40 | 3 | 8.3e-05 | 8.6e-16 |
| 40 | 4 | 3.3e-02 | 2.5e-15 |
| 40 | 5 | 6.6e+00 | 2.8e-15 |
| 40 | 6 | 7.3e+01 | 3.0e-15 |
| 40 | 8 | 3.5e+02 | 3.3e-15 |
| 40 | 10 | 3.5e+00 | 1.1e-15 |
| 80 | 0.1 | 6.7e-16 | 5.5e-14 |
| 80 | 0.5 | 7.8e-14 | 5.8e-15 |
| 80 | 1 | 3.5e-10 | 4.0e-15 |
| 80 | 2 | 4.9e-04 | 2.2e-15 |
| 80 | 3 | 8.3e+02 | 2.0e-15 |
| 80 | 4 | 1.6e+08 | 3.5e-15 |
| 80 | 5 | 1.9e+11 | 4.2e-15 |
| 80 | 6 | 1.6e+15 | 3.2e-15 |
| 80 | 8 | 1.7e+20 | 6.2e-15 |
| 80 | 10 | 2.9e+21 | 4.4e-15 |

## Proposed amendment A-7 (AN 7.4)

Replace the column recurrences of AN 7.4 by the methods above, with the published region and the fallback.
The AN 7.4 recurrences multiply by `cosh r`, `sinh r` and `|alpha|` at every column, so rounding errors grow
geometrically: at C = 80 they are wrong by more than 1 for r >= 0.5 and |alpha| >= 3. The padded eigen
fallback costs O(size^3) with a dense eigensolver (about a second at size 1,358); a tridiagonal eigensolver
that keeps only the first C/2 + 1 rows of the eigenvectors would make it O(size^2 C), a Phase 3 optimisation.

## Follow-ups

- GT-C02 (squeezing and displacement on the Gaussian side) in Phase 1.
- Cutoffs above 80 and |alpha| above 10 are untested.
