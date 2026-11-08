# PQCF and FORM's Forcer: four-loop massless propagators

FORM's own test suite checks Forcer, its program for four-loop massless propagator integrals, on fixed inputs with
exact expected results ([`check/extra/forcer.frm`](https://github.com/form-dev/form/blob/master/check/extra/forcer.frm)).
These are the integrals behind the five-loop QCD β-function. This directory reduces the same integrals with PQCF
and checks the answers against Forcer's and against the numbers FORM's checks assert.

- **Forcer:** B. Ruijl, T. Ueda, J.A.M. Vermaseren, *Forcer, a FORM program for the parametric reduction of
  four-loop massless propagator diagrams*, [arXiv:1704.06650](https://arxiv.org/abs/1704.06650); also
  [arXiv:1604.08767](https://arxiv.org/abs/1604.08767).
- **FORM 5:** [arXiv:2601.19982](https://arxiv.org/abs/2601.19982).
- **The master integrals:** R.N. Lee, A.V. Smirnov, V.A. Smirnov, [arXiv:1108.0732](https://arxiv.org/abs/1108.0732),
  and P.A. Baikov, K.G. Chetyrkin, [arXiv:1004.1153](https://arxiv.org/abs/1004.1153). Their values are taken as
  Forcer ships them, in Forcer's basis.
- **Where these integrals are used:** F. Herzog, B. Ruijl, T. Ueda, J.A.M. Vermaseren, A. Vogt, *The five-loop beta
  function of Yang-Mills theory with fermions*, [arXiv:1701.01404](https://arxiv.org/abs/1701.01404).
- **PQCF's method:** Laporta's algorithm ([arXiv:hep-ph/0102033](https://arxiv.org/abs/hep-ph/0102033)) with
  finite-field reconstruction (decision record R1).

## The checks

| FORM's check | integrals |
| --- | --- |
| `Forcer_example` | a three-loop propagator with two scalar-product numerators, and the two-loop kite with insertions (lines of power 1 + ε) on two lines |
| `Forcer_1` | the three four-loop non-planar top topologies no1, no2, no3 (11 lines, 14 invariants): every line once (F1), and the sum of their 14 variants with one squared line or one numerator (F2) |
| `Forcer_2` | no1 and no2 only: F1 as above, F2 the numerator variants only |

## How the comparison works

1. **Same integrals.** `forcer_suite.erl` holds each case once. PQCF reduces it and writes the same integrals in
   Forcer's notation (`vx(...)` vertices, line powers, numerators, `ex(...)` insertions), so both programs work on
   the same input.
2. **PQCF's side.** From the graph, PQCF derives:
   - a loop basis, from momentum conservation;
   - auxiliary invariants that complete the scalar products;
   - insertions, as lines whose power is index + ε (`pq_ibp_family`, `eps`).

   `pq_ibp:reduce/3` then reduces every target to PQCF's own masters, exactly in d.
3. **Check 1: against Forcer's results.** Forcer evaluates each of PQCF's masters. PQCF's coefficients times those
   values must equal Forcer's result for the target. Both are compared as ε-expansions with ζ-values (`epsexp.erl`),
   as deep as the master tables allow.
   - **What PQCF computes itself:** the G-scheme constants come from PQCF's own Gamma-function expansion, which
     reproduces Forcer's tables, and so does the MS-bar factor ε·G(1,1)·e^{γε}.
   - **What comes from Forcer:** the values of the genuinely four-loop master integrals.
4. **Check 2: against FORM's assertions.** The sums are converted to MS-bar as Forcer's `msbarexpand=4` does and
   compared with the expected results written in `forcer.frm`.

`run.sh` does all of it on the committed HEAD (images: `comparison/form/Dockerfile`, then
`comparison/forcer/Dockerfile`) and writes `results/<host>-<date>.json`. It exits 1 if any check fails.

## Results (8 Nov 2026, dgx-spark, 16 cores for PQCF, FORM on one core; commit c166d7e)

**Every answer agrees exactly:**
- 10 of 10 components equal Forcer's results, as ε-expansions to between ε³ and ε⁷;
- 5 of 5 MS-bar results equal FORM's assertions, every coefficient including ζ₃ to ζ₇.

For example, `Forcer_1` F2 = 15/2 ε⁻⁴ + 383/12 ε⁻³ − 2089/18 ε⁻² + (1466/9 − 728ζ₃ − 1245/8 ζ₅) ε⁻¹ − 2183/12 − 17198/9 ζ₃ − 3705/8 ζ₃² − 2169/2 ζ₄ − 57155/48 ζ₅ − 6225/16 ζ₆ + 441/2 ζ₇.

**Time:**

| check | FORM + Forcer (1 core) | PQCF (reduction only, 16 cores) |
| --- | --- | --- |
| `Forcer_example` | 0.51 s | 0.07 s |
| `Forcer_1` | 10.3 s | 918 s |
| `Forcer_2` | 5.2 s | 767 s |

PQCF's per-case times:

| case | seconds | equations | needed | masters |
| --- | --- | --- | --- | --- |
| no1 F1 | 6.7 | 1,864 | 0 (a master) | 1 |
| no1 F2 | 384 | 116,202 | 67,035 | 21 |
| no2 F1 | 15.8 | 1,417 | 0 (a master) | 1 |
| no2 F2 | 221 | 87,919 | 51,715 | 16 |
| no3 F1 | 66 | 44,053 | 31,283 | 10 |
| no3 F2 | 225 | 92,723 | 53,259 | 13 |
| no1 numerators | 16 | 13,244 | 10,076 | 21 |
| no2 numerators | 729 | 374,416 | 244,587 | 16 |

**Reading these honestly:**
- **On the four-loop checks Forcer is 75–150 times faster.** Forcer reduces with identities solved once,
  parametrically, for each of its 437 topologies. PQCF solves a Laporta system for each problem. This is the
  expected advantage of a program built for one class of integrals, and the result PQCF has to improve on.
- **The automatic seed choice is most of PQCF's cost here.** For no2's numerators it settles on 374,416 equations,
  where the system from the next smaller seeds is 42,217 equations and reduces in 9 s; the smaller seeds fail only
  the completeness check. Choosing seeds per sector is the next step.
- **PQCF's symmetry classes for a four-loop family took 75 s before this comparison;** they now take 0.7 s, by
  matching lines by invariants instead of trying all 11! permutations.
- **On the three-loop example PQCF is 7 times faster** (0.07 s against 0.51 s).
