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

## PQCF's symbolic results

The results themselves, not only the comparison, are in [`results/symbolic/`](results/symbolic), rewritten by every
run of `run.sh`:
- **`<case>.reduction.txt`:** PQCF's reduction of each check's integrals. It lists the integral (in Forcer's
  notation), then every master integral PQCF reduced it to, with its coefficient as an exact rational function of
  d = 4 − 2ε.
- **`<case>.expansion.txt`:** that reduction times the masters' values, expanded in ε, in Forcer's G-scheme
  normalisation.
- **`<check>.msbar.txt`:** each of FORM's checks in MS-bar, as PQCF computes it (every one equal to FORM's
  assertion).

## Results

**Every answer agrees exactly, in every run:**
- 10 of 10 components equal Forcer's results, as ε-expansions to between ε³ and ε⁷;
- 5 of 5 MS-bar results equal FORM's assertions, every coefficient including ζ₃ to ζ₇.

For example, `Forcer_1` F2 = 15/2 ε⁻⁴ + 383/12 ε⁻³ − 2089/18 ε⁻² + (1466/9 − 728ζ₃ − 1245/8 ζ₅) ε⁻¹ − 2183/12 − 17198/9 ζ₃ − 3705/8 ζ₃² − 2169/2 ζ₄ − 57155/48 ζ₅ − 6225/16 ζ₆ + 441/2 ζ₇.

**Time** (dgx-spark; FORM on one core, PQCF on 16; the files in `results/`):

| check | FORM + Forcer | PQCF, 8 Nov (c166d7e) | 9 Nov (f615b0c) | 9 Nov (9779719) |
| --- | --- | --- | --- | --- |
| `Forcer_example` | 0.51 s | 0.07 s | 0.06 s | **0.04 s** |
| `Forcer_1` | 10.3 s | 918 s | 272 s | **175 s** |
| `Forcer_2` | 5.2 s | 767 s | 174 s | **132 s** |

PQCF's per-case times in the latest run:

| case | seconds | equations | needed | masters |
| --- | --- | --- | --- | --- |
| no1 F1 | 6.8 | 1,864 | 0 (a master) | 1 |
| no1 F2 | 48 | 116,202 | 67,035 | 21 |
| no2 F1 | 16.6 | 1,417 | 0 (a master) | 1 |
| no2 F2 | 34 | 87,919 | 51,715 | 16 |
| no3 F1 | 47 | 33,454 | 22,098 | 12 |
| no3 F2 | 23 | 92,723 | 53,259 | 13 |
| no1 numerators | 2.9 | 13,244 | 10,076 | 20 |
| no2 numerators | 105 | 97,949 | 58,123 | 16 |

**What changed between the first two runs:**

1. **Seeds sector by sector.** The search starts from the targets' budget. Each master that is not itself a seed
   raises only its own sector's budget, and the system is built again. no2's numerator variants: 374,416 equations
   → 97,949. no3 F1: 44,053 → 33,454.
2. **One reconstruction for a sum** (`pq_ibp:reduce_sum/3`). Each check needs only the sum of its integrals, so one
   extra unknown x and one equation x − Σ c_T T = 0 are added, and only x is reconstructed. Before, every target's
   coefficient of every master was reconstructed separately: dozens of targets, after the numerators are written
   as integrals.
3. **Zero right-hand sides left out** (`pq_recon:linear/4`). The large systems carry one right-hand side per
   unreduced integral at their edge, 19,841 for no1, and only 21 are non-zero. They are now found first and
   dropped.

**What changed in the third run:**

4. **The completeness check solves once.** The seed search checks that every master with a coefficient is a seed.
   It now does this with one solve on a random combination of the targets, not one solve per target (no1 F2:
   37.7 s → 5.2 s).
5. **Fewer, cheaper reconstructions.**
   - Only the masters that check found non-zero are reconstructed.
   - The 1.36M-entry program is stored once per reconstruction, not once per column.
   - With one variable, the columns run concurrently.
6. **The seed search orders its budgets differently.** The global budgets come first, while their seeds stay
   within three times the smallest budget's; a larger budget hands over to the sector-by-sector search.
   - This fixes a regression the sector search had caused on the double box's harder target: 3.0 s instead of
     1.2 s, now 1.4 s.
   - The four-loop cases still reach the sector search. no3 F1 pays for the global attempts first: 47 s here,
     30 s in the run before.

**Tried and not kept:**
- **Every sector but the targets' own starting from no dots and no numerators.** The learning elimination filled
  in badly (9 GB, no result after 15 minutes on no1 F2).
- **Each target's own budget instead of the maximum over all targets:** {1 dot} ∪ {1 numerator} instead of {1 dot,
  1 numerator}. Its ragged budgets did the same: no result in 4 hours.

**Reading these honestly:**
- **On the four-loop checks Forcer is still 17–25 times faster** (75–150 in the first run). Forcer reduces with identities solved once,
  parametrically, for each of its 437 topologies. PQCF solves a Laporta system for each problem.
- **PQCF's remaining cost** is mostly the system: building it and learning its elimination (about 30 s of no1 F2's
  48 s), and the seeds that a dotted target needs in every sector.
- **On the three-loop example PQCF is 13 times faster** (0.04 s against 0.51 s).
- **PQCF's symmetry classes for a four-loop family took 75 s before this comparison;** they now take 0.7 s, by
  matching lines by invariants instead of trying all 11! permutations.
