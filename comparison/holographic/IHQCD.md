# Improved holographic QCD: scalar and tensor glueballs

Produced by `ihqcd.erl` (decision record SY2, phase 6), PQCF 3ab99f7.
The model is Gursoy, Kiritsis and Nitti, [arXiv:0707.1349](https://arxiv.org/abs/0707.1349), section 8.1
(background I; the framework is [arXiv:0707.1324](https://arxiv.org/abs/0707.1324)).

## Set-up

- **Model:** X(lambda) from (8.4) and W(lambda) from (8.6), with b0 = 4.2, b1/b0^2 = 51/121 and a = 5/16 (linear
  confinement, alpha = 2).
- **Background** (8.5): lambda' = -(4/3) X W lambda e^A and A' = -(4/9) W e^A (l = 1), from
  lambda(0.36) = lambda0, A(0.36) = 0.
- **Glueballs** (6.2): xi'' + 2 B' xi' + m^2 xi = 0, with B = (3/2) A + (1/2) log X^2 for 0++ and (3/2) A for
  2++.
- **Integration:** everything is integrated together by `pq_jet` (12 digits).
- **The UV boundary** r_b, where e^(-A) -> 0, is located by extrapolation: the solution with A(0.36) = 0 has it
  near r = -0.46, not 0. The normalisable solution xi ~ (r - r_b)^4 (psi ~ (r - r_b)^(5/2), (8.10)) starts
  10^-3 from it.
- **Eigenvalues:** where e^(2B) xi', the growing mode's coefficient, changes sign at r = 10.

## Results, b0 = 4.2 (masses in units of 1/l; ratios to the lowest 0++)

| state | m | m/m(0++), PQCF | the paper (8.12), (8.13) |
| --- | --- | --- | --- |
| 0++ 1 | 1.2921 | 1.000 | 1.000 |
| 0++ 2 | 2.1741 | 1.683 | 1.866 |
| 0++ 3 | 2.6819 | 2.076 | 2.414 |
| 0++ 4 | 3.1247 | 2.418 | 2.883 |
| 2++ 1 | 1.6915 | 1.309 | 1.393 |
| 2++ 2 | 2.3292 | 1.803 | 2.028 |
| 2++ 3 | 2.8130 | 2.177 | 2.535 |
| 2++ 4 | 3.2260 | 2.497 | 2.980 |

## Checks

- **Independence of the integration constants** (section 7.2 of the paper): with lambda0 = 0.25 the 0++ masses
  change (to [2.264,3.809,4.699,5.475]), but the ratios do not: [1.0,1.683,2.076,2.418].
  lambda0 = 0.08 gives the same ratios too.
- **Cutoff:** r = 14 gives the same values as r = 10.
- **The model's own consistency:** d log W/d lambda = -(4/3) X/lambda (2.22) holds for (8.4) and (8.6).

## Discrepancy with the paper (open)

- **The values differ:** PQCF gives R00 = 1.683 and R20 = 1.309 at b0 = 4.2. The paper fits b0 = 4.2 to
  R00 = 1.87 and finds R20 = 1.40. The excited 0++ ratios are lower too.
- **What the checks rule out:** the computation does not depend on the integration constants or the cutoff,
  as the paper's own argument requires, so the difference is not numerical noise.
- **Possible sources:** a convention not stated in the paper's formulas (the normalisation of lambda in X
  against W, or where the UV condition is imposed), or a difference in the paper's numerics. With A(r_in) = 0
  at r_in = 0.36, the paper's stated initial conditions put the boundary at r ~ -0.46, while (8.10) imposes
  the UV behaviour at r -> 0.
- **Not resolved:** a scan in b0 found R00 = 1.87 near b0 = 2, but with R20 = 1.09, so a single change of b0
  does not reproduce both ratios. That scan used a coarse grid and is not reliable enough to settle this.

Times: 0++ 181.9 s, 2++ 95.5 s, the lambda0 check 374.6 s.
