# Improved holographic QCD: scalar and tensor glueballs

Produced by `ihqcd.erl` (decision record SY2, phase 6), PQCF b91fba4.
The model is Gursoy, Kiritsis and Nitti, [arXiv:0707.1349](https://arxiv.org/abs/0707.1349), section 8.1
(background I; the framework is [arXiv:0707.1324](https://arxiv.org/abs/0707.1324)).

## Set-up

- **Model:** X(lambda) from (8.4) and W(lambda) from (8.6), with b0 = 4.2, b1/b0^2 = 51/121 and a = 3/8.
- **The value of a.** The paper uses the symbol a in two ways. Its IR asymptotics (8.2) are
  X = -1/2 - a/log(lambda), with a = (3/8)(alpha - 1)/alpha, so a = 3/16 for linear confinement (alpha = 2). Its
  interpolating X (8.4) tends instead to -1/2 - a/(2 log lambda), so through (8.4) and (8.6) linear confinement
  needs a = 3/8. The background confirms it: A(r) ~ -C r^alpha with alpha -> 2 (below).
- **Background** (8.5): lambda' = -(4/3) X W lambda e^A and A' = -(4/9) W e^A (l = 1), from
  lambda(0.36) = lambda0, A(0.36) = 0.
- **Glueballs** (6.2): xi'' + 2 B' xi' + m^2 xi = 0, with B = (3/2) A + (1/2) log X^2 for 0++ and (3/2) A for
  2++.
- **Integration:** everything is integrated together by `pq_jet` at 12 digits; eigenvalues by the Illinois
  method to a relative width of 1e-9. The ratios below are numerical results at that precision, not the 25-30
  digits of the soft-wall, hard-wall and Sakai-Sugimoto checks.
- **The UV boundary** r_b, where e^(-A) -> 0, is located by extrapolation (the solution with A(0.36) = 0 has it
  below r = 0). The solution starts 10^-3 from it as xi ~ (r - r_b)^4 (psi ~ (r - r_b)^(5/2), (8.10)).
- **Eigenvalues:** where e^(2B) xi', the growing mode's coefficient, changes sign at the cutoff.
- **Units of the runs:** the IR scale r1, where lambda = 1, moves by a factor 40 between lambda0 = 0.25 and
  0.05. The cutoff and the scanned range of m^2 are set in units of r1, so that every lambda0 covers the same
  states: cutoff 8 r1, and 11 r1 for the check and the long spectra (beyond about 12 r1, lambda exceeds
  10^150 and the integration overflows).

## Results, b0 = 4.2, lambda0 = 0.15 (masses in units of 1/l; ratios to the lowest 0++)

| state | m | m/m(0++), PQCF | the paper (8.12), (8.13) | difference |
| --- | --- | --- | --- | --- |
| 0++ 1 | 1.3465 | 1.000 | 1.000 | +0.0 % |
| 0++ 2 | 2.5251 | 1.875 | 1.866 | +0.5 % |
| 0++ 3 | 3.2699 | 2.428 | 2.414 | +0.6 % |
| 0++ 4 | 3.8996 | 2.896 | 2.883 | +0.4 % |
| 2++ 1 | 1.8765 | 1.394 | 1.393 | +0.0 % |
| 2++ 2 | 2.7419 | 2.036 | 2.028 | +0.4 % |
| 2++ 3 | 3.4330 | 2.550 | 2.535 | +0.6 % |
| 2++ 4 | 4.0369 | 2.998 | 2.980 | +0.6 % |

R00 = 1.8753 (the paper: 1.87 by its fit, 2753/1475 = 1.866), R20 = 1.3936 (the paper: 1.40).

## Checks

**Independence of the integration constants** (section 7.2 of the paper: lambda0 and A0 only set the scale).
lambda0 = 0.05 is the paper's own value (its Fig. 4 and fit).

| lambda0 | r1 | m(0++) | R00 | R20 |
| --- | --- | --- | --- | --- |
| 0.05 | 34.916 | 0.0832 | 1.8753 | 1.3936 |
| 0.08 | 6.651 | 0.4111 | 1.8753 | 1.3936 |
| 0.15 | 1.731 | 1.3465 | 1.8753 | 1.3936 |
| 0.25 | 0.879 | 2.3183 | 1.8753 | 1.3936 |

**Cutoff:** 11 r1 instead of 8 r1 leaves every mass of the table unchanged, to the last digit computed.

**Linear confinement** (alpha = 2). The local exponent d log(-A)/d log r of the background (lambda0 = 0.15):

| range | exponent |
| --- | --- |
| 2 r1 to 4 r1 | 2.094 |
| 4 r1 to 8 r1 | 2.096 |
| 8 r1 to 16 r1 | 2.053 |

It falls toward 2 slowly, through the logarithmic corrections of (8.2). The spectrum is close to linear in
m^2, as the paper's Fig. 7a shows: the gaps m_(n+1)^2 - m_n^2 of the first 12 0++ states are
4.56, 4.32, 4.51, 4.69, 4.82, 4.87, 4.95, 5.00, 5.04, 5.07, 5.10,
and of the first 8 2++ states 4.00, 4.27, 4.51, 4.68, 4.80, 4.88, 4.94:
nearly constant, and still growing slowly toward their asymptotic value, as the logarithmic corrections
predict.

**The fit of b0.** The paper fits b0 to R00 = 1.87 and finds b0 = 4.2. The secant method on R00(b0) gives
b0 = 4.147 (R00 - 1.87 = 6.5e-7). Our fit lands 1.3 % below the paper's b0,
and at its b0 = 4.2 our R00 is 1.875: the same small offset as in the table above, not a separate
difference.

**The model's own consistency:** d log W/d lambda = -(4/3) X/lambda (2.22) holds for (8.4) and (8.6).

## The earlier discrepancy (resolved, 15 Nov 2026)

This script first used a = 5/16, which is neither reading of the paper (it corresponds to alpha = 12/7 through
(8.4)). That run gave R00 = 1.683 and R20 = 1.309, 6-15 % below the paper, and the cause was not
found at the time. Rerun with a = 5/16 now: R00 = 1.683, R20 = 1.309.

With a = 3/8 the ratios agree with the paper to within 0.6 %. That remaining difference is larger than the
rounding of the paper's masses (to 1 MeV in about 1500) and is not yet explained. The next step is a
convergence study at higher precision and an independent solver.

Time: 2287 s, all runs in parallel on 16 cores (the longest, the first 12 0++ states).
