# Gubser-Nellore: a black hole with the equation of state of QCD

Produced by `gubser.erl` (decision record SY2, phase 6), PQCF 3acee2e. The model is Gubser and Nellore, [arXiv:0804.0434](https://arxiv.org/abs/0804.0434), with the potential (62)
and the parameters (64): V = -12 cosh(gamma phi) + b phi^2, gamma = 0.606, b = 2.06 (L = 1).

## Method

1. **Master equation:** the nonlinear master equation (34) for G = dA/dphi, written as
   G'' = (G P/R) [a G'/P + G'^2/G^2 + G'/(6G^2) + 4G' - G'(G' + R')/P^2], with P = G + R, R = V/(3V') and
   a = G'/G + 1/(6G) - 4G - G'/P.
2. **Integration:** from the horizon series (36) at phi_H - 10^-8 down to phi = 10^-6, by `pq_jet` (the
   Taylor-series method, 20 digits). The two integrals of (45) and (51) are integrated along with it.
3. **Speed of sound:** c_s^2 = d log T / d log s by central differences in phi_H (h = 10^-3).

## Checks

- **Chamblin-Reall:** for V = -12 e^(phi/2), the master equation keeps the exact solution G = -1/(3 gamma)
  to every digit (`test/pq_holo_tests.erl`).
- **Dimension:** Delta = 3.9270, as the paper's 3.93.
- **Infrared limit:** 1/3 - gamma^2/2 = 0.1497, the paper's 0.15. The table approaches it at large phi_H.
- **Ultraviolet limit:** c_s^2 -> 1/3 at small phi_H (high T).
- **The crossover:** c_s^2 dips near T_c, then rises toward its infrared value, as in the paper's figures 3 and 4.
  The figures give no table to compare digits with.
- **An independent double-precision solution** (SciPy DOP853, the same equations) agrees to 6 digits up to
  phi_H = 4. Beyond that it is unstable: at phi_H = 6 it swings between 0.0602 and 0.32 as its tolerance
  changes, because of the 0/0 at the horizon. PQCF's 0.06021565 does not change from 20 to 30 digits or
  with h halved.

## Results (T_c: the temperature of the minimum of c_s^2 on this grid)

| phi_H | T/T_c | c_s^2 |
| --- | --- | --- |
| 0.25 | 106470961.8858 | 0.333313 |
| 0.50 | 24235.4062 | 0.333225 |
| 1.00 | 46.9137 | 0.332265 |
| 1.50 | 5.5074 | 0.327208 |
| 2.00 | 2.2754 | 0.309131 |
| 2.50 | 1.5145 | 0.264719 |
| 3.00 | 1.2409 | 0.192827 |
| 3.50 | 1.1218 | 0.119136 |
| 4.00 | 1.0630 | 0.070171 |
| 4.50 | 1.0281 | 0.048662 |
| 5.00 | 1.0000 | 0.045018 |
| 6.00 | 0.9353 | 0.060216 |
| 7.00 | 0.8462 | 0.082245 |
| 8.00 | 0.7383 | 0.101764 |
| 10.00 | 0.5155 | 0.127537 |
| 12.00 | 0.3348 | 0.140095 |
| 16.00 | 0.1295 | 0.148126 |

Time: 245.5 s.
