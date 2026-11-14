# Holographic QCD: the soft and hard walls

Produced by `walls.erl` (decision record SY2, phase 3), PQCF ee18606, 2026-11-14. Each model goes from its action to its spectrum in PQCF alone:

1. **Equation of motion:** the Euler-Lagrange equation of the radial Lagrangian of a vector mode,
   (1/2) e^(-Phi) (1/z) (psi'^2 - m^2 psi^2).
2. **Schrödinger form:** the Liouville transformation gives -phi'' + V(z) phi = m^2 phi, symbolically
   (`pq_ode:schrodinger/5`).
3. **Spectrum:** the eigenvalues to 30 digits (`pq_ode:eigenvalues/3`): a Frobenius
   series at z = 0, Taylor-series integration in arbitrary precision, then the Illinois method.

## Soft wall

Karch, Katz, Son and Stephanov, [hep-ph/0602229](https://arxiv.org/abs/hep-ph/0602229). The dilaton is
Phi = z^2 (c = 1), and the decaying solution is imposed at z = 10.

- **Potential derived:** `{add,[{mul,[{rat,3,4},{pow,{sym,z},-2}]},{mul,[1,{pow,{sym,z},2}]}]}`, which is z^2 + 3/(4 z^2).
- **Prediction:** m_n^2 = 4 (n + 1).

| n | m_n^2, PQCF | 4 (n + 1) | |
| --- | --- | --- | --- |
| 0 | `4.00000000000000000000000000000e0` | `4.00000000000000000000000000000e0` | equal |
| 1 | `8.00000000000000000000000000000e0` | `8.00000000000000000000000000000e0` | equal |
| 2 | `1.20000000000000000000000000000e1` | `1.20000000000000000000000000000e1` | equal |
| 3 | `1.60000000000000000000000000000e1` | `1.60000000000000000000000000000e1` | equal |

Times: potential 0.01 s, spectrum 7.73 s.

## Hard wall

Erlich, Katz, Son and Stephanov, [hep-ph/0501128](https://arxiv.org/abs/hep-ph/0501128). There is no
dilaton, 0 < z =< z_m with z_m = 1, and psi'(z_m) = 0, which is phi'(1) + phi(1)/2 = 0.

- **Potential derived:** `{mul,[{rat,3,4},{pow,{sym,z},-2}]}`, which is 3/(4 z^2).
- **Prediction:** m_n z_m = j_(0,n), the zeros of J_0. Here they come from `pq_mp:besselj_zero/3`, which
  applies Newton's method to the Bessel series: a computation independent of the ODE solver.

| n | (m_n z_m)^2, PQCF | j_(0,n)^2 | |
| --- | --- | --- | --- |
| 1 | `5.78318596294678452117599575846e0` | `5.78318596294678452117599575846e0` | equal |
| 2 | `3.04712623436620863990781631750e1` | `3.04712623436620863990781631750e1` | equal |
| 3 | `7.48870067906951834448890413101e1` | `7.48870067906951834448890413101e1` | equal |

Times: potential 0.00 s, spectrum 3.22 s.
