# The hard-wall model: vector, axial and pion observables

Produced by `hardwall.erl` (decision record SY2, phase 3), PQCF 9fcdec0. The model is Erlich, Katz, Son and Stephanov, [hep-ph/0501128](https://arxiv.org/abs/hep-ph/0501128),
with its two parameter sets (Table II) and g5^2 = 12 pi^2/N_c = 4 pi^2. Every value is computed by PQCF in
arbitrary precision, in units of 1/z_m, then converted to MeV.

## Method

- **rho:** m_rho = j_(0,1)/z_m. F_rho^2 = psi''(0)^2/g5^2 for the mode normalised by int dz psi^2/z = 1.
  The norm `int phi^2` comes from `pq_ode:frobenius_solution/4` and equals the closed form 2 J_1(j)^2/j^2.
- **a_1:** the axial equation (eq. 16) in Schrödinger form, -phi'' + [3/(4u^2) + g5^2 (m_q z_m +
  sigma z_m^3 u^2)^2] phi = (m z_m)^2 phi, with psi'(z_m) = 0. Its first eigenvalue comes from
  `pq_ode:eigenvalues/3`, and F_a1 is found like F_rho.
- **f_pi:** eq. 20 in the chiral limit m_q = 0. A(0, z) combines the Frobenius solutions of both exponents,
  with A(0) = 1 and A'(z_m) = 0. Then f_pi^2 = -(1/g5^2) A'(z)/z at z -> 0.
- **m_pi:** the Gell-Mann-Oakes-Renner relation m_pi^2 f_pi^2 = 2 m_q sigma, as in the paper.

## Results

| model | observable (MeV) | PQCF | paper |
| --- | --- | --- | --- |
| A | m_rho | `7.7675866e2` | 775.8* |
| A | m_a1 | `1.3582432e3` | 1363 |
| A | f_pi (chiral limit) | `9.0903245e1` | 92.4* |
| A | m_pi (GOR) | `1.3921146e2` | 139.6* |
| A | F_rho^(1/2) | `3.2981203e2` | 329 |
| A | F_a1^(1/2) | `4.8583344e2` | 486 |
| B | m_rho | `8.3206964e2` | 832 |
| B | m_a1 | `1.2195034e3` | 1220 |
| B | f_pi (chiral limit) | `8.2254526e1` | 84.0 |
| B | m_pi (GOR) | `1.4094354e2` | 141 |
| B | F_rho^(1/2) | `3.5329710e2` | 353 |
| B | F_a1^(1/2) | `4.4095291e2` | 440 |

## Comparison with the paper

- **Agreement:** m_rho, F_rho, F_a1 and m_pi agree with Table II to the paper's rounding, and so does model B's m_a1.
- **f_pi** is 1.6% (A) and 2.1% (B) below the paper. The paper keeps m_q: an independent floating-point
  solution of the same boundary-value problem (SciPy, `solve_ivp`) gives the chiral-limit values above
  (90.90 and 82.26 MeV). With m_q included it gives 92.1-92.3 and 83.6-83.9 MeV, depending
  logarithmically on the cutoff, which is the O(m_q^2) contact term. So the difference is the quark-mass
  contribution, not an error.
- **Model A's m_a1** is 0.35% below the paper's 1363 MeV. The same independent SciPy solution gives 1358.2432
  (A) and 1219.5034 MeV (B), equal to PQCF's. The paper's 1363 probably comes from slightly different
  parameter values (sigma = (327 MeV)^3 is rounded); this is open.
- **Not computed:** g_(rho pi pi) (eq. 22).

Time: 6.4 s for both models.
