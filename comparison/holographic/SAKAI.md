# The Sakai-Sugimoto model: vector and axial-vector mesons

Produced by `sakai.erl` (decision record SY2, phase 6), PQCF 75dbdf0.
The model is Sakai and Sugimoto, [hep-th/0412141](https://arxiv.org/abs/hep-th/0412141), eq. (3.28):
-K^(1/3) d_z (K d_z psi) = lambda psi with K = 1 + z^2, and m_n^2 = lambda_n M_KK^2.

## Method

- **Coordinates:** in u = K^(1/3) the equation has polynomial coefficients,
  u (u^3 - 1) psi'' + (5u^3/2 - 1) psi' + (9u/4) lambda psi = 0, on [1, infinity).
- **Ends:** both are regular singular points. At u = 1 the exponent is 0 for the vectors (psi even in z) and
  1/2 for the axial vectors (odd). At infinity the exponent is 3/2 in w = 1/u, the decaying solution.
- **Solver:** `pq_ode:eigenvalues_sl/3` takes Frobenius series from both ends, integrates them by the
  Taylor-series method to u = 2, and finds where their Wronskian vanishes, to 25 digits.

## Results

| n | meson | psi(z) | lambda_n, PQCF | the paper (4.12) |
| --- | --- | --- | --- | --- |
| 1 | vector | even | `6.693144922442929586336821e-1` | 0.67 (rho) |
| 2 | axial | odd | `1.568765751391725077151082e0` | 1.6 (a_1(1260)) |
| 3 | vector | even | `2.874322526900746101374512e0` | 2.9 (rho(1450)) |
| 4 | axial | odd | `4.546102269155662318816202e0` | 4.5 |
| 5 | vector | even | `6.591184618626396633856336e0` | - |
| 6 | axial | odd | `9.007971092821969541329568e0` | - |
| 7 | vector | even | `1.179669410957417557764904e1` | - |

- **Mass ratios:** m_a1^2/m_rho^2 = lambda_2/lambda_1 = 2.3438 (paper ~2.4, experiment 2.51);
  m_rho(1450)^2/m_rho^2 = lambda_3/lambda_1 = 4.2944 (paper ~4.3, experiment 3.56).
- **Independent check:** a SciPy shooting in z (DOP853) with the decaying condition at a cutoff Z, extrapolated
  in 1/Z, gives lambda_1 = 0.6693144 and lambda_2 = 1.5687644. It agrees with the values above to its own
  accuracy (about 1e-7).
- **Self-check:** the result does not depend on the match point (`test/pq_ode_tests.erl`).

Times: vectors 5.4 s, axial vectors 9.3 s.
