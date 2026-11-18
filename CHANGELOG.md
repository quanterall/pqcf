# Changelog

PQCF follows semantic versioning from the API freeze (`docs/decisions/API-freeze.md`). Until 2.0, every change to
the public API is additive.

## Unreleased

### Fixed

- **`pq_ode`: Frobenius series with gaps were cut short.** The series stopped once two coefficients in a row were
  negligible; a series in z^k (k > 2) has k - 1 zeros between its terms. The hard wall's chiral-limit f_pi
  (`comparison/holographic/HARDWALL.md`) was 1e-6 low: 90.903245 and 82.254526 MeV, now 90.903337 and 82.254552,
  equal to the exact power series; m_pi (GOR) moves in its fourth decimal. Found by the new SciPy cross-check. No
  other published result changed: walls, Sakai-Sugimoto, the brane, the Weyl anomaly, Gubser-Nellore and the IHQCD
  quick check were rerun.

### Added

- **Feynman diagrams** (`pq_diagrams`, step F4, first part): every graph of a theory (vertex types as multisets of
  fields, the propagators that join them) with given external legs and loops, and its weight 1/S, computed exactly as
  the number of Wick pairings over the vertex normalisations (no automorphism group); legs of a vertex with the same
  field paired as classes; graphs identified by a canonical form from colour refinement. Filters as in QGRAF:
  connected, one_pi, no_tadpoles, no_snails.
  - **Checks** against zero-dimensional field theory (Cvitanovic, Lautrup and Pearson,
    [doi:10.1103/PhysRevD.18.1939](https://doi.org/10.1103/PhysRevD.18.1939)), its generating functions computed in the
    test in exact rationals: phi^4 vacuum graphs = ln Z to six loops (1/8, 1/12, 11/96, 17/72, 619/960; 1, 2, 4, 10,
    28 graphs), connected and one-particle-irreducible two- and four-point graphs = the cumulants of Z(J), the
    self-energy 1 - 1/G and the vertex G_c4 / G^4; phi^3 and a two-field theory (phi^2 chi) vacuum graphs = ln Z.
  - **Reach:** labelled pairings, so phi^4 to six loops in 0.3 s, phi^3 to four loops in 4 s; eight cubic vertices
    are out of reach. Next: orderly generation, fermions (signs, oriented lines), and the QED and QCD counts.
- **QCD's renormalisation group** (`pq_qcd`, step QCD1): beta_0 .. beta_4 and gamma_m 0 .. 4 exactly (rationals and
  zeta values, polynomials in n_f); alpha_s and MS-bar masses run at arbitrary precision (ln(mu1^2/mu0^2) =
  -int da/beta(a) solved by Newton's method, ln(m1/m0) = int gamma/beta, the integrals by pq_quad); alpha_s decoupled
  at an MS-bar threshold to three loops, upward by the inverse series. beta_0 .. beta_3 also for any gauge group
  (vRVL's formula in pq_lie's invariants), which for SU(3) equals the SU(3) coefficients exactly. Against mpmath (a
  separate transcription; its ODE solver at 45 digits): coefficients to 38 digits, running and decoupling to 28;
  against CRunDec 3.1.2: running to 1e-8 (its step tolerance), decoupling to 1e-14 (`comparison/qcd`).
- `pq_quad`: an integral with its limits reversed is minus the integral the other way (it failed before).
- **The one-loop four-point function D0** (`pq_oneloop:d0/11`, step F3, third part), at arbitrary precision for
  real kinematics. One analytic integration of 1/F^2 leaves singularities that are not integrable, so D0 is reduced by
  the dimension-shift relation (Bern, Dixon and Kosower, [hep-ph/9306240](https://arxiv.org/abs/hep-ph/9306240);
  Binoth, Guillet and Heinrich, [hep-ph/9911342](https://arxiv.org/abs/hep-ph/9911342)):
  D0 = -(1/2) sum_i b_i C0^(i) - (1/2) B E, b = S^-1 (1, 1, 1, 1), E the six-dimensional box (the integral of 1/F over
  the simplex, finite), computed like C0 one dimension higher: the innermost integral in closed form, the two others
  by pq_quad, split where the singular curves (quadratics in y with coefficients polynomial in x) cross and meet.
  - **Checks** (`test/oracle/gt_d0.sh`): three Euclidean configurations equal mpmath's direct 3D integral to 16
    digits; four above threshold equal LoopTools to 1e-15, two of them also QCDLoop (which returns nan for the other
    two).
  - **Speed:** 7-11 s for a Euclidean configuration at 18 digits, 45 s to 5 minutes above threshold (20 cores; the
    outer integral's nodes run in parallel). Refused: det S = 0, where the reduction does not apply.
- **Real roots of rational polynomials** (`pq_rroots`): Sturm sequences in exact rationals, bisection to any
  precision; for the singular points of loop integrands.
- `pq_quad`: `ends => 3` takes the nodes to 2^-(3P + 48) from the ends (singularities up to x^-2/3, met where two
  singular curves of the box meet; the default stays 2^-(2P + 48)), `parallel => true` for expensive integrands, and for an integral inside an integral
  `on_no_convergence => last` with `max_level` (it fails only at outer nodes whose weight makes any value
  negligible; the outer integral keeps its own test).
- **The one-loop three-point function C0** (`pq_oneloop:c0/7`, step F3, second part), at arbitrary precision for
  real kinematics, every region in one method: the y integral of the Feynman parameters in closed form (the roots of
  D as a quadratic in y with their -i eps prescription, the logarithms through the roots' products so nothing cancels
  where a root meets an end), the x integral by `pq_quad`, split where h(x) is singular. 't Hooft and Veltman's real
  substitution was tried first and dropped: it does not exist when the Kallen function lambda(p1^2, p2^2, p3^2) is
  negative (Euclidean configurations among them). Zero momenta in closed form (a divided difference of z log z);
  soft and collinear divergences refused.
  - **Checks** (`test/oracle/gt_c0.sh`, 40 configurations, 0.3-3 s each at 36 digits): mpmath's direct integral to
    28 digits where D > 0 on the whole triangle (9 configurations); mpmath's contour deformation to 1e-12 where two
    deformation parameters agree (28); elsewhere the value two of the deformation, QCDLoop 2.1.0 and LoopTools 2.16
    agree on; C0(1, 1, 1; 1, 1, 1) = -0.671253105748004686.
  - **The two external programs**, in our aarch64 builds: QCDLoop is wrong in 6 of the 40 configurations (and off
    by 7e-9 in one more), LoopTools in 7; each of them with lambda < 0, where PQCF equals the deformation to 1e-14 or
    better. Below threshold with lambda < 0 LoopTools is wrong as well (C0(1, 1, 1; 1, 1, 1) = 6.87 + 34.2 i).
    Details: `comparison/oneloop/README.md`.
- **Quadrature at arbitrary precision** (`pq_quad`): tanh-sinh (Takahasi and Mori), nodes placed from the ends so
  logarithmic and inverse-square-root end singularities cost nothing, breakpoints for interior singular points,
  complex integrands.
- **One-loop two-point functions** (`pq_oneloop`, step F3, first part): A0, B0, B1, B00, B11 and dB0/dp^2 at
  arbitrary precision for real kinematics, in LoopTools' conventions; B0 in closed form in every region (space-like,
  between the thresholds through an angle, above threshold with Im B0 = pi sqrt(lambda)/p^2, massless lines, p^2 = 0),
  the rest by Passarino-Veltman reduction. 60 configurations equal mpmath's Feynman-parameter integrals to 35 digits
  and LoopTools 2.16 to 1e-11 (`comparison/oneloop/`, `test/oracle/gt_oneloop.sh`). `pq_mp:atan/2` added (the header
  listed it already). (This entry said first that LoopTools' default C0 is wrong for a negative Kallen function and
  its Denner version right; corrected below under the three-point function: no version of our build is right there.)
- **Dirac algebra in d dimensions** (`pq_dirac`, step F2, first part): traces of gamma matrices and slashed momenta
  natively. Repeated indices are removed first by the identity gamma^mu X S gamma_mu = 2 S X - X (gamma^mu S gamma_mu)
  (gamma^mu gamma_mu = d), then the rest by the pairing recursion, memoised: a trace of 18 matrices with six
  contracted indices takes 0.2 s. 60 fully contracted traces of up to 12 matrices equal FORM's `tracen` exactly
  (`test/oracle/gt_dirac.sh`).
- **Spectra on the half-line** (step QM2, first part): `pq_ode:eigenvalues/3` takes a Neumann left end
  (`{neumann, A}`, psi'(A) = 0), for the even states of symmetric potentials. The quartic oscillator H = p^2 + x^4: its
  first four eigenvalues (1.0603620904841828996470460166926635455...) equal an independent mpmath shooting
  (`test/oracle/gt_quartic.py`) to 38 digits.
  The Frobenius start at 0 also takes a simple pole of V (a regular singular point: N and D are multiplied by z), so
  hydrogen's s-states are solved like the others: every level n =< 3 for l = 0, 1, 2 equals -1/n^2 to 30 digits.
- **Perturbation theory to high order** (`pq_perturb`, step QM1): the ground state of H = p^2 + x^2 + g V(x), V a
  polynomial, as an exact power series in g, by Bender and Wu's recursion (rational arithmetic; 100 orders in about
  5 s). The anharmonic oscillator's coefficients equal Bender and Wu's table
  ([Phys. Rev. 184, 1231](https://doi.org/10.1103/PhysRev.184.1231)) exactly, and their large-order formula holds to
  2e-4 at order 100 (with its first correction, the remainder falling like 1/n^2); V = x^2 gives the binomial series of
  sqrt(1 + g) exactly.
- **Lie algebras for physics** (`pq_lie`; plan "theoretical physics instruments", step F1).
  - **Wigner symbols exactly:** 3j, 6j and Clebsch-Gordan coefficients as signed square roots of rationals, 9j symbols
    as sums over 6j products (Racah's formulas, half-integers as doubled integers). Equal to SymPy's on every 3j and 6j
    with j =< 3, a sample up to j = 10 (5,402 in all, exactly) and 150 9j symbols (to 40 digits).
  - **SU(N) colour algebra with N symbolic:** products of traces of generators reduced by the Fierz identity to
    Laurent polynomials in N; f and d as traces. T_F, C_F, C_A, d33 and the quartic Casimirs of the four-loop beta
    function (d_F d_F, d_F d_A, d_A d_A over N_A) equal the closed forms of van Ritbergen, Schellekens and Vermaseren
    ([arXiv:hep-ph/9802376](https://arxiv.org/abs/hep-ph/9802376)) exactly; d_A d_A takes 0.5 s.
- **Input-output networks** (`pq_slh`, quantum optics Q6).
  - **What:** SLH triples (Gough and James, [arXiv:0707.0048](https://arxiv.org/abs/0707.0048); Combes, Kerckhoff,
    Sarovar, [arXiv:1611.00375](https://arxiv.org/abs/1611.00375)) with pq_qop operators: concatenation, series,
    feedback reduction, static elements and coherent inputs; the reduced network's master equation.
  - **Checks:** two cavities in cascade reduce exactly to Gardiner's and Carmichael's cascaded system; a cavity's output
    fed back through a phase gives L = (sqrt(k2) + e^(i phi) sqrt(k1)) a and the shift sqrt(k1 k2) sin phi (three
    phases, exactly); a three-mode mesh's S equals PQCF's own symbolic circuit unitary entry by entry; a coherent drive
    cascaded through two cavities relaxes to <a> = -2 alpha/sqrt(ka), <b> = 2 alpha/sqrt(kb).
- **Collective systems in the Dicke basis** (`pq_dicke`, quantum optics Q5).
  - **What:** N identical two-level atoms with collective and local emission, dephasing and pumping (PIQS's
    coefficients, Shammah et al., [arXiv:1805.05129](https://arxiv.org/abs/1805.05129)) and H = omega J_z +
    omega_x J_x, in arbitrary precision; only the populations when nothing mixes m - m' (about N^2/4 numbers).
  - **Checks:** the superradiant burst of 100 atoms equals the exact population equations (mpmath, 40 digits) to 34
    digits; its peak at ln N/(N gamma) as Gross and Haroche give; QuTiP's PIQS on three problems within 1e-8, PIQS's
    own precision (its rates are single-precision floats, `comparison/dicke`).
- **Cumulant expansions** (`pq_cumulant`, quantum optics Q4).
  - **Equations for averages** from a Hamiltonian and jump operators (`pq_qop:heisenberg/4`), closed at order K by
    setting the cumulants above K to zero (the sum over partitions with blocks of at most K operators), completed until
    the set closes; an average and its adjoint's share one variable. Solved by `pq_jet` in real and imaginary parts.
  - **Checks:** the driven damped cavity exactly (order 2); against QuantumCumulants.jl 0.3.7
    (`comparison/cumulants`): the single-atom laser (8 equations, the same <a^dagger a> and P_e to 1e-13) and the
    driven Jaynes-Cummings model at orders 2, 3, 4 (8, 15, 24 equations, the same <a^dagger a>(2) to 1e-15); against
    the exact master equation, order 4 within 6e-6.
- **Exact steady states and spectra** (`pq_qo_exact`, quantum optics Q3).
  - **Exact Liouvillian:** pq_qop operators with rational parameters, in Gaussian rationals. In the basis
    f_n = sqrt(n!) |n> the boson operators are integer matrices, so the master equation of a rational model is rational;
    the jump terms take the products of the coefficients symbolically, under the declared domains.
  - **Steady states** by exact sparse elimination; **spectra** by the quantum regression theorem, S(w) at rational w,
    exactly.
  - **Checks:** equal, as rationals, to SymPy's dense exact solutions in the Fock basis (`test/oracle/gt_qo_exact.py`)
    for resonance fluorescence (steady state and spectrum, five drives and detunings), a driven Kerr cavity and the
    driven, damped Jaynes-Cummings model; the Mollow triplet's strong-field limit (central peak three times the side
    peaks, at +-Omega); g2(tau) of resonance fluorescence against Kimble and Mandel's closed form to 35 digits, and
    g2(0) = 0.
  - **The driven Kerr cavity, exactly:** `pq_qo_exact:kerr_moment/3`, Drummond and Walls's complex-P solution
    ([J. Phys. A 13, 725](https://doi.org/10.1088/0305-4470/13/2/034)) with its Gamma functions as finite products,
    summed to any precision. Across the semiclassically bistable region (Delta = -2, chi = 1/2, kappa = 1, four
    drives), <a^dagger a> and <a^dagger^2 a^2> equal the master equation's steady state at 40 photons to 30 digits
    (to 60 at the first drive).
  - **Steady states as rational functions of the parameters:** `pq_qo_exact:steady_state_symbolic/6`. rho' in real
    coordinates (Hermitian), the master equation's real and imaginary parts with the trace replacing its dependent
    row, solved by finite-field reconstruction (`pq_recon:linear/4`, decision record R1): no symbolic elimination.
    Resonance fluorescence with Omega, Delta and gamma symbolic gives the optical Bloch steady state exactly
    (rho_ee = Omega^2/(4 Delta^2 + gamma^2 + 2 Omega^2), and <sigma_x>, <sigma_y>) in 0.05 s; the Kerr cavity with all
    four parameters symbolic equals the exact numeric steady state at every point tried.
  - **One- and two-photon driving and dissipation, exactly:** `pq_qo_exact:kerr2_moment/4`, Bartolo, Minganti,
    Casteels and Ciuti's complex-P solution ([arXiv:1607.06739](https://arxiv.org/abs/1607.06739)): the Kerr resonator
    with coherent and parametric drives and one- and two-photon loss, its terminating 2F1 functions summed in
    arbitrary precision. With U = 0 and F = 0 it is the degenerate parametric oscillator with its pump eliminated.
    <a^dagger a> and <a^2> equal the master equation at 40 photons to 30 digits in three cases, that one included.
  - `pq_lindblad`: zero to a positive real power is 0 (it took the logarithm of 0).
- **Master equations** (`pq_lindblad`, quantum optics Q2).
  - **Matrices:** pq_qop operators in a truncated basis (Fock states of each mode, levels of each atom), the
    coefficients evaluated at given parameters; the Liouvillian as a sparse superoperator.
  - **Time evolution:** the Schrödinger and master equations by the Taylor series to the working precision, in
    arbitrary precision (pq_mp); expectation values.
  - **Steady states** by sparse Gaussian elimination; **quantum trajectories** (Monte Carlo wave functions,
    Dalibard, Castin and Mølmer, [PRL 68, 580](https://doi.org/10.1103/PhysRevLett.68.580)).
  - **Checks, to 100-200 bits:** vacuum Rabi oscillations cos^2(g t); Jaynes-Cummings collapse and revival against
    the closed-form sum of Eberly, Narozhny and Sanchez-Mondragon ([PRL 44, 1323](https://doi.org/10.1103/PhysRevLett.44.1323));
    a damped coherent state, alpha e^(-kappa t/2); atomic decay e^(-gamma t); the steady state of resonance
    fluorescence (Mollow) and of the driven damped cavity; the average of 400 trajectories against the master
    equation, within four standard errors.
  - **Speed:** the Taylor terms in block fixed point (exact integers at one scale, 24 guard bits) and steps of
    8/|M|: the same 36 digits as pq_mp throughout, 40-100 times faster.
  - **Against QuTiP** (`comparison/quantum-optics`, a pinned QuTiP 5.3.1 image): the same problems on the same core,
    digits against the closed forms and times; `run.sh` writes RESULTS.md, `results/` and history.csv.
- **Operator algebra for quantum optics** (`pq_qop`, the first step of quantum optics in PQCF).
  - **Objects:** bosonic modes and N-level atoms (sigma^-, sigma^+, sigma_z for two levels), with exact coefficients.
  - **Canonical form:** normal-ordered within each mode, subsystems commuting, the ground projector eliminated by
    completeness; equal operators have equal forms.
  - **Operations:** products, powers, commutators, adjoints, and the Heisenberg-Lindblad equation of any operator for
    a Hamiltonian and jump operators. Symbol domains (positive, real) are declared with the spaces.
  - **Checks:** identical to SymPy's normal ordering on 60 random words of two modes (`test/oracle/gt_qop.py`);
    (a + a^dagger)^4, [a, (a^dagger)^n] = n (a^dagger)^(n-1), the BCH shift e^(-alpha a^dagger) a e^(alpha a^dagger)
    = a + alpha, the two- and three-level atom identities, and the Jaynes-Cummings equations with cavity loss.
- **Multivariate polynomial gcd** (decision record SY2, which lifts the scope fence of SY 1.2):
  - `pq_symbolic:gcd/2,3` on expressions; `pq_sym_gcd:gcd/2` and `divide/2` on polynomials over the rationals;
    `pq_form:gcd/2,3` for comparison.
  - The gcd is primitive, with integer coefficients and its leading term positive.
  - **The modular method:**
    - the gcd modulo primes near 2^63 in Rust (`mgcd.rs`): monomial contents first, a univariate coprimality
      proof, Zippel's sparse interpolation, and Brown's dense algorithm as a fallback, in Montgomery arithmetic;
    - Chinese remaindering in Erlang;
    - each candidate accepted only when it divides both inputs exactly.
  - Square-free decomposition, `pq_sym_gcd:squarefree/1` (Yun's algorithm, contents recursively).
- **Holographic renormalisation, first result** (SY2 phase 5):
  - **What:** the AdS_5 black brane's free energy from its renormalised on-shell action. Bulk, Gibbons-Hawking and
    counterterm, the eps^-4 divergences cancelling.
  - **Result:** F = -L^3/(16 pi G z_h^4) = -pi^2 N^2 T^4/8 (`comparison/holographic/BRANE.md`,
    `test/pq_holo_tests.erl`).
  - **New:** `pq_symbolic:integrate/2` for Laurent polynomials, the first case of integration.
- **The holographic Weyl anomaly** (SY2 phase 5): a = c = N^2/4 (`comparison/holographic/WEYL.md`,
  `test/pq_holo_tests.erl`).
  - **Method:** a Fefferman-Graham expansion on S^4 and on S^2 x S^2, its coefficients (g_(2), Tr g_(4)) derived
    from the 5-dimensional Einstein equations, then the volume coefficient a_4.
- **Abstract-index tensor monomials** (SY2 phase 4): `pq_tensor`.
  - **Canonical forms:** under slot symmetries (Riemann, Ricci, metric), the commutation of equal tensors and dummy
    renaming; zero detection.
  - **Checked against SymPy's Butler-Portugal `canon_bp`:** the same classes, zeros and signs on 160 random
    monomials (`test/oracle/gt_tensor.py`).
  - **Invariants:** the full contractions of k Riemann tensors, and their number of independent invariants by exact
    rank on random algebraic curvature tensors: 3 quadratic and 8 cubic (Fulling, King, Wybourne, Cummins 1992).
- **Curvature in components** (SY2 phase 4): `pq_geometry:curvature/2`.
  - **Output:** the inverse metric, Christoffel symbols, Riemann, Ricci, scalar curvature, Einstein, a Weyl test and
    the Kretschmann scalar, in any dimension.
  - **Components:** exact rational functions in normal form.
  - **Checks:**
    - AdS_5: R = -20/L^2, R_ab = -(4/L^2) g_ab, Weyl 0, Kretschmann 40/L^4;
    - the black brane solves G_ab + Lambda g_ab = 0 with Lambda = -6/L^2, has Kretschmann
      (40 + 72 z^8/z_h^8)/L^4, and T = 1/(pi z_h);
    - Schwarzschild: Ricci flat, Kretschmann 48 M^2/r^6.
- **Equations of motion and spectra** (SY2 phase 3): `pq_ode`.
  - `euler_lagrange/5` and `schrodinger/5`: from a mode's Lagrangian to -phi'' + V phi = lambda phi, symbolically.
  - `eigenvalues/3`: the spectrum to any number of digits, for V a rational function. It uses a Frobenius start at a
    singular end, Taylor-series integration on the polynomial form of the equation, Dirichlet or Robin conditions,
    and the Illinois method.
  - `eigenvalues_sl/3`: Sturm-Liouville problems with polynomial coefficients and two regular singular ends, one of
    them at infinity (Frobenius from both ends, Taylor integration, the Wronskian at a match point).
    - **Sakai-Sugimoto** (`comparison/holographic/SAKAI.md`): lambda_n = 0.66931449, 1.56876575, 2.87432253,
      4.54610227, ... to 25 digits.
    - **Check:** against the paper's 0.67, 1.6, 2.9, 4.5, and an independent SciPy solution to its 1e-7.
  - `frobenius_solution/4`: the Frobenius solution of either exponent at the right end, with int phi^2 (eigenfunction
    norms).
  - **The hard wall's axial and pion sector** (`comparison/holographic/HARDWALL.md`): m_a1, F_rho, F_a1, f_pi (chiral
    limit) and m_pi (GOR) for the paper's two parameter sets.
    - **Agreement:** with Table II of hep-ph/0501128 to its rounding, apart from f_pi, which is 1.6-2.1% lower: the
      quark-mass term, confirmed by an independent SciPy solution.
    - **Open:** model A's m_a1, 0.35% lower than the paper's but equal to the independent solution.
  - `eigenvalues_collocation/3`: the same spectrum by an independent method, Chebyshev collocation of phi / z^s
    with shifted inverse iteration.
  - **The holographic QCD walls from the action to the spectrum** (`comparison/holographic/RESULTS.md`), equal to
    30 digits:
    - soft wall: m^2 = 4(n+1);
    - hard wall: the zeros of J_0, computed independently.
- **Arbitrary-precision evaluation** (SY2 phase 2): `pq_symbolic:evaluate_digits/3`, `pq_mp`.
  - **Numbers:** binary floating point on Erlang's integers.
  - **Constants:** pi, log 2, Euler's gamma (Brent-McMillan), zeta(k) (Borwein).
  - **Functions:** exp, log, sin, cos, powers; Gamma and polygamma (Stirling with exact Bernoulli numbers); Bessel J
    and Y; Laguerre; 2F1.
  - **Bessel zeros:** `pq_mp:besselj_zero/3`.
  - **Checked against mpmath 1.3:** 37 values and 8 zeros to 48 digits (`test/oracle/gt_mp.py`).
  - **Integer relations:** `pq_mp:pslq/3` (Ferguson-Bailey PSLQ) finds exact forms of high-precision numbers. It
    recovers 6 zeta(2) = pi^2 and Y_0's series constant (2/pi)(gamma - log 2).
- **The Taylor-series method for nonlinear ODE systems** (SY2 phase 6): `pq_jet:integrate/6,7`.
  - **What it does:** solves y' = F(x, y), with F built from +, *, powers, exp, log, sin, cos, sinh, cosh, sqrt, by
    automatic differentiation on truncated series, in arbitrary precision.
  - **Improved holographic QCD** (arXiv:0707.1349, `comparison/holographic/IHQCD.md`): background I and the 0++ and
    2++ glueballs.
    - **Checks:** independent of the integration constants (the paper's section 7.2) and of the cutoff.
    - **Results:** R00 = 1.875 and R20 = 1.394 at b0 = 4.2, against the paper's 1.866 and 1.393; every ratio of the
      first four 0++ and 2++ states within 0.6 %. The first run used a = 5/16 in X and W, where the paper's (8.4) and
      (8.6) need a = 3/8 for linear confinement, and was 6-15 % low.
  - `pq_jet:diff/2`: derivatives of right-hand sides with log.
  - **Gubser-Nellore** (arXiv:0804.0434, `comparison/holographic/GUBSER.md`): the nonlinear master equation from the
    horizon series, then the entropy, the temperature and the speed of sound.
    - **Checks:** the exact Chamblin-Reall solution is kept to 25 digits; Delta = 3.927 (paper 3.93); c_s^2 goes from
      1/3 in the UV through a dip at the crossover toward 1/3 - gamma^2/2 = 0.1497 (paper 0.15).
    - **Against SciPy:** a double-precision solution agrees to 6 digits where it is stable.
- **Truncated series and limits** (SY2 phase 2): `pq_symbolic:series/3,4` and `limit/3,4`, `pq_series`.
  - **Terms:** C x^e log(x)^m, e rational, at 0, at a point or at infinity.
  - **Coefficients:** exact rational functions of the other symbols. Constants that are not rational (log 2, pi,
    Euler's gamma, zeta(k), Gamma(q), ...) are kept as opaque atoms.
  - **Functions:** sums, products, powers (integer, rational, symbolic), exp, log, sin, cos, sqrt; Bessel J_nu (nu
    rational) and Y_n (with its logarithms); Laguerre L_n^alpha; 2F1; Gamma at any rational point.
  - **Precision:** the working order rises until cancellations are covered.
  - **Checked against SymPy 1.14 and mpmath:** 32 expansions (`test/oracle/gt_series.py`,
    `test/fixtures/gt_series.json`), every coefficient to 1e-10.
- **Rational functions in normal form** (SY2): `pq_symbolic:cancel/1,2` (N D^-1, N and D coprime, D primitive with a
  positive leading coefficient: equal rational functions give equal results) and `numer_denom/1,2`;
  `pq_sym_ratfun`.
- **Changed:** `pq_symbolic:equiv/3` decides rational functions of symbols by their normal forms, true or false,
  where it answered `unknown` before (the scope fence of SY 10.3, lifted by SY2).
- **Factorisation over the rationals** (SY2): `pq_symbolic:factor/1,2` -> `{ok, {Constant, [{Factor, Multiplicity}]}}`,
  `pq_sym_factor`, `pq_form:factor/1,2` for comparison.
  - **One variable:** Zassenhaus's method. Factors modulo a small prime in Rust (`pq_nif:factor_mod/3`:
    distinct-degree by the Frobenius matrix, then Cantor-Zassenhaus), quadratic Hensel lifting on a factor tree,
    recombination by trial division.
  - **Several variables:** Wang's method.
    - Every variable is tried as the main one: a content splits the polynomial, and an irreducible image proves
      irreducibility.
    - The leading coefficient is imposed on each factor.
    - Lifting runs in Rust modulo 63-bit primes (`pq_nif:mfactor_lift/6`, in powers of (x_j - a_j)), followed by
      Chinese remaindering.
    - A two-factor subset search is the fallback.
  - Every factorisation is checked by multiplying back.
  - **FORM comparison:** new task T8, polybench's factorisation problems (5 each). Results equal to FORM's, factor
    for factor. On core 8: factor(a b) 1.02 s for PQCF against 0.13 s for FORM 5; factor(a b + c) 0.077 s against
    0.019 s.
  - `pq_simplify`'s gcd now uses it; the old pseudo-remainder gcd stays as `gcd_prs/2`, the reference in the tests.
  - **FORM comparison:** new task T7, polybench's defaults (T. Ueda; 5 variables, 10 problems each), results equal
    to FORM's. On core 8, 10 nontrivial gcds take PQCF 0.083 s and FORM 5 with FLINT 0.075 s; 10 with a monomial
    gcd take 0.052 s and 0.044 s.
- **The task behaviour** (`pq_task`, decision record P1): PQCF's executors run any application's work.
  - A task is `{pq_task, Module, Args}`; `Module` implements `run/2`, and optionally `resources/1` and `describe/1`.
  - A node runs only the modules it allows (PQCF's own, the `task_modules` environment, `pq_task:allow/1`).
  - `pq_rexec:map_traced/3` returns each result with its provenance record: task hash, module and code checksum,
    node, version, times.
  - The finite-field reconstruction's tasks are now `pq_ff_task`, with their bare tuples still accepted; results
    unchanged.
- Feynman-integral families with **massive propagators** and **any number of loops** (`pq_ibp_family`):
  - a propagator may be `{Sign, Momentum, MassSquared}`, D = Sign (q^2 - m^2), the mass squared an expression in
    the family's variables; the masses enter the scalar products' constants and the Symanzik polynomial F (zero
    sectors and symmetries see them);
  - U and F by determinant and adjugate for any number of loops (the two-loop special case is gone);
  - derivatives in a variable of the masses (the propagator's index raised), and vacuum families (no external
    momenta);
  - tests: the massive tadpole's I(2) = (d - 2)/(2 m^2) I(1) and its differential equation in m^2, a bubble with
    one massive line (zero sectors), and the three-loop massless banana (four lines, five numerators), whose
    differential equation in s is (3d/2 - 4)/s as dimensional analysis requires.
- **Sector symmetries with numerators** (`pq_ibp_family:expand/2`, `mapped/2`): for a sector equivalent to its
  representative, a change of the loop momenta k -> A k + B p (A integer, det A = +-1) carrying its propagators onto
  the representative's, found from the permutation (Pak) and an exact linear solve over the sign choices. Every
  integral of the sector, numerators included, then expands into integrals of the representative and its
  subsectors. Such sectors get no identities of their own, and a target in one is defined by its expansion. The
  kite's identities 2,274 -> 1,579; the double box 7,296 -> 6,807 and 18,366 -> 16,967, answers unchanged (equal to
  FIRE's). The maps may also move the external momenta (p1 <-> p2 and the like), when the move preserves every
  scalar product of the kinematics: the double box has 31 sector maps, its example target needs 3,018 equations
  (0.27 s on one core, FIRE 8.0 s). A sector whose dotted integrals expand into ones its representative's
  identities do not reduce is handled by the `sector_maps => false` option (every sector's identities, scalar
  symmetries only), which the automatic seeds try after `true` at each budget: I(1,1,1,1,1,1,1,-2,-1) takes
  1.28 s on one core this way (18,366 equations; FIRE 12.0 s).
- **Comparison with FORM's Forcer on FORM's own four-loop checks** (`comparison/forcer`; Forcer:
  [arXiv:1704.06650](https://arxiv.org/abs/1704.06650)): PQCF reduces the integrals of `forcer.frm`
  (`Forcer_example`, `Forcer_1`, `Forcer_2`: the four-loop non-planar topologies no1-no3 with dots and numerators,
  and a two-loop kite with insertions), and every answer is exact: equal to Forcer's results as epsilon expansions,
  and equal in MS-bar to the numbers FORM's checks assert (zeta_3 to zeta_7). Forcer is 17-25x faster on the
  four-loop checks (`Forcer_1`: 10.3 s against 175 s, from 918 s before the changes below); PQCF 13x faster on the
  three-loop example.
- **Reduction overheads:** the seed search's completeness check is one solve on a random combination of the
  targets (four loops: 37.7 s -> 5.2 s); `reduce_sum/3` reconstructs only the masters found non-zero;
  `pq_recon:linear/4` stores its program once per call and runs univariate columns concurrently. The seed search
  takes the global budgets first while they stay within 3x the smallest, then the sector-by-sector search: the
  double box's harder target is back to 1.4 s (3.0 s with the sector search first).
- **Seeds chosen sector by sector** (`pq_ibp:auto_seeds/3`): each master that is not a seed raises only its own
  sector's budget (`pq_ibp_family:sector_seeds/4`); the global budgets stay as the fallback, `sector_seeds =>
  false` restores them. Four loops: 374,416 equations -> 97,949 for one case.
- **`pq_ibp:reduce_sum/3`**: a combination of targets reduced as one result (one extra unknown and equation, one
  reconstruction); **`pq_recon:linear/4`** leaves out the zero right-hand sides (19,841 -> 21 on a four-loop
  system).
- **`pq_ibp_family`: insertions and numerators** -- `eps => #{Line => M}` (a line's power index + M epsilon, as a
  one-loop insertion leaves it) and `numerator/3` (an integral times scalar products, as integrals of the family).
- **Four-loop symmetry classes 75 s -> 0.7 s**: sector maps found by backtracking over lines with equal invariants
  instead of every permutation (11! for a four-loop top sector).
- **A reduction benchmark suite** (`bench/ibp.sh`, `bench/ibp_suite.erl`, `bench/fire_run.sh`): PQCF and FIRE
  6.5.2 on FIRE's own examples (the one-loop box, the two-loop vertex v2, the double box) and harder double-box
  targets, on the same core back to back. Every answer is checked equal to FIRE's (`fire_check:compare/5`,
  `tables/1`), and each run is compared with the previous one. PQCF-only cases: the kite, the massive tadpole and
  the three-loop banana. The first run: PQCF at 0.03–0.78 of FIRE's time on one core, every answer equal.
- **Smaller results** (decision record R1, "Smaller results"):
  - `pq_ibp:choose_masters/3`: a greedy search over master bases (each master swapped for a neighbour in its
    sector, kept when the simplified coefficients shrink), every candidate an exact reduction on the default
    reduction's seeds; `reduce/3` takes `masters => [Integral]` (preferred masters) and `pq_ibp:size_of/2` measures
    a reduction. The double box: 18 -> 12, 56 -> 32 and 102 -> 69 additions (12.7 s, 53 s and 295 s of search on
    16 cores); the expanded coefficients 25 -> 16, 79 -> 48 and 148 -> 103.
  - `pq_simplify`: multivariate partial fractions (Leinartas's decomposition for affine denominator factors, then
    the numerators reduced modulo the denominators; nested: in d, then in the kinematics) as a fourth form, chosen
    when smaller or forced by `form => partial`; `apart => [[Var]]` the variable sets tried.
  - `pq_simplify:discover/2` and `gcd_poly/2`: factors outside the alphabet, repeated in a coefficient or shared
    by several, found by exact multivariate gcds (primitive pseudo-remainder sequences) and a coprime refinement,
    to pass as `extra`.
- **Large systems: Laporta's elimination equation by equation** (`pq_nif:ff_pivots/4`): the equations from the
  simplest to the most complex, each reduced only while its leading integral already has a pivot, its other entries
  left as they are, the redundant ones dropped at once. The pivot columns (the masters) are the same; mark and
  sweep is unchanged. The column-by-column elimination it replaces carried every row along: on a double-box system
  of 194,739 equations (I(1,1,1,1,1,1,2,-2,-2): seeds with a dot and four numerators), 8.8 billion merged entries
  against a final 654,082, 31 s; now 3.2 s, and the target needs 18,543 equations instead of 67,663. That reduction:
  10.7 s on one core (FIRE 6.5.2: 12.8 s on one thread), 9.6 s on 20 (FIRE: 6.9 s on 20 threads), peak memory
  1.9-3.7 GB; FIRE's example targets 0.31 s and 0.79 s on one core (FIRE 8.0 s and 12.2 s), answers equal to FIRE's.
- The systems' coefficients are evaluated mod p in the processes that generate them (one copy of the equations
  instead of two): the large double-box system 6.45 -> 5.28 s on 20 cores, its whole reduction about 7.1 s (FIRE: 6.9
  s on 20 threads). A batched parallel learning elimination (equations reduced in parallel against earlier pivots,
  then inserted in order: exact, the same result as sequential) was measured and dropped: no faster, since
  neighbouring equations reduce mostly against each other's pivots.
- Fixed: `pq_ibp_family:canon/2` mapped integrals even in a family with `symmetries => false`; the automatic seeds
  could then never be complete for such a family. Found by the sector-map test, which compares mapped reductions
  with unmapped ones, master by master.

### Documentation

- After the external reviews of 15 Nov 2026:
  - `REPRODUCE.md`: every published result (PQCF's, GenHub2's and quantum_computing_simulations'), its step-by-step
    page, the commit it was tested at, its time and what it needs;
  - every step-by-step page has two paths: the published run, at the commit it was tested at, and the latest.
  - every Dockerfile names its base image by digest; the digests are in REPRODUCE.md.
  - the SciPy cross-checks quoted in SAKAI.md, HARDWALL.md and GUBSER.md as committed scripts with tolerances
    (`comparison/holographic/scipy_*.py`, in the oracle image);
  - `REPRODUCE.md`: how the answers of every comparison are compared, exactly, to how many digits or statistically,
    and why.
  - CI in four groups on both targets (`tools/ci.sh`: static, fast, integration, scientific), each its own job; the
    README shows its status. The x86 failures since 14 Nov were one test past EUnit's default 5 s timeout.
  - the README opens with PQCF's two capabilities, the scientific engine and distributed execution, and an
    architecture diagram (`docs/architecture/`); the mathematical physics and the task layer are in "What is done".
- After an external review of the repository (8 Nov 2026):
  - the user guide's introduction said PQCF was not ready and its API not frozen (stale since 1 Nov); it now
    states the freeze and the release;
  - a new guide section, 17, on reducing results: a photonic quantity (the permanent of a symbolic interferometer
    block) reconstructed over finite fields, a Feynman-integral reduction, a differential equation; every example
    is a test;
  - the README says that PQCF is not post-quantum cryptography;
  - "One computational foundation": the four layers (the symbolic engine, compiled programs, reconstruction over
    finite fields, execution on the BEAM) that the photonic simulator and the symbolic work share, what each serves
    on either side, and how each is validated;
  - "Next" brought up to date.

### Changed

- **The symbolic engine, faster where FORM was ahead** (found by profiling the FORM comparison, 14 Nov 2026):
  - a result kept in Rust comes back as its handle, its node count taken in Rust (`pq_nif:rpoly_nodes/2`, exactly
    `pq_symbolic:size/1`, so `max_expr_nodes` refuses at the same count, before building), and is built in the caller
    with a faster `rpoly_to_ast`;
  - powers of a base of up to 64 terms multiply by merging instead of hashing;
  - the permanent of a matrix of symbols with integer coefficients runs in Rust (`pq_sym_poly:kernel_permanent/3`).

  Results unchanged. Old against new on one core, same load: powers 1.6-2.2x faster, the sparse product 1.6-1.9x,
  permanents 2.9-5.3x (now faster than FORM), the 8-mode mesh 1.5x.
- Programs are evaluated over finite fields in Montgomery arithmetic (`Mont` in `ffslp.rs`: three 64-bit products per
  multiplication, no 128-bit division, which is a library call on aarch64): a double-box sample 2.86 -> 2.14 ms,
  the same values.
- `pq_ibp:system/4` builds its systems from the identities' coefficients as exact polynomials
  (`pq_ibp_family:identities_poly/2`, the zero test and symmetry name of each integral computed once per call),
  evaluated mod p directly, only the needed rows turned into ASTs; Shoup multiplication in the search-based
  eliminations: the harder double-box target's seeds and system 1.12 -> 0.63 s, the whole reduction on one worker
  2.04 -> 1.63 s, the same answers (equal to FIRE's).
- Reductions use the cores: the systems' identities and their values are generated in parallel chunks (the
  family's caches computed once and handed to the chunks, `pq_ibp_family:warm/2`, `load/1`); the search-based
  eliminations update their column index only where a row changes (fill-in, cancellations), which made the
  learning elimination about 4x faster; a system's solved points are cells that concurrent requests wait on, and
  its learned elimination is shared by its black boxes; the guessing slices run as parallel tasks. Two double-box
  targets with four numerators each: 9.4 s on one core before, now 6.3 s on one and 4.4 s on 20 (the harder single
  target 1.14 s on 20); answers unchanged (equal to FIRE's). Reconstructing the masters' coefficients concurrently
  was measured slower than one after another (each then repeats the others' setup, and points are shared anyway),
  so they stay sequential.

## 1.4.0 (8 November 2026)

Reducing large symbolic results, additive under the API freeze (decision record `R1`). Results are reconstructed
exactly from evaluations modulo primes, spread over local workers or Erlang nodes with the same answer on any of
them. Feynman integrals reduce to master integrals through `pq_ibp`; on FIRE's own two-loop double box, PQCF is
3-7x faster than FIRE 6.5.2 on one core (1.21 s against 7.99, 3.73 s against 11.99), with answers checked equal
master by master. Results are simplified by verified rewriting (`pq_simplify`), and master-integral differential
equations are derived exactly.

### Added

- **Feynman-integral reduction as part of the library** (`pq_ibp`, `pq_ibp_family`; R1):
  - `pq_ibp_family:new/1`: a family as data (loop and external momenta, propagators with numerators, external
    scalar products as a map, variables, top sector), validated; `example/1` gives the benchmarks' box, kite and
    FIRE's double box; identities, scaleless sectors (Lee-Pak), sector symmetries (Pak), seeds and derivatives with
    respect to the invariants, cached per process by the family's hash, random points from local states;
  - `pq_ibp:reduce/3`: several targets at once, seeds automatic (`seeds => auto`, the default) or given
    (`{dots, D, S}`, Kira's `{R, S}`, a list), the masters' order (`order => numerators | dots`), pq_recon's options
    (executors, denominator guessing); a zero target maps to `[]`, a master to itself; `system/4`, `seeds/3` and
    `auto_seeds/3` for the system itself; `differential_equation/4`: the masters' equation in an invariant;
  - the double box's two targets in one call: 2.7 s on one core; the benchmark drivers (`comparison/form/ibp.erl`,
    `comparison/pinn/box_de.erl`, `comparison/fire/fire_check.erl`) now use the library, with the same results;

- **Reducing results by finite-field reconstruction, phase 1** (`R1`; the project owner's decision, 7 Nov 2026):
  - `pq_symbolic:compile/3` with `exact_constants => true`: constants known exactly stay exact in the program
    (float evaluation unchanged);
  - `pq_ff`: compiled programs evaluated modulo primes below 2^63 in Rust (`ffslp.rs`), a batch of points at a
    time, a division by zero marking a point as bad; deterministic primes (Miller-Rabin); inverses, Chinese
    remaindering and Wang's rational reconstruction;
  - `pq_recon:polynomial/3`: a polynomial reconstructed from its program's values (Zippel over the first prime,
    the support reused in the next, Chinese remaindering and rational reconstruction, a final check); the resultants
    res(7,4) and res(7,5) exactly, in 12.6 s and 248 s on one core;
  - the store and the executors (phase 3): `pq_rstore` (content-addressed blobs, fetched once by key, programs
    kept prepared) and `pq_rexec` (`local`, the default, or `{nodes, Nodes}` over erpc, a failed node's tasks run
    elsewhere); a reconstruction's columns run in parallel with the same result for any executor (res(7,4) in
    2.5 s, res(7,5) in 42 s on 16 local workers);
  - `pq_recon:rational/3` (phase 4): a rational function N / D reconstructed from its program, along lines with a
    shift, the coefficients by one multi-output Zippel, the shift taken out by a Taylor recursion; checked against
    FORM 5's `PolyRatFun` on sums of fractions (`comparison/form/fractions.erl`): the same results exactly, FORM
    300-12,000x faster on that task (an explicit sum, where adding fractions is the right tool);
  - `pq_recon:linear/4` (phase 5): the components of the solution of a linear system with polynomial entries, each
    reconstructed from eliminations mod p at the sample points (`pq_nif:ff_linsys/3`), never solved symbolically;
    against FORM 5's symbolic elimination under `PolyRatFun` (`comparison/form/linsys.erl`), the same answers, FORM
    12-18x faster when the answer is about as large as the elimination, PQCF 25 to 35,000x faster when it is
    compact (as in integration-by-parts reductions; n = 24: FORM 697 s, PQCF 0.02 s);
  - several right-hand sides at once (one per master integral): sparse elimination per point in Laporta order, all
    right-hand sides tested for 0 together; on the integration-by-parts identities of the one-loop massless box
    (`comparison/form/ibp.erl`, 86 to 497 equations), the same reductions as FORM 5 and, core for core, 1.9-4.5x
    faster than its symbolic elimination from about 260 equations up; on a two-loop family, the massless propagator "kite" (222 to 2,102
    equations, the kite reduced to the bubble squared and the sunset), 3-4.8x faster; `linear/4` passes its
    options' expression limits to the compiler, and takes A's rows sparse (maps Column => entry), only the non-zero
    entries compiled; B's rows may be sparse too (`rhs => M`), the right-hand sides carried as sparse columns of the
    elimination (`pq_nif:ff_linsys/7` with both patterns, a blob in the store), so thousands of masters cost their
    non-zero entries: the kite 11-38x faster than FORM on one core; `pq_nif:ff_pivots/3` (the pivots of a sparse
    matrix mod p, shortest row first); `comparison/form/ibp_family.erl`, IBP identities for any family from its
    propagators (scalar products solved, zero sectors by the Lee-Pak criterion, Kira-style seeds; it found and fixed
    a gap in the kite's zero rule); against FIRE 6.5.2 (FLINT) on its double box, FIRE reduces in 16 s what costs
    PQCF 56 s per sample point -- the gap is the Laporta front end (no symmetries, no equation selection), recorded
    in R1 with the steps that would close it; then closed: mark and sweep (`pq_nif:ff_pivots/4` records which
    rows each row was combined with, and only the equations the target needs are solved: 0.19 s per point instead
    of 56 s), sector symmetries by Pak's criterion (`ibp_family:canon/2`: the double box reduces to its 8 masters), a
    cache of solved points shared by a system's right-hand sides (prefix-stable sample points), and dot-free seeds
    (`{dots, D}`); the double box checked equal to FIRE's at a random point mod p, and reduced in 3.5 s on one core
    against FIRE's 7.7-8.4 s; seeds chosen automatically (`{auto, Family[, Target]}`: the smallest seeds whose
    reduction is complete, every non-zero master itself a seed): the double box in 3.95 s (FIRE 7.7-8.4 s), and
    I(1,1,1,1,1,1,1,-2,-1) in 32.2 s (FIRE 11.9-12.2 s), both checked equal to FIRE's;
  - `pq_simplify`: results simplified by verified rewriting (exact division by an alphabet of candidate factors,
    collecting by monomials, the smallest of the expanded, factored and collected forms by the additions written,
    checked at random points; `format/1` prints the structure); the double box's coefficients 25 -> 18 additions
    and 79 -> 58, in the textbook form (factors in d, monomials in s, t); the remaining numerators are irreducible
    (checked with FORM), so a network proposing factors has nothing to add on these results;
  - master-integral differential equations derived by PQCF (`ibp_family:derivative/3`, external-momentum
    operators reduced by IBP; `comparison/pinn/box_de.erl`): the one-loop box's textbook equation, checked to 1e-8
    against direct Feynman-parameter integration; a physics-informed neural network on it
    (`comparison/pinn/pinn_box.py`, after TechUPGRADE's PINN-PONN, in its JAX image): 1e-4 accuracy after 52 s of
    training against an ODE solver's 1e-13 in 19 ms, and 0.3-5.7% as a surrogate over unseen dimensions; recorded in
    R1 with what was positive and what was negative; `ibp:system/4` with `order => numerators | dots`;
  - denominator guessing in `pq_recon:linear/4` (default on, `guess => false` to turn off): univariate denominators
    on one line per variable matched against the factor alphabet, F * D reconstructed as a polynomial through a
    scaled black box (`pq_nif:ff_linsys/8`), falling back to the rational reconstruction when the guess fails its
    check; bit-identical results; the double box I(1,1,1,1,1,1,1,-2,-1) 32 -> 9 s reconstruction, whole runs on one
    core 3.16 s and 11.87 s against FIRE's 8.10 s and 12.21 s;
  - a learned elimination for linear systems (`Plan` in `ffslp.rs`: recorded at the first point, replayed on a
    flat array at the others, the wanted component through the transposed triangular system), `invmod` by the
    extended Euclidean algorithm instead of Fermat (every inversion in PQCF), Shoup multiplication in the replay:
    a double-box sample point 31 -> 2.9 ms, bit-identical results; whole runs on one core 1.21 s and 3.73 s against
    FIRE's 7.99 s and 11.99 s (6.6x and 3.2x); three executors
    behind one task interface (this node, Erlang nodes, FVLLMP federation peers), from a single node to 100+.

## 1.3.0 (7 November 2026)

The symbolic engine against FORM 5. Everything is additive; the default programs and canonical forms are the same as
1.2.0's, term for term. Expansion: dense products about 4 times faster than FORM, sparse products, powers and the
symbolic mesh 1.5–4.3 times FORM's time (from 1.3–32 times on 2 November). Compiled programs: on FORM's resultant
benchmarks, fewer operations than FORM 5's `O4` on all three (decision records K1, K2, H1–H4).

### Added

- **Horner schemes for compiled expressions** (`H1`):
  - `pq_symbolic:compile/3` with `#{optimize => horner}`: every polynomial is compiled as a multivariate Horner
    scheme, in a fixed occurrence order, with identical sub-polynomials shared;
  - on the resultant benchmarks of FORM's expression-simplification papers: 2.6–3.0× fewer operations than the
    default program, 4–6 % fewer than FORM 5's `Format O1`, and evaluation 3.4–3.6× faster;
  - the default is unchanged, bit for bit.
- **Common subexpressions across compiled programs** (`H2`):
  - `optimize => horner_cse` (Horner, then the pass) and `optimize => cse` (the default program, then the pass);
  - sums and products are treated as n-ary nodes, and the operand pair that most nodes share becomes one
    operation, repeatedly (FORM's `O2` step);
  - on the resultants: 4,285 / 18,166 / 61,404 operations, against Horner's 4,678 / 19,490 / 67,928 and FORM 5's
    `O2` 3,969 / 16,398 / 55,685.
- **A search over Horner variable orders** (`H3`):
  - `order_search => #{trials, seed, batch, nodes}` with `optimize => horner | horner_cse`;
  - batches of candidate orders are scored in parallel, one process each, on this node or spread over nodes, with
    the same program for any number of cores or nodes;
  - on the resultants, 160 trials give 4,273 / 16,408 / 54,643 operations, below FORM 5's `O2` on the two larger
    ones; 16 cores are 6.3 times faster than one.

- **Toward FORM's `O3`** (`H4`):
  - Horner nodes are computed as their primitive parts (content and sign taken out), so negated and scaled copies of a
    sub-polynomial share one slot;
  - a subtraction operation `{sub, Slot, A, B}` in Horner programs;
  - `order_search => #{strategy => mcts}`: a parallel, deterministic Monte Carlo tree search over order prefixes,
    with `exploration` and `direction`;
  - the search scores an order by counting its Horner operations in Rust (`pq_nif:horner_prepare/2`,
    `horner_count/2`), equal to the built program's count;
  - on the resultants: below FORM 5's `O2` at the occurrence order (3,771 / 15,875 / 52,423); with the tree search
    2,741 / 10,243 / 34,799 (5,000 / 5,000 / 10,000 trials), below FORM's `O4` (2,936 / 10,722 / 35,419) on all
    three.

### Changed

- **Polynomials that stay in Rust** (`K2`):
  - inside `canonicalise`, large plain products and powers stay in Rust between operations, and the result's
    expression tree is built there in print order;
  - the canonical form is the same, term for term; `poly_kernel => false` turns it off with the kernel;
  - against FORM 5: powers 4.3×, sparse products 2.4× (from 25–32× on 2 Nov), dense products about 4× faster.
  - exponentials exp(k i x) of symbols count as kernel variables, and `circuit_unitary` keeps its entries in Rust:
    the 8-mode symbolic mesh in 0.6 s (1.5× FORM 5, from 5×).
- **Products of polynomials in Rust** (`K1`):
  - `pq_sym_poly` multiplies plain polynomials in a NIF: symbols, conjugates and π with integer exponents, and
    exact coefficients;
  - exponents are packed into 64 bits and coefficients are exact i128, with a fallback to the Erlang product on
    overflow, on exponents that don't pack, or when no NIF is loaded;
  - powers of plain polynomials too (`pq_nif:sym_poly_pow/5`), with the multinomial enumeration as the fallback;
  - the result is the same polynomial, term for term (D1); `poly_kernel => false` turns the kernel off;
  - the canonical print order sorts once instead of twice, and by binary keys when every monomial is a product of
    symbols with integer exponents: 2.4 times faster on 591k terms, with the same output.

## 1.2.0 (5 November 2026)

Beyond the specification's roadmap: PQCF from Python, and faster density matrices. Everything is additive and
results are bit for bit the same as 1.1.0's. Upgrade every node of a cluster together, as for 1.1.0.

### Added

- **A JSON service and a Python client** (`C1`):
  - `pq_service`: version, validate, backend_for, run, sample, probabilities, evolve and gradient as JSON
    requests;
  - two transports: stdio (`-s pq_service stdio`), and HTTP (`-s pq_service http`, `start_http/1,2`, 127.0.0.1 by
    default);
  - atom-safe: option keys, enum values and symbol names are matched against fixed tables or the circuit's own
    symbols.
- **The `pqcf` package** (`python/`, standard library only):
  - `Circuit`, a builder of the JSON circuit format;
  - `connect()`, the client;
  - observables;
  - `pqcf.convert`, to and from Perceval and Piquasso;
  - `pqcf.linear`, a passive circuit's unitary.
- `pq_json:predicate_from_json/1` (internal).

### Changed

- **Faster density-matrix gates** (`B2`): each gate is prepared once (index groups, sector matrices, the Clements
  mesh, the truncated matrices) instead of for every column and row. On photon-number blocks, number-conserving
  gates touch only the stored sectors.
  - Dense ρ: 5–23× per gate.
  - Blocks: 25× per beamsplitter.
  - Cutoff 20 on 4 modes, end to end: 12 s → 1.2 s.

  Results are identical bit for bit.

## 1.1.0 (5 November 2026)

The first release after 1.0.0-rc.1: Phase 4 of the specification's roadmap (FS 15). Everything is additive. Code
written for 1.0.0-rc.1 runs unchanged, and with the defaults its results are bit for bit the same, including the
sampling streams. Each item has a decision record in `docs/decisions/`.

### Added

- **Threshold detection of displaced Gaussian states** (the loop torontonian, `L1`). Threshold probabilities and
  the threshold sampler now accept non-zero means instead of returning `nonzero_mean_threshold`.
- **Feed-forward** (`F1`):
  - the instructions `{feed_forward, Gate, #{Symbol => Mode}}` and `{conditional, Predicate, [Gate]}`;
  - gate parameters and conditions may depend on earlier measurement outcomes in the same shot;
  - the JSON ops `feed_forward` and `conditional`.
- **Gradients** (`AD1`):
  - `mode => adjoint`: one forward and one backward pass over the Fock state, on density matrices when there is
    loss, with the loss transmissivity differentiated too;
  - `mode => parameter_shift`: exact shift rules with frequencies bounded by the cutoff, also with loss;
  - the strategy modules `pq_gradient_adjoint` and `pq_gradient_shift`, and the type `pq_gradient:observable/0`.
- **Noise models** (`N1`):
  - `visibility` and `gram`: partially distinguishable photons in boson sampling;
  - `multi_photon`: multi-photon source terms;
  - `dark_counts`: detector dark counts, also in exact probabilities;
  - the instruction `{thermal_loss, Mode, Eta, Nth}`, exact on the Gaussian backend and truncated on density
    matrices, with the JSON op `thermal_loss`.
- **GPU offload** (`G1`): `device => cpu | gpu | auto` for permanents, boson probabilities and hafnians, on NVIDIA
  GPUs through the CUDA driver.
  - The GPU computes every term bit-identically to the CPU kernel, and the results are reproducible bit for bit.
  - `auto` uses the GPU for permanents from n = 24 and keeps hafnians on the CPU (spike G1 on the GB10).
  - Building needs no CUDA.
- **Density matrices by photon-number blocks** (`B1`): `dm_storage => blocks` and `max_dm_entries`. Results are the
  same, with far less memory for number-conserving circuits: cutoff 20 on 4 modes takes 0.18 GB instead of 1.81 GB.
- **New metadata and info keys:** `dark_count_truncated_mass`, and `device` in the permanent and hafnian kernels'
  info maps.

### Fixed

- `pq_stats:mean_photon_number/2` had no clause for density matrices, so the numeric gradient of a mean photon
  number through loss crashed with `function_clause`.

### Changed

- The NIF contract is 0.2: new kernels only, and `pq_nif` remains internal.
- **Upgrade every node of a cluster together.** The node handshake (unchanged) shares work only between nodes with
  the same `pq` version and NIF minor version. A 1.1.0 node and a 1.0.0-rc.1 node do not split jobs: each runs its
  own jobs locally, with the same results.
- `pq_api:backend_for/2` chooses `boson` when `visibility`, `gram` or `multi_photon` is given, and `fock_dm` for
  circuits with `thermal_loss`. These are the only backends that model them.

## 1.0.0-rc.1 (1 November 2026)

- The API frozen (`docs/decisions/API-freeze.md`).
- Every functional requirement of the specification implemented and tested.
- Benchmarks on a provisional reference machine (`bench/HARDWARE.md`).
