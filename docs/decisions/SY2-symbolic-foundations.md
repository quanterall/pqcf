# SY2: the symbolic engine beyond FORM's domain (gcd, factorisation, then theoretical physics)

**Status: accepted (the project owner's decision, 14 Nov 2026: "yes, allow gcd and factorisation … start
implementation in your order"). Phase 1 in progress.**

## Decision

The scope fence of the symbolic specification (SY 1.2) excluded polynomial gcd and factorisation, and with them a
canonical form for rational functions. **That exclusion is lifted.** The engine gains, in this order:

1. exact multivariate gcd, square-free decomposition and factorisation, then rational functions in normal form;
2. truncated series (fractional powers and logarithms) and special functions with known derivatives and series;
3. variational calculus, ODE reduction and high-precision numerical ODE solvers;
4. tensors: components for a given metric (curvature), and abstract indices with canonicalisation;
5. holographic renormalisation;
6. holographic QCD models, from the action to the numbers.

Holographic QCD is the first target area because it uses nearly every one of these steps. Amplitudes
(polylogarithms) and gravitational-wave tooling reuse the same foundation and come after.

**What does not change:**
- results stay exact (numbers to a stated precision);
- every computation stays bounded (time, heap, size);
- the answer is the same on any number of cores or machines;
- PQCF stays Apache 2.0. GPL and LGPL systems (GiNaC, Cadabra2, FLINT, Singular, PARI) are used only as independent
  oracles in the comparison images, never linked. Arbitrary precision in Rust uses permissive code only.

## Why

FORM serves one kind of computation: very large polynomial expressions of perturbative quantum field theory. Other
areas of theoretical physics need a different set of tools. Holographic QCD needs:

| step | example | engine needs |
| --- | --- | --- |
| geometry | the metric of AdS₅, a black brane, a dilaton background | component tensors: Christoffel symbols, Riemann, Ricci, Einstein tensors |
| action | Einstein-dilaton, Yang-Mills in AdS, the DBI action of a D8 brane | functions of fields, square roots, expansion in fields |
| equations of motion | variation; Kaluza-Klein reduction to radial modes | variational derivatives, separation of variables |
| reduction | to Schrödinger form | changes of variable, normal forms of second-order ODEs |
| near-boundary analysis | the Fefferman-Graham expansion, with logarithms | series in fractional powers and logarithms, Frobenius solutions |
| renormalisation | counterterms, one-point functions, the Weyl anomaly | series algebra, local covariant expressions |
| spectra, thermodynamics | meson and glueball masses; entropy, speed of sound | high-precision eigenvalue solvers, root finding, quadrature |
| closed forms | Laguerre spectra, exact constants | special functions; integer relations (PSLQ) |

Every one of these steps needs reduced rational functions, which is why the gcd comes first.

## Phases and acceptance

Each phase ends with a checked result, tests in the suite and its results as files under `comparison/`.

### Phase 1: algebraic foundations

**Build:**
- **Multivariate gcd over ℚ by the modular method:**
  - `pq_sym_gcd` in Erlang: integer contents, the leading-coefficient correction, Chinese remaindering with
    Erlang's integers, trial division;
  - Rust (`native/pq_nif/src/mgcd.rs`): the gcd modulo a 63-bit prime, by Brown's dense interpolation.
- **Then:** square-free decomposition; univariate and multivariate factorisation (Berlekamp-Zassenhaus, Wang's
  Hensel lifting); rational functions in normal form; exact algebraic numbers.

**Acceptance:**
- results equal to the existing exact (but slow) Erlang gcd, on random cases;
- gcd(a·g, b·g) = g·gcd(a, b);
- the multivariate gcd and factorisation tasks of T. Ueda's polybench give FORM's and FLINT's results, and join the
  weekly FORM comparison.

### Phase 2: series and function atoms

**Build:**
- truncated series in x^(1/q) and log x;
- the function atoms exp, log, x^a, Bessel, Laguerre, ₂F₁, Gamma, digamma, each with its derivative, its series and
  high-precision evaluation;
- limits by series.

**Acceptance:** agreement with mpmath to 50 digits; the zeros of Bessel functions; the Frobenius solutions of
Bessel's equation.

### Phase 3: variational calculus and ODE solvers

**Build:**
- Euler-Lagrange equations from a Lagrangian;
- Liouville reduction to Schrödinger form;
- two independent arbitrary-precision eigenvalue solvers: shooting from Frobenius starts, and Chebyshev
  collocation;
- PSLQ.

**Acceptance:**
- the soft wall's m²ₙ = 4c(n + 1) (Karch, Katz, Son, Stephanov,
  [hep-ph/0602229](https://arxiv.org/abs/hep-ph/0602229)) to 30 digits;
- the hard wall's zeros of J₀ (Erlich, Katz, Son, Stephanov, [hep-ph/0501128](https://arxiv.org/abs/hep-ph/0501128)).

### Phase 4: tensors and geometry

**Build:**
- component curvature in any dimension;
- abstract-index canonicalisation by the Butler-Portugal algorithm (as xPerm, Martín-García,
  [0803.0862](https://arxiv.org/abs/0803.0862)).

**Acceptance:**
- AdS₅: R = −20/L²;
- the black brane solves Einstein's equations, with T = 1/(π z_h);
- Schwarzschild's Kretschmann scalar 48M²/r⁶;
- canonical forms equal to Cadabra2's (Peeters, [hep-th/0701238](https://arxiv.org/abs/hep-th/0701238)).

### Phase 5: holographic renormalisation

**Build:** Fefferman-Graham expansions with logarithms; counterterms.

**Acceptance:**
- the Weyl anomaly a = c = N²/4 (Henningson, Skenderis, [hep-th/9806087](https://arxiv.org/abs/hep-th/9806087);
  de Haro, Skenderis, Solodukhin, [hep-th/0002230](https://arxiv.org/abs/hep-th/0002230));
- the black brane's free energy −π²N²T⁴/8.

### Phase 6: holographic QCD models

**Acceptance:**
- Sakai-Sugimoto meson spectra ([hep-th/0412141](https://arxiv.org/abs/hep-th/0412141));
- improved holographic QCD glueball ratios ([0707.1324](https://arxiv.org/abs/0707.1324));
- Gubser-Nellore thermodynamics ([0804.0434](https://arxiv.org/abs/0804.0434)), its temperature scan distributed
  as PQCF tasks.

## Risks

| risk | what we do |
| --- | --- |
| multivariate factorisation is a large piece of work | gcd and square-free decomposition first, since most simplification needs only these |
| series with logarithms are subtle | two checks each: oracle series and numerical truncation errors |
| precision loss near singular points | Frobenius starts; two methods that must agree; precision raised until stable |
| no Mathematica or xAct to compare with | open-source oracles and published numbers |

## Progress

- **14 Nov 2026: record written; the modular gcd done.**
  - **API:** `pq_symbolic:gcd/2,3`, `pq_sym_gcd`, Rust `mgcd.rs`.
  - **How it got level with FORM (polybench defaults, core 8, 10 problems each, results equal to FORM's):**

    | version | nontrivial | coprime |
    | --- | --- | --- |
    | Zippel's sparse interpolation, Montgomery arithmetic, candidates tried before a confirming prime | 0.29 s | 0.28 s |
    | monomial contents first; univariate coprimality proof; contents left out of the images | **0.083 s** | **0.052 s** |
    | *FORM 5 with FLINT* | *0.075 s* | *0.044 s* |

  - **Before that,** on similar problems: Brown's dense method alone took 0.75–1.65 s per problem, and Zippel's
    sparse interpolation cut that about 20-fold.
  - **Why the monomial contents mattered:** random polybench factors have every variable in nearly every term, so
    even the "coprime" pairs share a monomial. Left in, it drove about 10,000 recursive gcds per problem.
  - **Square-free decomposition** (`pq_sym_gcd:squarefree/1`, Yun's algorithm) done the same day.
- **14 Nov 2026: factorisation** (`pq_symbolic:factor/1,2`, `pq_sym_factor`, Rust `mfactor.rs`).
  - **How it works:**
    - one variable by Zassenhaus's method;
    - several variables by Wang's method, with the multivariate Hensel lifting in Rust modulo 63-bit primes and
      Chinese remaindering;
    - every result checked by multiplying back.
  - **FORM comparison** (task T8, polybench's factorisation problems, 5 each, core 8): results equal to FORM's,
    factor for factor.

    | version | factor(a b) | factor(a b + c) |
    | --- | --- | --- |
    | first, all in Erlang modulo one large prime | 59 s | 0.49 s |
    | lifting in Rust modulo 63-bit primes | 2.0 s | 0.49 s |
    | monomial content first in the square-free step; univariate factoring mod p in Rust; every variable tried as the main one | **1.02 s** | **0.077 s** |
    | *FORM 5 with FLINT (`#factdollar`)* | *0.13 s* | *0.019 s* |

  - **Notes:**
    - FORM's `Factorize` statement takes 7.7 s for the same five nontrivial problems; its `#factdollar` route is
      the one compared here.
    - The remaining gap (8×) is in the lifting: each factor carries the whole leading coefficient, and its powers of
      (x_j - a_j) are dense. Wang's distribution of the leading coefficient's factors would remove the first.
- **14 Nov 2026: rational functions in normal form** (`pq_symbolic:cancel/1,2`, `numer_denom/1,2`, `pq_sym_ratfun`).
  - `equiv/3` now decides rational functions of symbols: `(x^2 - 1)/(x - 1)` against `x + 1` is true. It was
    `unknown` at the old fence (SY 10.3).
  - **Phase 1 done**, apart from exact algebraic numbers, which come with phase 2's needs.
- **14 Nov 2026: phase 2, first part: series and limits** (`pq_symbolic:series/3,4`, `limit/3,4`, `pq_series`).
  - **Series:** generalised power series with logarithms. Coefficients are rational functions of the parameters,
    with opaque atoms for the transcendental constants.
  - **Functions:** elementary ones, Bessel J_nu and Y_n, Laguerre, 2F1 and Gamma.
  - **Check:** the 32 expansions of `test/fixtures/gt_series.json` (SymPy 1.14, mpmath) agree in every coefficient.
  - **Convention:** the terms returned have exponents below N. SymPy also keeps x^N log(x)^m.
- **14 Nov 2026: phase 2, second part: high-precision evaluation** (`pq_symbolic:evaluate_digits/3`, `pq_mp`).
  - **Coverage:** the constants, the elementary and special functions, every atom the series produce, and the zeros
    of J_nu.
  - **Check:** against mpmath 1.3 to 48 digits (37 values, 8 zeros: `test/fixtures/gt_mp.json`), among them the
    phase 3 hard-wall values, the zeros of J_0.
  - **Phase 2 done.**
- **14 Nov 2026: phase 3, first holographic QCD results** (`pq_ode`, `comparison/holographic/walls.erl`,
  `RESULTS.md`).
  - **Soft wall:** from the action, PQCF derives V = z^2 + 3/(4 z^2) and the spectrum 4, 8, 12, 16 to 30 digits:
    m_n^2 = 4(n+1) exactly (hep-ph/0602229).
  - **Hard wall:** from the action, PQCF derives V = 3/(4 z^2), and the eigenvalues equal j_(0,n)^2 to 30 digits
    (hep-ph/0501128). The zeros come from `pq_mp:besselj_zero/3`, independent of the ODE solver.
  - **Second solver:** Chebyshev collocation of phi / z^s with inverse iteration (`pq_ode:eigenvalues_collocation/3`).
    It agrees with the closed forms to 25 digits on both walls (64 and 32 points), so the plan's "two methods that
    must agree" holds.
  - **PSLQ** (`pq_mp:pslq/3`): integer relations among high-precision numbers, tested on known identities.
- **14 Nov 2026: phase 3, the hard wall's axial and pion sector** (`comparison/holographic/hardwall.erl`,
  `HARDWALL.md`).
  - **Computed:** m_rho, m_a1, F_rho^(1/2), F_a1^(1/2), f_pi (chiral limit) and m_pi (GOR), for models A and B.
  - **Agreement:** Table II's values to their rounding, and F_rho to its closed form (int phi^2 = 2 J_1(j)^2/j^2).
  - **Two differences, both checked against an independent SciPy solution:**
    - f_pi is 1.6-2.1% low: the chiral limit against the paper's m_q /= 0, whose O(m_q^2) contact term depends on
      the cutoff;
    - model A's m_a1 is 1358.24 MeV against the paper's 1363; SciPy gives 1358.2432 too. The parameter values behind
      the paper's number are open.
  - **Not computed:** g_(rho pi pi).
  - **Phase 3 done.**
- **14 Nov 2026: phase 6, Sakai-Sugimoto** (`comparison/holographic/sakai.erl`, `SAKAI.md`).
  - **Method:** eq. (3.28) in u = K^(1/3) has polynomial coefficients and two regular singular ends; the new
    `pq_ode:eigenvalues_sl/3` solves it.
  - **Results:** lambda_1..4 = 0.669314492, 1.568765751, 2.874322527, 4.546102269 (25 digits). The paper gives 0.67,
    1.6, 2.9, 4.5, and m_a1^2/m_rho^2 = 2.344 against its ~2.4.
  - **Checks:** an independent SciPy shooting agrees as its cutoff grows (0.66931458 at Z = 16000), and the result
    does not depend on the match point.
- **14 Nov 2026: phase 6, Gubser-Nellore** (`comparison/holographic/gubser.erl`, `GUBSER.md`; `pq_jet`, a new
  arbitrary-precision Taylor-series integrator for nonlinear systems).
  - **Method:** the master equation (34) is integrated from the horizon series (36), with the integrals of (45) and
    (51) alongside; then c_s^2 = d log T / d log s.
  - **Checks:**
    - Chamblin-Reall (G = -1/(3 gamma)) is kept exactly;
    - Delta = 3.927 (paper 3.93), and the infrared limit is 0.1497 (paper 0.15);
    - c_s^2 dips to about 0.06 near T_c, as the paper's figures 3 and 4 show.
  - **Against SciPy:** a double-precision solution of the same equations agrees to 6 digits up to phi_H = 4 and
    becomes unstable beyond (the 0/0 at the horizon). PQCF's values do not change from 20 to 30 digits.
  - **Limit of the comparison:** the paper gives figures, not tables, so the dip agrees qualitatively.
- **14 Nov 2026: phase 4, component tensors** (`pq_geometry:curvature/2`). The acceptance checks of phase 4 that
  concern components pass exactly:
  - AdS_5: R = -20/L^2, R_ab = -(4/L^2) g_ab, Weyl tensor zero (Kretschmann 40/L^4);
  - the AdS_5 black brane solves Einstein's equations with Lambda = -6/L^2, and its Hawking temperature is
    T = 1/(pi z_h); its Kretschmann scalar is (40 + 72 z^8/z_h^8)/L^4;
  - Schwarzschild, in x = cos theta: Ricci flat, Kretschmann 48 M^2/r^6.

  Each takes a fraction of a second.
- **14 Nov 2026: phase 4, abstract indices** (`pq_tensor`).
  - **Canonical forms of tensor monomials:** the least index string over the slot-symmetry group and the orders of
    equal tensors, the dummies renamed by first appearance. This is the canonical representative of the double coset
    that Butler-Portugal computes, found by exhaustive search; it is exact and takes 0.3 s for every cubic Riemann
    contraction.
  - **Oracle:** Cadabra2 is not installed. SymPy 1.14's `canon_bp`, an independent Butler-Portugal implementation,
    serves instead: the same classes, zeros and relative signs on 160 random monomials.
  - **Literature:** there are 4 quadratic and 13 cubic classes under the slot symmetries. With the cyclic identity
    (exact rank on random algebraic curvature tensors), 3 and 8 are independent, the numbers of Fulling et al. (1992,
    Class. Quantum Grav. 9 1151, https://doi.org/10.1088/0264-9381/9/5/003).
  - **Phase 4 done.** Butler-Portugal's group-theoretic algorithm (Schreier-Sims) remains the route for products of
    more than about four Riemann tensors.
- **14 Nov 2026: phase 5, the black brane's free energy** (`comparison/holographic/brane.erl`, `BRANE.md`,
  `test/pq_holo_tests.erl`).
  - **Computed in PQCF, from the metric up:**
    - the curvature;
    - sqrt(g) and the bulk integral;
    - the extrinsic curvature;
    - the Gibbons-Hawking term and the counterterm (3/L) sqrt(gamma);
    - the series in the cutoff, where the eps^-4 terms (1/8 - 1/2 + 3/8) cancel;
    - the finite limit -L^3/(16 pi G z_h^4), which is -pi^2 N^2 T^4/8.
  - **Two limitations found:**
    - the engine does not cancel f f^-1 inside an expanded product (det g goes through `cancel` first);
    - `pq_series` keeps sqrt(L^8) as an atom, since it does not know the sign of L (sqrt(gamma) is factored by hand
      as sqrt(A^4) sqrt(f)).
- **14 Nov 2026: phase 5, the Weyl anomaly** (`comparison/holographic/weyl.erl`, `WEYL.md`).
  - **Method:** the Fefferman-Graham expansion on two Einstein boundaries, S^4 (W = 0) and S^2 x S^2 (W^2 /= 0).
    The 5-dimensional Einstein equations, solved order by order, give:
    - g_(2) = -(1/2r^2) g0 and -(1/6r^2) g0, the Henningson-Skenderis -(R_ij - R g_ij/6)/2;
    - Tr g_(4) = (1/4) Tr g_(2)^2.
  - **Result:** the log divergence gives a = c = N^2/4.
  - **Phase 5 done.**
