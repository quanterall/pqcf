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
  - **Next:** factorisation, then rational functions in normal form.
