# The three-loop QCD β-function

The β-function of QCD to three loops, computed from the gluon and ghost self-energies in Landau gauge, as a test
of PQCF on a real physics calculation. The result is the known one: Tarasov, Vladimirov, Zharkov (1980), and
Larin, Vermaseren, [arXiv:hep-ph/9302208](https://arxiv.org/abs/hep-ph/9302208).

**Result** (a = g²/16π², β(a) = −Σ β_k a^(k+2)):

- β₀ = 11/3 C_A − 4/3 T_F n_f
- β₁ = 34/3 C_A² − 20/3 C_A T_F n_f − 4 C_F T_F n_f
- β₂ = 2857/54 C_A³ − 1415/27 C_A² T_F n_f − 205/9 C_F C_A T_F n_f + 2 C_F² T_F n_f + 158/27 C_A T_F² n_f²
  + 44/9 C_F T_F² n_f²

For SU(3), that is 11 − 2n_f/3, 102 − 38n_f/3, and 2857/2 − 5033n_f/18 + 325n_f²/54. There is no ζ₃ in β₂,
and the higher poles of Z_g² agree with what the lower orders predict (for example −β₀³ at a³/ε³).

## How

1. **Diagrams:** FORM 5's own generator (`diagrams_`). The gluon and ghost self-energies are generated one-particle
   irreducible, and only the scaleless diagrams are dropped (a line from a vertex to itself, or both external legs
   on one vertex): 3/1, 18/6 and 257/81 diagrams (gluon/ghost) at one, two and three loops.
   - FORM's `NOSNAIL_` option is not used. At three loops it also removes non-zero diagrams, 3 gluon and 3 ghost
     diagrams; with it, β₂'s C_A³ and C_A² n_f terms come out wrong.
2. **Feynman rules, colour, Dirac algebra** (`se.frm`, run by FORM):
   - QCD's rules in Peskin–Schroeder's conventions, Landau gauge;
   - the colour factors by FORM's colour package (van Ritbergen, Schellekens, Vermaseren,
     [arXiv:hep-ph/9802376](https://arxiv.org/abs/hep-ph/9802376));
   - Dirac traces in d dimensions, then a Wick rotation to Euclidean scalar integrals.

   **Checked:** at one and two loops the gluon self-energy is exactly transverse; at one loop its gauge dependence
   is the textbook (13/6 − ξ/2) C_A.
3. **The integrals,** two ways:
   - **FORM + Forcer** ([arXiv:1704.06650](https://arxiv.org/abs/1704.06650)), as the reference, at one to three
     loops.
   - **PQCF**, at one and two loops (`beta_pqcf.erl`):
     - every graph is a PQCF integral family;
     - the terms are split by colour factor and reduced on one system (`pq_ibp:reduce_sums/3`);
     - the values of PQCF's master integrals come from Forcer;
     - the expansion in ε is PQCF's own (`../forcer/epsexp.erl`).
4. **Renormalisation** (`renorm.erl`, PQCF):
   - Z₃ and Z̃₃ make the two self-energies finite, order by order, with the bare coupling a_B = Z_a a;
   - Z_a = Z_g² = 1/(Z̃₃² Z₃): the ghost-gluon vertex is finite in Landau gauge (Taylor);
   - β_(L−1) = −L [a^L ε⁻¹] Z_a.

## Where PQCF stands

- **One and two loops:** PQCF reduces every integral itself, and its β₀ and β₁ are exact. Its bare self-energies
  equal Forcer's.
- **Three loops: not yet.** Landau gauge gives the integrals many squared propagators and high numerator powers.
  PQCF's choice of equations, one budget of dots and numerators for every sector, then grows past a million
  equations for a single graph: one ghost graph took 10 hours. Forcer reduces these with its precomputed rules: the
  three-loop gluon self-energy in 704 s on 8 cores (FORM + Forcer, everything included), the ghost in 67 s on one.
- **Next:** equations chosen for such numerators, or graphs split into simpler integrals first, as Forcer does.

## Files

- **`results/renormalisation_form_forcer_3loops.txt`:** Z₃, Z̃₃, Z_g² and β₀..β₂ at three loops, from FORM + Forcer's
  bare self-energies.
- **`results/renormalisation_pqcf_2loops.txt`:** the same at two loops, with every integral reduced by PQCF.
- **`results/form_forcer/bare_<field>_<loops>.txt`:** the bare self-energies (MS-bar, the gluon's times D − 1).
- **`results/pqcf/`:** PQCF's reductions (`reduction_<field>_<loops>.txt`: colour factor, master number, exact
  coefficient in d), its master integrals in Forcer's notation (`masters_*.frm`), and its bare self-energies.
