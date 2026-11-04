# N1: noise models

**Decision: noise models as run options, not instructions, wherever the noise belongs to the photons, the sources
or the detectors rather than to a point in the circuit:**

- **`visibility` and `gram`:** partially distinguishable photons in boson sampling;
- **`multi_photon`:** multi-photon source terms;
- **`dark_counts`:** detector dark counts.

**Noise at a point in the circuit is an instruction: `{thermal_loss, Mode, Eta, Nth}`.** Phase 4, item 4 (FS 15 item 4). Released after 1.0.0-rc.1,
purely additive under the freeze:

- new options;
- new functions in `pq_noise`;
- `pq_sampling:gram_distribution/4`;
- new kernels: `pq_nif:boson_partial_probability/4`, `gaussian_apply_thermal_loss/6` and
  `fock_dm_apply_thermal_loss/6`;
- optional `shared` and `extra` entries for `sample_boson`;
- the new instruction `thermal_loss`.

## Partial distinguishability

**The model behind `visibility`.** Photon a of input mode i has the internal state √p_i |0⟩ + √(1 − p_i) |a⟩,
with p_i = √M_i and every |a⟩ orthogonal to the others and to |0⟩. Then ⟨φ_a|φ_b⟩ = √(p_a p_b), and two photons
with the same M show the Hong–Ou–Mandel visibility |⟨φ_a|φ_b⟩|² = M, the number experiments quote. The unitary
acts on the spatial modes only, so the branches for different sets of photons in |0⟩ stay orthogonal through the
circuit and in the output.

The output distribution is therefore exactly a mixture:

1. draw for each photon whether it is in |0⟩;
2. sample those photons together by Clifford–Clifford;
3. sample each other photon alone from |U[j][i]|².

The sampler costs no more than Clifford–Clifford and is exact for this model (`boson::noisy_shot`).

**Stream order (D1),** per shot:

1. input loss, as before;
2. one uniform per surviving photon in mode order;
3. the Clifford–Clifford shot;
4. one uniform per private photon;
5. output loss.

Without the option nothing changes, so earlier samples are reproduced bit for bit.

**General overlaps (`gram => Rows`).** For any Gram matrix S of the internal states (photons in input-mode order):

P(T) = Σ_π (Π_k S[π(k)][k]) perm(M ∘ conj(M_π)) / (Π T! · Π_i perm(S_i))

- M = U[rows of T; the photons' columns], and M_π[r][k] = M[r][π(k)].
- S_i is the block of S for the photons of input mode i, so the normaliser is the input state's squared norm:
  Π S! for identical photons, 1 for distinguishable ones.

Derivation: write the ket and the bra as sums over assignments of photons to output rows; the overlap of the
internal states depends only on π = ρ⁻¹σ, and the remaining sum over σ is the permanent of M ∘ conj(M_π).

The formula costs n! permanents per pattern, so it serves as:

- the exact distribution behind `pq_api:probabilities` with `visibility` or `gram`;
- the sampler for `gram` (one uniform per shot from the full distribution, like `brute_force`);
- the oracle for the mixture sampler.

The matrix is checked to be Hermitian with a unit diagonal, and positive semidefinite by a Cholesky factorisation
with zero pivots allowed. The probabilities always sum to 1, so they cannot reveal a matrix that is not a Gram
matrix: the first test matrix, with overlaps 0.99 of mixed signs, passed the sum check with every probability
positive.

**Scope:** boson sampling only. Distinguishability options are refused on the other backends, and with them the
backend chosen by default is `boson`. Not with `brute_force` or `metropolis`.

## Dark counts

`dark_counts => ν` (or per mode): each detector adds Poisson(ν) counts per shot.

- A PNR column reports n + k, drawn by inversion.
- A threshold column clicks if it saw light, or otherwise with probability 1 − e^(−ν).
- Homodyne columns are untouched.

The uniforms come from stream (dark_seed(Seed), k) for shot k, where dark_seed is the seed XOR a fixed constant:
one uniform per counting column. So the noise does not depend on chunking, on the backend's own stream order, or
on the sampler.

The noise is applied inside the post-selection loop before the predicate, so post-selection sees what a detector
would report. Mid-circuit measurements with dark counts are refused: feed-forward would have to act on the noisy
outcome, which is a different model.

Exact probabilities convolve the noiseless distribution:

- **Threshold:** exactly.
- **PNR:** with at most K extra counts per detector, where K is the smallest with a Poisson tail below 1e-13. The
  remainder is reported as `dark_count_truncated_mass`.

## Multi-photon source terms

`multi_photon => q`, or a list of q per input mode. An occupied input mode emits one more photon with probability
q_i, in a private internal state: overlap 0 with every other photon. This is the noise photon of a quantum-dot
source (Ollivier et al., PRL 126, 063602 (2021)). A single-photon source then has g2 = 2q/(1 + q)², whether or not
the extra photon is distinguishable.

- **Sampling:** the shot draws one uniform per occupied input mode before input loss. The extra photons pass input
  loss (after the nominal ones) and are sampled alone.
- **Exact probabilities:** the mixture over which modes emitted, each branch by the Gram formula, with the extra
  photons' rows and columns zero off the diagonal. It combines with `visibility`, not with `gram`.

## Thermal loss

`{thermal_loss, Mode, Eta, Nth}` mixes the mode with a thermal state of mean occupation Nth at a beamsplitter of
transmissivity η:

- **Gaussian, exact:** μ → √η μ and Σ → ηΣ + (1 − η)(2Nth + 1)ħ/2 on the mode.
- **Density matrices:** the quantum-limited amplifier of gain G = 1 + (1 − η)Nth after pure loss of transmissivity
  η/G. The output mean photon number is then G(η/G)n + (G − 1) = ηn + (1 − η)Nth. The amplifier's Kraus operators
  have ⟨n+k|A_k|n⟩ = √C(n+k, k) G^(−(n+1)/2) ((G − 1)/G)^(k/2). The amplifier raises photon numbers, so components
  beyond the cutoff are dropped and reported as a trace deficit, like displacement (`truncated`; a warning above
  `truncation_tol`).

On the pure-state Fock backend thermal loss needs `promote_dm`, as plain loss does, and a Fock input with it goes to
the density-matrix backend by default. Its parameters are not periodic and the adjoint is for pure states, so both
gate-level gradient modes refuse them.

## Tests (`test/phase4_noise_tests.erl`, Rust `boson::tests`)

- **The Gram formula's limits (Rust):** identical photons equal `bs_probability`, distinguishable ones
  perm(|U|²)/T!, with a bunched input [2, 1, 0] within 1e-14; Hong–Ou–Mandel P(1,1) = (1 − x²)/2.
- **HOM dip:** P(1,1) = (1 − M)/2 and P(2,0) = (1 + M)/4 exactly for M = 0, 0.3, 0.6, 1. 400,000 shots at M = 0.6:
  fraction 0.1991, chi-square p = 0.356.
- **Mixture sampler against the Gram formula:** three photons, visibilities 0.9/0.5/0.7, a four-mode mesh,
  1,000,000 shots, p = 0.482. Power: the same samples against indistinguishable photons, p < 1e-10. M = 1 equals
  the boson distribution within 1e-13.
- **Mutation check:** p = M instead of √M in the shot is caught by both sampler tests (p = 0.000).
- **General Gram matrix** (complex overlaps, two photons in one mode): 300,000 shots, p = 0.670.
- **Refusals:** a non-Gram matrix (`not_positive_semidefinite`) and a non-Hermitian one.
- **Dark counts, threshold:** the closed forms, empty detectors clicking with 1 − e^(−ν), per-mode means;
  400,000 shots, p = 0.672.
- **Dark counts, PNR:** P(0, 3) = P₀(0, 3) e^(−2ν) and the mean rising by exactly ν; 400,000 shots on the Fock and
  boson backends, p = 0.305.
- **Post-selection with dark counts:** collision-free acceptance falls, and no accepted row has a collision.
- **D1:** identical samples for 1 or 8 workers and chunk 37, with both options.
- **Refusals:** other backends, visibility outside [0, 1], Metropolis, negative means, mid-circuit measurements.
- **Thermal loss, Rust:** vacuum gives the thermal distribution of mean (1 − η)Nth within 1e-15; one photon gives
  the mean η + (1 − η)Nth; Nth = 0 equals pure loss.
- **Thermal loss, Gaussian:** ⟨n⟩ and Σ_xx against the closed forms within 1e-13; Nth = 0 bit-identical to loss.
- **Thermal loss, density matrix against Gaussian:** a squeezed, displaced, split two-mode state with thermal loss
  on both modes, max difference 3.3e-16 at cutoff 18. Mutation check: the amplifier's G exponent off by one, caught
  (4.8e-2).
- **Thermal loss, other:** a Fock input's mean photon number; the JSON round trip; domain violations; the pure-state
  backend without `promote_dm`; both gradient modes refusing.
- **Multi-photon, Hanbury Brown–Twiss:** P(1,1) = q/2 and g2 = 2q/(1 + q)² exactly (q = 0.08: 0.13717); 300,000
  shots fit the distribution.
- **Multi-photon with partial distinguishability,** three photons on a mesh: the distribution sums to 1, and
  500,000 shots give p = 0.380.
- **Multi-photon, other:** identical samples for any chunking; other backends refused.
