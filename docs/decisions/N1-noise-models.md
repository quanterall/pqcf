# N1: noise models

**Decision: noise models as run options, not instructions, wherever the noise belongs to the photons or the
detectors rather than to a point in the circuit: `visibility` and `gram` for partially distinguishable photons in
boson sampling, and `dark_counts` for the detectors.** Phase 4, item 4 (FS 15 item 4). Released after 1.0.0-rc.1,
purely additive under the freeze:

- new options;
- new functions in `pq_noise`;
- `pq_sampling:gram_distribution/4`;
- one new kernel, `pq_nif:boson_partial_probability/4`;
- an optional `shared` entry for `sample_boson`.

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
