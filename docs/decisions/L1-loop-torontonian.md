# L1: the loop torontonian (displaced threshold detection)

**Decision: threshold probabilities and the threshold chain-rule sampler support displaced Gaussian states through
the loop torontonian, in thewalrus 0.22.0's conventions, computed incrementally inside the torontonian's existing
Cholesky walk.** Phase 4, item 1 (FS 15 item 2; FR-GB02 had deferred it). Released after 1.0.0-rc.1: an addition
under the API freeze, since displaced threshold detection used to return `{unsupported_input, nonzero_mean_threshold}`
and now returns probabilities and samples.

## The formula

P(B) = P(vacuum) · Re ltor(O_B, γ_B), where:

- ltor(O, γ) = Σ over the subsets Z of the clicked modes of (−1)^(c−|Z|) · exp(½ γ_Zᵀ (I − O_Z)⁻¹ conj(γ_Z)) /
  √det(I − O_Z);
- O = I − conj(Q⁻¹), and γ = Q⁻¹ conj(β) with β = (α, conj α), α = (x + ip)/√(2ħ);
- P(vacuum) = 1/Norm, the same normaliser as the photon-counting probabilities (`gaussian_to_amat`).

## The kernel

`toronto.rs` already walks the subsets depth-first, growing a Cholesky factor one mode (two indices) at a time.
The loop term needs a quadratic form for every subset. The factor is of conj(I − O_Z) = L L†, so the form equals
‖w‖² with w = L⁻¹ γ_Z, and w grows by one forward-substitution entry per index pushed. The work per term stays
O(|Z|), the total O(c² 2^c), and every term is real. The subtree blocks of D1 are unchanged, so results do not
depend on the thread count.

## The trap

The factor is of the *conjugate* of I − O. For the plain torontonian that is invisible, since only the real
determinant is used. The loop term then needs w = L⁻¹ γ, not L⁻¹ conj(γ). The first version used conj(γ) and gave
9.6·10⁶¹ for a probability of 0.97. A brute-force test of the definition built the same conjugated matrix and agreed
with it. The golden test against thewalrus caught it. The brute-force test now uses I − O itself.

## Tests

- `toronto.rs`, the loop torontonian by its definition (explicit inverses and determinants of I − O_Z), c = 1 to 7,
  within 1e-12 of the scale; bit-identical for 1 and 4 threads; with γ = 0 it equals the torontonian bit for bit.
- GT-H03 (`native/pq_nif/tests/gt_h03.rs`, fixtures from `test/oracle/gt_h03.py`): 18 random mixed, displaced
  states against thewalrus's `threshold_detection_prob`, 0 to 12 clicks, within 1e-10 relative.
- The coherent state's closed form: no click with probability exp(−|α|²), in Rust and through `pq_api`.
- The sampler: 1,000,000 shots of a displaced three-mode state against the exact probabilities: chi-square p = 0.851
  (FS 12.4).

## Limits

The exponent ½‖w‖² grows with the displacement. P(vacuum) shrinks it back, but for very large displacements the
terms could overflow. `max_clicks` (default 20) still bounds the cost.
