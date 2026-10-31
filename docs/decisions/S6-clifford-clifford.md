# S6: the Clifford-Clifford sampler

**Decision: `clifford_clifford` as specified in AN 5.1, with the weights corrected** (amendment A-8 below).
Closed 30 October 2026.

## The correction

AN v0.1, section 5.1, gives the weights of the k-th photon's output mode as
`w_j = sum_l |U[j][c_l]|^2 * |perm(minor_l)|^2`. That is the distribution of distinguishable photons: for the
Hong-Ou-Mandel case (50:50 beamsplitter, one photon in each input) it gives p(1,1) = 1/2 instead of 0.
Clifford and Clifford's algorithm B uses the Laplace expansion of the k x k permanent along the new row:

    w_j = | sum_l U[j][c_l] * perm(U[y_1..y_(k-1); c_1..c_k without c_l]) |^2

with the complex permanents of the minors. The draw is normalised by the sum of the weights. AN's normaliser,
`sum_l |perm(minor_l)|^2`, equals that sum only when the input columns are orthonormal, so it fails for bunched
inputs (repeated columns); the sum of the weights covers both.

## Exit test (AN 10): chi-square p above 1e-4 against the exact distribution, n <= 6, m <= 8, 10^6 shots

`native/pq_nif/examples/s6.rs` samples each case 10^6 times (seed 2026), and 2 x 10^5 times with the AN v0.1
weights as a mutant; `test/oracle/s6_chi2.py` computes the chi-square statistic with bins of expected count below
5 pooled (FS 12.4). Unitaries are Haar-random (QR of a complex Gaussian matrix, phases fixed).

| m | input | patterns | chi-square / dof | p | AN v0.1 mutant p | 10^6 shots |
| --- | --- | --- | --- | --- | --- | --- |
| 2 | 1,1 | 3 | 0.3 / 2 | 0.86 | 0 | 0.6 s |
| 4 | 1,1,1,0 | 20 | 25.6 / 19 | 0.14 | 0 | 1.4 s |
| 5 | 1,1,1,1,0 | 70 | 88.7 / 69 | 0.055 | 0 | 2.7 s |
| 6 | 1,1,1,1,1,1 | 462 | 499.8 / 460 | 0.097 | 0 | 8.6 s |
| 8 | 1,1,1,1,1,1,0,0 | 1716 | 1749.9 / 1690 | 0.15 | 0 | 9.4 s |
| 4 | 2,1,0,0 | 20 | 23.7 / 19 | 0.21 | 0 | 1.4 s |
| 5 | 3,0,1,0,0 | 70 | 72.2 / 69 | 0.37 | 0 | 2.7 s |
| 8 | 2,2,1,1,0,0,0,0 | 1716 | 1766.1 / 1673 | 0.056 | 0 | 8.1 s |

All pass; the mutant is rejected in every case. Also: no (1,1) in 20,000 Hong-Ou-Mandel shots (unit test).

## Time per sample (one core, aarch64, shared machine)

| n | m | seconds per sample |
| --- | --- | --- |
| 12 | 24 | 0.0006 |
| 16 | 32 | 0.017 |
| 20 | 40 | 0.50 |

The growth is the O(n^2 2^n) of AN 5.1. **Phase 3 (31 October 2026):** all k minors' permanents now come from one
Glynn walk over the rows with prefix and suffix products of the column sums (`boson::minor_permanents`, O(n 2^n)).
The samples are unchanged (the same chi-square statistics in all 8 cases; GT-S7 still bit for bit), and the time
per sample fell to 0.0002 / 0.0037 / 0.078 s at n = 12 / 16 / 20 (3x, 4.6x, 6.4x; not on a pinned core, so the
ratios are approximate).

## Proposed amendment A-8 (AN 5.1)

Replace the weight formula as above; normalise by the sum of the weights; the self-test "weights sum to
`sum_l r_l`" applies only to collision-free inputs.
