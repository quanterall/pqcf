# AD1: adjoint gradients and parameter-shift rules

**Decision: two gate-level gradient modes. The first, `mode => adjoint`, is reverse mode over the pure Fock state
in one NIF call. The second, `mode => parameter_shift`, uses exact shift rules with frequencies bounded by the
cutoff. Both turn gate-parameter derivatives into symbol derivatives through one shared chain rule.** Phase 4,
item 3 (FS 15 item 3). Released after 1.0.0-rc.1, purely additive under the freeze:

- two new values of `mode`;
- two strategy modules implementing the `pq_gradient` behaviour (`pq_gradient_adjoint`, `pq_gradient_shift`);
- `pq_gradient:observable/0` exported as a type;
- one new kernel, `pq_nif:fock_adjoint_gradient/5` (NIF contract 0.2).

## Adjoint mode

The observable is f = ⟨ψ_K|O|ψ_K⟩, where ψ_k = G_k ψ_(k-1) and O is real diagonal in the basis of FS 3.4 (a
weighted sum of marginal PNR projectors and photon numbers). The backward pass starts from λ_K = O ψ_K and steps
with λ_(k-1) = G_k† λ_k. The derivative for each parameter p of gate k is

df/dp = 2 Re ⟨λ_k| (∂G_k/∂p) ψ_(k-1)⟩.

That is one forward pass and one backward pass, whatever the number of parameters
(`native/pq_nif/src/fockgrad.rs`).

**The derivatives of the gates the backend applies.** The result is exact for the truncated model, not
approximately exact for the untruncated one:

| Gate | ∂/∂(first parameter) | ∂/∂(second parameter) |
| --- | --- | --- |
| phase shift e^(iφn) | i n ψ_k | — |
| Kerr e^(iκn²) | i n² ψ_k | — |
| beamsplitter | the sector generator with the phase, G_φ ψ_k (it commutes with B) | i(B n_a − n_a B): U_N(θ,φ)[j][k] = e^(iφ(k−j)) U_N(θ,0)[j][k] |
| displacement | (e^(iφ) a† − e^(−iφ) a) D: one row of D beyond the cutoff | i(m − n) D[m][n] |
| squeezing | ½(e^(−iφ) a² − e^(iφ) a†²) S: two rows beyond | i(m − n)/2 S[m][n] |

The generator of a truncated gate does not commute with the cutoff's projector. So P (A D) P is computed from the
matrix one or two levels larger, not as (P A P)(P D P), which would be wrong at the cutoff.

**Memory:** earlier states come back by inverting each exact gate: B(−θ, φ), the phase and Kerr with the opposite
sign, and the Clements mesh of U†. Only the state before each truncated gate is kept, so memory is
(1 + number of truncated gates) × D amplitudes, not K × D.

**Scope:** pure states only. Loss (the density-matrix backend), mid-circuit conditioning and feed-forward are
refused with `differentiate_unsupported`. Reverse mode through Kraus channels needs the adjoint channel kernel and is
left for later; until then, use `parameter_shift` or `numeric` there.

## Parameter-shift mode

When the photon number is bounded by P, a gate parameter enters as a trigonometric polynomial with frequencies
k·ω, k ≤ R. The equidistant rule of Wierichs, Izaac, Wang and Lin (Quantum 6, 677 (2022)) is then exact with 2R
evaluations:

f′(p) = ω Σ_(μ=1..2R) f(p + x_μ/ω) (−1)^(μ−1) / (4R sin²(x_μ/2)), where x_μ = (2μ − 1)π/(2R).

| Parameter | ω | R |
| --- | --- | --- |
| phase shift φ | 1 | P |
| beamsplitter θ | 1 | 2P |
| beamsplitter φ | 1 | 2P |
| Kerr κ | 1 | P² |
| displacement φ | 1 | 2P |
| squeezing φ | 1 | P |

P is the backend's cutoff, which for a Fock input without a cutoff is its photon number. Displacement and squeezing
magnitudes and the loss transmissivity are not periodic and are refused, and so is the Gaussian backend, whose photon
number is unbounded. Loss is allowed elsewhere in the circuit (the density-matrix backend).

**A bound found wrong by the tests: beamsplitter θ.** Within one sector N of the two modes, the eigenvalues of −iG
are N, N−2, …, −N. This suggests even frequencies only (ω = 2, R = P). But a later gate that mixes one of the two
modes with a third couples sectors N and N′ of different parity, so λ − λ′ takes odd values too. On a three-mode
circuit the spectrum of f(θ) had components at 1, 2, 3 and 4 (P = 2), and the ω = 2 rule was off by 0.28 for one
observable. The rule is now ω = 1, R = 2P. The same coupling argument was checked for the other rows: every
exponent there is a signed sum of occupation numbers, each at most P.

## The chain rule (`pq_gradient_chain`)

Every parameter expression is expanded the way lowering resolves it (FR-C04). A symbol other than the
differentiated ones takes its value from the binding, otherwise from its default, which may itself use the
differentiated symbols. Then d(expression)/ds is evaluated at the effective binding.

A weighted sum is flattened into terms, and a weight that depends on s adds dW/ds times its observable's value
(the kernel returns every term's value). Custom instructions, feed-forward and symbolic unitaries in the
differentiated symbols are refused.

Found on the way: `numeric` and `symbolic` mode take only *numeric* defaults for a weight (`numbers/1`), so a
weight whose default is an exact rational is `unbound_symbols` there. The gate-level modes resolve it. Left as it
is in the frozen modes, since an error there is not a wrong answer.

## Tests (`test/phase4_gradients_tests.erl`, Rust `fockgrad::tests`)

- **Rust:**
  - every derivative of every gate kind, the truncated ones included, against a fourth-order central difference
    within 1e-9;
  - mutation check: the sign of the squeezing a†² term flipped, caught.
- **Passive circuit, three modes, two photons, five symbols** (one used twice, through 2a): adjoint against the
  symbolic gradient within 1e-12 and parameter shift within 1e-10, for probabilities, photon numbers and a weighted
  sum.
- **Kerr, displacement and squeezing at cutoff 9** (a parameter given as p·t): adjoint against central differences
  within 1e-8.
- **Parameter shift against adjoint:**
  - Kerr (R = P²) and the phases of displacement and squeezing, within 1e-10;
  - mutation check: Kerr's R lowered to P, caught.
- **Parameter shift with loss,** against central differences within 1e-9.
- **Squeezed vacuum:** d⟨n⟩/dr = sinh 2r within 1e-10 at cutoff 30.
- **Chain rule:**
  - a default q = 3s;
  - the same default bound explicitly (no longer following s);
  - a weight t² in a weighted sum, against the closed form.
- **Cost:** 24 gate parameters, 8 symbols: adjoint 0.4–0.8 ms against 4–12 ms for central differences (48
  evaluations), with equal results within 1e-8.
- **Determinism:** bit-identical on repetition.
- **Refusals:** loss in adjoint mode, a backend other than Fock, a displacement magnitude, the Gaussian backend in
  shift mode, symbolic unitaries, feed-forward, a malformed gate list.
