# PQCF against the photonic simulators, head to head

The photonic simulators benchmark each other on four tasks: the torontonian, the permanent, boson-sampling
throughput, and Fock-space simulation of a 10-mode continuous-variable neural network. Their papers show these as
graphs, each on its own hardware:
- **Piquasso**, [arXiv:2403.04006](https://arxiv.org/abs/2403.04006), Figs. 7, 10 and 11;
- **DeepQuantum**, [arXiv:2512.18995](https://arxiv.org/abs/2512.18995), Fig. 3.

None of them includes PQCF. This directory runs the same tasks with PQCF and its peers, on the same inputs, the same
machine and the same core. It checks that every program gives the same answers, and publishes the times as tables.
`run.sh` does it all on the committed HEAD, and a weekly job runs it (Sundays, 03:30).

**Latest run:** [RESULTS.md](RESULTS.md). **Every run:** [`results/`](results) and [history.csv](history.csv).

## The programs

| program | version | what it is | tasks here |
| --- | --- | --- | --- |
| **PQCF** | the committed HEAD | Erlang/OTP with Rust kernels | all four |
| **The Walrus** | 0.22.0 | Xanadu's hafnian and torontonian library, Python and Numba ([JOSS 2019](https://doi.org/10.21105/joss.01705)) | permanent, torontonian |
| **Piquasso** | 8.0.2 | Python with C++ kernels ([arXiv:2403.04006](https://arxiv.org/abs/2403.04006)) | all four |
| **Perceval** | 1.3.1 | Quandela's platform; permanents from its C++ library exqalibur ([arXiv:2204.00602](https://arxiv.org/abs/2204.00602)) | permanent, boson sampling |
| **DeepQuantum** | 4.5.0 | PyTorch, CPU build ([arXiv:2512.18995](https://arxiv.org/abs/2512.18995)) | all four |
| **Strawberry Fields** | 0.23.0, the last release (archived) | Xanadu's continuous-variable platform ([arXiv:1804.03159](https://arxiv.org/abs/1804.03159)) | Fock simulation |

**PiquassoBoost**, which Piquasso's paper also benchmarks, is not on PyPI. Its kernels are now part of Piquasso's
own C++ extension.

## The tasks

The inputs come from fixed seeds (`gen_inputs.py`), so every program gets exactly the same matrices and circuits.

1. **Permanent** of the top-left n × n block of a random 2n × 2n unitary, n = 12 to 28. This is the amplitude of
   n photons through an interferometer.
2. **Torontonian**, for threshold-detector Gaussian boson sampling ([arXiv:1807.01639](https://arxiv.org/abs/1807.01639)).
   The state is a pure Gaussian state of k modes: every mode squeezed (r between 1 and 1.5, as in Gaussian boson sampling experiments), then a random
   interferometer, with every mode clicking. The matrix is O = I − Q⁻¹, of size 2k × 2k, for k = 8 to 24. Piquasso
   takes the same matrix written in quadratures, mode by mode. That change of basis is done before timing and
   leaves the value unchanged.
3. **Boson sampling** ([arXiv:1011.3245](https://arxiv.org/abs/1011.3245)): n single photons in the first n of n²
   modes, through a random unitary, n = 4 to 24. Every program draws exact samples, and its rate is reported:
   - PQCF and Perceval use Clifford and Clifford's algorithm ([arXiv:1706.01260](https://arxiv.org/abs/1706.01260));
   - Piquasso uses its own sampler;
   - DeepQuantum first computes the probability of every output pattern, so it runs only where there are at most
     200,000 patterns.
4. **Fock simulation**: one layer of a continuous-variable quantum neural network
   ([arXiv:1806.06871](https://arxiv.org/abs/1806.06871)) on 10 modes, starting from vacuum. The layer is an
   interferometer, squeezing, another interferometer, displacement, then Kerr: 66 gates with small parameters.
   The programs cut off the infinite state space in two ways:
   - PQCF and Piquasso keep every state with at most c photons in total;
   - Strawberry Fields and DeepQuantum keep at most c photons in each mode.

   Their tables therefore report the number of states kept, and how far each result is from the converged one.

## How the programs are run

- **Same core.** Every program runs alone, in its own container, pinned to the same core (a Cortex-X925 of the
  DGX Spark), with one thread (`OMP_NUM_THREADS`, `NUMBA_NUM_THREADS`, `torch.set_num_threads(1)`; PQCF's kernels
  `threads => 1`).
- **No compile time.** Each task first runs once on its smallest case without timing, so just-in-time compilation
  (Numba, PyTorch) is not counted.
- **Medians.** A time is the median of up to 7 runs: one run if it takes a second or more, else as many as fit in
  about two seconds. Boson sampling times one call drawing all the shots.
- **Time limit.** A size that takes more than 120 s is the last one run for that program and task.
- **PQCF** is called through its public Erlang interface: the kernels in `pq_nif`, and the Fock backend through
  `pq_api:evolve/3`. The peers are called through their Python interfaces (`peers.py`).

## How the answers are checked

`compare.py` checks every answer:
- **Permanents and torontonians:** every program against The Walrus, the reference both PQCF's and Piquasso's
  tests use. The limit is a relative difference of 1e-8 for permanents and 1e-7 for torontonians.
  - At n = 28, a permanent of 5e-12 is a sum of terms of order 1, and the programs differ among themselves by a few
    1e-9 in double precision. That is why the permanent limit is 1e-8.
  - The torontonian's state is squeezed strongly enough that the answer is well-conditioned. With weak squeezing
    (r between 0.2 and 0.6, our first run), the probability that every mode clicks is tiny: 1e-13 at k = 22. The
    torontonian's 2^k terms then cancel, and every program's answer at k = 22 to 24 was noise, differing by factors
    of up to a thousand.
- **Boson sampling:** each program draws 20,000 samples of 3 photons in 6 modes. These are tested against the exact
  output distribution (its permanents from The Walrus) by a chi-square test over all 56 patterns, and must give
  p > 0.001.
- **Fock simulation:**
  - Piquasso must give PQCF's mean photon number on every mode within 1e-9 at the same total cutoff, since it keeps
    the same states.
  - DeepQuantum must give Strawberry Fields' within 1e-6 at the same per-mode cutoff.

`run.sh` exits 1 if any check fails.

## Results

From the first run (12 November 2026, PQCF at ef79196). The full tables are in [RESULTS.md](RESULTS.md).

**Every answer agrees: 51 of 51 checks.**
- Permanents and torontonians agree within 1e-8 and 1e-7.
- Every program's boson samples pass the chi-square test against the exact distribution.
- Piquasso's Fock states match PQCF's to 5e-16 on every mode.
- DeepQuantum's match Strawberry Fields' at the same cutoff.

**Time.** In the table, "faster" or "slower" says how the program compares with PQCF, at the largest size every
program in that row ran:

| task | PQCF's place | details |
| --- | --- | --- |
| **Permanent** (n = 28) | **slowest of the four compiled programs** | Perceval 1.7× faster (6.4 s against 11.1 s), Piquasso 1.4×, The Walrus 1.3×. DeepQuantum about 100× slower (n = 26). |
| **Torontonian** (k = 24) | **second** | Piquasso 2.7× faster (3.0 s against 8.0 s). The Walrus 3.4× slower, DeepQuantum 68× slower. |
| **Boson sampling** | **fastest at every size** | Piquasso 4–360× slower, Perceval 13–66× slower, DeepQuantum 39× slower (4 photons). At 20 photons in 400 modes, PQCF draws 12.8 samples a second and Piquasso 3.2; Perceval holds at most 256 modes. |
| **Fock simulation** (total cutoff 10, 184,756 states) | **second** | Piquasso 3.2× faster (0.35 s against 1.10 s), with the same states and answers. At cutoff 4 PQCF is slightly faster. |

**Fock simulation with a per-mode cutoff.** Strawberry Fields and DeepQuantum keep a different set of states, so
their times are not directly comparable with PQCF's:
- At 4 photons per mode (1,048,576 states), DeepQuantum takes 0.15 s and Strawberry Fields 10.9 s.
- Their truncation error is still 5e-4.
- PQCF reaches 3e-7 with 43,758 states in 0.23 s.

**What this says about PQCF:**
- **Its sampler is the strongest part.**
- **Its permanent and torontonian kernels are behind the best C++ ones**, by 1.4–1.7× and 2.7×. Perceval's Glynn
  permanent and Piquasso's recursive torontonian ([arXiv:2109.04528](https://arxiv.org/abs/2109.04528)) show what
  is reachable on one core.
- **Its Fock backend is 3× behind Piquasso's** at large cutoffs.

These three are the targets for PQCF's kernels.

## A bug found on the way: Piquasso's permanent on ARM

**Piquasso 8.0.2, built on aarch64, returns wrong permanents.** The cause is in its `src/permanent.cpp` (and
`permanent_laplace.cpp`): `char parity = (minus_signs_all % 2 == 0) ? 1 : -1;`.
- `char` is signed on x86 but unsigned on ARM, so −1 becomes 255 there.
- For [[1, 2], [3, 4]] the permanent is (24 − 4)/2 = 10, but this build returns (24 + 255·4)/2 = 522.
- Random matrices come out wrong by factors of hundreds to thousands.

Piquasso has no aarch64 wheel, so on ARM it is always built from source. Its higher-level probabilities that we
checked were right, so the effect is limited to code that calls these kernels.

`peers/Dockerfile` builds Piquasso with `-fsigned-char`, as on x86. The code is unchanged and the permanents are
then exact (6e-12 relative at n = 20).

## Running it

```sh
PQ_BENCH_HOST=dgx-spark sh comparison/photonic/run.sh [core]   # default core 8
```

The images are `pqcf-dev`, `peers/Dockerfile` (The Walrus, Piquasso, Perceval, DeepQuantum) and `sf/Dockerfile`
(Strawberry Fields, which needs NumPy 1.x).
