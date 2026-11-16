# Reproducing the quantum optics results

This page takes you from nothing to the numbers in [RESULTS.md](RESULTS.md): PQCF's master equations against QuTiP
on the same tasks, on your machine, and then to the two other quantum optics checks behind the published results:
cumulant expansions against QuantumCumulants.jl and collective atoms against QuTiP's PIQS. What the tasks are and why:
[README.md](README.md).

## What you need

- **Docker**, and about 8 GB of disk for four images: PQCF's development image (2.5 GB), the QuTiP image (QuTiP
  5.3.1, about 1 GB), the Julia image for QuantumCumulants.jl (about 3 GB) and the oracle image with mpmath.
- **git**, and a network connection for the first builds.
- Any Linux machine, x86_64 or aarch64. One core is used for the timings.

The images' base images are named by digest; the digests are in [REPRODUCE.md](../../REPRODUCE.md#base-images).

## 1. PQCF against QuTiP

```sh
git clone https://github.com/quanterall/pqcf.git
cd pqcf
git checkout aef6143
PQ_BENCH_HOST=my-laptop sh comparison/quantum-optics/run.sh 0
```

**Run the latest** version: the same commands without `git checkout aef6143`. The answers must still be the
published ones; what changed since is in [CHANGELOG.md](../../CHANGELOG.md).

- `PQ_BENCH_HOST` is a label for your machine; it goes into the results.
- The last argument is the core the tasks run on (default 8). Use a core with little else on it.

The script:

1. builds `pqcf-dev` from the repository's [Dockerfile](../../Dockerfile) (Erlang/OTP 27, Rust, PQCF) and
   `pqcf-qutip` from [qutip/Dockerfile](qutip/Dockerfile) (QuTiP 5.3.1, NumPy 2.4.6, SciPy 1.17.1);
2. runs PQCF's side ([qo_pqcf.erl](qo_pqcf.erl), 128 bits) and prints `PQCF done`;
3. runs QuTiP's side ([qutip_side.py](qutip_side.py)) with its default tolerances and with the tightest ones;
4. checks every answer against its closed form ([compare.py](compare.py)), writes `RESULTS.md` and
   `results/<label>-<date>.json`, and appends a line to `history.csv`. It exits 1 if any answer is wrong.

With both images built, the run takes about 75 seconds on our machine; the first time add the two image builds
(several minutes, plus Docker's base images).

## What you should see

`PQCF done`, then the script ends with exit status 0 (`echo $?` prints 0). `RESULTS.md` then has two tables:

```
| task | PQCF (128 bits) | QuTiP (default) | QuTiP (tight) |
| rabi | 36 digits, ... | 6 digits, ... | 11 digits, ... |
| revival | 36 digits, ... | 5 digits, ... | 9 digits, ... |
| cavity_decay | 36 digits, ... | 6 digits, ... | 12 digits, ... |
| fluorescence | 36 digits, ... | 16 digits, ... | 16 digits, ... |
```

and, for the driven, damped Jaynes-Cummings model at photon cutoffs 4 to 32, PQCF's time and QuTiP's time with the
difference between their answers (QuTiP default about 1e-7, tight about 1e-13 to 1e-14).

## How to compare with the published results

- **The digits** must be the published ones: 36 for PQCF on every task (the working precision, 128 bits), and the
  same counts for QuTiP, which are a property of its integrators, not of the machine.
- **The differences** in the Jaynes-Cummings table must be the published ones to the digits printed: they are
  QuTiP's errors, measured against PQCF's 36-digit answer.
- **The times** depend on your machine and core. Compare the ratios (PQCF's time over QuTiP's), not the seconds: the
  published ratios are in RESULTS.md and every run in [history.csv](history.csv).

`git diff comparison/quantum-optics/RESULTS.md` shows exactly what your run changed: only times should differ.

## 2. Cumulant expansions against QuantumCumulants.jl

```sh
sh comparison/cumulants/run.sh
docker run --rm pqcf-dev rebar3 eunit --module=pq_cumulant_tests
```

`run.sh` builds the Julia image ([../cumulants/julia/Dockerfile](../cumulants/julia/Dockerfile), QuantumCumulants
0.3.7; about 5 minutes the first time), solves the single-atom laser and the driven Jaynes-Cummings model at orders 2,
3 and 4, and rewrites `test/fixtures/qc_cumulants.json`. The tests then compare PQCF's equations and solutions with
it: the laser to 1e-13, the Jaynes-Cummings model to 1e-15. Expected: `N tests, 0 failures`. What is compared and
why: [../cumulants/README.md](../cumulants/README.md).

## 3. Collective atoms against QuTiP's PIQS, and superradiance against mpmath

```sh
docker run --rm -v "$PWD":/w -w /w pqcf-qutip python comparison/dicke/piqs_check.py         # about 8 minutes
docker run --rm -v "$PWD":/w -w /w pqcf-oracle python comparison/dicke/superradiance_mp.py  # about 1 minute
docker run --rm pqcf-dev rebar3 eunit --module=pq_dicke_tests
```

Build the oracle image first if you do not have it: `docker build -t pqcf-oracle test/oracle`. The two scripts
rewrite `test/fixtures/piqs_dicke.json` and `test/fixtures/superradiance_mp.json`; the tests require PQCF to agree with
PIQS within 1e-8 (PIQS computes its rates in single precision, so that is its own accuracy) and with the exact
superradiance equations (mpmath, 40 digits) to 34 digits. Expected: `N tests, 0 failures`. Details:
[../dicke/README.md](../dicke/README.md).

`docker run --rm pqcf-dev` uses the image built in step 1, which holds the repository as it was then: rebuild it
(`docker build -t pqcf-dev .`) if you changed anything since, or the tests read the old fixtures.

## The exact results (no external program)

The exact steady states and spectra (`pq_qo_exact`: SymPy's exact solutions, the Drummond-Walls and Bartolo et al.
solutions of the Kerr cavity) and the input-output networks (`pq_slh`) are checked by PQCF's own test suite:

```sh
docker run --rm pqcf-dev rebar3 eunit --module=pq_qo_exact_tests,pq_slh_tests,pq_lindblad_tests,pq_qop_tests
```

## Tested

From a fresh clone of the public repository on 16 November 2026, at aef6143, on an aarch64 machine (Cortex-X925 and
Cortex-A725): step 1 in 75 seconds with the images built, every digit count and every difference equal to the
published ones; step 2 (QuantumCumulants.jl image already built) and step 3 (PIQS in 6.6 minutes, mpmath in 1 minute)
rewrote the three fixtures byte for byte (`git diff` empty), and `pq_cumulant_tests` and `pq_dicke_tests` passed
(8 tests, 0 failures).
