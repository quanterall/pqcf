# Reproducing the published results

Every result we publish (here and on [quanterall.com](https://quanterall.com)) has a step-by-step page that takes you
from a fresh clone to the published numbers. This is the list. Each page has two paths: **reproduce the published
run** (check out the commit below after cloning) and **run the latest** version (no checkout); both must give the
published answers. Each page ends with a "Tested" section: when we last ran it from a fresh clone, where, and how long
it took.

Everything runs in Docker on any Linux machine, x86_64 or aarch64; the tests were run on aarch64 (Cortex-X925 and
Cortex-A725). The times are ours, with the images built; the first build adds the minutes given on each page.

## PQCF

| result | page | commit | time | needs |
| --- | --- | --- | --- | --- |
| Symbolic computation against FORM 5 | [comparison/form/REPRODUCE.md](comparison/form/REPRODUCE.md) | 92e2aeb | 1 min (4 the first time) | 5 GB, 1 core |
| Feynman-integral reduction against FIRE 6 | [comparison/fire/REPRODUCE.md](comparison/fire/REPRODUCE.md) | 92e2aeb | 1.5 min (+6 to build) | 5 GB, 1 core, python3 |
| Four-loop integrals against Forcer | [comparison/forcer/REPRODUCE.md](comparison/forcer/REPRODUCE.md) | 92e2aeb | 8 min (+4 to build) | 8 GB, 8 cores |
| Photonic head-to-head against five simulators | [comparison/photonic/REPRODUCE.md](comparison/photonic/REPRODUCE.md) | 92e2aeb | 14 min (+4.5 to build) | 5 GB, 1 core |
| The three-loop QCD β-function | [comparison/beta/REPRODUCE.md](comparison/beta/REPRODUCE.md) | eef0e6f | 15 min, images included | 8 GB, 8 cores |
| Holographic QCD: seven computations | [comparison/holographic/README.md](comparison/holographic/README.md) | d77ead6 | seconds to 4 min each; IHQCD 40 min on 16 cores (quick check 4 min) | 2.5 GB |
| Quantum optics against QuTiP, QuantumCumulants.jl and PIQS | [comparison/quantum-optics/REPRODUCE.md](comparison/quantum-optics/REPRODUCE.md) | aef6143 | 75 s against QuTiP; 7 min PIQS (+ images) | 8 GB, 1 core |
| QFT instruments: one-loop integrals against LoopTools, QCDLoop and mpmath; Dirac traces against FORM | [comparison/oneloop/REPRODUCE.md](comparison/oneloop/REPRODUCE.md) | d7d95e1 | tests 1 min; fixtures 30 min (C0) | 10 GB, 1 core |
| Genome research: the worked example and G0-G4 (branch `research/genome`) | [research/genome/REPRODUCE.md](https://github.com/quanterall/pqcf/blob/research/genome/research/genome/REPRODUCE.md) | bb8a5dd | 6 min for the example on 8 cores; G3 1.5 min, G4 4.5 min | 4 GB, network (NCBI) |
| The test suite (EUnit, Dialyzer, xref, the Python tests, cargo test, clippy) | `docker build -t pqcf-dev . && docker run --rm pqcf-dev` | any | about 10 min | 2.5 GB |

The quantum-optics comparison with QuTiP ([comparison/quantum-optics](comparison/quantum-optics/README.md)) runs with
its `run.sh`; it has no tested step-by-step page yet.

## GenHub2

| result | page | commits | time | needs |
| --- | --- | --- | --- | --- |
| Phase A's acceptance test (20 checks) | [acceptance/REPRODUCE.md](https://github.com/quanterall/genhub2/blob/main/acceptance/REPRODUCE.md) | GenHub2 d6a1733, PQCF d77ead6, `research/variants` 98fcfda | 20 min on 10 cores, plus 10 min to download GRCh38 | 25 GB, 16 GB of memory, network (NCBI) |

## Quantum computing simulations on the BEAM

All at commit 3c31766 of [quantum_computing_simulations](https://github.com/quanterall/quantum_computing_simulations).

| result | page | time | needs |
| --- | --- | --- | --- |
| Logical against physical error rate (the sweep) | [reproduce/sweeps.md](https://github.com/quanterall/quantum_computing_simulations/blob/main/reproduce/sweeps.md) | seconds | Python and NumPy |
| The surface code under circuit noise, against Stim | [reproduce/noise.md](https://github.com/quanterall/quantum_computing_simulations/blob/main/reproduce/noise.md) | 7 min (+7 to build) | 4 GB |
| Benchmarks against Qiskit Aer and Stim | [reproduce/benchmarks.md](https://github.com/quanterall/quantum_computing_simulations/blob/main/reproduce/benchmarks.md) | 7 min | 4 GB |
| Photons, the KLM gates, fusion, detectors and loss | [reproduce/photonic.md](https://github.com/quanterall/quantum_computing_simulations/blob/main/reproduce/photonic.md) | 7 min (+5 to build the tools) | the tools image, the service image |
| The ten demo circuits (39 checks) | [reproduce/demos.md](https://github.com/quanterall/quantum_computing_simulations/blob/main/reproduce/demos.md) | 2 min to build, 11 s to run | 2 GB |

## If a result does not reproduce

Please open an issue with your machine and the output: <https://github.com/quanterall/pqcf/issues> (or the
repository of the result).

## How the answers are compared

Times depend on the machine; answers do not. Each comparison is exact, to a stated number of digits, or statistical,
and says so:

| result | the answers must be | why this tolerance |
| --- | --- | --- |
| FORM | identical: PQCF's expression and FORM's, term for term, and the term counts | exact rational arithmetic on both sides |
| FIRE | equal, master by master, at a random point modulo a large prime | exact: unequal rational functions differ at a random point with probability below 1e-15 |
| Forcer | equal: PQCF's reduction times Forcer's masters, as expansions in ε, and the MS-bar sums equal to FORM's assertions | exact rationals and ζ-values |
| The β-function | equal to the known coefficients, line by line | exact rationals |
| Photonic head-to-head | permanents within 1e-8, torontonians within 1e-7, Fock states within 1e-9 (1e-6 between two of the peers); boson samples against the exact distribution by a chi-square test, p > 1e-3 | double precision on the peers' side: a permanent of 5e-12 at n = 28 is a sum of terms of order 1 |
| Quantum optics (QuTiP) | within 1e-4 of the closed form, or the answer is an error; the digits each tool gets are reported | QuTiP's integrators give about 6 digits at their defaults |
| Holographic QCD | the result files equal to the committed ones (only the commit line, date and times change); against the papers, to the papers' rounding; the SciPy checks to 1e-9, 1e-7 and 2e-6 ([the table](comparison/holographic/README.md#independent-checks-in-double-precision)) | PQCF computes to 20-30 digits; the papers and double precision limit the comparison |
| IHQCD quick check | the masses equal to the committed raw results to 1e-7 | 12 digits are computed |
| Genome research | every result equal to the published one (assemblies, phasing, modules, optima) | seeded and deterministic |
| GenHub2 acceptance | all 20 checks pass, every result equal to the published run, both trees' hashes identical | seeded and deterministic |
| Sweeps | Monte Carlo counts equal; the exact and closed-form logical rates within 1e-9 relative (3.9e-11 found) | the same seeded draws; the published rates carry 11 significant digits |
| Surface code under noise | every count equal; Stim's own counts within their statistical errors when its version differs | seeded on both sides |
| Benchmarks | the circuits (qubits, gates) equal; times compared as ratios | times are the measurement |
| Photonic schemes | the output files identical (`cmp`), except two files whose floating-point sums differ by about 1e-14 relative | seeded and deterministic |
| Demos | 39 of 39 checks pass | |

## Base images

Every Dockerfile names its base image by digest, so a build today starts from the same image as the published runs
(the digests are of multi-architecture indexes: the same line works on x86_64 and aarch64). Everything else in the
images is pinned in the Dockerfiles: versions, commits, checksums.

| base image | digest | used by |
| --- | --- | --- |
| `erlang:27` | `sha256:b944db88c674b941248f67562c6e7d56972d6ad9663876ff5ddeafebd8f28e35` | PQCF's [Dockerfile](Dockerfile) (`pqcf-dev`, and on it the FORM, Forcer, β, genome, variants and GenHub2 images); FORM's build stage |
| `debian:bookworm` | `sha256:268fb24c39bca92d326dd83ea3abebe4a3f7751c804ffef23bb6d857dd343ff8` | FIRE ([comparison/fire](comparison/fire/Dockerfile)) |
| `python:3.11-slim` | `sha256:e88e9763f943ec1834f992a4b51e0f24500486803e8bc534e5767af9ea65f6ce` | the photonic peers, QuTiP, the test oracles, the genome image |
| `python:3.10-slim` | `sha256:6ff506466f8b1e981719b468cc0023b8ef000b983cdec0eafe119d01f379e31d` | Strawberry Fields ([comparison/photonic/sf](comparison/photonic/sf/Dockerfile)) |
| `python:3.12-slim` | `sha256:78387bc3881b8273120a12ebe6c1ab22b018ccc2c9adf565ae1ac9b536e184ea` | quantum_computing_simulations: the tools and `qsim` |
| `elixir:1.18` | `sha256:16382709b3765b946748714575d8748ea7dcb7cb9a4849730bcdc56533503691` | quantum_computing_simulations: the service |

`erlang:27` and the two newer Python images are the ones the published runs were built from; for the other three the
images of the published runs were no longer on our machine, and the digests are those of the tags on 15 November
2026.
