# Reproducing the Feynman diagrams and Feynman rules

This page takes you from nothing to the published checks of `pq_diagrams` (every Feynman graph of a theory with its
symmetry factor and fermion sign) and `pq_lagrangian` (Feynman rules from a Lagrangian), against FORM's diagram
generator, zero-dimensional field theory and Peskin and Schroeder. What is compared and why: [README.md](README.md).

The comparison with FORM has two halves. A script runs FORM's generator and writes a fixture
(`test/fixtures/gt_diagrams.json`); PQCF's tests then compare PQCF's diagrams with it. The other checks (zero-dimensional
field theory, the Feynman rules) are computed inside the tests and need no fixture.

## What you need

- **Docker**, and about 5 GB of disk for two images: PQCF's development image (2.5 GB) and the FORM image (2.5 GB,
  FORM 5 compiled from its sources at a fixed commit, on top of PQCF's image).
- **git**, **python3** (only the standard library), and a network connection for the first builds.
- Any Linux machine, x86_64 or aarch64.

The images' base images are named by digest; the digests are in [REPRODUCE.md](../../REPRODUCE.md#base-images).

## 1. Get the code and build the images

**Reproduce the published run**, at the commit it was tested at (see "Tested" below):

```sh
git clone https://github.com/quanterall/pqcf.git
cd pqcf
git checkout 7570e62
docker build -t pqcf-dev .                                 # Erlang/OTP 27, Rust, PQCF
docker build -t pqcf-form -f comparison/form/Dockerfile .  # FORM 5
```

**Run the latest** version: the same commands without `git checkout 7570e62`. The fixture must still come out the
same; what changed since is in [CHANGELOG.md](../../CHANGELOG.md).

## 2. Only the tests (seconds)

```sh
docker run --rm pqcf-dev rebar3 eunit --module=pq_diagrams_tests,pq_lagrangian_tests
```

Expected: `N tests, 0 failures`. What they require ([README.md](README.md) has the table):

- for 24 cases of QED and QCD (self-energies to three loops, vertices to two, QED vacuum graphs to four), the same
  number of diagrams as FORM's generator and the same multiset of weights sign/S;
- phi^4 and phi^3 graphs equal to zero-dimensional field theory, exactly, and QED vacuum graphs with N fermion
  flavours (the fermion-loop signs);
- the vertices of phi^4, phi^3, Yukawa, QED, scalar QED and QCD equal to Peskin and Schroeder's, term by term, and the
  propagators equal to the inverse of their quadratic Lagrangians.

## 3. Regenerate the fixture and compare

```sh
sh test/oracle/gt_diagrams.sh
git diff --stat test/fixtures/gt_diagrams.json
```

The script prints the 24 cases with their loops and numbers of diagrams, for example:

```
qcd_glu_se 1 4
qcd_glu_se 2 23
qcd_glu_se 3 324
```

**How to compare:** `git diff` on the fixture should be empty: FORM's generator is deterministic. Then rerun the tests
on your fixture (step 2).

## Tested

From a fresh clone of the public repository on 19 November 2026, at 7570e62, on an aarch64 machine (Cortex-X925 and
Cortex-A725): the two images built; `gt_diagrams.sh` listed the 24 cases (the three-loop gluon self-energy with 324
diagrams, the ghost and quark self-energies with 97 each, the two-loop three-gluon vertex with 112) in under a second
and rewrote the fixture byte for byte (`git diff` empty); the tests passed (23 tests, 0 failures, 15 seconds).
