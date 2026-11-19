# Reproducing the running of QCD

This page takes you from nothing to the published checks of `pq_qcd`: the beta function and the quark mass anomalous
dimension of QCD to five loops, exactly, and the running and decoupling of alpha_s and the MS-bar quark masses at
arbitrary precision, against mpmath and CRunDec. What is compared and why: [README.md](README.md).

Each check has two halves. A script runs the reference programs and writes a fixture
(`test/fixtures/gt_qcd.json`); PQCF's tests then compare PQCF's results with it. You can run only the tests (they use
the committed fixture), or regenerate the fixture first and check that it comes out the same.

## What you need

- **Docker**, and about 4 GB of disk for three images: PQCF's development image (2.5 GB), the oracle image with
  mpmath (about 1 GB) and CRunDec 3.1.2 (small, built from its authors' repository at a fixed commit).
- **git**, **python3** (only the standard library), and a network connection for the first builds.
- Any Linux machine, x86_64 or aarch64.

The images' base images are named by digest; the digests are in [REPRODUCE.md](../../REPRODUCE.md#base-images).

## 1. Get the code and build the images

**Reproduce the published run**, at the commit it was tested at (see "Tested" below):

```sh
git clone https://github.com/quanterall/pqcf.git
cd pqcf
git checkout 7570e62
docker build -t pqcf-dev .                               # Erlang/OTP 27, Rust, PQCF
docker build -t pqcf-oracle test/oracle                  # mpmath
docker build -t pqcf-crundec comparison/qcd/crundec      # CRunDec 3.1.2
```

**Run the latest** version: the same commands without `git checkout 7570e62`. The fixture must still come out the
same; what changed since is in [CHANGELOG.md](../../CHANGELOG.md).

## 2. Only the tests (seconds)

```sh
docker run --rm pqcf-dev rebar3 eunit --module=pq_qcd_tests
```

Expected: `N tests, 0 failures`. What they require ([README.md](README.md) has the table):

- beta_0 .. beta_3 equal, term by term, to van Ritbergen, Vermaseren and Larin's formula for any gauge group, with
  SU(3)'s invariants computed by `pq_lie`;
- beta_0 .. beta_4 and gamma_0 .. gamma_4, n_f = 3..6, equal to mpmath's separate transcription to 38 digits;
- alpha_s(mu) and m(mu) in 10 cases (1 to 5 loops, n_f = 3..6) equal to mpmath's ODE solver to 28 digits, and to
  CRunDec's `AsmMSrunexact` to 1e-8;
- the decoupling of alpha_s, down and up, 8 cases: mpmath to 28 digits, CRunDec's `DecAsDownMS` and `DecAsUpMS` to
  1e-14;
- one loop against the closed form to 55 digits.

## 3. Regenerate the fixture and compare

```sh
sh test/oracle/gt_qcd.sh
git diff --stat test/fixtures/gt_qcd.json
```

The script runs mpmath, then CRunDec on the same cases, checks the two references against each other, and prints
`... coefficient sets, 10 running cases, 8 decoupling cases; mpmath and CRunDec agree (1e-8 running, 1e-14
decoupling)`. It stops with an error if they do not agree.

**How to compare:** `git diff` on the fixture should be empty: the cases are fixed and both references are
deterministic. Then rerun the tests on your fixture (step 2).

**About the tolerances.** PQCF is compared with mpmath to 28-38 digits. CRunDec is compared only to 1e-8 for the
running, because its `AsmMSrunexact` integrates with adaptive steps at a tolerance of 1e-10 per step: it is accurate
to about 1e-9 by design, not wrong ([README.md](README.md), "About CRunDec's numbers").

## Tested

From a fresh clone of the public repository on 19 November 2026, at 7570e62, on an aarch64 machine (Cortex-X925 and
Cortex-A725): the three images built; `gt_qcd.sh` printed `4 coefficient sets, 10 running cases, 8 decoupling cases;
mpmath and CRunDec agree (1e-8 running, 1e-14 decoupling)` in 8 seconds and rewrote the fixture byte for byte
(`git diff` empty); the tests passed (5 tests, 0 failures, 7 seconds).
