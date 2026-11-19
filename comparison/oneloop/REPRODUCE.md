# Reproducing the QFT instruments' results

This page takes you from nothing to the published checks of PQCF's first instruments for quantum field theory: the
one-loop integrals (two-point functions against LoopTools and mpmath, the three-point function C0 against mpmath,
QCDLoop and LoopTools, the four-point function D0 against mpmath, LoopTools and QCDLoop) and the Dirac traces in d dimensions (against FORM). What is compared and why:
[README.md](README.md); the results: the CHANGELOG entries for `pq_oneloop`, `pq_quad` and `pq_dirac`.

Each check has two halves. A script runs the reference programs and writes a fixture (a JSON file in
`test/fixtures/`); PQCF's tests then compare PQCF's results with that fixture. You can run only the tests (they use
the committed fixtures), or regenerate the fixtures first and check that they come out the same.

## What you need

- **Docker**, and about 10 GB of disk for five images: PQCF's development image (2.5 GB), the oracle image with
  mpmath and SymPy (about 1 GB), LoopTools 2.16 and QCDLoop 2.1.0 (small, built from their release tarballs, checked
  by SHA-256), and the FORM image (2.5 GB, FORM 5 compiled from its sources).
- **git**, **python3** (only the standard library), and a network connection for the first builds.
- Any Linux machine, x86_64 or aarch64. Note that our results about LoopTools and QCDLoop come from builds on aarch64;
  on x86_64 those two programs may agree with PQCF in more configurations (that is worth knowing: please tell us).

The images' base images are named by digest; the digests are in [REPRODUCE.md](../../REPRODUCE.md#base-images).

## 1. Get the code and build the images

```sh
git clone https://github.com/quanterall/pqcf.git
cd pqcf
git checkout 7570e62
docker build -t pqcf-dev .                                 # Erlang/OTP 27, Rust, PQCF
docker build -t pqcf-oracle test/oracle                    # mpmath, SymPy
docker build -t pqcf-looptools comparison/oneloop/looptools  # LoopTools 2.16
docker build -t pqcf-qcdloop comparison/oneloop/qcdloop      # QCDLoop 2.1.0
docker build -t pqcf-form -f comparison/form/Dockerfile .  # FORM 5 (for the Dirac traces)
```

**Run the latest** version: the same commands without `git checkout 7570e62`; the fixtures must still come out the
same, and what changed since is in [CHANGELOG.md](../../CHANGELOG.md).

## 2. Only the tests (a few minutes)

```sh
docker run --rm pqcf-dev rebar3 eunit --module=pq_oneloop_tests,pq_dirac_tests,pq_quad_tests
```

Expected: `N tests, 0 failures`. What they require:

- `pq_oneloop_tests`: the two-point functions B0, B1, B11 equal to mpmath to 35 digits and all five (B0, B1, B00,
  B11, dB0) to LoopTools to 1e-11, in 60 configurations; C0 in 40 configurations equal to mpmath's direct integral to
  28 digits where it exists, to mpmath's contour deformation to 1e-12 where that is stable, and elsewhere to the value
  two of the references agree on (README.md has the table); the refusals of divergent configurations.
- `pq_oneloop_tests` also requires D0 equal to mpmath's direct integral to 16 digits in three Euclidean
  configurations, and to LoopTools to about 1e-15 in three of the four configurations above threshold (the fourth,
  slow, is left out of the suite). D0 above threshold takes 45 s to 5 minutes at 18 digits on 20 cores, so this module
  takes a few minutes.
- `pq_dirac_tests`: 60 traces of up to 12 gamma matrices equal to FORM's `tracen` exactly.
- `pq_quad_tests`: the quadrature on integrals with closed forms, at 200 bits.

## 3. Regenerate the fixtures and compare

```sh
sh test/oracle/gt_oneloop.sh   # LoopTools and mpmath, two-point functions: seconds
sh test/oracle/gt_dirac.sh     # FORM, Dirac traces: seconds
sh test/oracle/gt_c0.sh        # mpmath (contour deformation, direct integral), QCDLoop, LoopTools: about 40 minutes
sh test/oracle/gt_d0.sh        # mpmath (D0's direct integral), LoopTools, QCDLoop: about 2 minutes
git status --short test/fixtures
git diff --stat test/fixtures
```

- `gt_oneloop.sh` prints `60 configurations; LoopTools and mpmath agree`: it checks the two references against each
  other before it writes the fixture.
- `gt_dirac.sh` prints `60 traces; longest 12`.
- `gt_d0.sh` prints the three Euclidean values, then `3 Euclidean, 4 physical`. It records LoopTools' and QCDLoop's
  values as they are: in our aarch64 build QCDLoop returns `nan` for two of the four configurations above threshold
  (D0(2, 3, 1, 4, 20, -5; 1, 2, 1, 3) and D0(1, 2, 3, 4, 30, -10; 2, 1, 2, 1)).
- `gt_c0.sh` prints one line per configuration (the deformed value, `stable` or `UNSTABLE`, and `direct` where the
  direct integral was done), then a summary: how many configurations are stable, how many have the direct integral,
  and how many QCDLoop and LoopTools disagree on.

**How to compare:** `git diff` on the fixtures should be empty: the configurations come from fixed seeds and the
references are deterministic. On another machine or architecture LoopTools and QCDLoop may give different values in
the configurations where they disagree with mpmath (their entries and the `qcdloop_agrees` / `looptools_agrees` flags
in `gt_c0.json`); mpmath's values must not change. Then rebuild the image and rerun the tests on your fixtures:

```sh
docker build -t pqcf-dev . && docker run --rm pqcf-dev rebar3 eunit --module=pq_oneloop_tests,pq_dirac_tests
```

## The finding about LoopTools and QCDLoop, by hand

```sh
printf 'C 1 1 1 1 1 1\nC 9 3 9 1 3 1\n' | docker run --rm -i pqcf-looptools lt_eval 2>/dev/null | grep -E '^-?[0-9]'
printf 'C 1 1 1 1 1 1\nC 9 3 9 1 3 1\n' | docker run --rm -i pqcf-qcdloop ql_eval 2>/dev/null | grep -E '^-?[0-9n]'
```

On our aarch64 machine they print (LoopTools, then QCDLoop; real and imaginary parts):

```
6.8721263867742088 34.18931254658434
2.6592102185045721 -1.2475576149154497
-0.67125310574800379 0
0.45122204525346032 0.23994830838944026
```

The first configuration, C0(1, 1, 1; 1, 1, 1), is real: D > 0 on the whole Feynman-parameter triangle, so
C0 = -0.671253105748004686 (mpmath's direct integral). The second, C0(9, 3, 9; 1, 3, 1), is above threshold with a
negative Kallen function; mpmath's contour deformation and PQCF give the same value. Compare with PQCF:

```sh
docker run --rm pqcf-dev sh -c 'erl -noshell -pa _build/default/lib/*/ebin -eval "
  [io:format(\"~p~n\", [pq_oneloop:value(c0, A, 20)]) || A <- [[1,1,1,1,1,1], [9,3,9,1,3,1]]], halt()."'
```

PQCF prints `{"-6.7125310574800468622e-1","0"}` and `{"2.5190306869311200980e-1","-5.5256131155165702499e-1"}`:
LoopTools is wrong in both configurations, QCDLoop in the second.

## Tested

From a fresh clone of the public repository on 19 November 2026, at 7570e62, on an aarch64 machine (Cortex-X925 and
Cortex-A725): the images built; `gt_oneloop.sh` (3 s), `gt_dirac.sh` (1 s) and `gt_d0.sh` (124 s) rewrote their
fixtures byte for byte (`git diff` empty); `gt_d0.sh` printed the three Euclidean values (0.035995465803132427694,
0.080700137002475906277, 0.032168270315717678003) and `3 Euclidean, 4 physical`; the tests passed:
`pq_oneloop_tests` 7 tests, 0 failures, in 205 s, and `pq_dirac_tests` with `pq_quad_tests` 5 tests, 0 failures.
`gt_c0.sh` was not rerun at this commit; its run follows.

Earlier, from a fresh clone of the public repository on 16 November 2026, at d7d95e1, on an aarch64 machine: the five images
built; `gt_oneloop.sh` (2 s), `gt_dirac.sh` (1 s) and `gt_c0.sh` (29 minutes) rewrote the three fixtures byte for byte
(`git diff` empty); `gt_c0.sh` printed "40 configurations; 28 stable under the deformation; 9 with the direct
integral; QCDLoop disagrees on 12; LoopTools disagrees or gives nothing on 12" (its flags compare each program with
the deformation, including where the deformation itself is unstable); the tests passed (10 tests, 0 failures).
