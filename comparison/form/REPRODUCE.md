# Reproducing the comparison with FORM

This page takes you from nothing to the numbers in [RESULTS.md](RESULTS.md): PQCF's symbolic engine against FORM 5 on
the same tasks, on your machine. What the tasks are and why: [README.md](README.md).

## What you need

- **Docker**, and about 5 GB of disk for two images: PQCF's development image (2.5 GB) and the FORM image (2.5 GB,
  FORM 5 and FLINT 3.2.2 compiled from their sources).
- **git**, and a network connection for the first build.
- Any Linux machine, x86_64 or aarch64. One core is used for the timings.

The images' base images are named by digest; the digests are in [REPRODUCE.md](../../REPRODUCE.md#base-images).

## Steps

**Reproduce the published run**, at the commit it was tested at (see "Tested" below):

```sh
git clone https://github.com/quanterall/pqcf.git
cd pqcf
git checkout 92e2aeb
COMPARISON_MACHINE=my-laptop CPUSET=0 sh comparison/form/run.sh
```

**Run the latest** version: the same commands without `git checkout 92e2aeb`. The answers must still be
the published ones; what changed since is in [CHANGELOG.md](../../CHANGELOG.md).

- `COMPARISON_MACHINE` is a label for your machine; it goes into the results.
- `CPUSET` is the core the tasks run on (leave it out to let Docker choose). Use a core with little else on it.

The script does four things:

1. builds the image `pqcf-dev` from the repository's [Dockerfile](../../Dockerfile) (Erlang/OTP 27, Rust, PQCF);
2. builds `pqcf-form` from [comparison/form/Dockerfile](Dockerfile): FORM 5 and FLINT from their sources, on top of
   `pqcf-dev`;
3. runs PQCF's FORM bridge tests (`pq_form_tests`, 15 tests): they must all pass;
4. runs every reference task three ways, one after another on the chosen core, and checks the answers.

On our machine the first run takes about 4 minutes (45 s for `pqcf-dev`, 90 s for FORM and FLINT, then the tasks),
plus downloading Docker's base images the first time; later runs take about 80 seconds.

## What you should see

The tests, then one line per task:

```
15 tests, 0 failures
T1 power n=12             0.014      0.210      0.008      18564 ok
T1 power n=16             0.055      1.089      0.030      74613 ok
...
T8 factor nontrivial      1.060      0.214      0.147         34 ok
T8 factor irreducibl      0.097      0.149      0.027         28 ok
T4 mesh k=4               0.007          -      0.007        400 ok
T4 mesh k=6               0.032          -      0.027      11616 ok
T4 mesh k=8               1.149          -      0.655     276800 ok
```

The columns are: the task; PQCF's time (native, seconds); the time through PQCF's FORM bridge; FORM's own time (best
of three); the number of terms of the answer; the check.

## How to compare with the published results

- **The answers** must be exactly the published ones: every line ends in `ok`, and the term counts (fifth column)
  are the same as in [RESULTS.md](RESULTS.md). `ok` means PQCF's own result and the result through FORM are
  identical, and all three term counts agree. A run with any other answer stops with an error.
- **The times** depend on your machine. Compare the ratio of PQCF's time to FORM's, not the seconds: the published
  ratios are in RESULTS.md and, for every run since 2 November 2026, in [history.csv](history.csv).

The script rewrites `RESULTS.md` and the charts in `charts/`, and appends your run to `history.csv`.
`git diff comparison/form` shows the differences; `git checkout comparison/form` puts the published files back.

## A quicker check

Only the bridge tests, which compare PQCF's results with FORM's on a set of expressions (a minute after the images
are built):

```sh
docker build -t pqcf-dev . && docker build -t pqcf-form -f comparison/form/Dockerfile .
docker run --rm pqcf-form rebar3 eunit --module=pq_form_tests
```

## If something goes wrong

- **The FORM build fails:** it needs the network (FORM and FLINT are downloaded at fixed versions: FORM at commit fb12f2e, the 5.0.2 of the published runs, and FLINT 3.2.2). Retry; the error
  is in Docker's output.
- **A check is not `ok`:** that is a real disagreement. Please open an issue with your machine and the output:
  <https://github.com/quanterall/pqcf/issues>.

## Tested

From a fresh clone on 15 November 2026 at commit 92e2aeb (the published run above), aarch64 (Cortex-X925
and Cortex-A725): every check `ok`, 48 seconds with the images built.

Earlier, from a fresh clone of the public repository on 15 November 2026 (commit d77ead6), on an aarch64 machine (Cortex-X925
and Cortex-A725): every check `ok`, 81 seconds with the images built; the images built without cache in 45 s and 87 s.
