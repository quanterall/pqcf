# Reproducing the four-loop comparison with Forcer

This page takes you from nothing to the results in [README.md](README.md): PQCF reduces the four-loop integrals of
FORM's own Forcer checks, and its answers are compared with Forcer's and with FORM's expected results, exactly.

## What you need

- **Docker**, and about 8 GB of disk for three images: PQCF's development image (2.5 GB), the FORM image
  (2.5 GB, FORM 5 and FLINT compiled from their sources) and the Forcer image (FORM plus Forcer, 2.5 GB, sharing most
  of its layers with the FORM image).
- **git**, and a network connection for the first build.
- Any Linux machine, x86_64 or aarch64. More cores make PQCF's reductions faster; 8 is plenty.

## Steps

**Reproduce the published run**, at the commit it was tested at (see "Tested" below):

```sh
git clone https://github.com/quanterall/pqcf.git
cd pqcf
git checkout 92e2aeb
docker build -t pqcf-dev .
docker build -t pqcf-form -f comparison/form/Dockerfile .
PQ_BENCH_HOST=my-laptop sh comparison/forcer/run.sh 0-7
```

**Run the latest** version: the same commands without `git checkout 92e2aeb`. The answers must still be
the published ones; what changed since is in [CHANGELOG.md](../../CHANGELOG.md).

- The two `docker build` lines are needed once: `run.sh` builds the Forcer image itself, on top of `pqcf-form`.
- `0-7` are the cores PQCF uses (a Docker cpuset: `0-7`, `0,2,4`, ...); Forcer runs on the first of them.
- `PQ_BENCH_HOST` is a label for your machine in the results file.

`run.sh` takes the committed version of the repository (`git archive HEAD`, so local edits are not measured) and:

1. builds `pqcf-forcer` from [Dockerfile](Dockerfile): the FORM image plus Forcer at a fixed commit and FORM's own
   checks `forcer.frm` and `mincer.frm` at a fixed commit;
2. **PQCF reduces** every case: each integral of FORM's checks to master integrals, with exact coefficients;
3. **Forcer** evaluates PQCF's master integrals and each target integral, and FORM runs its own checks;
4. **compares** exactly: PQCF's reduction times Forcer's master values must equal Forcer's result for each target, as
   expansions in ε; and the MS-bar sums must equal the results FORM's checks assert.

On our machine the images build in about 4 minutes and the run takes about 7 minutes on 8 cores.

## What you should see

PQCF's reductions, one line per case:

```
example_3loop        0.04 s      181 equations      79 needed    1 masters  seeds {dots,0,1}
...
zno1_f2             48.81 s   116202 equations   67035 needed   21 masters  seeds {dots,1,1}
```

then the checks, every one `equal`:

```
example_3loop    Forcer's result, G-scheme: equal (compared to O(ep^7))
...
zno2_f2n         Forcer's result, G-scheme: equal (compared to O(ep^3))
Forcer_example F, MS-bar, against FORM's assertion: equal
...
Forcer_1 F2, MS-bar, against FORM's assertion: equal
  ep^-4: 15/2
  ep^-3: 383/12
  ep^-2: -2089/18
  ep^-1: 1466/9 + -728*z3 + -1245/8*z5
  ...
```

and finally `written comparison/forcer/results/<host>-<date>.json`. The script exits with 1 if any check differs.

## How to compare with the published results

- **The answers:** all 10 components `equal` to Forcer's results and all 5 MS-bar results `equal` to FORM's
  assertions, as in [README.md](README.md). The coefficients printed (for example `Forcer_1 F2` above) are the
  published ones, with ζ-values as `z3`, `z5`, ...
- **PQCF's own results:** the run rewrites [`results/symbolic/`](results/symbolic), every reduction, expansion and
  MS-bar result. `git status` must show them unchanged: only your new `results/<host>-<date>.json` is new.
- **The times** are in your results file, next to the published runs' in [`results/`](results); they depend on the
  machine.

## A quicker check

There is none shorter than the full run: every case is part of one of FORM's checks, and a check needs all of its
cases. The full run takes about 7 minutes on 8 cores once the images are built.

## If something goes wrong

- **`pqcf-form` not found:** build it first (the second `docker build` above).
- **A check says `DIFFER`:** that is a real disagreement; please open an issue with the output:
  <https://github.com/quanterall/pqcf/issues>.

## Tested

From a fresh clone on 15 November 2026 at commit 92e2aeb (the published run above), aarch64 (Cortex-X925
and Cortex-A725), with nothing else using the comparison images: 15 of 15 checks `equal`, `results/symbolic/`
unchanged, 465 seconds.

Earlier, from a fresh clone of the public repository on 15 November 2026 (commit d77ead6), aarch64 (Cortex-X925 and
Cortex-A725), cores 0-7: every check `equal`, `results/symbolic/` unchanged byte for byte, 400 seconds.
