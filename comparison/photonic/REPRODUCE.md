# Reproducing the photonic head-to-head

This page takes you from nothing to the tables in [RESULTS.md](RESULTS.md): PQCF against five photonic simulators
(The Walrus, Piquasso, Perceval, DeepQuantum, Strawberry Fields) on the same inputs, on one core of your machine.
What the tasks are, and what the results mean: [README.md](README.md).

## What you need

- **Docker**, and about 5 GB of disk for three images: PQCF's development image (2.5 GB); the peers' image (1.9 GB:
  The Walrus 0.22.0, Piquasso 8.0.2, Perceval 1.3.1, DeepQuantum 4.5.0 with PyTorch 2.9.1 for the CPU); and
  Strawberry Fields' (0.6 GB: 0.23.0, which needs NumPy 1.x, so it has an image of its own).
- **git**, and a network connection for the first build.
- Any Linux machine, x86_64 or aarch64. Every program runs on one core, with one thread.

The images' base images are named by digest; the digests are in [REPRODUCE.md](../../REPRODUCE.md#base-images).

## Steps

**Reproduce the published run**, at the commit it was tested at (see "Tested" below):

```sh
git clone https://github.com/quanterall/pqcf.git
cd pqcf
git checkout 92e2aeb
PQ_BENCH_HOST=my-laptop sh comparison/photonic/run.sh 0
```

**Run the latest** version: the same commands without `git checkout 92e2aeb`. The answers must still be
the published ones; what changed since is in [CHANGELOG.md](../../CHANGELOG.md).

- `0` is the core every program runs on, one after another; pick one with little else on it (default 8).
- `PQ_BENCH_HOST` is a label for your machine in the results.

`run.sh` builds the three images (pinned versions: [peers/Dockerfile](peers/Dockerfile), [sf/Dockerfile](sf/Dockerfile)),
takes the committed version of the repository (`git archive HEAD`), and:

1. **writes the inputs** from fixed seeds ([gen_inputs.py](gen_inputs.py)): random unitaries, Gaussian states, the
   boson-sampling circuits and the 10-mode continuous-variable network;
2. **runs PQCF** ([photonic_pqcf.erl](photonic_pqcf.erl)), then **each peer** ([peers.py](peers.py)), each pinned to
   the core with one thread; a size that takes more than 120 s is the last one run for that program;
3. **checks the answers against each other** and writes the tables ([compare.py](compare.py)): permanents and
   torontonians within 1e-8 and 1e-7, every program's boson samples against the exact distribution by a
   chi-square test, the Fock states against each other.

On our machine the run takes 13 minutes once the images are built; building them the first time adds about 4.5
minutes (3 of them for the peers' image, where Piquasso is compiled from source).

## What you should see

```
PQCF done
thewalrus done
piquasso done
perceval done
deepquantum done
strawberryfields done
50 of 50 checks agree
```

(with some deprecation warnings from the peers in between). The script exits with 1 if any check fails.

## How to compare with the published results

- **The answers:** every check agrees: `50 of 50 checks agree` in the published [RESULTS.md](RESULTS.md) (the run of
  15 November). The number of checks depends on how many sizes each program reaches within the 120 s limit on your
  machine; the summary names every size that was not run, for example
  `not run (time limit or capacity): ..., DeepQuantum torontonian 24, ...`. The first published run (12 November) had
  51 checks: there DeepQuantum's torontonian at k = 22 took less than 120 s, so k = 24 ran too.
- **Incomplete runs fail:** if a program produces no results, or no answer for a task it covers, `compare.py` prints
  `INCOMPLETE: ...` and the script exits with 1, whatever the other checks say.
- **The times:** the run rewrites `RESULTS.md` with your tables, appends to [history.csv](history.csv) and writes
  `results/<host>-<date>.json`. The seconds depend on your machine; compare each program's time relative to PQCF's
  (the factor in brackets in every table), and PQCF's place in each task, with [README.md](README.md)'s summary.
  `git diff comparison/photonic/RESULTS.md` shows your tables against the published ones; `git checkout
  comparison/photonic` puts them back.
- **On aarch64,** Piquasso is built with signed `char`: its released source gives wrong permanents on ARM otherwise
  ([README.md](README.md), "A bug found on the way"). On x86 the flag changes nothing.

## A quicker check

One program and one task, for example the permanents of PQCF and The Walrus, after the images are built (the full
run is the comparison; this only shows that the pieces work). It takes about 2 minutes and writes `inputs.json` and
`thewalrus.json`, The Walrus's times and answers for every size:

```sh
docker run --rm -v "$PWD/comparison/photonic":/p -v /tmp/ph:/out pqcf-photonic-peers \
  sh -c 'python /p/gen_inputs.py > /out/inputs.json && python /p/peers.py thewalrus /out && ls /out'
```

## If something goes wrong

- **A program fails:** the run goes on and says `<program> failed: the run is incomplete`; the others are still
  compared and tabled, but the run fails (exit 1). The versions are pinned, so this should not happen; please open an
  issue.
- **A check fails:** please open an issue with the output: <https://github.com/quanterall/pqcf/issues>.

## Tested

From a fresh clone on 15 November 2026 at commit 92e2aeb (the published run above), aarch64 (Cortex-X925
and Cortex-A725): 50 of 50 checks agree, the sizes not run named in the summary, 835 seconds.

Earlier, from a fresh clone of the public repository on 15 November 2026 (commit d77ead6), aarch64 (Cortex-X925 and
Cortex-A725), core 19: 50 of 50 checks agree, 788 seconds.
