# Reproducing the reduction comparison with FIRE

This page takes you from nothing to the comparison of PQCF's Feynman-integral reduction with FIRE 6 (A. V. Smirnov
and F. S. Chuharev, [arXiv:1901.07808](https://arxiv.org/abs/1901.07808)), on your machine: the same integrals reduced
by both, on the same core, the answers compared exactly. The method is decision record
[R1](../../docs/decisions/R1-finite-field-reconstruction.md).

## What you need

- **Docker**, and about 5 GB of disk: PQCF's development image (2.5 GB) and the FIRE image (2.1 GB, FIRE 6 with the
  FLINT simplifier, compiled from its sources).
- **git**, **python3** (it writes the summary, standard library only) and a network connection for the first build.
- Any Linux machine, x86_64 or aarch64.

## Steps

**Reproduce the published run**, at the commit it was tested at (see "Tested" below):

```sh
git clone https://github.com/quanterall/pqcf.git
cd pqcf
git checkout 92e2aeb
docker build -t pqcf-dev .
docker build -t pqcf-fire2 comparison/fire
PQ_BENCH_HOST=my-laptop bench/ibp.sh 0
```

**Run the latest** version: the same commands without `git checkout 92e2aeb`. The answers must still be
the published ones; what changed since is in [CHANGELOG.md](../../CHANGELOG.md).

- `0` is the core the timed runs use; pick one with little else on it.
- `PQ_BENCH_HOST` is a label for your machine in the results file.
- The FIRE image builds FIRE at the commit of our published runs, with FIRE's own pins for its submodules; on
  aarch64 it is compiled with `-fsigned-char` (FIRE assumes a signed `char`, as on x86).

`bench/ibp.sh` takes the committed version of the repository (`git archive HEAD`) and:

1. **PQCF** reduces every case of the suite ([bench/ibp_suite.erl](../../bench/ibp_suite.erl)) on the chosen core,
   then on all cores;
2. **FIRE** reduces the cases it ships as its own examples (the box, the two-loop vertex, the double box at three
   numerator settings), one thread, on the same core, timed; then it reduces PQCF's master integrals into its own
   basis, untimed;
3. the answers are **compared exactly** at a random point modulo a large prime, master by master in FIRE's basis
   ([fire_check.erl](fire_check.erl)): the programs may choose different masters, and FIRE's examples do not apply
   every symmetry (the double box has 12 masters in FIRE, 8 up to symmetries, as in PQCF).

On our machine the images build in about 6 minutes (4.5 for FIRE) and the suite runs in 90 seconds.

## What you should see

```
fire_box         FIRE 3 masters (2 up to symmetries), PQCF 2: equal
fire_v2          FIRE 5 masters (3 up to symmetries), PQCF 2: equal
doublebox_1      FIRE 12 masters (8 up to symmetries), PQCF 8: equal
doublebox_2      FIRE 12 masters (8 up to symmetries), PQCF 8: equal
doublebox_3      FIRE 12 masters (8 up to symmetries), PQCF 8: equal
case                1 core      all     FIRE  ratio  adds  check   change
fire_box            0.133s   0.075s   0.374s   0.35    24  equal
fire_v2             0.039s   0.037s   0.505s   0.08     4  equal
doublebox_1         0.317s   0.167s   7.568s   0.04    18  equal
doublebox_2         2.293s   1.051s  11.658s   0.20    56  equal
doublebox_3        20.876s   6.287s  12.865s   1.62   102  equal
kite                0.004s   0.004s        -      -     5  -
massive_tadpole     0.002s   0.002s        -      -     3  -
banana3             0.025s   0.021s        -      -     2  -
written bench/results/ibp-my-laptop-<date>.json
```

The columns: PQCF's time on one core and on all cores, FIRE's time on the same core, their ratio (below 1: PQCF is
faster), the number of master-integral coefficients PQCF needed (`adds`), and the check. The last three cases are
PQCF's own (FIRE has no example for them).

## How to compare with the published results

- **The answers:** every case with a FIRE example says `equal`. The script exits with 1 otherwise.
- **The times** depend on the machine; compare the ratio column with the published runs in
  [bench/results/](../../bench/results) (`ibp-<host>-<date>.json`), and with the table on quanterall.com. A second
  run on the same machine also prints the change against the previous one.
- The machine's load average is recorded in the results file: the timings are only as good as the core is quiet.

## If something goes wrong

- **`pqcf-fire2` not found:** build it first (the second `docker build` above).
- **A check says `DIFFER`:** a real disagreement; please open an issue with the output:
  <https://github.com/quanterall/pqcf/issues>.

## Tested

From a fresh clone on 15 November 2026 at commit 92e2aeb (the published run above), aarch64 (Cortex-X925
and Cortex-A725): every check `equal`, 97 seconds.

Earlier, from a fresh clone of the public repository on 15 November 2026 (commit d77ead6, with the FIRE image built from the
pinned Dockerfile of this page's commit), aarch64 (Cortex-X925 and Cortex-A725), core 18: every check `equal`,
90 seconds.
