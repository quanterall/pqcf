# Reproducing the genome research

This page takes you from nothing to the results of the genome research on this branch: the worked example on public
genomes ([example/EXAMPLE.md](example/EXAMPLE.md)) and the five studies G0-G4 ([README.md](README.md)), on your
machine.

## What you need

- **Docker**, and about 4 GB of disk: PQCF's development image (2.5 GB) and the genome image
  ([docker/Dockerfile](docker/Dockerfile): Python with the solvers and samplers, pinned).
- **git**, **curl** and a network connection: the large genomes (E. coli K-12, NC_000913.3, and human chromosome 21,
  NC_000021.9, about 50 MB) come from NCBI and are checked by SHA-256 ([data/SHA256SUMS](data/SHA256SUMS)).
- Any Linux machine, x86_64 or aarch64. The example and G1 use 8 cores; G0, G2 and G3 one core each.

## Steps

The research lives on the `research/genome` branch:

```sh
git clone -b research/genome https://github.com/quanterall/pqcf.git
cd pqcf
git checkout bb8a5dd
docker build -t pqcf-dev .
sh research/genome/data/fetch.sh
```

`git checkout bb8a5dd` **reproduces the published runs**: the commit of this page, with the code they were tested
at (see "Tested"). To **run the latest** version of the branch, leave it out; the results must still be the published
ones.

`fetch.sh` downloads the two genomes into `research/genome/data/` and checks them; it prints `OK` twice. The smaller
genomes (SARS-CoV-2 and the others) are in the repository already.

### The worked example (start here)

```sh
sh research/genome/example/run.sh "0 1 2 3 4 5 6 7"
```

The argument is the list of cores (8 of them here; the default suits our machine only). In about 6 minutes on 8
cores it runs the four steps of [example/README.md](example/README.md) on the public genomes:

1. **assembly** of SARS-CoV-2 and an E. coli region from simulated reads (overlap graph, exact layout);
2. **phasing** of a simulated individual on a region of human chromosome 21, from short and from long reads;
3. **dense modules** of E. coli's DNA replication and repair network, by PQCF's simulated Gaussian boson sampling
   device against a quantum-inspired sampler and uniform sampling;
4. the genomes **as quantum states** (matrix product states).

It prints one line per step and ends with:

```
6 of 6 checks hold
written example.json and EXAMPLE.md
```

### The studies G0-G4

Each writes `research/genome/results/<study>-<host>-<date>.json` and its results page. Give `PQ_BENCH_HOST=<label>`
to name your machine in the files.

| study | command | what it does | time on our machine |
| --- | --- | --- | --- |
| G0 | `sh research/genome/run.sh 0` | the classical baselines on the genomic instances, one core | about 45 minutes |
| G1 | `sh research/genome/g1_run.sh "0 1 2 3 4 5 6 7"` | Gaussian boson sampling (PQCF) against Oh et al.'s quantum-inspired sampler and uniform sampling | about 30 minutes on 8 cores |
| G2 | `sh research/genome/g2_run.sh "0 1 2 3"` | coherent-Ising-machine solvers against simulated annealing and parallel tempering, by time to solution | about 40 minutes |
| G3 | `sh research/genome/g3_run.sh 0` | the GBS graph kernel against classical features, telling organisms apart | about 1.5 minutes |
| G4 | `sh research/genome/g4_run.sh "0-3"` | genomes as matrix product states: the bond dimension they need | about 4.5 minutes on 4 cores |

G0 writes only the JSON; its page is made by
`python3 research/genome/report.py research/genome/results/<run>.json > research/genome/RESULTS.md`.

## How to compare with the published results

- **The worked example:** `git diff research/genome/example/EXAMPLE.md`. Every result must be the same: the
  assemblies equal the reference, the same phasing errors and corrections, the same modules found and the same
  success rates. Only the times in seconds change.
- **The studies:** compare your `RESULTS*.md` with the published ones in the same way: the instances, the solutions
  and the checks are fixed by seeds; the times (and so TTS\*, the time to solution) depend on the machine.
- **The conclusion** of the research, that no photonic method has a practical advantage on these problems, rests on
  the comparisons, not on one machine's seconds: a slower or faster machine moves every program together.

## If something goes wrong

- **`fetch.sh` reports a checksum mismatch:** NCBI served a different file; delete it and run `fetch.sh` again.
- **A check does not hold:** please open an issue with the output: <https://github.com/quanterall/pqcf/issues>.

## Tested

From a fresh clone of the public `research/genome` branch on 15 November 2026 (commit d4c3e4f), aarch64 (Cortex-X925
and Cortex-A725): the worked example on cores 0-7, 6 of 6 checks, every result equal to the published one (only the
times changed), 336 seconds; G3 on one core (11 of 11 checks, 80 s) and G4 on four (2 of 2 checks, 265 s), each with only its
times changed. G0, G1 and G2 were not rerun for this page; their times are those of the published runs.
