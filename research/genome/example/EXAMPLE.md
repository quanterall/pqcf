# Worked example: a virus, a bacterium and a human chromosome

Results of `run.sh`, on public genomes from NCBI. What each step means is in [README.md](README.md).

**Checks: 6 of 6 hold.**

## 1. Assembling genomes from reads

Simulated error-free reads, laid out by the exact solver on the overlap graph, then joined.

| genome | bases | reads | overlaps | layout (exact) | assembly = reference | QUBO for an annealer: n² encoding / edge encoding |
| --- | --- | --- | --- | --- | --- | --- |
| SARS-CoV-2 (NC_045512.2), whole genome | 29,903 | 599 × 300 bp | 2,862 | 0.05 s, optimal | yes | 358,801 / 2,862 variables |
| E. coli K-12 (NC_000913.3), 340,001-390,000 (the lac operon region) | 50,000 | 1,000 × 300 bp | 4,777 | 0.12 s, optimal | yes | 1,000,000 / 4,777 variables |

## 2. Phasing a simulated individual on human chromosome 21

chr21:30,000,001-31,000,000 (GRCh38): 648 heterozygous sites placed on the real sequence (about one per 1,500 bases), in a random phase; reads from either haplotype.

| reads | fragments joining sites | phase blocks | largest block (sites) | sites phased |
| --- | --- | --- | --- | --- |
| paired 2 x 150 bp, insert 500 bp, coverage 30, 0.1% errors | 1,715 | 110 | 6 | 268 |
| 15 kb, coverage 10, 1% errors | 666 | 1 | 644 | 644 |

The long-read phasing as an Ising problem:

| method | time to the optimum | phasing |
| --- | --- | --- |
| exact (CP-SAT, 8 workers) | 0.17 s (optimal) | 0 switch errors, 63 allele corrections |
| parallel tempering, G2's tuned settings | never reached | — |
| CIM-derived AHC, G2's tuned settings | never reached | — |

## 3. Dense modules in E. coli's DNA replication and repair network

STRING 12.5, 45 genes, 244 links with score ≥ 0.9; GBS simulated by PQCF.

| k | densest module | links | reached within 200 samples: GBS / Oh et al. / uniform | PQCF CPU-seconds per GBS sample |
| --- | --- | --- | --- | --- |
| 4 | recA, polA, dnaE, dnaB | 6 of 6 | 100% / 100% / 100% | 0.006 |
| 6 | dnaN, dnaE, dnaC, dnaA, dnaG, dnaX | 15 of 15 | 100% / 100% / 0% | 0.021 |
| 8 | recA, recQ, recR, ruvA, ruvB, recJ, recF, recO | 28 of 28 | 100% / 30% / 0% | 0.241 |

## 4. The genomes as quantum states

From G4 (`g4-dgx-spark-2026-11-13T1129Z.json`).

| genome | bases | qubits | Schmidt values needed (error ≤ 1e-6) | a random sequence of the same length |
| --- | --- | --- | --- | --- |
| SARS-CoV-2 | 29,903 | 17 | largest bond 234 | 234 |
| E. coli K-12 (middle cut) | 4,641,652 | 25 | 2,267 | (G4: the same) |
| human chromosome 21 (middle cut) | 40,088,619 | 28 | 9,787 | (every block of the chromosome) |

## Every check

| check | holds | detail |
| --- | --- | --- |
| SARS-CoV-2 (NC_045512.2), whole genome: the assembly equals the reference | yes | 29,903 bases |
| E. coli K-12 (NC_000913.3), 340,001-390,000 (the lac operon region): the assembly equals the reference | yes | 50,000 bases |
| chr21 region, long reads: the exact phasing is proven | yes | 0.2 s |
| E. coli network, k = 4: densest module proven, PQCF's truncation negligible | yes | optimal, truncated 0.0e+00 |
| E. coli network, k = 6: densest module proven, PQCF's truncation negligible | yes | optimal, truncated 2.0e-05 |
| E. coli network, k = 8: densest module proven, PQCF's truncation negligible | yes | optimal, truncated 2.0e-05 |
