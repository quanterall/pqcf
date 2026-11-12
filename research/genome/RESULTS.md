# G0: classical baselines on genomic instances, latest run

From [`results/g0-dgx-spark-2026-11-12T2021Z.json`](results/g0-dgx-spark-2026-11-12T2021Z.json), written by `run.sh` (one core, one thread). The method is in [README.md](README.md). Times in seconds. TTS\* is the optimal time to solution: the shortest expected time to reach the exact optimum with 99% probability, over run lengths of 10 to 10,000 sweeps.

**Checks: 19 of 19 hold.**

## Haplotype phasing (Ising, from simulated reads: coverage 10, 2–6 sites per read, 5% errors)

| sites | Ising couplings | exact (CP-SAT) | switch errors of the optimum | SA TTS\* | PT TTS\* |
| --- | --- | --- | --- | --- | --- |
| 16 | 65 | 2.07e-03 (optimal) | 1 | 1.77e-05 (10 sweeps) | 6.92e-05 (10 sweeps) |
| 32 | 137 | 4.18e-03 (optimal) | 0 | 4.72e-04 (10 sweeps) | 1.68e-04 (10 sweeps) |
| 64 | 313 | 0.054 (optimal) | 2 | 2.53e-03 (100 sweeps) | 2.23e-03 (10 sweeps) |
| 128 | 611 | 1.009 (optimal) | 2 | 0.096 (1000 sweeps) | 0.063 (10 sweeps) |
| 256 | 1266 | 44.2 (optimal) | 0 | 14.1 (10000 sweeps) | 53.6 (10000 sweeps) |

## φX174 assembly: the layout of reads as a Hamiltonian path of largest overlap

Reads of 100 bp, coverage 4, from a random segment of the genome. *Native*: CP-SAT's circuit constraint on the overlap graph. *Positions*: Boev et al.'s QUBO with n² variables. *Edges*: one variable per overlap edge.

| reads | overlap edges | native exact | positions: variables | CP-SAT on it | SA TTS\* | PT TTS\* | edges: variables | CP-SAT on it | SA TTS\* | PT TTS\* |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 6 | 15 | 1.21e-03 | 36 | optimum, feasible at the 60.0 s limit | 1.29e-03 (316 sweeps) | 1.76e-03 (100 sweeps) | 15 | optimum, proven, 0.112 s | 2.07e-04 (316 sweeps) | 2.72e-04 (32 sweeps) |
| 8 | 25 | 4.99e-04 | 64 | **not reached**, feasible at the 60.0 s limit | 0.015 (1000 sweeps) | 0.015 (1000 sweeps) | 25 | optimum, feasible at the 60.0 s limit | 4.72e-03 (10000 sweeps) | 3.56e-03 (316 sweeps) |
| 12 | 37 | 5.62e-04 | 144 | **not reached**, feasible at the 60.0 s limit | 0.526 (3162 sweeps) | 0.473 (10000 sweeps) | 37 | optimum, feasible at the 60.0 s limit | 0.361 (10000 sweeps) | 0.109 (10000 sweeps) |
| 16 | 46 | 6.62e-04 | 256 | **not reached**, feasible at the 60.0 s limit | 20.3 (10000 sweeps) | 16.7 (3162 sweeps) | 46 | **not reached**, feasible at the 60.0 s limit | not reached | 0.619 (1000 sweeps) |
| 24 | 68 | 7.11e-04 | 576 | **not reached**, feasible at the 60.0 s limit | not reached | not reached | 68 | **not reached**, feasible at the 60.0 s limit | not reached | not reached |

## Dense modules in the STRING DNA-repair network (45 genes, links with score ≥ 0.9)

| k | links in the densest module | the module | exact (CP-SAT) | CP-SAT on the QUBO | SA TTS\* | PT TTS\* |
| --- | --- | --- | --- | --- | --- | --- |
| 3 | 3 of 3 | BRCA2, H2AX, MDC1 | 0.020 (optimal) | optimum, feasible at the 60.0 s limit | 9.03e-05 (100 sweeps) | 1.59e-04 (10 sweeps) |
| 5 | 10 of 10 | BRCA1, MSH2, RAD51, ATM, MSH6 | 0.035 (optimal) | **not reached**, feasible at the 60.0 s limit | 1.12e-03 (1000 sweeps) | 8.30e-04 (10 sweeps) |
| 8 | 28 of 28 | RAD51, RAD50, ATM, TP53, MSH6, ATR, XRCC6, H2AX | 2.454 (optimal) | optimum, feasible at the 60.0 s limit | 0.256 (10000 sweeps) | 0.118 (316 sweeps) |
| 10 | 45 of 45 | BRCA1, RAD51, MRE11, RAD50, ATM, TP53, CHEK2, ATR, H2AX, BARD1 | 2.437 (optimal) | optimum, feasible at the 60.0 s limit | not reached | not reached |
| 12 | 64 of 66 | BRCA1, RAD51, MRE11, RAD50, ATM, TP53, MSH6, CHEK2, ATR, XRCC6, H2AX, BARD1 | 250.4 (optimal; beyond the 60 s limit of the other rows, so not timed against the heuristics) | — | — | — |

## Every check

| check | holds | detail |
| --- | --- | --- |
| haplotype n=16: CP-SAT proves the optimum | yes | optimal |
| haplotype n=16: the optimum's phasing is close to the truth | yes | 1 switch errors, MEC 9 |
| haplotype n=32: CP-SAT proves the optimum | yes | optimal |
| haplotype n=32: the optimum's phasing is close to the truth | yes | 0 switch errors, MEC 17 |
| haplotype n=64: CP-SAT proves the optimum | yes | optimal |
| haplotype n=64: the optimum's phasing is close to the truth | yes | 2 switch errors, MEC 32 |
| haplotype n=128: CP-SAT proves the optimum | yes | optimal |
| haplotype n=128: the optimum's phasing is close to the truth | yes | 2 switch errors, MEC 68 |
| haplotype n=256: CP-SAT proves the optimum | yes | optimal |
| haplotype n=256: the optimum's phasing is close to the truth | yes | 0 switch errors, MEC 132 |
| assembly n=6: the native optimum is the reads' true order | yes | optimal, overlap 450 |
| assembly n=8: the native optimum is the reads' true order | yes | optimal, overlap 600 |
| assembly n=12: the native optimum is the reads' true order | yes | optimal, overlap 900 |
| assembly n=16: the native optimum is the reads' true order | yes | optimal, overlap 1200 |
| assembly n=24: the native optimum is the reads' true order | yes | optimal, overlap 1800 |
| dense 3-subgraph: CP-SAT proves the optimum | yes | optimal |
| dense 5-subgraph: CP-SAT proves the optimum | yes | optimal |
| dense 8-subgraph: CP-SAT proves the optimum | yes | optimal |
| dense 10-subgraph: CP-SAT proves the optimum | yes | optimal |
