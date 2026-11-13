# G1: Gaussian boson sampling for dense modules and cliques, latest run

From [`results/g1-dgx-spark-2026-11-13T0620Z.json`](results/g1-dgx-spark-2026-11-13T0620Z.json), written by `g1_run.sh`. The method is in [README.md](README.md). Every sampler gives post-selected k-vertex samples. *After s samples*: over 10 disjoint blocks of 200 samples, the fraction of blocks whose best sample so far is optimal.

**Checks: 28 of 28 hold.**

## Planted graph: 30 vertices, a dense 10-vertex subgraph of low-degree vertices (as Arrazola & Bromley's Fig. 3)

### Densest k-subgraph by random search: reached the optimum after s samples

| k | optimum (links) | sampler | mean density of a sample | s = 1 | 10 | 50 | 200 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 4 | 6 of 6 | GBS (PQCF) | 0.60 | 10% | 90% | 100% | 100% |
| 4 | 6 of 6 | Oh et al. | 0.51 | 0% | 40% | 100% | 100% |
| 4 | 6 of 6 | uniform | 0.27 | 0% | 10% | 40% | 80% |
| 6 | 15 of 15 | GBS (PQCF) | 0.51 | 10% | 60% | 100% | 100% |
| 6 | 15 of 15 | Oh et al. | 0.41 | 0% | 0% | 20% | 80% |
| 6 | 15 of 15 | uniform | 0.27 | 0% | 0% | 0% | 10% |
| 8 | 28 of 28 | GBS (PQCF) | 0.47 | 0% | 10% | 60% | 90% |
| 8 | 28 of 28 | Oh et al. | 0.37 | 0% | 0% | 0% | 0% |
| 8 | 28 of 28 | uniform | 0.27 | 0% | 0% | 0% | 0% |
| 10 | 44 of 45 | GBS (PQCF) | 0.42 | 0% | 0% | 10% | 30% |
| 10 | 44 of 45 | Oh et al. | 0.34 | 0% | 0% | 0% | 0% |
| 10 | 44 of 45 | uniform | 0.27 | 0% | 0% | 0% | 0% |

### Maximum clique (size 9) by shrink-and-search: found after s samples

| k | sampler | s = 1 | 10 | 50 | 200 |
| --- | --- | --- | --- | --- | --- |
| 4 | GBS (PQCF) | 10% | 90% | 100% | 100% |
| 4 | Oh et al. | 50% | 100% | 100% | 100% |
| 4 | uniform | 20% | 100% | 100% | 100% |
| 6 | GBS (PQCF) | 50% | 100% | 100% | 100% |
| 6 | Oh et al. | 40% | 80% | 100% | 100% |
| 6 | uniform | 30% | 100% | 100% | 100% |
| 8 | GBS (PQCF) | 60% | 100% | 100% | 100% |
| 8 | Oh et al. | 60% | 100% | 100% | 100% |
| 8 | uniform | 40% | 100% | 100% | 100% |
| 10 | GBS (PQCF) | 70% | 100% | 100% | 100% |
| 10 | Oh et al. | 70% | 100% | 100% | 100% |
| 10 | uniform | 70% | 100% | 100% | 100% |

**The planted subgraph itself, per sample:** GBS (PQCF) 0.20%, Oh et al. 0.00%, uniform 0.00%.

### What a sample costs

| k | GBS: simulated shots per kept sample | GBS: PQCF CPU-seconds per kept sample | GBS: mean photon number, cap | Oh et al.: seconds per sample | uniform: seconds per sample |
| --- | --- | --- | --- | --- | --- |
| 4 | 400 | 0.005 | 0.46, 10 | 1.2e-06 | 3.1e-06 |
| 6 | 400 | 0.011 | 0.63, 12 | 1.8e-06 | 3.2e-06 |
| 8 | 1,700 | 0.173 | 0.83, 14 | 2.9e-06 | 3.3e-06 |
| 10 | 10,500 | 2.929 | 0.98, 16 | 5.4e-06 | 3.4e-06 |

## STRING DNA-repair network: 45 genes, links with score ≥ 0.9

### Densest k-subgraph by random search: reached the optimum after s samples

| k | optimum (links) | sampler | mean density of a sample | s = 1 | 10 | 50 | 200 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 4 | 6 of 6 | GBS (PQCF) | 0.73 | 30% | 90% | 100% | 100% |
| 4 | 6 of 6 | Oh et al. | 0.65 | 10% | 80% | 100% | 100% |
| 4 | 6 of 6 | uniform | 0.33 | 0% | 10% | 50% | 90% |
| 6 | 15 of 15 | GBS (PQCF) | 0.69 | 0% | 30% | 80% | 100% |
| 6 | 15 of 15 | Oh et al. | 0.56 | 0% | 10% | 10% | 70% |
| 6 | 15 of 15 | uniform | 0.33 | 0% | 0% | 0% | 20% |
| 8 | 28 of 28 | GBS (PQCF) | 0.67 | 0% | 0% | 20% | 40% |
| 8 | 28 of 28 | Oh et al. | 0.52 | 0% | 0% | 0% | 10% |
| 8 | 28 of 28 | uniform | 0.33 | 0% | 0% | 0% | 0% |
| 10 | 45 of 45 | GBS (PQCF) | 0.65 | 0% | 0% | 0% | 0% |
| 10 | 45 of 45 | Oh et al. | 0.50 | 0% | 0% | 0% | 0% |
| 10 | 45 of 45 | uniform | 0.33 | 0% | 0% | 0% | 0% |

### Maximum clique (size 11) by shrink-and-search: found after s samples

| k | sampler | s = 1 | 10 | 50 | 200 |
| --- | --- | --- | --- | --- | --- |
| 4 | GBS (PQCF) | 10% | 50% | 100% | 100% |
| 4 | Oh et al. | 0% | 70% | 100% | 100% |
| 4 | uniform | 30% | 100% | 100% | 100% |
| 6 | GBS (PQCF) | 0% | 60% | 100% | 100% |
| 6 | Oh et al. | 10% | 90% | 100% | 100% |
| 6 | uniform | 10% | 60% | 100% | 100% |
| 8 | GBS (PQCF) | 0% | 90% | 100% | 100% |
| 8 | Oh et al. | 10% | 70% | 100% | 100% |
| 8 | uniform | 10% | 100% | 100% | 100% |
| 10 | GBS (PQCF) | 0% | 40% | 90% | 100% |
| 10 | Oh et al. | 0% | 70% | 100% | 100% |
| 10 | uniform | 10% | 90% | 100% | 100% |

### What a sample costs

| k | GBS: simulated shots per kept sample | GBS: PQCF CPU-seconds per kept sample | GBS: mean photon number, cap | Oh et al.: seconds per sample | uniform: seconds per sample |
| --- | --- | --- | --- | --- | --- |
| 4 | 400 | 0.007 | 0.37, 10 | 1.2e-06 | 3.1e-06 |
| 6 | 400 | 0.016 | 0.50, 12 | 1.6e-06 | 3.2e-06 |
| 8 | 1,600 | 0.152 | 0.62, 14 | 2.2e-06 | 3.2e-06 |
| 10 | 6,550 | 2.094 | 0.76, 16 | 3.9e-06 | 3.3e-06 |

## At equal cost: samples and time to reach the densest k-subgraph

GBS: the samples for a 50% chance, read from the curves above (so at their resolution: 1, 2, 5, 10, 20, 50, 100, 200), and PQCF's CPU time to simulate them. Oh et al. and uniform: the median number of samples until the optimum first appears, over 20 runs of up to 2,000,000 samples, and its time on one core. *Break-even rate*: the post-selected samples per second a GBS device would need to match Oh et al.'s sampler on one core.

| graph | k | GBS samples | Oh et al. samples | uniform samples | simulated GBS, CPU s | Oh et al., s | uniform, s | break-even rate |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| string | 4 | 2 | 5 | 62 | 0.0149 | 6e-06 | 0.0002 | 333,655 per s |
| string | 6 | 20 | 80 | 2,458 | 0.318 | 0.00013 | 0.0081 | 155,716 per s |
| string | 8 | > 200 (under 50% within the run's 200) | 1,760 | 321,208 | — | 0.0041 | 1.1 | — |
| string | 10 | > 200 (under 50% within the run's 200) | 359,189 | not reached | — | 1.3 | not reached | — |
| planted | 4 | 5 | 10 | 64 | 0.0231 | 1.1e-05 | 0.00021 | 450,275 per s |
| planted | 6 | 10 | 78 | 3,290 | 0.111 | 0.00014 | 0.011 | 74,007 per s |
| planted | 8 | 50 | 1,040 | 221,086 | 8.63 | 0.0029 | 0.76 | 17,531 per s |
| planted | 10 | > 200 (under 50% within the run's 200) | 20,254 | not reached | — | 0.11 | not reached | — |

## The densest modules of the DNA-repair network

| k | links | genes |
| --- | --- | --- |
| 4 | 6 of 6 | ATM, XRCC6, PRKDC, NHEJ1 |
| 6 | 15 of 15 | RAD51, RAD50, TP53, CHEK2, ATR, XRCC6 |
| 8 | 28 of 28 | BRCA1, BRCA2, ATM, CHEK2, XRCC6, H2AX, MDC1, PRKDC |
| 10 | 45 of 45 | BRCA1, RAD51, MRE11, RAD50, ATM, NBN, MSH6, ATR, H2AX, BARD1 |

## Every check

| check | holds | detail |
| --- | --- | --- |
| PQCF's GBS on a 10-vertex graph, k = 4: samples follow |haf(A_S)|^2 | yes | chi-square p = 0.401, 143 degrees of freedom, 2000 samples |
| Oh et al.'s sampler on a 10-vertex graph, k = 4: samples follow haf(A_S) | yes | chi-square p = 0.554, 143 degrees of freedom |
| PQCF's GBS on a 10-vertex graph, k = 6: samples follow |haf(A_S)|^2 | yes | chi-square p = 0.736, 138 degrees of freedom, 2000 samples |
| Oh et al.'s sampler on a 10-vertex graph, k = 6: samples follow haf(A_S) | yes | chi-square p = 0.127, 186 degrees of freedom |
| string, k = 4: the densest 4-subgraph is proven | yes | optimal |
| string, k = 4: PQCF's truncated probability is negligible | yes | largest per chunk 6.0e-05 |
| string, k = 4: enough GBS samples | yes | 2000 of 2000 |
| string, k = 6: the densest 6-subgraph is proven | yes | optimal |
| string, k = 6: PQCF's truncated probability is negligible | yes | largest per chunk 5.0e-05 |
| string, k = 6: enough GBS samples | yes | 2000 of 2000 |
| string, k = 8: the densest 8-subgraph is proven | yes | optimal |
| string, k = 8: PQCF's truncated probability is negligible | yes | largest per chunk 4.0e-05 |
| string, k = 8: enough GBS samples | yes | 2000 of 2000 |
| string, k = 10: the densest 10-subgraph is proven | yes | optimal |
| string, k = 10: PQCF's truncated probability is negligible | yes | largest per chunk 9.0e-05 |
| string, k = 10: enough GBS samples | yes | 2000 of 2000 |
| planted, k = 4: the densest 4-subgraph is proven | yes | optimal |
| planted, k = 4: PQCF's truncated probability is negligible | yes | largest per chunk 4.0e-05 |
| planted, k = 4: enough GBS samples | yes | 2000 of 2000 |
| planted, k = 6: the densest 6-subgraph is proven | yes | optimal |
| planted, k = 6: PQCF's truncated probability is negligible | yes | largest per chunk 3.0e-05 |
| planted, k = 6: enough GBS samples | yes | 2000 of 2000 |
| planted, k = 8: the densest 8-subgraph is proven | yes | optimal |
| planted, k = 8: PQCF's truncated probability is negligible | yes | largest per chunk 6.0e-05 |
| planted, k = 8: enough GBS samples | yes | 2000 of 2000 |
| planted, k = 10: the densest 10-subgraph is proven | yes | optimal |
| planted, k = 10: PQCF's truncated probability is negligible | yes | largest per chunk 7.0e-05 |
| planted, k = 10: enough GBS samples | yes | 2000 of 2000 |
