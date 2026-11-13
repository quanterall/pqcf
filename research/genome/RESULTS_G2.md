# G2: coherent-Ising-machine solvers against SA and PT, latest run

From [`results/g2-dgx-spark-2026-11-13T0744Z.json`](results/g2-dgx-spark-2026-11-13T0744Z.json), written by `g2_run.sh`. The method is in [README.md](README.md). TTS\* in seconds on one core: the expected time to reach the exact optimum with 99% probability, at the best run length. \* marks a best run length at the edge of the grid (an upper bound).

**Checks: 46 of 46 hold.**

## Haplotype phasing from short reads (2–6 sites per read; 3 instances per size)

The couplings form a band 5 sites wide, so a dynamic programme over the band solves it exactly in linear time (the *DP* column, one core).

| sites | instance | exact: CP-SAT, 8 workers | exact: DP, 1 core | SA | PT | SimCIM | AHC (CIM) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 16 | 16 | 4.07e-03 | 7.82e-04 | 3.79e-06 * | 3.35e-05 * | 9.68e-05 * | 4.64e-05 |
| 16 | 17 | 2.92e-03 | 6.82e-04 | 9.81e-06 * | 3.23e-05 * | 1.21e-04 * | 4.92e-05 |
| 16 | 18 | 2.84e-03 | 6.64e-04 | 9.68e-06 * | 3.17e-05 * | 1.71e-04 * | 4.75e-05 |
| 32 | 32 | 3.97e-03 | 1.32e-03 | 9.50e-05 * | 1.72e-04 * | 2.20e-03 * | 3.57e-04 |
| 32 | 33 | 4.39e-03 | 1.35e-03 | 6.90e-05 * | 6.97e-05 * | 3.31e-03 * | 1.03e-04 |
| 32 | 34 | 4.05e-03 | 1.42e-03 | 9.58e-05 | 8.61e-05 * | 1.70e-03 * | 4.18e-04 * |
| 64 | 64 | 0.012 | 2.75e-03 | 8.45e-04 * | 1.94e-03 * | 0.011 * | 1.08e-03 |
| 64 | 65 | 7.64e-03 | 2.73e-03 | 8.46e-04 * | 1.38e-03 * | 0.026 | 1.08e-03 |
| 64 | 66 | 5.75e-03 | 2.83e-03 | 5.52e-04 * | 8.29e-04 * | 6.37e-03 * | 2.37e-04 |
| 128 | 128 | 0.022 | 5.51e-03 | 0.046 | 0.270 | 0.068 * | 3.24e-03 |
| 128 | 129 | 0.019 | 5.53e-03 | 0.035 | 0.360 | 6.611 | 8.26e-03 |
| 128 | 130 | 0.024 | 5.56e-03 | 0.046 | 0.120 | 0.214 | 5.30e-03 |
| 256 | 256 | 0.088 | 0.011 | not reached | not reached | not reached | 0.209 |
| 256 | 257 | 0.073 | 0.011 | 6.608 | not reached | not reached | 0.319 |
| 256 | 258 | 0.107 | 0.011 | 6.552 | not reached | not reached | 0.068 |

## Haplotype phasing from long reads (10–40 sites per read; 3 instances per size)

The band is up to 39 sites wide: too wide for the dynamic programme above.

| sites | instance | exact: CP-SAT, 8 workers | switch errors of the optimum | SA | PT | SimCIM | AHC (CIM) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 64 | 64 | 0.051 | 3 | 1.24e-05 * | 9.29e-05 * | 3.90e-04 * | 1.93e-04 * |
| 64 | 65 | 0.022 | 0 | 1.23e-05 * | 9.18e-05 * | 3.67e-04 * | 1.70e-04 * |
| 64 | 66 | 0.036 | 3 | 1.22e-05 * | 9.44e-05 * | 3.72e-04 * | 1.75e-04 * |
| 128 | 128 | 0.272 | 6 | 7.19e-05 * | 2.99e-04 * | 2.45e-03 * | 3.77e-04 * |
| 128 | 129 | 0.312 | 7 | 1.33e-04 * | 3.07e-04 * | 1.98e-03 * | 3.61e-04 * |
| 128 | 130 | 0.238 | 1 | 5.75e-05 * | 3.07e-04 * | 3.63e-03 * | 3.95e-04 * |
| 256 | 256 | 1.903 | 2 | 3.80e-03 * | 1.16e-03 * | 0.052 * | 0.011 * |
| 256 | 257 | 1.661 | 3 | 1.23e-03 * | 1.21e-03 * | 0.024 * | 5.70e-03 * |
| 256 | 258 | 1.441 | 5 | 2.40e-03 * | 1.19e-03 * | 0.011 * | 4.65e-03 * |
| 512 | 512 | 23.5 | 0 | 0.072 | 8.15e-03 * | 0.470 | 0.030 * |
| 512 | 513 | 9.250 | 0 | 0.161 | 0.012 * | 0.337 | 0.059 |
| 512 | 514 | 23.5 | 4 | 0.036 * | 8.53e-03 * | 0.464 | 0.037 * |

## How the time grows (median over the 3 instances per size)

| family | method | median time at the three largest sizes | factor per doubling of the sites |
| --- | --- | --- | --- |
| short reads, 64 / 128 / 256 sites | exact (CP-SAT, 8 workers) | 7.64e-03 / 0.022 / 0.088 | 3.4 |
| short reads, 64 / 128 / 256 sites | SA TTS\* | 8.45e-04 / 0.046 / 6.608 | 88.4 |
| short reads, 64 / 128 / 256 sites | PT TTS\* | 1.38e-03 / 0.270 / not reached | — |
| short reads, 64 / 128 / 256 sites | SimCIM TTS\* | 0.011 / 0.214 / not reached | — |
| short reads, 64 / 128 / 256 sites | AHC (CIM) TTS\* | 1.08e-03 / 5.30e-03 / 0.209 | 13.9 |
| long reads, 128 / 256 / 512 sites | exact (CP-SAT, 8 workers) | 0.272 / 1.661 / 23.5 | 9.3 |
| long reads, 128 / 256 / 512 sites | SA TTS\* | 7.19e-05 / 2.40e-03 / 0.072 | 31.6 |
| long reads, 128 / 256 / 512 sites | PT TTS\* | 3.07e-04 / 1.19e-03 / 8.53e-03 | 5.3 |
| long reads, 128 / 256 / 512 sites | SimCIM TTS\* | 2.45e-03 / 0.024 / 0.464 | 13.8 |
| long reads, 128 / 256 / 512 sites | AHC (CIM) TTS\* | 3.77e-04 / 5.70e-03 / 0.037 | 9.9 |

## φX174 assembly, Boev et al.'s QUBO with one variable per overlap edge

| reads | native exact | SA | PT | SimCIM | AHC (CIM) |
| --- | --- | --- | --- | --- | --- |
| 6 | 8.11e-04 | 9.30e-05 | 2.01e-04 | 7.04e-04 | 7.41e-05 |
| 8 | 4.38e-04 | 1.00e-03 | 5.79e-04 | 0.016 | 1.53e-04 |
| 12 | 4.22e-04 | 7.86e-03 | 1.58e-03 | 1.324 | 2.89e-04 |
| 16 | 3.78e-04 | 0.033 * | 3.48e-03 | not reached | not reached |

## Tuning: the parameters chosen (best geometric-mean TTS\* on separate tuning instances)

| solver | short reads | long reads | assembly | grid |
| --- | --- | --- | --- | --- |
| SA | {'b0': 1.0, 'b1': 1.0} | {'b0': 1.0, 'b1': 0.5} | {'b0': 1.0, 'b1': 2.0} | 8 settings |
| PT | {'replicas': 8, 'b1': 2.0} | {'replicas': 8, 'b1': 4.0} | {'replicas': 8, 'b1': 4.0} | 8 settings |
| SimCIM | {'amp': 0.3, 'zeta': 1.0, 'sigma': 0.03} | {'amp': 0.3, 'zeta': 1.0, 'sigma': 0.03} | {'amp': 0.1, 'zeta': 1.0, 'sigma': 0.03} | 8 settings |
| AHC (CIM) | {'p': 0.0, 'eps': 0.3, 'beta': 0.1} | {'p': 0.0, 'eps': 0.3, 'beta': 0.4} | {'p': 0.5, 'eps': 0.3, 'beta': 0.4} | 8 settings |

## Every check

| check | holds | detail |
| --- | --- | --- |
| hap-16-16: the exact optimum is proven | yes | optimal, 0.00 s |
| hap-16-16: the banded dynamic programme gives the same optimum | yes | -164.0 in 0.001 s, band 5 |
| hap-16-17: the exact optimum is proven | yes | optimal, 0.00 s |
| hap-16-17: the banded dynamic programme gives the same optimum | yes | -258.0 in 0.001 s, band 5 |
| hap-16-18: the exact optimum is proven | yes | optimal, 0.00 s |
| hap-16-18: the banded dynamic programme gives the same optimum | yes | -221.0 in 0.001 s, band 5 |
| hap-32-32: the exact optimum is proven | yes | optimal, 0.00 s |
| hap-32-32: the banded dynamic programme gives the same optimum | yes | -380.0 in 0.001 s, band 5 |
| hap-32-33: the exact optimum is proven | yes | optimal, 0.00 s |
| hap-32-33: the banded dynamic programme gives the same optimum | yes | -433.0 in 0.001 s, band 5 |
| hap-32-34: the exact optimum is proven | yes | optimal, 0.00 s |
| hap-32-34: the banded dynamic programme gives the same optimum | yes | -533.0 in 0.001 s, band 5 |
| hap-64-64: the exact optimum is proven | yes | optimal, 0.01 s |
| hap-64-64: the banded dynamic programme gives the same optimum | yes | -921.0 in 0.003 s, band 5 |
| hap-64-65: the exact optimum is proven | yes | optimal, 0.01 s |
| hap-64-65: the banded dynamic programme gives the same optimum | yes | -872.0 in 0.003 s, band 5 |
| hap-64-66: the exact optimum is proven | yes | optimal, 0.01 s |
| hap-64-66: the banded dynamic programme gives the same optimum | yes | -1007.0 in 0.003 s, band 5 |
| hap-128-128: the exact optimum is proven | yes | optimal, 0.02 s |
| hap-128-128: the banded dynamic programme gives the same optimum | yes | -1794.0 in 0.006 s, band 5 |
| hap-128-129: the exact optimum is proven | yes | optimal, 0.02 s |
| hap-128-129: the banded dynamic programme gives the same optimum | yes | -1860.0 in 0.006 s, band 5 |
| hap-128-130: the exact optimum is proven | yes | optimal, 0.02 s |
| hap-128-130: the banded dynamic programme gives the same optimum | yes | -1682.0 in 0.006 s, band 5 |
| hap-256-256: the exact optimum is proven | yes | optimal, 0.09 s |
| hap-256-256: the banded dynamic programme gives the same optimum | yes | -3601.0 in 0.011 s, band 5 |
| hap-256-257: the exact optimum is proven | yes | optimal, 0.07 s |
| hap-256-257: the banded dynamic programme gives the same optimum | yes | -3486.0 in 0.011 s, band 5 |
| hap-256-258: the exact optimum is proven | yes | optimal, 0.11 s |
| hap-256-258: the banded dynamic programme gives the same optimum | yes | -3684.0 in 0.011 s, band 5 |
| long-64-64: the exact optimum is proven | yes | optimal, 0.05 s |
| long-64-65: the exact optimum is proven | yes | optimal, 0.02 s |
| long-64-66: the exact optimum is proven | yes | optimal, 0.04 s |
| long-128-128: the exact optimum is proven | yes | optimal, 0.27 s |
| long-128-129: the exact optimum is proven | yes | optimal, 0.31 s |
| long-128-130: the exact optimum is proven | yes | optimal, 0.24 s |
| long-256-256: the exact optimum is proven | yes | optimal, 1.90 s |
| long-256-257: the exact optimum is proven | yes | optimal, 1.66 s |
| long-256-258: the exact optimum is proven | yes | optimal, 1.44 s |
| long-512-512: the exact optimum is proven | yes | optimal, 23.49 s |
| long-512-513: the exact optimum is proven | yes | optimal, 9.25 s |
| long-512-514: the exact optimum is proven | yes | optimal, 23.48 s |
| asm-6-6: the exact optimum is proven | yes | optimal, 0.00 s |
| asm-8-8: the exact optimum is proven | yes | optimal, 0.00 s |
| asm-12-12: the exact optimum is proven | yes | optimal, 0.00 s |
| asm-16-16: the exact optimum is proven | yes | optimal, 0.00 s |
