# Quantum optics head-to-head: latest run

**Machine:** dgx-spark, Cortex-A725+Cortex-X925, one core (8). **Date:** 2026-11-14T2045Z. **PQCF:** c0131ad. **QuTiP:** 5.3.1.

How the tasks are set and checked: [README.md](README.md). Every run: [`results/`](results), [history.csv](history.csv).

## Against the closed forms

Correct decimal digits (the worst point of each task) and seconds (the median of repeated runs).

| task | PQCF (128 bits) | QuTiP (default) | QuTiP (tight) |
| --- | --- | --- | --- |
| rabi | 36 digits, 0.00015 s | 6 digits, 0.000114 s | 11 digits, 0.000284 s |
| revival | 36 digits, 0.0117 s | 5 digits, 0.00114 s | 9 digits, 0.00996 s |
| cavity_decay | 36 digits, 0.0137 s | 6 digits, 0.000387 s | 12 digits, 0.000955 s |
| fluorescence | 36 digits, 0.00051 s | 16 digits, 0.0373 s | 16 digits, 0.0375 s |

## Driven, damped Jaynes-Cummings model

No closed form: QuTiP's answers against PQCF's at 128 bits (the largest difference of <a^dagger a> and P_e at t = 5), and the seconds of one run. Sizes are the photon cutoff N; the Hilbert space has 2(N + 1) states and the Liouvillian 4(N + 1)^2.

| N | dimension | PQCF (128 bits) | QuTiP (default) | QuTiP (tight) |
| --- | --- | --- | --- | --- |
| 4 | 10 | 0.0244 s | 0.000714 s (differs by 1.4e-07) | 0.00204 s (differs by 1.4e-13) |
| 8 | 18 | 0.0936 s | 0.00145 s (differs by 8.0e-08) | 0.00478 s (differs by 3.7e-14) |
| 12 | 26 | 0.346 s | 0.00277 s (differs by 5.6e-08) | 0.00889 s (differs by 1.3e-14) |
| 16 | 34 | 1.08 s | 0.00456 s (differs by 7.6e-08) | 0.0145 s (differs by 6.3e-15) |
| 24 | 50 | 2.69 s | 0.01 s (differs by 1.1e-07) | 0.0277 s (differs by 7.8e-14) |
| 32 | 66 | 4.86 s | 0.0174 s (differs by 1.7e-07) | 0.0518 s (differs by 6.1e-14) |
