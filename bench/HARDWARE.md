# Reference hardware (NFR-01, amendment A-6)

The NFR-01 baseline is recorded per target triple on a dedicated machine, with `bench/run.sh`, which writes the
CPU model, core count, pinned core, frequency governor and kernel version next to the timings in
`bench/results/`. Until the reference machines are in place, results come from a shared aarch64 development
machine and are marked as preliminary.

| Target triple | Machine | CPU | Cores | Governor | Kernel | Recorded |
| --- | --- | --- | --- | --- | --- | --- |
| aarch64-unknown-linux-gnu | to be set up (dedicated) | | | | | |
| x86_64-unknown-linux-gnu | to be decided | | | | | |

## Preliminary results

| Date | Machine | CPU | Core | n = 16 | n = 18 | n = 20 |
| --- | --- | --- | --- | --- | --- | --- |
| 2026-10-30 | shared development machine (aarch64, 20 cores) | Cortex-X925, governor performance | 6 | 0.98 | 1.05 | 1.09 |

Ratios are our median time over thewalrus 0.22.0's fastest method (bbfg, glynn or ryser) on the same core;
NFR-01 requires at most 1.5.

## NFR-02 (preliminary)

`bench/nfr02.erl`: 100,000 boson-sampling shots, 6 photons in 12 modes, clifford_clifford, seed 1, through the job
server with p workers; efficiency `T1 / (p Tp)`, NFR-02 requires at least 0.8 up to the dirty CPU schedulers.

| Date | Machine | Cores | p = 1 | p = 2 | p = 4 | samples identical |
| --- | --- | --- | --- | --- | --- | --- |
| 2026-10-30 | shared development machine, 4 pinned cores | 4 | 1.17 s (1.00) | 0.61 s (0.95) | 0.33 s (0.88) | yes |
