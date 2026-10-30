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
