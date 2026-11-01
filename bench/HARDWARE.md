# Reference hardware (NFR-01, amendment A-6)

The NFR-01 baseline is recorded per target triple on a dedicated machine, with `bench/run.sh`, which writes the
CPU model, core count, pinned core, frequency governor and kernel version next to the timings in
`bench/results/`. Until the reference machines are in place, results come from a shared aarch64 development
machine and are marked as preliminary.

| Target triple | Machine | CPU | Cores | Governor | Kernel | Recorded |
| --- | --- | --- | --- | --- | --- | --- |
| aarch64-unknown-linux-gnu | dev-aarch64-shared (provisional, project owner's decision 31 Oct 2026; shared, not dedicated) | 10 Cortex-X925 + 10 Cortex-A725 | 20 | performance | 6.17 | 2026-10-31 |
| x86_64-unknown-linux-gnu | to be decided | | | | | |

## Preliminary results

| Date | Machine | CPU | Core | n = 16 | n = 18 | n = 20 |
| --- | --- | --- | --- | --- | --- | --- |
| 2026-10-30 | shared development machine (aarch64, 20 cores) | Cortex-X925, governor performance | 6 | 0.98 | 1.05 | 1.09 |
| 2026-10-31 | the same, now the provisional reference; three runs, identical to two decimals | Cortex-X925, governor performance | 8 | 0.98 | 1.06 | 1.11 |

Ratios are our median time over thewalrus 0.22.0's fastest method (bbfg, glynn or ryser) on the same core;
NFR-01 requires at most 1.5.

## NFR-02 (preliminary)

`bench/nfr02.erl`: 100,000 boson-sampling shots, 6 photons in 12 modes, clifford_clifford, seed 1, through the job
server with p workers; efficiency `T1 / (p Tp)`, NFR-02 requires at least 0.8 up to the dirty CPU schedulers.

| Date | Machine | Cores | p = 1 | p = 2 | p = 4 | samples identical |
| --- | --- | --- | --- | --- | --- | --- |
| 2026-10-30 | shared development machine, 4 pinned cores | 4 | 1.17 s (1.00) | 0.61 s (0.95) | 0.33 s (0.88) | yes |
| 2026-10-31 | the provisional reference, the 10 X925 cores (cpuset), best of 5 per p | 10 | 0.19 s (1.00) | 0.10 s (0.94) | 0.06 s (0.83) | yes |

On 31 Oct the configuration above **fails NFR-02 beyond p = 4** (p = 8: 0.68, p = 10: 0.58 in the best runs; all 20
mixed cores: 0.19 at p = 20). Since Phase 1 the Clifford-Clifford kernel became about 4x faster (1.17 s -> 0.18 s for
the 100,000 shots), and what is left is serial work in the calling process: turning the samples into Erlang lists,
1.2 million cells (about 14-18 ms, half of it garbage collection; the kernels themselves scale: 10 bare processes
18.8 ms against 17.8 ms for one). Rows now cross processes as packed binaries and become lists once, in the
caller (no list copied between processes), which brought T1 from 0.27 s to 0.19 s and p = 4 from 0.54 to 0.83.

The same benchmark with a heavier instance -- 10 photons in 20 modes, the rest unchanged -- where the kernel's work
per shot grows faster than the output: 3.43 s at p = 1, efficiency 0.97 (p = 2), 0.92 (4), 0.83 (8), 0.77 (10), on
cores another workload keeps 20-60 % busy. Open for the project owner: keep the Phase 1 configuration (and accept
the serial floor of list results) or fix the heavier one as the NFR-02 benchmark.

