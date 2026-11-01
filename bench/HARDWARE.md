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


## NFR-02 across nodes (provisional)

`bench/nfr02_nodes.erl`: the same 100,000 shots on 1 and on 4 peer nodes on the provisional reference (one host,
the 10 X925 cores; 2 dirty schedulers per node; the calling node only dispatches), best of 5; NFR-02 requires at
least 0.7 across 4 nodes. On one host it measures the distribution layer (erpc, the spec cache, merging), not a
network.

| Date | Instance | 1 node | 4 nodes | efficiency | samples identical |
| --- | --- | --- | --- | --- | --- |
| 2026-11-01 | 6 photons in 12 modes (Phase 1) | 0.11 s | 0.06 s | 0.45 | yes |
| 2026-11-01 | 10 photons in 20 modes | 1.76 s | 0.50 s | 0.88 | yes |

The same picture as on one node: the heavier instance meets the requirement, the Phase 1 one stops at the serial
cost of the list results in the caller.

## NFR-02: the benchmark instance (1 Nov 2026)

The project owner fixed **10 photons in 20 modes** as the NFR-02 benchmark (`bench/nfr02.erl`, `bench/nfr02_nodes.erl`
measures both). On the provisional reference: one node 0.97 / 0.92 / 0.83 at p = 2 / 4 / 8 (stable across runs);
across 4 nodes 0.88. At p = 10 -- all ten fast cores -- the result follows the shared machine's other load: 0.77 at
a load average near 3, 0.65 near 14. That last step, "up to the dirty CPU scheduler count", needs a machine without
other work to confirm.
