# API review before the freeze (Phase 3)

The public API, checked function by function against the interface summary and entry points of the functional
specification (sections 7, 7.1 and 9.4–9.6) and the symbolic engine's signatures (SY 8.1).

## Public modules

These modules are the API that will be frozen:

- `pq_api`
- `pq_circuit`
- `pq_symbolic`
- `pq_state`
- `pq_sampling`
- `pq_stats`
- `pq_noise`
- `pq_gradient`
- `pq_convert`
- `pq_parallel`: `run_job/2`, `cancel/1`, `status/1`, `register_node/0`, `job_opts/1`
- `pq_core`: configuration, the error taxonomy, `pack`/`unpack`

The records in `include/pq.hrl` are public as well: `#result{}` and the states. Build circuits through
`pq_circuit`.

Everything else is internal and may change: the `pq_sym_*` modules, `pq_job`, `pq_worker`, the supervisors,
`pq_node_registry`, `pq_spec_cache`, `pq_json`, `pq_telemetry` and the backend modules. `pq_simulator` and
`pq_backend` are public only to backend authors.

## Return shapes

Every data-dependent call returns `{ok, Value}` or `{error, Reason}`. The exceptions follow the specification:

- `pq_api:evolve` returns `{ok, State, Warnings}`.
- `probabilities` returns `{ok, Map, Info}`.

Plain values come back only where nothing can fail:

- `pq_circuit:add/2` and `to_binary/1`
- the `pq_state` constructors
- the gate helpers
- the term builders of `pq_symbolic`

## Found and fixed in the review

- **`pq_symbolic:evaluate/2,3` walked the expression tree.** SY 8.1 defines it as
  `eval_slp(compile([E], free_symbols(E)), B)`, so that evaluation and batch evaluation share one operation
  order. It now goes through the program. The tree evaluator remains as the test oracle (`tree => true`).
  `compile` of an identity matrix returned a list of entries instead of a matrix (fixed).
- **Terminal sampling ran in the caller's process.** This covered GBS, Fock photon counting and homodyne
  sampling without mid-circuit measurements, so these runs were never shot-parallel or distributed, contrary to
  FR-P01 and FR-P02. They now go through the job server in chunks, each chunk from its `stream_start`, so the
  samples are unchanged. Chunks return their metadata. `truncated_mass` (the sampler's mean over shots) is merged
  as a mean weighted by each chunk's shots, and `strict` is judged on the merged value. That value equals the
  one-call value up to the order of the floating-point sums.
- **The simulator passed only five job options to the job server.** `distributed`, `nodes`, `isolate_kernels`
  and `max_job_memory` were dropped for every circuit except boson sampling. There is now one list,
  `pq_parallel:job_opts/1`, for every caller.
- **Without the `pq` application started,** a sampling job runs inline as one chunk (the same samples) instead of
  failing with `noproc`.
- **Missing entry points** (FS 7, 7.1):
  - `pq_api:run/3`, `sample/3` and `evolve/3` (the backend first, as in `pq_simulator` and `pq_sampling`)
  - `pq_sampling:probabilities/3` and `postselect/2` (after the run, so no stream moves)
  - `pq_parallel:cancel/1`, `status/1` and `register_node/0`
  - `pq_circuit:two_mode_squeezing/3`
  - `pq_state:from_xpxp/1`
- **Job control** (FS 9.4):
  - `notify => Pid` sends `{pq_job, Job}` when the job starts; this is the handle for `cancel/1` and `status/1`.
  - `partial_results => true` turns a timeout into `{error, {timeout, #{partial => #{samples, shot_indices}}}}`.
    The rows are the full run's rows at those shot indices.
  - `auto_throttle => true` (FR-P08) lowers `max_inflight` to fit `max_job_memory` instead of refusing the job.
- **Telemetry** (FS 9.6): the events `[pq, job, start | stop | exception]` and `[pq, chunk, start | stop | retry]`,
  sent through `telemetry` when the application has it. It is a soft dependency, used only in tests here.

## Settled after the review

- **Cancelling chunks on other nodes** (FS 9.4): done. A proxy worker traps exits. When its job stops
  (cancelled, timed out or failed), it sets the chunk's cancel token on the remote node; each node keeps the
  tokens of its running chunks in the table `pq_remote_cancel`. Test: one chunk of 2,000,000 shots on a peer
  stops 11 ms after the cancel. Without the cancel, the chunk was still running 4 s later, which is the mutation
  check.
- **`kernel_hard_timeout`:** satisfied by construction. A cancelled or timed-out job replies at once and kills
  its local workers; it never waits for a kernel, so no kernel can hold it.

## Deviations, accepted (31 October 2026)

- **`ordered => false` (FR-P03, arrival-order merge) is not offered.** The samples of a run depend only on the
  circuit, the options and the seed, never on scheduling. An arrival-order merge would break that for no gain
  in throughput, because merging is not the bottleneck. Callers who want early samples have `notify`/`status`
  and `partial_results`.
- **`[pq, nif, call]` telemetry is not emitted.** The Fock backends call small kernels in loops, sometimes
  millions of calls per run, and an event per call would cost more than many of the calls themselves. The chunk
  events carry each chunk's duration, and a chunk is one kernel call for boson sampling and terminal sampling.

Both deviations were accepted by the project owner on 31 October 2026. The API is frozen with them.
