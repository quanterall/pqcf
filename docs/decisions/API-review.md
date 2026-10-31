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

## Still open

- `ordered => false` (FR-P03): merging in arrival order, with `shot_indices`. Samples are always merged in shot
  order.
- `[pq, nif, call]` telemetry for each kernel call.
- Cancellation of chunks already running on another node. The job stops waiting and kills the local proxies, but
  the remote kernel runs to the end of its chunk, because cancel tokens are node-local resources.
- `kernel_hard_timeout`. A cancelled or timed-out job already returns without waiting for its kernels, so the
  option is not used yet.
