# P1: the task behaviour, PQCF as a platform for applications

Status: accepted, 13 November 2026. Related: R1 (finite-field reconstruction: the executors), D1 (results
independent of how work is spread).

## Context

PQCF's executors (`pq_rexec`: this node, or one owner's Erlang nodes; the FVLLMP federation to come) were written for
one kind of work, the finite-field reconstruction's tasks. GenHub2, a virus-genomics application planned on the
FVLLMP federation, is to be PQCF's first application, and needs the same executors for its own work: aligning reads,
classifying genomes, building trees. Its work must not be written into PQCF's core, so that the next application
needs nothing from GenHub2.

## Decision

1. **A task is `{pq_task, Module, Args}`,** where `Module` implements the `pq_task` behaviour:
   - `run(Args, #{home := Node}) -> Result`, required. The result is a pure function of the task, so it is the same
     on any node and however the work is spread (D1). `home` is the coordinating node, whose store (`pq_rstore`)
     holds the tasks' data.
   - `resources(Args) -> #{cpu, gpu, memory_mb, data}`, optional: what the task needs, for scheduling (to come).
   - `describe(Args) -> map()`, optional: what the task is, for its provenance record.
2. **The reconstruction's tuples become one such module,** `pq_ff_task`, moved unchanged from `pq_rexec`. The bare
   tuples are still accepted and mean `{pq_task, pq_ff_task, Tuple}`, so `pq_recon` and its results are unchanged.
3. **A node runs only the task modules it allows:**
   - PQCF's own (`pq_ff_task`);
   - those in the application environment (`{pq, [{task_modules, [...]}]}`);
   - those registered with `pq_task:allow/1`.

   Anything else fails with `{task_not_allowed, Module}` on that node. An executor never runs code it was not told
   to, which the federation, where nodes belong to different owners, requires.
4. **Provenance:** `pq_rexec:map_traced/3` returns each result with a record:
   - the task's SHA-256 (of its deterministic external term);
   - the module, and the MD5 of its loaded code;
   - its own description;
   - the node that ran it, PQCF's version;
   - the start time and the wall time.

   `pq_rexec:map/3` is unchanged and returns the results alone.

## Consequences

- **An application adds work by writing task modules and allowing them** on its nodes. The executors, and the
  federation executor when it comes, run them as they run PQCF's own.
- **Results do not depend on tracing:** `map_traced/3`'s results are `map/3`'s. Only the records carry times.
- **Allowing a module is per node.** A peer that has not allowed a module refuses its tasks, and the error is raised
  at the coordinator; it is not retried elsewhere.
- **To come:** scheduling by `resources/1` and by where the data is, and the federation executor (GA4GH TES through
  the FVLLMP gate).

## Tests

`test/pq_task_tests.erl`:
- a test module refused until allowed, then run in task order with 1 and 7 workers;
- the reconstruction's tuples giving the same results bare, wrapped and traced;
- the provenance record's fields, with the task hash stable for one task and different for another;
- on two peer nodes, refused until each peer allows the module, then equal to this node's results, with every peer
  named in the records.

`test/pq_rexec_tests.erl` (unchanged) passes, including the reconstruction across peers with a node lost mid-run.
