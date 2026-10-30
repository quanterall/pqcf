# Mid-circuit measurements: stream order and modes

**Decisions** (30 October 2026), for FR-G08, FR-M05-M07 and FS 9.3.

1. **One stream per shot, consumed in instruction order.** A shot k that measures mid-circuit runs its
   instructions from the state before its first measurement, drawing from stream (Seed, k) held in a NIF
   resource (`stream_new/2`, passed to kernels as `stream_ref`): each mid-circuit homodyne takes one normal pair
   (Box-Muller, the cosine normal used), each mid-circuit PNR one uniform, then the terminal sampler continues on
   the same stream. Without mid-circuit measurements the terminal sampler draws shot k from (Seed, k) as before.
   Samples are identical for any number of workers and chunk sizes (tested).
2. **Measured modes keep their place.** FS A.7 deletes the measured mode and has the backend track a live-mode
   map. Instead: in the Fock backends the measured modes are projected onto the outcome and stay in the state; in
   the Gaussian backend the conditioned state has the measured mode put back as vacuum, uncorrelated. No later
   instruction may act on a measured mode (FR-M06) and the measured columns come from the measurement, so this
   changes no outcome and needs no renumbering.
3. **Post-selection** (FR-M05): `attempts` returns the accepted samples among `shots` attempts; `accepted` runs
   further attempts (indices continuing, so attempt i always uses stream (Seed, i)) until `shots` are accepted or
   `max_attempts` (default 1000 x shots) is reached (`postselection_exhausted`). The accepted samples are the first
   accepted attempts, so they equal the attempts-mode samples truncated.
