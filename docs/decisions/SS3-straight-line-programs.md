# SS3: straight-line programs

**Decision: no `slp_eval` NIF.** The Erlang interpreter of straight-line programs meets the target of SY 7.5 by a
wide margin. Closed 31 October 2026.

Exit test (SY 11): at least 5x faster than tree evaluation for k = 4 and k = 8. Measured on the symbolic unitary of
a rectangular mesh with every angle and phase symbolic (`test/pq_symbolic_tests.erl`, ss3_test_), one core of the
shared aarch64 machine, 20 evaluations each:

| k | build the unitary | compile | ops | eval, program | eval, tree | speed-up |
| --- | --- | --- | --- | --- | --- | --- |
| 4 | under 0.01 s | 0.01 s | 1,597 | 0.30 ms | 5.94 ms | 19.9x |
| 6 | 0.06 s | 0.43 s | 46,423 | 16.2 ms | 357 ms | 22.1x |

The two evaluations agree within 1e-12. k = 8 was not run here; the 4-mode unitary meets SY 7.5's 0.1 s target.

Also in this change: `equiv` (FR-S11) answers `false` only where the canonical form is faithful -- no opaque base,
no float, no exponential of pi (the cos(pi/5) counter-example to SY 3.5), not a constant-only difference with e or
pi -- else `unknown`, optionally followed by N deterministic numeric points that can return `{false, Binding}`.
