# F1: feed-forward

**Decision: feed-forward is two new instructions, `{feed_forward, Gate, #{Symbol => Mode}}` and
`{conditional, Predicate, [Gate]}`, with conditions and parameters as data (symbolic expressions over named
outcomes and the post-selection predicates), resolved shot by shot in the simulator.** Phase 4, item 2 (FS 15
item 1). Released after 1.0.0-rc.1, purely additive under the freeze.

## Why explicit instructions, not implicit outcome symbols

Implicit names such as "`m0` is the outcome on mode 0" would change the meaning of an existing circuit that uses a
symbol `m0` after a measurement on mode 0. That is not additive. With the wrapper, outcome symbols exist only where
a `feed_forward` declares them. They are excluded from `free_symbols/1`, from the caller's bindings and from
renaming.

## Why data, not closures

Conditions and parameters are expressions and predicates, not funs, so a circuit with feed-forward:

- serialises (JSON, the safe binary format);
- is sent to other nodes (the run specification);
- hashes deterministically (the spec cache);
- validates up front (every outcome comes from a measurement earlier in the circuit).

`local_fun` predicates are refused inside `conditional`.

## Semantics

- **Validation:** only gates may go inside (no measurements, no custom instructions, no nesting), so every shot has
  the same measured columns, and the gates' kinds must be supported by the backend.
- **Lowering:** binds everything except the outcome symbols. In each shot, `pq_circuit:resolve_feed_forward/2`
  substitutes that shot's outcomes and evaluates the parameters to floats.
- **The shot loop:** records outcomes by mode as measurements happen (photon counts, clicks 0/1, homodyne values).
- **Determinism (D1):** the shot's stream is consumed in instruction order, as for any mid-circuit measurement, so
  samples are identical for any chunking.
- **Not defined:** exact probabilities and `evolve`, because the state branches on the outcomes. Both return
  `unsupported_input` with `feed_forward => true`.

## Tests (`test/phase4_feed_forward_tests.erl`)

- **Fock, conditional:** a photon sent on to mode 2 only when mode 0 was found empty. The exact branch
  distribution cos²θ₀, sin²θ₀cos²θ₁, sin²θ₀sin²θ₁ fits with chi-square p = 0.346 over 200,000 shots.
- **Fock, feed-forward:** the same branching expressed as an outcome-dependent angle θ₁(1 − k) gives the same
  distribution.
- **Gaussian:**
  - Setup: homodyne x on half of a two-mode squeezed state, then a displacement of the other half by g·x₀.
  - Closed form: the final x₁'s variance is Σ₁₁ + 2gΣ₀₁ + g²Σ₀₀; at the optimal g it is the conditional variance.
  - Result: measured 0.3859 against 0.3880, against 2.58 uncorrected.
- **Determinism:** identical samples for 1 or 8 workers and chunk 37.
- **Refusals:** unmeasured outcome modes, measurements and funs inside, exact probabilities and `evolve`, kinds
  the backend lacks.
- **JSON round trip,** with outcome symbols checked like every symbol name before they become atoms.
