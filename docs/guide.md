# PQCF user guide

This guide covers the Erlang API of the `pq` application. Every example here is also a test,
`test/guide_examples_tests.erl`, so the examples run as written. If you change one, change the test as well.

PQCF is not ready for use yet and the API may still change. It will be frozen at the end of Phase 3.

## 1. Setup

Build and test in the container (see the README), or with OTP 27+, rebar3 and Rust 1.90:

```sh
rebar3 compile
rebar3 shell
```

Start the application before you sample. The job server, the worker pool and the node registry run under it:

```erlang
{ok, _} = application:ensure_all_started(pq).
```

The examples use the records in `include/pq.hrl` (`#circuit{}`, `#result{}` and the states). In the shell, load
them with `rr("include/pq.hrl").`

## 2. Conventions

| | |
| --- | --- |
| Modes | numbered from 0 |
| ħ | 2 by default (option `hbar`), so the vacuum covariance is the identity |
| Phase space | means and covariances are ordered `x_0 .. x_{m-1}, p_0 .. p_{m-1}` |
| Beamsplitter | `{beamsplitter, M1, M2, Theta, Phi}` acts on the mode operators as `[[cos θ, -e^{-iφ} sin θ], [e^{iφ} sin θ, cos θ]]`, so θ = π/4 is 50:50 |
| Loss | `{loss, M, Eta}`: `Eta` is the transmissivity, 1 means no loss |
| Fock basis | patterns are lists of photon numbers per mode, e.g. `[2, 0]` |
| Results | `{ok, Value}` or `{error, Reason}`. No exceptions, except for programming errors such as a bad record |

## 3. A first circuit: Hong–Ou–Mandel

Two photons enter a 50:50 beamsplitter. They always leave on the same side:

```erlang
{ok, C0} = pq_api:circuit(2, {fock, [1, 1]}),
C1 = pq_circuit:add(C0, {beamsplitter, 0, 1, math:pi() / 4, 0}),
C  = pq_circuit:add(C1, {measure_pnr, [0, 1]}),

fock = pq_api:backend_for(C, #{}),
{ok, P, _Meta} = pq_api:probabilities(C, #{}),
%% P: [2,0] => 0.5, [0,2] => 0.5, [1,1] => ~0

{ok, #result{kind = pnr, measured = [0, 1], samples = S}} =
    pq_api:sample(C, #{shots => 1000, seed => 42}).
%% S has 1000 rows, none of them [1, 1]
```

Circuits are plain values, and every edit returns a new circuit: `add/2`, `insert/3`, `replace/3`, `delete/2`,
`compose/2,3` (with a mode map). `pq_circuit:validate/1` checks modes, parameters, unitarity and measured-mode
reuse. The functions in `pq_api` call it for you.

### Instructions

| Instruction | |
| --- | --- |
| `{beamsplitter, M1, M2, Theta, Phi}` | see Conventions |
| `{phaseshift, M, Phi}` | |
| `{displacement, M, R, Phi}` | α = r e^{iφ} |
| `{squeezing, M, R, Phi}` | |
| `{kerr, M, Kappa}` | Fock backends only |
| `{loss, M, Eta}` | pure-loss channel |
| `{thermal_loss, M, Eta, Nth}` | loss to a thermal environment of mean occupation Nth (1.1) |
| `{unitary, Modes, U}` | a passive unitary on the listed modes (numeric or symbolic matrix) |
| `{measure_pnr, Modes}` | photon counting |
| `{measure_threshold, Modes}` | click / no click |
| `{measure_homodyne, M, Phi}` | quadrature x cos φ + p sin φ |
| `{custom, Name, Params}` | a registered custom instruction (option `custom_instructions`) |

Measurements may sit in the middle of a circuit. A measured mode can't be used again.
`pq_circuit:two_mode_squeezing/4` and `pq_circuit:mzi/4` return lists of instructions for common blocks.
`pq_circuit:mesh/1` turns every `unitary` into a Clements mesh of beamsplitters and phases.

Inputs: `vacuum`, `{fock, Pattern}`, `{gaussian, #gaussian_state{}}` or `{fock_state, #fock_state{}}`.

`pq_state` builds input states:

- **Gaussian:** `vacuum/2`, `coherent/3`, `squeezed_vacuum/3`, `thermal/3`, `from_mean_cov/3`.
- **Fock:** `fock_state(Pattern, Cutoff)`, and `fock_superposition([{Coeff, Pattern}], Cutoff)`, which is
  normalised for you.
- **Density matrices:** `fock_dm(PureState)`, or `fock_dm([{Weight, PureState}])` for a mixture.
- **Gaussian to Fock:** `gaussian_to_fock_dm(GaussianState, Cutoff)` converts a Gaussian state into a density
  matrix within the cutoff, for cross-validation. It returns `{ok, Dm, #{truncated_mass := M}}`; the matrix is not
  renormalised.

## 4. Backends

`pq_api:backend_for/2` picks a backend unless you pass `backend => B`:

| Backend | Chosen when | State |
| --- | --- | --- |
| `gaussian` | the input is vacuum or Gaussian, and every instruction is Gaussian (no Kerr) | mean and covariance |
| `fock` | otherwise, without loss | amplitudes under a cutoff |
| `fock_dm` | otherwise, with loss | density matrix under a cutoff |
| `boson` | only when asked for: passive circuits with a Fock input | none; exact boson sampling (Clifford–Clifford) |

On the Fock backends, the cutoff defaults to the input's photon number. Pass `cutoff => N` when displacement or
squeezing creates photons; a vacuum input without a cutoff is `{unsupported_input, cutoff_required}`. Truncation is reported as a warning (above `truncation_tol`), or as an error under
`strict => true`.

Gaussian circuits with photon-counting measurements are Gaussian boson sampling. The samplers use the hafnian
(PNR) or the torontonian (threshold; for a displaced state, the loop torontonian). Probabilities lost to the `max_photons` bound are reported as
`truncated_mass`:

```erlang
R = 0.5,
{ok, C0} = pq_api:circuit(1),
C = pq_circuit:add(pq_circuit:add(C0, {squeezing, 0, R, 0}), {measure_pnr, [0]}),
gaussian = pq_api:backend_for(C, #{}),
{ok, P, #{truncated_mass := T}} = pq_api:probabilities(C, #{}),
%% P: [0] => 1 / cosh(R), [1] => 0.0, ...; T < 1e-6

{ok, S, _Warnings} = pq_api:evolve(C0#circuit{instrs = [{squeezing, 0, R, 0}]}, #{}),
{ok, N} = pq_stats:mean_photon_number(S, [0]).
%% N = sinh(R)^2
```

## 5. Running

| Function | Returns |
| --- | --- |
| `pq_api:run(C, Opts)` | `{ok, #result{}}`; `shots` defaults to 1. `run(Backend, C, Opts)`, `sample/3` and `evolve/3` name the backend |
| `pq_api:sample(C, Opts)` | the same, with `shots` required |
| `pq_api:evolve(C, Opts)` | `{ok, State, Warnings}`: the state before any measurement |
| `pq_api:probabilities(C, Opts)` | `{ok, #{Pattern => P}, Meta}`: exact probabilities of the measured modes' outcomes |
| `pq_api:run_batch(C, [Binding], Opts)` | `{ok, [#result{}]}`: one run per binding, in parallel, in binding order |

`#result{}` has `kind` (`pnr | threshold | homodyne | mixed`), `measured` (the column order, i.e. the modes in
measurement order), `samples` (one list per accepted shot) and `metadata` (`seed`, `shots_requested`,
`shots_accepted`, warnings, ...).

Common options:

| Option | Default | |
| --- | --- | --- |
| `shots` | 1 | |
| `seed` | random, reported in metadata | see Determinism |
| `bindings` | `#{}` | values for the circuit's symbols |
| `backend` | chosen | `gaussian`, `fock`, `fock_dm`, `boson` |
| `cutoff` | the input's photon number | Fock backends |
| `noise` | `none` | `trajectory`: sample loss as quantum jumps on pure states instead of a density matrix |
| `postselect` | none | a predicate, see Post-selection |
| `shots_mode` | `attempts` | `accepted`: keep going until `shots` samples are accepted (at most `max_attempts`) |
| `max_photons` | 16 | Gaussian enumeration and GBS bound |
| `max_enumeration` | 1 000 000 | the largest outcome set `probabilities/2` will list |
| `max_inflight`, `chunk` | derived | job server parallelism; they never change the samples |
| `timeout` | `infinity` | ms; `{error, timeout}`, or with `partial_results => true` the samples finished so far |
| `notify` | none | a pid that receives `{pq_job, Job}` when the sampling job starts, see Jobs |
| `strict` | false | turn truncation warnings into errors |
| `collision_free` | false | boson sampling: keep only outcomes with at most one photon per mode. This is post-selection, so `attempts` and `acceptance_rate` are reported |
| `check_physical` | false | Gaussian backend: validate the state after every instruction, not only the input |
| `visibility`, `gram`, `multi_photon`, `dark_counts` | none | noise models (1.1), see Loss and noise |
| `dm_storage` | `dense` | `blocks`: the density-matrix backend stores ρ by photon-number blocks inside, skipping zero blocks (number-conserving circuits keep only the diagonal ones), within `max_dm_entries` instead of `max_dm_dim`; results are the same, and `evolve` still returns the full matrix (1.1) |
| `device` | `cpu` | permanents (boson probabilities, `pq_nif:permanent/2`): `gpu` runs them on an NVIDIA GPU, `auto` from n = 24 when one is present (1.1) |

`pq_core:defaults/0` lists every option. Application environment values override the defaults, and call options
override both.

## 6. Determinism

A run is a function of the circuit, the bindings, the options and the seed. Shot k always uses the random stream
(seed, k), whatever the number of workers, the chunk size or the nodes:

```erlang
{ok, A} = pq_api:sample(C, #{backend => boson, shots => 2000, seed => 7, max_inflight => 1}),
{ok, B} = pq_api:sample(C, #{backend => boson, shots => 2000, seed => 7, max_inflight => 8, chunk => 37}),
A#result.samples =:= B#result.samples.   %% true
```

With mid-circuit measurements, a shot draws from its stream in instruction order. With post-selection, the
attempt index is the stream index. Without a seed, a random seed is drawn and recorded in
`metadata`, so any run can be repeated.

## 7. Parameters

Parameters can be symbols, written as expressions. Bind them per run:

```erlang
{ok, Th} = pq_symbolic:parse(<<"theta">>),
{ok, C0} = pq_api:circuit(2, {fock, [1, 0]}),
C = pq_circuit:add(pq_circuit:add(C0, {beamsplitter, 0, 1, Th, 0}), {measure_pnr, [0, 1]}),
[theta] = pq_circuit:free_symbols(C),
{error, {unbound_symbols, [theta]}} = pq_api:probabilities(C, #{}),
{ok, P, _} = pq_api:probabilities(C, #{bindings => #{theta => math:pi() / 3}}),
%% P: [1,0] => 0.25

{ok, [R0, R1]} = pq_api:run_batch(C, [#{theta => 0}, #{theta => math:pi() / 2}], #{shots => 10, seed => 1}).
%% R0: ten [1,0]; R1: ten [0,1]
```

A circuit can carry default bindings (`#circuit.params`) and symbol domains (`real`, `positive`, `nonnegative`,
`complex`). `pq_circuit:effective_binding/2` shows what a run will use.

## 8. Gradients

`pq_gradient:gradient(C, Observable, Binding, Opts)` differentiates an observable with respect to the
symbols in `Binding`. The observables are:

- `{probability, Modes, Pattern}`
- `{mean_photon_number, Modes}`
- `{weighted_sum, [{Weight, Obs}]}`

```erlang
C = pq_circuit:add(C0, {beamsplitter, 0, 1, Th, 0}),
Obs = {probability, [0, 1], [1, 0]},
{ok, #{theta := Gn}} = pq_gradient:gradient(C, Obs, #{theta => 0.3}, #{}),                    %% central differences
{ok, #{theta := Gs}} = pq_gradient:gradient(C, Obs, #{theta => 0.3}, #{mode => symbolic}),  %% exact
{ok, #{theta := E}}  = pq_gradient:symbolic_gradient(C, Obs, [theta]).                      %% the expression, -sin(2 theta)
```

Symbolic gradients are available in two cases:

- **Passive circuits with a Fock input:** computed from the symbolic unitary and the permanent.
- **Gaussian circuits with a vacuum or Gaussian input:** computed from the symbolic mean and covariance. Only
  `mean_photon_number` (and weighted sums of it) is supported here. In every other case the call returns `differentiate_unsupported`; use
numeric mode instead. A custom rule can be plugged in with `strategy => Module`, a module implementing the
`pq_gradient` behaviour.

Two gate-level modes (1.1) work on the Fock backends, non-Gaussian gates included:

- **`mode => adjoint`:** one forward and one backward pass over the pure state, so the cost does not grow with the
  number of parameters. It covers beamsplitters, phase shifts, Kerr gates, numeric unitaries, displacement and
  squeezing. The derivatives are exact for the truncated gates the backend applies, so the result is the gradient
  of the model under the cutoff. With loss or thermal loss in the circuit it runs on density matrices (within
  `max_dm_dim`), and the loss transmissivity is differentiated too. Feed-forward and the parameters of thermal loss
  are refused.
- **`mode => parameter_shift`:** exact shift rules, two evaluations per frequency of each gate parameter. The
  frequencies are bounded by the cutoff P: P for a phase shift, 2P for a beamsplitter or a displacement angle, P² for
  Kerr. It also runs with loss (the density-matrix backend). The magnitudes of displacement and squeezing and the
  loss transmissivity are not periodic and are refused, and so are circuits on the Gaussian backend.

```erlang
{ok, K} = pq_symbolic:parse(<<"kappa">>),
{ok, C0} = pq_api:circuit(1, vacuum),
C = lists:foldl(fun(I, Acc) -> pq_circuit:add(Acc, I) end, C0,
                [{displacement, 0, 0.6, 0}, {kerr, 0, K}, {displacement, 0, 0.6, math:pi()}]),
Obs = {probability, [0], [0]},
{ok, #{kappa := Ga}} = pq_gradient:gradient(C, Obs, #{kappa => 0.4}, #{mode => adjoint, cutoff => 12}),
{ok, #{kappa := Gp}} = pq_gradient:gradient(C, Obs, #{kappa => 0.4}, #{mode => parameter_shift, cutoff => 12}).
%% Ga and Gp agree to 1e-10; central differences to 1e-6
```

Both modes resolve defaults the way lowering does: a gate parameter given by a default that uses a differentiated
symbol is differentiated through it, and a weight in a `weighted_sum` that depends on one adds dW/ds times its
observable.

## 9. Loss and noise

Loss moves the circuit to the density-matrix backend, which applies the exact Kraus sum. With
`noise => trajectory`, the pure-state backend samples the same distribution instead:

```erlang
{ok, C0} = pq_api:circuit(1, {fock, [2]}),
C = pq_circuit:add(pq_circuit:add(C0, {loss, 0, 0.5}), {measure_pnr, [0]}),
fock_dm = pq_api:backend_for(C, #{}),
{ok, P, _} = pq_api:probabilities(C, #{}),          %% [0] => 0.25, [1] => 0.5, [2] => 0.25
{ok, _} = pq_api:sample(C, #{shots => 500, seed => 3, noise => trajectory}).
```

`pq_noise:insert_loss/3` and `pq_noise:detector_efficiency/2` add loss to an existing circuit. Boson sampling
takes loss at the inputs or the outputs, and has a Metropolis sampler for large instances.

Two further noise options (1.1) work on sampling and on exact probabilities:

- **`visibility => M`, or a list of M per input mode:** partially distinguishable photons in boson sampling. Each
  photon shares an internal state with the others with probability √M, so two photons with the same M show the
  Hong–Ou–Mandel visibility M. The sampler stays as fast as Clifford–Clifford, and it is exact for this model. For
  arbitrary overlaps, give `gram => Rows` instead: the n × n Gram matrix of the photons' internal states, with the
  photons listed in input-mode order. The probabilities then cost n! permanents per pattern, so keep n small.
- **`dark_counts => Nu`, or `#{Mode => Nu}`:** every detector adds Poisson-distributed counts with mean ν per shot.
  A photon-counting detector reports n + k; a threshold detector clicks when it saw light, or otherwise with
  probability 1 − e^(−ν). The noise comes from its own stream per shot, so samples still do not depend on chunking,
  and post-selection sees the noisy outcomes. It applies to terminal measurements only.

```erlang
{ok, H0} = pq_api:circuit(2, {fock, [1, 1]}),
H = lists:foldl(fun(I, Acc) -> pq_circuit:add(Acc, I) end, H0, [{beamsplitter, 0, 1, math:pi() / 4, 0}, {measure_pnr, [0, 1]}]),
{ok, Pv, _} = pq_api:probabilities(H, #{visibility => 0.9}),      %% [1, 1] => 0.05: the dip is (1 - M) / 2
{ok, Pd, _} = pq_api:probabilities(H, #{dark_counts => 0.01}),    %% [2, 1] => 0.5 * 0.01 * exp(-0.02); [1, 1] stays 0
{ok, _} = pq_api:sample(H, #{shots => 1000, seed => 4, visibility => 0.9, dark_counts => 0.01}).
```

Two more models:

- **`multi_photon => q`, or a list per input mode:** an occupied input mode emits one extra, distinguishable photon
  with probability q, so a single-photon source has g2 = 2q/(1 + q)². It works in boson sampling and its exact
  probabilities, and it combines with `visibility`.
- **The instruction `{thermal_loss, Mode, Eta, Nth}`:** loss to a thermal environment of mean occupation Nth. It is
  exact on the Gaussian backend. On density matrices it is truncated at the cutoff, and the trace deficit is
  reported like that of displacement.

```erlang
{ok, T0} = pq_api:circuit(1, {fock, [1]}),
T = pq_circuit:add(pq_circuit:add(T0, {thermal_loss, 0, 0.6, 0.8}), {measure_pnr, [0]}),
fock_dm = pq_api:backend_for(T, #{}),
{ok, Pt, _} = pq_api:probabilities(T, #{cutoff => 40}).     %% mean photon number 0.6 + 0.4 * 0.8
```

`pq_stats` provides:

- `mean_photon_number/2`
- `g2/2`
- `tv_distance/2`
- `histogram/2`
- `fidelity/2`

## 10. Post-selection

`postselect => Predicate` keeps only the samples that satisfy the predicate. Modes in the predicate are circuit
modes, not column numbers:

```erlang
{ok, C0} = pq_api:circuit(2, {fock, [1, 0]}),
C = pq_circuit:add(pq_circuit:add(C0, {beamsplitter, 0, 1, math:pi() / 4, 0}), {measure_pnr, [0, 1]}),
{ok, #result{samples = S}} =
    pq_api:sample(C, #{shots => 100, seed => 5, postselect => {eq, 1, 0}, shots_mode => accepted}).
%% S: a hundred [1, 0]
```

The predicates are:

- `{eq | ge | le, Mode, N}`
- `{in, Mode, Ns}`
- `{total_eq | total_ge | total_le, Modes, N}`
- `{'and', Ps}`, `{'or', Ps}`, `{'not', P}`
- `{local_fun, F}`: an Erlang fun. It can't be sent to other nodes, so it gives `predicate_not_distributable` when
  the run is distributed.

If nothing is accepted within `max_attempts`, the run fails with `postselection_exhausted`.
`pq_sampling:postselect(Result, Predicate)` applies a predicate after the run.

### Feed-forward

From 1.1 on, gates can depend on earlier measurements in the same shot. There are two instructions for this:

- **`{feed_forward, Gate, #{Symbol => Mode}}`:** each symbol is bound, per shot, to the outcome measured earlier
  on its mode. Outcomes are photon counts, clicks (0 or 1) or homodyne values.
- **`{conditional, Predicate, [Gate]}`:** the gates are applied only when the predicate holds. It uses the
  post-selection predicates (data only, no `local_fun`).

```erlang
{ok, C0} = pq_api:circuit(3, {fock, [1, 0, 0]}),
C = lists:foldl(fun(I, Acc) -> pq_circuit:add(Acc, I) end, C0,
                [{beamsplitter, 0, 1, 0.6, 0}, {measure_pnr, [0]},
                 %% only when mode 0 was found empty: send the photon on from mode 1
                 {conditional, {eq, 0, 0}, [{beamsplitter, 1, 2, 1.1, 0}]},
                 %% a parameter from an outcome: the phase is pi times the count on mode 0
                 {feed_forward, {phaseshift, 2, {mul, [{const, pi}, {sym, k}]}}, #{k => 0}},
                 {measure_pnr, [1, 2]}]),
{ok, #result{measured = [0, 1, 2]}} = pq_api:sample(C, #{shots => 1000, seed => 1}).
```

The outcome must come from a measurement earlier in the circuit. Only gates go inside, so every shot has the same
measured columns. The outcome symbols are not the caller's to bind. A shot draws from its stream in instruction
order as before, so the samples don't depend on chunking. Exact probabilities and `evolve` return
`unsupported_input` for such circuits, because the state depends on the outcomes: sample it instead. JSON carries
both instructions (`"op": "feed_forward"` with `"gate"` and `"outcomes"`, and `"op": "conditional"` with
`"predicate"` and `"gates"`).

## 11. The symbolic engine

`pq_symbolic` keeps expressions in an exact canonical form: Gaussian rationals, radicals, π, e, and
trigonometric functions written as exponentials:

```erlang
{ok, E} = pq_symbolic:parse(<<"cos(x)^2 + sin(x)^2">>),
{ok, One} = pq_symbolic:simplify(E),
<<"1">> = pq_symbolic:pretty(One),

{ok, F} = pq_symbolic:parse(<<"x^3 * exp(i*y)">>),
{ok, Dx} = pq_symbolic:differentiate(F, x),
{ok, Prog} = pq_symbolic:compile([Dx], [x, y]),          %% a straight-line program
{ok, [V]} = pq_symbolic:eval_slp(Prog, #{x => 2.0, y => 0.0}),   %% 12

{ok, A} = pq_symbolic:parse(<<"sqrt(x^2)">>),
{ok, X} = pq_symbolic:parse(<<"x">>),
{ok, {false, _Witness}} = pq_symbolic:equiv(A, X, #{numeric => 8}).
```

The engine also provides:

- **Substitution and evaluation:** `subs`, `evaluate` (complex doubles, or `exact`).
- **Derivatives:** `gradient`, `jacobian`.
- **Matrices:** `matmul`, `kron`, `dagger`, ...
- **Printing:** `to_latex`.
- **Circuits:** `circuit_unitary/1` (the symbolic unitary of a passive circuit), `is_unitary`, `permanent`,
  `bs_probability`.

Gaussian circuits evolve symbolically too: the mean and covariance become canonical matrices. This covers
squeezing, displacement, phase shifts, beamsplitters, passive unitaries and loss, for up to
`max_symbolic_modes` modes. Declare loss symbols `nonnegative`. To sample, lower the state first:

```erlang
{ok, R} = pq_symbolic:parse(<<"r">>),
{ok, C0} = pq_api:circuit(1),
C = pq_circuit:add(C0, {squeezing, 0, R, 0}),
{ok, S} = pq_gaussian:symbolic_evolve(C, #{}),
{matrix, [[Vx, _], [_, _]]} = S#gaussian_state.cov,
<<"exp(-2*r)">> = pq_symbolic:pretty(Vx),
{ok, L} = pq_gaussian:lower(S, #{r => 0.5}),                 %% a numeric state again
{ok, #{r := G}} = pq_gradient:gradient(C, {mean_photon_number, [0]}, #{r => 0.5}, #{mode => symbolic}).
%% G = sinh(1)
```

`pq_gaussian:to_symbolic/1` turns a numeric state into a symbolic one.

`equiv/3` returns `true` or `false` only when the canonical form decides the question. Otherwise it returns
`unknown`, or, with `numeric => N`, it may find a binding where the two expressions differ.

Every call is bounded by `max_expr_nodes`, `sym_timeout_ms` and `sym_max_heap_words`. When a limit is hit, the
call returns `symbolic_limit` instead of running away. Large products of plain polynomials (symbols, `conj`
and π with integer exponents, exact coefficients) are multiplied in Rust. The result is the same, and
`poly_kernel => false` turns this off. The parser creates atoms only for allowed names
(`allowed_symbols`, or at most `max_symbols`).

Use `compile/2` and `eval_slp/2` (or `eval_batch/2`) when one expression is evaluated many times. They share
common subexpressions and run about 20× faster than evaluating the expression tree.

With `compile/3` and `#{optimize => horner}` (1.3), every polynomial becomes a multivariate Horner scheme. It needs
2.6 to 3 times fewer operations on large polynomials, and evaluation is about 3.5 times faster. The operation order
differs from the default's, so values agree within rounding, not bit for bit. The default stays as it was.

```erlang
{ok, Fast} = pq_symbolic:compile([Dx], [x, y], #{optimize => horner}).
```

`#{optimize => horner_cse}` (1.3) also shares common subexpressions across the Horner scheme. On the resultant
benchmarks it needs 7–9 % fewer operations than `horner`, within 8–10 % of FORM 5's `Format O2`.
`#{optimize => cse}` applies the same pass to the default program. Values again agree within rounding.

With Horner, `order_search => #{trials => N}` (1.3) also searches N variable orders and keeps the one with the
fewest operations. Each batch of candidates (`batch`, 16 by default) is scored in parallel, one process per
candidate, on this node or spread over `nodes => [...]`. The program is the same whatever the cores or nodes, and
`seed` picks another search. On the larger resultants, 160 trials go below FORM 5's `O2`.

```erlang
{ok, P} = pq_symbolic:compile([E], Syms, #{optimize => horner_cse,
                                          order_search => #{trials => 160, nodes => [node() | nodes()]}}).
```

## 12. Interchange

- **JSON:** `pq_circuit:to_json/1` and `pq_circuit:from_json/1,2`. Symbols survive the round trip, and
  `from_json/2` takes options such as `allowed_symbols`.
- **Erlang binary:** `pq_circuit:to_binary/1` and `pq_circuit:from_binary/1`, a safe format for Erlang peers.
- **Perceval:** `pq_convert:to_perceval/2` and `from_perceval/1`.
- **Piquasso:** `pq_convert:to_piquasso/2` and `from_piquasso/1`. Their covariances use ħ = 1 against our 2; the
  converter rescales.

The conventions of both converters are checked against the libraries themselves in
`test/oracle/converter_conventions.py`.

```erlang
{ok, Json} = pq_circuit:to_json(C),
{ok, Back} = pq_circuit:from_json(Json).
```

## 13. Jobs

Every run with at least `parallel_threshold` shots (64 by default) is a job. The job server splits the shots
into chunks, runs at most `max_inflight` chunks at a time in worker processes, retries a crashed chunk up to
`max_retries` times and merges the results in shot order. If the `pq` application isn't started, the job runs
in the calling process instead, with the same samples.

```erlang
Self = self(),
spawn(fun() -> Self ! {done, pq_api:run(C, #{shots => 1000000, seed => 1, notify => Self})} end),
receive {pq_job, Job} -> ok end,
{ok, #{state := running, done := Done, chunks := Total}} = pq_parallel:status(Job),
ok = pq_parallel:cancel(Job).                 %% the run returns {error, cancelled}
```

With `timeout => Ms, partial_results => true`, a timeout returns
`{error, {timeout, #{partial => #{samples => Rows, shot_indices => Is}}}}`. These are the rows the full run would
have had at those shot indices, as the worker produced them, before any projection to columns.

When the `telemetry` library is loaded, jobs emit `[pq, job, start | stop | exception]` and
`[pq, chunk, start | stop | retry]`. PQCF doesn't depend on it: add it to your application to get the events.

## 14. Several nodes

Connect the nodes and start `pq` on each. Nodes find each other through the `pq` process group and a version
handshake (same NIF, same architecture unless `allow_heterogeneous_arch`). Then:

- `distributed => true` spreads the chunks of a sampling job over the registered nodes. To choose the nodes
  yourself, pass `distributed => nodes, nodes => [...]`.
- If a node goes down, its chunks are retried on another node. The samples are the same as on one node.
- `isolate_kernels => true` runs the kernels on peer nodes, so a crash in native code takes down only the
  peer.
- A job is refused up front if its estimated memory is more than `max_job_memory` (default: half of the
  machine's memory). With `auto_throttle => true`, it runs with fewer chunks in flight instead.

## 15. Errors

Every error is `{error, {Category, Detail}}` or `{error, Category}`. `pq_core:categories/0` lists the
categories. The common ones are:

| Category | Meaning |
| --- | --- |
| `bad_shape`, `mode_out_of_range`, `duplicate_modes`, `measured_mode_reuse` | the circuit or the input is malformed |
| `unbound_symbols`, `param_conflict`, `domain_violation` | parameters |
| `unsupported_instruction`, `unsupported_input` | the backend can't run this, e.g. Kerr on `gaussian`, or a missing cutoff |
| `dimension_limit`, `symbolic_limit` | a configured bound was hit, and the detail names it |
| `postselection_exhausted`, `cancelled`, `timeout`, `node_down` | jobs |
| `sym_parse`, `json_*` | parsing, with the byte position |

Warnings (`truncation`, `gbs_truncation`, `unused_bindings`) never fail a run. They are returned in the
metadata, or as the third element of `evolve/2`.

## 16. From Python (1.2)

The `pqcf` Python package (`python/`) builds circuits in the JSON format of section 12 and runs them through
PQCF's JSON service, `pq_service`. Simulations run in PQCF itself: Python only describes circuits and reads
results.

```python
import pqcf
c = pqcf.Circuit(2, fock=[1, 1]).beamsplitter(0, 1, "pi/4").measure_pnr([0, 1])
with pqcf.connect() as q:                       # starts `erl ... -s pq_service stdio`; no network
    q.probabilities(c)                          # {(0, 2): 0.5, (1, 1): ~0, (2, 0): 0.5, ...}
    q.sample(c, shots=1000, seed=1).samples     # the same samples as pq_api:sample/2 with that seed
    q.probabilities(c, visibility=0.9)          # every run option of section 5, by name
    q.gradient(pqcf.Circuit(2, fock=[1, 0]).beamsplitter(0, 1, "t"),
               pqcf.probability([0], [1]), {"t": 0.3}, mode="adjoint")
```

- **Transports:**
  - **stdio**, the default: `pqcf.connect()` starts the service as a subprocess, using `PQCF_HOME`'s build or
    `PQCF_COMMAND`, and talks over its standard input and output.
  - **HTTP:** start `erl ... -s pq_service http`, which listens on 127.0.0.1 with `PQ_HTTP_PORT`, or `start_http/2`
    from Erlang, then use `pqcf.connect("http://127.0.0.1:PORT")`. The service has no authentication: bind it to
    another address only behind something that adds it.
- **Errors** raise `pqcf.PqcfError`, carrying PQCF's error category (section 15).
- **Converters:** `pqcf.convert` converts to and from Perceval and Piquasso, with numpy, and `pqcf.linear.unitary`
  gives a passive circuit's unitary.
- **Tests:** the package's tests run with the container's test suite, and the converters' tests in the oracle image
  (`python/README.md`).
