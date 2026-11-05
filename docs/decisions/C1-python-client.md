# C1: a JSON service and a Python client

**Decision: PQCF gets a transport-independent JSON request handler, `pq_service`, and a Python package, `pqcf`, that
talks to it. The default transport is a subprocess over standard input and output; HTTP is opt-in. Python describes
circuits and reads results, and every simulation runs in PQCF.** The project owner chose this as the first item
beyond the specification's roadmap (5 Nov 2026). It is additive: a new module and a new directory, with nothing
frozen changed.

## Why

The photonics community works in Python: Perceval, Piquasso, Strawberry Fields, thewalrus. Without a Python
interface, PQCF can only be used from Erlang.

## Why a service rather than bindings

- **No second engine:** the kernels, the job server, the determinism rules (D1) and the error taxonomy stay where
  they are. Python gets PQCF's results, including the same samples for the same seed, not a reimplementation.
- **A boundary that can be checked:** JSON in, JSON out, with the circuit in the format of FR-C09 (Appendix B),
  which already has an atom-safe reader.
- **Two transports, one handler:**
  - **stdio:** the client starts `erl ... -s pq_service stdio` and exchanges one JSON document per line. Nothing
    listens on the network, and the service's lifetime is the client's.
  - **HTTP:** `POST /api` on OTP's `inets` httpd, for a long-running service shared by several clients. It binds
    127.0.0.1 unless the caller says otherwise. It has no authentication of its own, so binding it elsewhere is a
    deployment decision that needs something in front of it.

## The protocol

A request is `{"id", "op", "circuit", "options", ...}`, with these ops:

- `version`, `validate`, `backend_for`;
- `run`, `sample`, `probabilities`, `evolve`;
- `gradient`, which also takes `observable` and `binding`.

The response is `{"id", "ok": result}` or `{"id", "error": {"category", "reason", "detail"}}`, with PQCF's error
category.

- **Results:**
  - samples as arrays;
  - probabilities as `[pattern, p]` pairs, sorted;
  - states as numbers: Gaussian mean and covariance (xxpp), Fock amplitudes, or a density matrix as `[re, im]`
    pairs, limited to 4,000,000 numbers;
  - metadata as JSON (atoms as strings, tuples as arrays).
- **Atom safety.** Nothing in a request becomes an atom:
  - option keys come from a fixed table, and enum values (`backend`, `mode`, `device`, ...) from fixed lists;
  - names in bindings, `symbols` and gradient weights are matched against the circuit's own symbols, which
    `pq_json` checked and created;
  - predicates go through `pq_json`'s reader.

  The tests send six unknown names in every position, then check that none exists as an atom.
- **Robustness:** a request that crashes the handler is answered with `pq_bug`, and the crash is logged on
  standard error. Requests over 64 MB are refused. In stdio mode the logger writes to standard error, so standard
  output carries only the protocol.

## The Python package (`python/`)

- `pqcf.Circuit`: a builder of the JSON format. Parameters are numbers or expressions in PQCF's syntax, and every
  gate and measurement of the guide is there, `thermal_loss` included.
- `pqcf.connect()`: the client. `run`, `sample`, `probabilities`, `evolve`, `gradient`, `validate` and `version`,
  with typed results (`Result`) and `PqcfError`.
- `pqcf.probability`, `mean_photon_number`, `weighted_sum`: the observables.
- `pqcf.convert`: to and from Perceval and Piquasso, in the conventions measured by
  `test/oracle/converter_conventions.py`:
  - Perceval: `BS.Rx(2θ, φ_bl = π/2 − φ, φ_br = φ − π/2)` on adjacent modes, `PS`; anything else carried exactly as
    a unitary.
  - Piquasso: its gates one for one, with its covariance twice ours.
- `pqcf.linear.unitary`: a passive circuit's single-photon unitary, the reference for the converters.

The client needs only the Python standard library; numpy is needed only by `convert` and `linear`.

## Tests

- **`python/tests/test_client.py`** (15 tests, in the container's suite, so in CI):
  - GT-P01's boson-sampling probabilities (thewalrus 0.22.0) through the client, within 1e-12 + 1e-9 relative;
  - Hong–Ou–Mandel;
  - reproducible samples;
  - all four gradient modes against −sin 2θ, and a weighted sum with a weight t²;
  - bindings and defaults;
  - a Gaussian state's covariance;
  - visibility, Gram matrix, dark counts (both forms) and block storage;
  - post-selection;
  - error categories;
  - the HTTP transport: identical probabilities and samples to stdio, 405 for GET, 404 elsewhere;
  - no ResourceWarnings.
- **`python/tests/test_convert.py`** (6 tests, in the oracle image):
  - converted circuits' unitaries against Perceval's `compute_unitary` within 1e-12, both directions, with
    non-adjacent and reversed beamsplitters;
  - a Piquasso round trip;
  - a Gaussian circuit's mean and covariance against Piquasso's simulator;
  - mutation check: the sign of `φ_bl` flipped, caught by two tests.
- **`test/pq_service_tests.erl`:** requests, the same samples as `pq_api`, atom safety, error categories with the
  id kept, limits.
