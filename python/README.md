# pqcf: the Python client of PQCF

PQCF, the Photonic Quantum Computing Framework, simulates photonic quantum computers in Erlang/OTP with Rust kernels
(https://github.com/quanterall/pqcf). This package lets Python programs use it:

- build circuits in its JSON format;
- run them through its JSON service: samples, exact probabilities, output states and gradients;
- convert circuits to and from Perceval and Piquasso.

```python
import pqcf
c = pqcf.Circuit(2, fock=[1, 1]).beamsplitter(0, 1, "pi/4").measure_pnr([0, 1])
with pqcf.connect() as q:
    print(q.probabilities(c))       # Hong-Ou-Mandel: {(0, 2): 0.5, (2, 0): 0.5, ...}
```

The client needs only the Python standard library. Converting needs numpy, plus Perceval or Piquasso.

## Requirements

A built PQCF: OTP 27 and `rebar3 compile` in a checkout, or the `pqcf-dev` image. `pqcf.connect()` starts
`erl -noshell -pa <PQCF_HOME>/_build/default/lib/*/ebin -s pq_service stdio`:

- `PQCF_HOME` defaults to the checkout this package sits in;
- `PQCF_COMMAND` replaces the whole command;
- `PQCF_QUIET=1` silences the service's log on standard error.

For HTTP, start `erl ... -s pq_service http` (`PQ_HTTP_PORT`, default a free port; it prints `{"listening": port}`)
and call `pqcf.connect("http://127.0.0.1:PORT")`.

## What it offers

| Call | Returns |
| --- | --- |
| `q.version()` | PQCF's version, the NIF contract, the target |
| `q.validate(c)` | `{"valid": True, "free_symbols": [...]}`, or `PqcfError` |
| `q.run(c, **opts)`, `q.sample(c, shots, **opts)` | `Result(kind, measured, samples, metadata)` |
| `q.probabilities(c, **opts)` | `{pattern tuple: probability}` |
| `q.evolve(c, **opts)` | the output state: Gaussian mean and covariance, Fock amplitudes, or a density matrix |
| `q.gradient(c, observable, binding, mode=...)` | `{symbol: derivative}`; observables from `pqcf.probability`, `pqcf.mean_photon_number`, `pqcf.weighted_sum` |

The options are PQCF's run options by name (guide section 5), for example `shots`, `seed`, `backend`, `cutoff`,
`visibility`, `gram`, `multi_photon`, `dark_counts`, `device`, `dm_storage`, `bindings` and `postselect` (as PQCF's
JSON predicates). An unknown option is an error, not ignored.

## Tests

```sh
docker run --rm pqcf-dev        # the whole suite, including python -m unittest discover -s python/tests
docker run --rm -v "$PWD":/src -w /src pqcf-oracle python -m unittest discover -s python/tests   # the converters
```

The client tests check GT-P01's boson-sampling probabilities (thewalrus 0.22.0) through the client, closed forms,
all four gradient modes, the noise options, errors, and identical samples over both transports.
