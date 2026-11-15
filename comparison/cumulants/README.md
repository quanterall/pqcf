# PQCF's cumulant expansions against QuantumCumulants.jl

`pq_cumulant` (quantum optics plan, Q4) derives the equations of motion for averages of operators from a Hamiltonian
and jump operators, closes them by setting the cumulants above a chosen order to zero, adds every average they need
until the set closes, and solves them. QuantumCumulants.jl (Plankensteiner, Hotter, Ritsch,
[Quantum 6, 617 (2022)](https://doi.org/10.22331/q-2022-01-04-617), [arXiv:2105.01657](https://arxiv.org/abs/2105.01657))
does the same in Julia; this directory runs the same models with it, for PQCF's tests.

| model | order | equations (both) | agreement |
| --- | --- | --- | --- |
| single-atom laser (QuantumCumulants' first documented example) | 2 | 8 | <a^dagger a> and P_e at t = 1, 2, 5, 10 within 1e-13 |
| driven, damped Jaynes-Cummings model | 2, 3, 4 | 8, 15, 24 | <a^dagger a>(2) within 1e-15 at each order |

Against the exact master equation (PQCF's `pq_lindblad`, 25 photons), the Jaynes-Cummings <a^dagger a>(2) is off by
1.3e-5 at order 2, 1.9e-4 at order 3 and 5.8e-6 at order 4: the expansion converges, but not monotonically
(QuantumCumulants gives the same numbers).

## Running it

```sh
sh comparison/cumulants/run.sh      # builds the Julia image (about 5 minutes), runs both models, writes the fixture
docker build -t pqcf-dev . && docker run --rm pqcf-dev rebar3 eunit --module=pq_cumulant_tests
```

`run.sh` writes `test/fixtures/qc_cumulants.json`, which `test/pq_cumulant_tests.erl` reads. The image
([julia/Dockerfile](julia/Dockerfile)) pins QuantumCumulants 0.3.7; its build prints the versions of the rest
(ModelingToolkit 9.84.0, OrdinaryDiffEq 6.107.0, Symbolics 6.58.0 when the fixture was written).
