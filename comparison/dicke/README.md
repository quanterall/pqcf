# Collective emission: PQCF's Dicke-basis master equation against PIQS and an exact solution

`pq_dicke` (quantum optics plan, Q5) evolves N identical two-level atoms in the permutation-invariant (Dicke) basis,
with collective and local emission, dephasing and pumping and a Hamiltonian omega J_z + omega_x J_x, in arbitrary
precision. The method and its coefficients are those of PIQS (Shammah, Ahmed, Lambert, De Liberato, Nori,
[PRA 98, 063815 (2018)](https://doi.org/10.1103/PhysRevA.98.063815), [arXiv:1805.05129](https://arxiv.org/abs/1805.05129)),
QuTiP's `qutip.piqs`.

| check | result |
| --- | --- |
| the superradiant burst of 100 atoms against the exact population equations (mpmath, 40 digits) | equal to 34 digits at six times (1.4 s) |
| the burst's peak against Gross and Haroche's t_D = ln N/(N gamma) = 0.046 | at 0.048 on a grid of 0.002, as in PIQS; height 1,972 (N^2/4 = 2,500 only asymptotically) |
| PIQS, N = 6 with every process; N = 5 driven (half-integer ladders); N = 100 burst | within 1e-8 to 1e-9 |
| the trace, every process at once | conserved to 1e-25 |

**PIQS computes in single precision.** Its Cython source (`qutip/piqs/_piqs.pyx`, QuTiP 5.3.1) declares the rates
and the Lindbladian's coefficients as `float`, which is 32-bit; its results are therefore good to about 1e-7
relative, whatever tolerance the solver is given. PQCF agrees with the exact solution to 34 digits and with PIQS to
1e-8, which locates the difference in PIQS.

## Running it

```sh
docker run --rm -v "$PWD":/w -w /w pqcf-qutip python comparison/dicke/piqs_check.py         # 8 minutes
docker run --rm -v "$PWD":/w -w /w pqcf-oracle python comparison/dicke/superradiance_mp.py  # 1 minute
docker run --rm pqcf-dev rebar3 eunit --module=pq_dicke_tests
```

The images: `pqcf-qutip` from [comparison/quantum-optics/qutip](../quantum-optics/qutip), `pqcf-oracle` from
[test/oracle](../../test/oracle). The two scripts write `test/fixtures/piqs_dicke.json` and
`test/fixtures/superradiance_mp.json`.
