# PQCF against QuTiP, quantum optics head to head

QuTiP ([Lambert et al., arXiv:2412.04705](https://arxiv.org/abs/2412.04705); QuTiP 2,
[arXiv:1211.6518](https://arxiv.org/abs/1211.6518)) is the standard open-source library for open quantum systems:
Python, sparse matrices, and SciPy's ODE integrators in double precision. PQCF's master equations (`pq_lindblad`)
work in arbitrary precision: the operators come exact from `pq_qop`, and the linear equations are summed by their
Taylor series to the working precision (here 128 bits, about 38 digits).

This directory runs the same tasks with both on the same machine and the same core, checks every answer against its
closed form where one exists, and publishes the digits and the times. `run.sh` does it all on the committed HEAD.

**Latest run:** [RESULTS.md](RESULTS.md). **Every run:** [`results/`](results) and [history.csv](history.csv).

## The programs

| program | version | how it is run |
| --- | --- | --- |
| **PQCF** | the committed HEAD | `pq_qop` and `pq_lindblad` at 128 bits ([qo_pqcf.erl](qo_pqcf.erl)) |
| **QuTiP (default)** | 5.3.1 (NumPy 2.4.6, SciPy 1.17.1) | `sesolve`, `mesolve` and `steadystate` with their default tolerances ([qutip_side.py](qutip_side.py)) |
| **QuTiP (tight)** | the same | atol 1e-14, rtol 1e-12: the tightest that its default integrator accepts on these problems |

## The tasks ([tasks.json](tasks.json))

| task | what | closed form |
| --- | --- | --- |
| `rabi` | the resonant Jaynes-Cummings model, |e, 0> | P_e(t) = cos^2(g t) |
| `revival` | the same, the field coherent (alpha = 2, 20 photons): collapse and revival up to g t = 30 | P_e(t) = sum_n p_n cos^2(g sqrt(n + 1) t), Eberly, Narozhny and Sanchez-Mondragon, [PRL 44, 1323 (1980)](https://doi.org/10.1103/PhysRevLett.44.1323) |
| `cavity_decay` | a damped cavity, L = sqrt(kappa) a, from a coherent state | <a^dagger a>(t) = <a^dagger a>(0) e^(-kappa t) |
| `fluorescence` | the steady state of a driven two-level atom (Rabi frequency Omega, detuning Delta, decay gamma) | rho_ee = (Omega^2/4) / (Delta^2 + gamma^2/4 + Omega^2/2), Mollow, [Phys. Rev. 188, 1969 (1969)](https://doi.org/10.1103/PhysRev.188.1969) |
| `damped_jc` | the driven Jaynes-Cummings model with cavity and atomic loss, photon cutoffs 4 to 32, to t = 5 | none: QuTiP against PQCF |

The truncations are chosen so that the closed forms hold exactly for the truncated problem: the coherent field of
`revival` stops one photon below the cutoff (every |e, n> keeps its partner |g, n + 1>), and the decay of `cavity_decay`
never needs more photons than it starts with. The parameters are binary fractions, the same numbers on both sides.

Times are the median of repeated runs after one untimed run (`damped_jc`: one run), as in the photonic comparison.
A size of `damped_jc` slower than 120 s is the last one run.

## What to expect

Correctness, not speed, is PQCF's case here. QuTiP's integrators are double precision and adaptive: at their defaults
they give about 6 digits, and at the tightest tolerances 10-14. PQCF's Taylor series gives every digit of its working
precision, and the precision is a parameter. QuTiP is faster on the larger master equations: compiled sparse linear
algebra in double precision against PQCF's exact integer arithmetic in Erlang (the Taylor terms in block fixed point,
P + 24 bits).
