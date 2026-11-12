"""The inputs of the photonic head-to-head (README.md), generated once per run from fixed seeds so every program
gets exactly the same matrices and circuits. Complex numbers are written as [re, im].
    python gen_inputs.py > inputs.json
"""
import json
import sys

import numpy as np
from scipy.stats import unitary_group
from thewalrus.quantum import Qmat
from thewalrus.symplectic import interferometer, squeezing


def cx(a):
    return [[[float(z.real), float(z.imag)] for z in row] for row in np.atleast_2d(a)]


def haar(n, seed):
    return unitary_group.rvs(n, random_state=np.random.default_rng(seed))


def torontonian_input(k, seed):
    # a pure Gaussian state of k modes, every mode squeezed, through a random interferometer; all k modes click. The
    # squeezing of Gaussian boson sampling experiments (r of 1 to 1.5): with weak squeezing the probability that every
    # mode clicks is tiny (1e-13 at k = 22 for r of 0.2 to 0.6), the 2^k terms of the torontonian cancel, and no
    # program's answer is accurate in double precision.
    rng = np.random.default_rng(seed)
    r = rng.uniform(1.0, 1.5, k)
    S = interferometer(haar(k, seed + 1)) @ squeezing(r, np.zeros(k))
    cov = S @ S.T  # hbar = 2, vacuum covariance I
    O = np.eye(2 * k) - np.linalg.inv(Qmat(cov, hbar=2))
    return O


def cvqnn_layer(modes, seed):
    # one layer of the continuous-variable neural network of Killoran et al. (arXiv:1806.06871): interferometer,
    # squeezing, interferometer, displacement, Kerr; each interferometer a phase on every mode, then beamsplitters
    # on the even, then the odd neighbouring pairs. Small parameters keep the truncation error small.
    rng = np.random.default_rng(seed)
    gates = []

    def interf():
        gates.extend(["phaseshift", m, float(rng.uniform(0, 2 * np.pi))] for m in range(modes))
        for start in (0, 1):
            for m in range(start, modes - 1, 2):
                gates.append(["beamsplitter", m, m + 1, float(rng.uniform(0, np.pi / 2)),
                              float(rng.uniform(0, 2 * np.pi))])

    interf()
    gates.extend(["squeezing", m, float(rng.uniform(0.05, 0.15)), float(rng.uniform(0, 2 * np.pi))] for m in range(modes))
    interf()
    gates.extend(["displacement", m, float(rng.uniform(0.05, 0.15)), float(rng.uniform(0, 2 * np.pi))]
                 for m in range(modes))
    gates.extend(["kerr", m, float(rng.uniform(0.05, 0.15))] for m in range(modes))
    return gates


inputs = {
    "permanent": [{"n": n, "matrix": cx(haar(2 * n, 100 + n)[:n, :n])} for n in (12, 16, 20, 24, 26, 28)],
    "torontonian": [{"k": k, "matrix": cx(torontonian_input(k, 200 + k))} for k in (8, 12, 16, 20, 22, 24)],
    "boson_sampling": [{"photons": n, "modes": n * n, "shots": shots, "unitary": cx(haar(n * n, 300 + n))}
                       for n, shots in ((4, 20000), (8, 2000), (12, 200), (16, 20), (20, 5), (24, 2))],
    # the statistics check: every program's samples against the exact output distribution
    "boson_check": {"photons": 3, "modes": 6, "shots": 20000, "unitary": cx(haar(6, 399))},
    "fock": {"modes": 10, "gates": cvqnn_layer(10, 400),
             "total_cutoffs": [4, 6, 8, 10], "per_mode_cutoffs": [3, 4]},
}
json.dump(inputs, sys.stdout)
