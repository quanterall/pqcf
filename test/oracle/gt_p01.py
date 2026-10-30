"""GT-P01 (FS 12.2): permanent values from thewalrus, random and structured, with repeated rows and columns, and
boson-sampling probabilities with bunching. Repeated multiplicities are expanded into the full matrix here
(np.repeat), so the oracle shares nothing with the kernel's repeated-row formula. Writes test/fixtures/gt_p01.json.

    docker build -t pqcf-oracle test/oracle
    docker run --rm -v "$PWD":/src -w /src pqcf-oracle python test/oracle/gt_p01.py
"""
import json
import math
import pathlib

import numpy as np
import scipy
import thewalrus
from scipy.stats import unitary_group

rng = np.random.default_rng(20261030)


def cm(a):
    return [[[float(z.real), float(z.imag)] for z in row] for row in a]


def perm(a):
    return complex(thewalrus.perm(np.asarray(a, dtype=complex), method="ryser")) if len(a) else 1 + 0j


def rand(r, c):
    return rng.normal(size=(r, c)) + 1j * rng.normal(size=(r, c))


cases = []
for n in list(range(1, 13)) + [14, 16]:
    a = rand(n, n)
    cases.append({"kind": "random", "matrix": cm(a), "perm": [perm(a).real, perm(a).imag]})
for n in (3, 6, 9):
    cases.append({"kind": "all_ones", "matrix": cm(np.ones((n, n))), "perm": [float(math.factorial(n)), 0.0]})
    cases.append({"kind": "identity", "matrix": cm(np.eye(n)), "perm": [1.0, 0.0]})
for n in (4, 8, 12):
    u = unitary_group.rvs(n, random_state=int(rng.integers(1 << 31)))
    cases.append({"kind": "unitary", "matrix": cm(u), "perm": [perm(u).real, perm(u).imag]})

repeated = []
for rows, cols in [([2, 0, 1], [1, 1, 1]), ([3, 1], [2, 2]), ([2, 2, 2], [3, 0, 3]), ([4, 1, 0, 2], [1, 3, 2, 1]),
                   ([5], [5]), ([1, 1, 1, 1, 1, 1], [2, 0, 2, 0, 2, 0]), ([3, 3, 2], [4, 4, 0])]:
    u = rand(len(rows), len(cols))
    p = perm(np.repeat(np.repeat(u, rows, axis=0), cols, axis=1))
    repeated.append({"matrix": cm(u), "rows": rows, "cols": cols, "perm": [p.real, p.imag]})

probabilities = []
for m, inp, outs in [(4, [1, 1, 0, 0], [[2, 0, 0, 0], [1, 1, 0, 0], [0, 1, 0, 1]]),
                     (5, [1, 1, 1, 0, 0], [[3, 0, 0, 0, 0], [1, 1, 1, 0, 0], [0, 2, 0, 0, 1]]),
                     (6, [2, 1, 1, 1, 1, 0], [[6, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1], [0, 0, 3, 0, 3, 0]])]:
    u = unitary_group.rvs(m, random_state=int(rng.integers(1 << 31)))
    for out in outs:
        sub = np.repeat(np.repeat(u, out, axis=0), inp, axis=1)
        p = abs(perm(sub)) ** 2 / math.prod(math.factorial(k) for k in inp + out)
        probabilities.append({"unitary": cm(u), "input": inp, "output": out, "probability": p})

out = {"test": "GT-P01", "generator": "test/oracle/gt_p01.py",
       "versions": {"thewalrus": thewalrus.__version__, "numpy": np.__version__, "scipy": scipy.__version__},
       "permanents": cases, "repeated": repeated, "probabilities": probabilities}
path = pathlib.Path(__file__).resolve().parents[1] / "fixtures" / "gt_p01.json"
path.write_text(json.dumps(out) + "\n")
print(f"wrote {path}: {len(cases)} permanents, {len(repeated)} repeated, {len(probabilities)} probabilities")
