"""GT-H01 (part) and spike S3: hafnian, loop hafnian and repeated-index hafnians from thewalrus 0.22.0, for M up
to 40, with thewalrus's own timings. Writes test/fixtures/gt_h01.json.

    docker run --rm -u $(id -u):$(id -g) -v "$PWD":/w -w /w pqcf-oracle python test/oracle/gt_h01.py
"""
import json
import pathlib
import time

import numpy as np
import thewalrus
from thewalrus import hafnian, hafnian_repeated

rng = np.random.default_rng(20261031)


def sym(m):
    a = rng.normal(size=(m, m)) + 1j * rng.normal(size=(m, m))
    return (a + a.T) / 2


def cm(a):
    return [[[float(z.real), float(z.imag)] for z in row] for row in a]


def c(z):
    return [float(complex(z).real), float(complex(z).imag)]


cases = []
for m in list(range(2, 41, 2)):
    a = sym(m)
    t0 = time.perf_counter()
    h = hafnian(a)
    dt = time.perf_counter() - t0
    cases.append({"kind": "haf", "matrix": cm(a), "value": c(h), "thewalrus_seconds": dt})
for m in list(range(1, 31)):
    a = sym(m)
    t0 = time.perf_counter()
    h = hafnian(a, loop=True)
    dt = time.perf_counter() - t0
    cases.append({"kind": "lhaf", "matrix": cm(a), "loops": [c(a[i, i]) for i in range(m)], "value": c(h),
                  "thewalrus_seconds": dt})
for mult in ([2, 1, 1, 0], [3, 3, 2], [1, 2, 3, 4], [5, 1], [2, 2, 2, 2, 2, 2]):
    a = sym(len(mult))
    cases.append({"kind": "haf_rep", "matrix": cm(a), "mult": mult, "value": c(hafnian_repeated(a, mult))})
    mu = rng.normal(size=len(mult)) + 1j * rng.normal(size=len(mult))
    cases.append({"kind": "lhaf_rep", "matrix": cm(a), "loops": [c(x) for x in mu], "mult": mult,
                  "value": c(hafnian_repeated(a, mult, mu=mu, loop=True))})
path = pathlib.Path(__file__).resolve().parents[1] / "fixtures" / "gt_h01.json"
path.write_text(json.dumps({"test": "GT-H01 (hafnians)", "generator": "test/oracle/gt_h01.py",
                            "versions": {"thewalrus": thewalrus.__version__, "numpy": np.__version__},
                            "cases": cases}) + "\n")
print("wrote", path, len(cases), "cases")
