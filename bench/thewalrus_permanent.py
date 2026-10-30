"""NFR-01 reference timings: thewalrus permanents of the matrices of bench_permanent.rs, every method, median."""
import json
import math
import sys
import time

import numpy as np
import thewalrus

out = []
for n, reps in ((16, 31), (18, 15), (20, 7)):
    a = np.array([[math.sin(i * n + j) + 1j * math.cos(i * j + 0.5) for j in range(1, n + 1)] for i in range(1, n + 1)])
    row = {"n": n}
    for method in ("bbfg", "glynn", "ryser"):
        thewalrus.perm(a, method=method)
        ts = []
        for _ in range(reps):
            t0 = time.perf_counter()
            v = thewalrus.perm(a, method=method)
            ts.append(time.perf_counter() - t0)
        row[method] = sorted(ts)[reps // 2]
        row["value"] = [float(np.real(v)), float(np.imag(v))]
    out.append(row)
json.dump({"thewalrus": thewalrus.__version__, "numpy": np.__version__, "results": out}, sys.stdout)
