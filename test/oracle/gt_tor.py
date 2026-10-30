"""GT-H01 (torontonian) and spike S4: tor(O_B) from thewalrus 0.22.0 for random mixed Gaussian states, click
patterns of 1 to 20 modes, and the threshold probability P(B) = tor(O_B) / sqrt(det Q) against thewalrus's
threshold_detection_prob. O_B is laid out as our kernel expects: the clicked modes j, then their conjugates j + N.
Writes test/fixtures/gt_tor.json.
"""
import json
import pathlib
import time

import numpy as np
import thewalrus
from thewalrus import tor
from thewalrus.quantum import Qmat
from thewalrus.random import random_covariance

try:
    from thewalrus.quantum import threshold_detection_prob
except ImportError:
    from thewalrus import threshold_detection_prob

rng = np.random.default_rng(20261032)
cases = []
for N, clicks in [(4, [1, 2, 3, 4]), (8, [5, 8]), (12, [10, 12]), (16, [14, 16]), (20, [18, 20])]:
    cov = random_covariance(N, hbar=2, pure=False)
    Q = Qmat(cov, hbar=2)
    O = np.eye(2 * N) - np.linalg.inv(Q)
    for c in clicks:
        modes = sorted(rng.choice(N, size=c, replace=False).tolist())
        idx = modes + [j + N for j in modes]
        OB = O[np.ix_(idx, idx)]
        t0 = time.perf_counter()
        t = tor(OB)
        dt = time.perf_counter() - t0
        pattern = [1 if j in modes else 0 for j in range(N)]
        p = threshold_detection_prob(np.zeros(2 * N), cov, np.array(pattern), hbar=2) if N <= 12 else None
        cases.append({"n_modes": N, "modes": modes, "o": [[[float(z.real), float(z.imag)] for z in row] for row in OB],
                      "tor": [float(np.real(t)), float(np.imag(t))], "sqrt_det_q": float(np.sqrt(np.linalg.det(Q)).real),
                      "threshold_prob": None if p is None else float(np.real(p)), "thewalrus_seconds": dt})
        print(N, c, dt, flush=True)
path = pathlib.Path(__file__).resolve().parents[1] / "fixtures" / "gt_tor.json"
path.write_text(json.dumps({"test": "GT-H01 (torontonian)", "generator": "test/oracle/gt_tor.py",
                            "versions": {"thewalrus": thewalrus.__version__}, "cases": cases}) + "\n")
print("wrote", path)
