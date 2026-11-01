"""GT-H03 (FS 15 item 2): threshold detection of displaced Gaussian states -- the loop torontonian -- against
thewalrus 0.22.0's threshold_detection_prob (Bulmer, Quesada and Paesani's formula: P(vacuum) Re ltor(O_B, gamma_B)),
for random mixed states with random displacements, click patterns of 1 to 12 modes. Mean and covariance are xxpp,
hbar = 2. Writes test/fixtures/gt_h03.json.

    docker run --rm -v "$PWD":/w -w /w -e NUMBA_CACHE_DIR=/tmp pqcf-oracle python test/oracle/gt_h03.py
"""
import json
import pathlib

import numpy as np
import thewalrus
from thewalrus.random import random_covariance

try:
    from thewalrus.quantum import threshold_detection_prob
except ImportError:
    from thewalrus import threshold_detection_prob

rng = np.random.default_rng(20261101)
cases = []
for N, clicks, scale in [(1, [0, 1], 0.8), (2, [0, 1, 2], 0.6), (4, [1, 2, 3, 4], 0.5), (6, [2, 4, 6], 0.4),
                         (8, [3, 8], 0.3), (12, [6, 12], 0.25)]:
    cov = random_covariance(N, hbar=2, pure=False)
    mu = rng.normal(scale=scale, size=2 * N)
    for c in clicks:
        modes = sorted(rng.choice(N, size=c, replace=False).tolist())
        pattern = [1 if j in modes else 0 for j in range(N)]
        p = float(np.real(threshold_detection_prob(mu, cov, np.array(pattern), hbar=2)))
        cases.append({"n_modes": N, "mean": mu.tolist(), "cov": cov.tolist(), "pattern": pattern, "prob": p})
        print(N, c, p, flush=True)
path = pathlib.Path(__file__).resolve().parents[1] / "fixtures" / "gt_h03.json"
path.write_text(json.dumps({"test": "GT-H03 (displaced threshold detection)", "generator": "test/oracle/gt_h03.py",
                            "versions": {"thewalrus": thewalrus.__version__}, "cases": cases}) + "\n")
print("wrote", path)
