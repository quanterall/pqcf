"""GT-G01 (FS 12.2): GBS probabilities on small mode counts from thewalrus 0.22.0 -- photon counting for random
mixed states with and without displacement (thewalrus.quantum.probabilities), and threshold detection for
zero-mean states (threshold_detection_prob). hbar = 2, xxpp. Writes test/fixtures/gt_g01.json.
"""
import itertools
import json
import pathlib

import numpy as np
import thewalrus
from thewalrus.quantum import probabilities
from thewalrus.random import random_covariance

try:
    from thewalrus.quantum import threshold_detection_prob
except ImportError:
    from thewalrus import threshold_detection_prob

rng = np.random.default_rng(20261033)
cases = []
for n_modes, pure, displaced, cutoff in [(1, True, False, 6), (2, False, False, 5), (2, True, True, 5),
                                         (3, False, True, 4), (3, False, False, 4), (4, True, False, 3)]:
    cov = random_covariance(n_modes, hbar=2, pure=pure)
    mu = rng.normal(scale=0.5, size=2 * n_modes) if displaced else np.zeros(2 * n_modes)
    probs = probabilities(mu, cov, cutoff, hbar=2)
    pnr = [{"pattern": list(map(int, idx)), "p": float(probs[idx])}
           for idx in itertools.product(range(cutoff), repeat=n_modes)]
    thr = []
    if not displaced:
        for b in itertools.product([0, 1], repeat=n_modes):
            thr.append({"pattern": list(b), "p": float(np.real(threshold_detection_prob(mu, cov, np.array(b), hbar=2)))})
    cases.append({"n_modes": n_modes, "mean": mu.tolist(), "cov": cov.tolist(), "pnr": pnr, "threshold": thr})
path = pathlib.Path(__file__).resolve().parents[1] / "fixtures" / "gt_g01.json"
path.write_text(json.dumps({"test": "GT-G01", "generator": "test/oracle/gt_g01.py",
                            "versions": {"thewalrus": thewalrus.__version__}, "cases": cases}) + "\n")
print("wrote", path, sum(len(c["pnr"]) + len(c["threshold"]) for c in cases), "probabilities")
