"""Spike S6: chi-square goodness of fit of the Clifford-Clifford samples against the exact distribution (FS 12.4:
bins with expected count below 5 pooled; pass if p > 1e-4), and of the AN v0.1 mutant (must fail).

    docker run --rm -v "$PWD":/w -w /w pqcf-oracle python test/oracle/s6_chi2.py s6.json
"""
import json
import sys

import numpy as np
from scipy.stats import chi2


def test(probs, counts):
    probs, counts = np.array(probs), np.array(counts, dtype=float)
    n = counts.sum()
    exp = probs * n
    order = np.argsort(exp)
    pooled_e, pooled_o, keep_e, keep_o = 0.0, 0.0, [], []
    for i in order:
        if exp[i] < 5:
            pooled_e += exp[i]
            pooled_o += counts[i]
        else:
            keep_e.append(exp[i])
            keep_o.append(counts[i])
    if pooled_e > 0:
        keep_e.append(pooled_e)
        keep_o.append(pooled_o)
    e, o = np.array(keep_e), np.array(keep_o)
    stat = ((o - e) ** 2 / e).sum()
    return stat, len(e) - 1, chi2.sf(stat, len(e) - 1)


data = json.load(open(sys.argv[1]))
ok = True
for c in data["cases"]:
    s, dof, p = test(c["probs"], c["counts"])
    sm, dm, pm = test(c["probs"], c["mutant_counts"])
    passed = p > 1e-4
    caught = pm < 1e-4
    ok &= passed and caught
    print(f"m={c['m']} input={c['input']}: {len(c['probs'])} patterns; Clifford-Clifford chi2 {s:.1f} on {dof} dof, "
          f"p = {p:.3g} ({'pass' if passed else 'FAIL'}); AN v0.1 mutant p = {pm:.2g} ({'caught' if caught else 'NOT caught'}); "
          f"{c['seconds']:.1f} s for 10^6 shots")
for t in data["timing"]:
    print(f"n = {t['n']}, m = {t['m']}: {t['seconds_per_sample']:.4f} s per sample")
print("S6:", "PASS" if ok else "FAIL")
