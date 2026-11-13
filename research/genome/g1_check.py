"""G1's check of the GBS encoding: PQCF's post-selected samples against the exact distribution |haf(A_S)|^2 over
all k-subsets of a small graph (chi-square, with subsets of small expected count pooled), and Oh et al.'s sampler
against haf(A_S).
    python g1_check.py <dir with small.npy, gbs.json>"""
import json
import os
import sys

import numpy as np
from scipy.stats import chisquare

import gbs_graph as G


def chi2(subsets, probs, samples):
    counts = {S: 0 for S in subsets}
    for S in samples:
        counts[tuple(S)] += 1
    n = len(samples)
    exp = probs * n
    obs = np.array([counts[S] for S in subsets])
    big = exp >= 5
    o = np.append(obs[big], obs[~big].sum())
    e = np.append(exp[big], exp[~big].sum())
    if e[-1] == 0:
        o, e = o[:-1], e[:-1]
    return chisquare(o, e).pvalue, len(o) - 1


def main(d):
    A = np.load(os.path.join(d, "small.npy"))
    gbs = {g["name"]: g for g in json.load(open(os.path.join(d, "gbs.json")))}
    out = []
    rng = np.random.default_rng(11)
    for k in (4, 6):
        g = gbs[f"small-k{k}"]
        subsets, p2 = G.exact_distribution(A.tolist(), k, 2)
        pv, dof = chi2(subsets, p2, g["samples"])
        out.append({"what": f"PQCF GBS, 10-vertex graph, k={k}: samples ~ |haf|^2", "p": pv, "dof": dof,
                    "kept": len(g["samples"]), "shots": g["shots"], "ok": bool(pv > 1e-3)})
        subsets, p1 = G.exact_distribution(A.tolist(), k, 1)
        s, _ = G.oh_sampler(A, k, 20000, rng)
        pv, dof = chi2(subsets, p1, s)
        out.append({"what": f"Oh et al. sampler, 10-vertex graph, k={k}: samples ~ haf", "p": pv, "dof": dof,
                    "kept": len(s), "ok": bool(pv > 1e-3)})
    for c in out:
        print(c)
    return out


if __name__ == "__main__":
    res = main(sys.argv[1])
    sys.exit(0 if all(c["ok"] for c in res) else 1)
