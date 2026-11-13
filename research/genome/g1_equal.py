"""G1 at equal cost: how many samples, and how much time on one core, each sampler needs to reach the densest
k-subgraph. GBS's figure comes from its measured curve (samples for a 50% chance, from the run's blocks), the
others' from sampling until the optimum is hit, 20 times, up to CAP samples.
    python g1_equal.py results/g1-<host>-<date>.json      # adds "equal_cost" to every row, rewrites the file"""
import json
import sys
import time

import numpy as np

import gbs_graph as G
import g1

CAP = 2_000_000
RUNS = 20
BATCH = 20_000


def densities(A, S):
    S = np.asarray(S)
    return A[S[:, :, None], S[:, None, :]].sum(axis=(1, 2)) // 2


def first_hit(A, k, opt, draw, rng):
    n = 0
    while n < CAP:
        S = draw(BATCH, rng)
        hit = np.flatnonzero(densities(A, S) >= opt)
        if hit.size:
            return n + int(hit[0]) + 1
        n += len(S)
    return None


def s50(curve, checkpoints):
    """The smallest checkpoint at which the chance of having reached the optimum is at least 50%."""
    for p, s in zip(curve, checkpoints):
        if p >= 0.5:
            return s
    return None


def main(path):
    rec = json.load(open(path))
    graphs = g1.graphs()
    rng = np.random.default_rng(3)
    for row in rec["results"]:
        A = graphs[row["graph"]][0]
        n, k, opt = len(A), row["k"], row["densest"]["edges"]
        samplers = {
            "oh": lambda b, r: np.array(G.oh_sampler(A, k, b, r)[0]),
            "uniform": lambda b, r: np.array(G.uniform_sampler(n, k, b, r)),
        }
        out = {}
        for lab, draw in samplers.items():
            t0 = time.perf_counter()
            draw(BATCH, rng)
            per = (time.perf_counter() - t0) / BATCH
            hits = [first_hit(A, k, opt, draw, rng) for _ in range(RUNS)]
            found = [h for h in hits if h is not None]
            med = float(np.median(found)) if len(found) * 2 > RUNS else None
            out[lab] = {"runs": RUNS, "reached": len(found), "median_samples": med,
                        "seconds_per_sample": per, "median_seconds": med * per if med else None}
        s = s50(row["gbs"]["dks_p_optimum"], row["gbs"]["checkpoints"])
        per = row["cost"]["pqcf_gbs"]["cpu_seconds_per_sample"]
        out["gbs"] = {"samples_for_50pct": s, "cpu_seconds_per_sample": per,
                      "cpu_seconds_for_50pct": s * per if s else None}
        row["equal_cost"] = out
        print(row["graph"], k, {l: (v.get("median_samples", v.get("samples_for_50pct")),
                                   v.get("median_seconds", v.get("cpu_seconds_for_50pct"))) for l, v in out.items()},
              flush=True)
    json.dump(rec, open(path, "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1])
