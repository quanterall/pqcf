"""G1: Gaussian boson sampling (simulated by PQCF) against Oh et al.'s quantum-inspired classical sampler and uniform
sampling, for the densest k-subgraph and the maximum clique, on a planted-subgraph graph like Arrazola & Bromley's
and the STRING DNA-repair network. See README.md.
    python g1.py specs <dir> <chunks>     # writes <dir>/chunk-<i>/specs.json for PQCF (g1_pqcf.erl)
    python g1.py analyse <dir>            # merges PQCF's samples, runs the other samplers and the searches -> g1.json
"""
import json
import math
import os
import sys
import time

import networkx as nx
import numpy as np

import gbs_graph as G
import instances as I
import solvers as Sv

HERE = os.path.dirname(os.path.abspath(__file__))
KS = (4, 6, 8, 10)
TARGET = 2000            # post-selected samples per graph and k, for every sampler
CHECKPOINTS = [1, 2, 5, 10, 20, 50, 100, 200]
REPS = TARGET // CHECKPOINTS[-1]
# PQCF shots per 100,000 that are k-photon collision-free, measured (12 Nov), to size the runs
RATE = {("string", 4): 1584, ("string", 6): 143, ("string", 8): 57, ("string", 10): 7,
        ("planted", 4): 1900, ("planted", 6): 367, ("planted", 8): 54, ("planted", 10): 5}


def graphs():
    net = I.string_network(os.path.join(HERE, "data", "string_dna_repair.tsv"), 0.9)
    n = len(net.genes)
    A = np.zeros((n, n), int)
    for (i, j) in net.edges:
        A[i, j] = A[j, i] = 1
    P, planted = G.planted_graph()
    return {"string": (A, net.genes, None), "planted": (P, [f"v{i}" for i in range(len(P))], planted)}


def small_graph():
    rng = np.random.default_rng(5)
    A = (rng.random((10, 10)) < 0.5).astype(int)
    A = np.triu(A, 1)
    return A + A.T


def cx(U):
    return [[[float(z.real), float(z.imag)] for z in row] for row in U]


def specs(out, chunks):
    jobs = []
    A = small_graph()
    for k in (4, 6):
        pt = G.sampling_point(A, k)
        e = G.encode(A, pt["mean_photons"])
        jobs.append({"name": f"small-k{k}", "modes": 10, "r": e["r"], "U": cx(e["U"]), "k": k, "target": 2000,
                     "max_shots": int(20 * 2000 / pt["p_k"]), "seed": 100 + k, "max_photons": pt["cap"], "point": pt})
    for name, (M, _, _) in graphs().items():
        for k in KS:
            pt = G.sampling_point(M, k)
            e = G.encode(M, pt["mean_photons"])
            jobs.append({"name": f"{name}-k{k}", "modes": len(M), "r": e["r"], "U": cx(e["U"]), "k": k,
                         "target": TARGET, "max_shots": int(3 * TARGET * 100000 / RATE[(name, k)]),
                         "seed": 1000 * k, "max_photons": pt["cap"], "point": pt})
    # split every job's shots over the chunks, each with its own seed (an independent stream)
    for c in range(chunks):
        d = os.path.join(out, f"chunk-{c}")
        os.makedirs(d, exist_ok=True)
        part = [{**j, "target": -(-j["target"] // chunks), "max_shots": -(-j["max_shots"] // chunks),
                 "seed": j["seed"] + 100 * c} for j in jobs]
        json.dump(part, open(os.path.join(d, "specs.json"), "w"))
    json.dump({j["name"]: {k: v for k, v in j.items() if k not in ("U", "r")} for j in jobs},
              open(os.path.join(out, "plan.json"), "w"), indent=1)
    print(len(jobs), "jobs over", chunks, "chunks")


def merge(out):
    merged = {}
    for d in sorted(os.listdir(out)):
        if not d.startswith("chunk-"):
            continue
        for g in json.load(open(os.path.join(out, d, "gbs.json"))):
            m = merged.setdefault(g["name"], {"samples": [], "shots": 0, "seconds": 0.0, "truncated": []})
            m["samples"] += [tuple(s) for s in g["samples"]]
            m["shots"] += g["shots"]
            m["seconds"] += g["seconds"]
            m["truncated"].append(g["truncated_mass"])
    return merged


def search_stats(A, samples, opt_edges, max_clique, rng):
    """Over REPS disjoint blocks of CHECKPOINTS[-1] samples: the mean best density and the probability of having
    reached the densest k-subgraph, after s samples; and the same for the maximum clique after shrink-and-search."""
    blocks = [samples[r * CHECKPOINTS[-1]:(r + 1) * CHECKPOINTS[-1]] for r in range(REPS)]
    blocks = [b for b in blocks if len(b) == CHECKPOINTS[-1]]
    curves = np.array([G.random_search_curve(A, b, CHECKPOINTS) for b in blocks])
    cliques = []
    for b in blocks:
        best, row = 0, []
        it = iter(b)
        done = 0
        for s in CHECKPOINTS:
            while done < s:
                best = max(best, len(G.shrink_and_search(A, next(it), rng)))
                done += 1
            row.append(best)
        cliques.append(row)
    cliques = np.array(cliques)
    dens = [G.edges_in(A, S) / (len(S) * (len(S) - 1) / 2) for S in samples]
    return {"blocks": len(blocks), "checkpoints": CHECKPOINTS,
            "mean_density_per_sample": float(np.mean(dens)),
            "dks_mean_best_edges": curves.mean(axis=0).tolist(),
            "dks_p_optimum": (curves >= opt_edges).mean(axis=0).tolist(),
            "clique_p_maximum": (cliques >= max_clique).mean(axis=0).tolist()}


def analyse(out):
    gbs = merge(out)
    checks, rng = [], np.random.default_rng(7)

    def check(what, ok, detail):
        checks.append({"what": what, "ok": bool(ok), "detail": detail})

    # the encoding, on the small graph: PQCF's GBS ~ |haf|^2, Oh et al.'s sampler ~ haf (chi-square)
    import g1_check
    A = small_graph()
    for k in (4, 6):
        g = gbs[f"small-k{k}"]
        subsets, p2 = G.exact_distribution(A.tolist(), k, 2)
        pv, dof = g1_check.chi2(subsets, p2, g["samples"])
        check(f"PQCF's GBS on a 10-vertex graph, k = {k}: samples follow |haf(A_S)|^2", pv > 1e-3,
              f"chi-square p = {pv:.3f}, {dof} degrees of freedom, {len(g['samples'])} samples")
        subsets, p1 = G.exact_distribution(A.tolist(), k, 1)
        s, _ = G.oh_sampler(A, k, 20000, rng)
        pv, dof = g1_check.chi2(subsets, p1, s)
        check(f"Oh et al.'s sampler on a 10-vertex graph, k = {k}: samples follow haf(A_S)", pv > 1e-3,
              f"chi-square p = {pv:.3f}, {dof} degrees of freedom")

    plan = json.load(open(os.path.join(out, "plan.json")))
    results = []
    for name, (M, labels, planted) in graphs().items():
        n = len(M)
        net = I.Network(labels, {(i, j): 1.0 for i in range(n) for j in range(i + 1, n) if M[i, j]})
        omega = max(len(c) for c in nx.find_cliques(nx.from_numpy_array(M)))
        for k in KS:
            ex = Sv.exact_dks(net, k, 120.0)
            check(f"{name}, k = {k}: the densest {k}-subgraph is proven", ex["status"] == "optimal", ex["status"])
            g = gbs[f"{name}-k{k}"]
            check(f"{name}, k = {k}: PQCF's truncated probability is negligible", max(g["truncated"]) < 1e-4,
                  f"largest per chunk {max(g['truncated']):.1e}")
            check(f"{name}, k = {k}: enough GBS samples", len(g["samples"]) >= TARGET,
                  f"{len(g['samples'])} of {TARGET}")
            t0 = time.perf_counter()
            oh, draws = G.oh_sampler(M, k, TARGET, rng)
            oh_secs = time.perf_counter() - t0
            t0 = time.perf_counter()
            un = G.uniform_sampler(n, k, TARGET, rng)
            un_secs = time.perf_counter() - t0
            row = {"graph": name, "k": k, "vertices": n, "edges": int(M.sum() // 2),
                   "densest": {"edges": ex["edges"], "of": k * (k - 1) // 2, "seconds": ex["seconds"],
                               "module": [labels[i] for i in ex["chosen"]]},
                   "max_clique": omega,
                   "sampling_point": plan[f"{name}-k{k}"]["point"],
                   "cost": {"pqcf_gbs": {"shots": g["shots"], "kept": len(g["samples"]),
                                         "cpu_seconds_per_sample": g["seconds"] / max(1, len(g["samples"]))},
                            "oh": {"draws": draws, "seconds_per_sample": oh_secs / TARGET},
                            "uniform": {"seconds_per_sample": un_secs / TARGET}}}
            row["gbs_samples"] = [list(map(int, S)) for S in g["samples"][:TARGET]]
            for label, samp in (("gbs", g["samples"][:TARGET]), ("oh", oh), ("uniform", un)):
                row[label] = search_stats(M, samp, ex["edges"], omega, rng)
            if planted is not None and k == len(planted):
                row["planted_found"] = {lab: float(np.mean([set(S) == set(planted) for S in samp]))
                                        for lab, samp in (("gbs", g["samples"]), ("oh", oh), ("uniform", un))}
            results.append(row)
            print(name, k, "GBS p(opt@200) %.2f, Oh %.2f, uniform %.2f" % (
                row["gbs"]["dks_p_optimum"][-1], row["oh"]["dks_p_optimum"][-1], row["uniform"]["dks_p_optimum"][-1]),
                flush=True)
    json.dump({"results": results, "checks": checks}, open(os.path.join(out, "g1.json"), "w"), indent=1)
    bad = [c for c in checks if not c["ok"]]
    print(f"{len(checks) - len(bad)} of {len(checks)} checks hold")
    for c in bad:
        print("FAILS:", c["what"], c["detail"])
    return not bad


if __name__ == "__main__":
    if sys.argv[1] == "specs":
        specs(sys.argv[2], int(sys.argv[3]))
    else:
        sys.exit(0 if analyse(sys.argv[2]) else 1)
