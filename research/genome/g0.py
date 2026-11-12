"""G0: the classical baselines on the genomic instances, measured as Maurizio & Mazzola ask (arXiv:2507.04111):
the optimal time to solution TTS* = min over run lengths t of t * ln(0.01) / ln(1 - p(t)), with p(t) the
probability that one run of length t reaches the exact optimum, at several sizes. The exact optimum comes from
CP-SAT on the native formulation; CP-SAT is also run on each QUBO, to measure what the translation costs.
    python g0.py <out dir>
Writes g0.json and prints a summary; exits 1 if any check fails."""
import json
import math
import os
import sys
import time

import numpy as np

import instances as I
import solvers as S

HERE = os.path.dirname(os.path.abspath(__file__))
SWEEPS = [10, 32, 100, 316, 1000, 3162, 10000]
RUNS = {"sa": 40, "pt": 20}
POINT_BUDGET = 60.0  # seconds for all runs at one run length; a longer point ends the grid
EXACT_LIMIT = 60.0
checks = []


def check(what, ok, detail):
    checks.append({"what": what, "ok": bool(ok), "detail": detail})


def tts(solver, d, Sm, c, target, decode_ok=lambda x: True):
    """TTS* over the run-length grid; a run succeeds if it reaches the target value (and decodes validly)."""
    fn = S.simulated_annealing if solver == "sa" else S.parallel_tempering
    fn(d, Sm, c, 2, 0)  # compile
    points = []
    for t in SWEEPS:
        ok, secs = 0, []
        t0 = time.perf_counter()
        for r in range(RUNS[solver]):
            a = time.perf_counter()
            v, x = fn(d, Sm, c, t, 1000 * t + r)
            secs.append(time.perf_counter() - a)
            ok += (v <= target + 1e-9) and decode_ok(x)
        p = ok / RUNS[solver]
        run = float(np.median(secs))
        ttsv = run if p == 1 else (math.inf if p == 0 else run * math.log(0.01) / math.log(1 - p))
        points.append({"sweeps": t, "p": p, "run_seconds": run, "tts": ttsv})
        if time.perf_counter() - t0 > POINT_BUDGET or p == 1:
            break
    best = min(points, key=lambda q: q["tts"])
    return {"tts_star": best["tts"], "at_sweeps": best["sweeps"], "points": points}


def solve_qubo_all(name, Q, c, n, target, decode_ok=lambda x: True):
    d, Sm = S.dense(Q, n)
    q = S.exact_qubo(Q, c, n, EXACT_LIMIT)
    out = {"variables": n, "terms": len(Q),
           "cpsat_on_qubo": {"status": q["status"], "seconds": q["seconds"], "value": q["value"],
                             "reached_optimum": q["value"] is not None and q["value"] <= target + 1e-6}}
    for s in ("sa", "pt"):
        out[s] = tts(s, d, Sm, c, target, decode_ok)
    return out


def haplotypes():
    rows = []
    for n in (16, 32, 64, 128, 256):
        h = I.haplotype(n, seed=n)
        e = S.exact_haplotype(h, EXACT_LIMIT)
        check(f"haplotype n={n}: CP-SAT proves the optimum", e["status"] == "optimal", e["status"])
        check(f"haplotype n={n}: the optimum's phasing is close to the truth", h.switch_errors(e["spins"]) <= 2,
              f"{h.switch_errors(e['spins'])} switch errors, MEC {h.mec(e['spins'])}")
        Q, c = I.ising_to_qubo(h.J)
        target = e["energy"]
        r = {"size": n, "exact": {"status": e["status"], "seconds": e["seconds"], "energy": int(target),
                                  "switch_errors": h.switch_errors(e["spins"]), "mec": h.mec(e["spins"])}}
        r.update(solve_qubo_all("haplotype", Q, c, n, target))
        rows.append(r)
        print("haplotype", n, "exact", e["status"], f"{e['seconds']:.3f}s", "SA TTS*", f"{r['sa']['tts_star']:.3g}",
              "PT TTS*", f"{r['pt']['tts_star']:.3g}", flush=True)
    return rows


def path_from_positions(x, n):
    X = x.reshape(n, n)  # X[v, p]
    if not (np.all(X.sum(axis=0) == 1) and np.all(X.sum(axis=1) == 1)):
        return None
    return [int(np.argmax(X[:, p])) for p in range(n)]


def path_from_edges(x, E, n):
    nxt = {}
    for k, e in enumerate(E):
        if x[k] > 0.5:
            if e[0] in nxt:
                return None
            nxt[e[0]] = e[1]
    starts = set(range(n)) - set(nxt.values())
    if len(starts) != 1:
        return None
    path = [starts.pop()]
    while path[-1] in nxt and len(path) <= n:
        path.append(nxt[path[-1]])
    return path


def assemblies():
    g = I.read_fasta(os.path.join(HERE, "data", "phix174.fasta"))
    rows = []
    for n in (6, 8, 12, 16, 24):
        a = I.assembly(g, n, seed=n)
        e = S.exact_layout(a, EXACT_LIMIT)
        check(f"assembly n={n}: the native optimum is the reads' true order", e["path"] == a.order,
              f"{e['status']}, overlap {e['overlap']}")
        L = max(len(r) for r in a.reads)
        target = -e["overlap"] / L
        r = {"size": n, "overlap_edges": len(a.edges),
             "exact_native": {"status": e["status"], "seconds": e["seconds"], "overlap": e["overlap"]}}
        Q, c = I.layout_qubo_positions(a)
        r["qubo_positions"] = solve_qubo_all("positions", Q, c, n * n, target,
                                             lambda x: path_from_positions(x, n) == a.order)
        Q, c, E = I.layout_qubo_edges(a)
        r["qubo_edges"] = solve_qubo_all("edges", Q, c, len(E), target,
                                         lambda x: path_from_edges(x, E, n) == a.order)
        rows.append(r)
        print("assembly", n, "native", e["status"], f"{e['seconds']:.3f}s",
              "positions SA", f"{r['qubo_positions']['sa']['tts_star']:.3g}",
              "edges SA", f"{r['qubo_edges']['sa']['tts_star']:.3g}", flush=True)
    return rows


def modules():
    net = I.string_network(os.path.join(HERE, "data", "string_dna_repair.tsv"), 0.9)
    rows = []
    for k in (3, 5, 8, 10):
        e = S.exact_dks(net, k, EXACT_LIMIT)
        check(f"dense {k}-subgraph: CP-SAT proves the optimum", e["status"] == "optimal", e["status"])
        Q, c = I.dks_qubo(net, k)
        r = {"size": k, "genes": len(net.genes), "network_edges": len(net.edges),
             "exact": {"status": e["status"], "seconds": e["seconds"], "edges": e["edges"],
                       "module": [net.genes[i] for i in e["chosen"]]}}
        r.update(solve_qubo_all("dks", Q, c, len(net.genes), -e["edges"]))
        rows.append(r)
        print("dks", k, e["edges"], r["exact"]["module"], "SA TTS*", f"{r['sa']['tts_star']:.3g}", flush=True)
    # an open instance, reported without a check: CP-SAT does not prove k = 12 within the limit
    e = S.exact_dks(net, 12, 300.0)
    rows.append({"size": 12, "open": True, "exact": {"status": e["status"], "seconds": e["seconds"],
                 "edges": e["edges"], "module": [net.genes[i] for i in e["chosen"]]}})
    print("dks 12 (open)", e["status"], e["edges"], flush=True)
    return rows


if __name__ == "__main__":
    out = sys.argv[1]
    rec = {"haplotype": haplotypes(), "assembly": assemblies(), "modules": modules(), "checks": checks}
    json.dump(rec, open(os.path.join(out, "g0.json"), "w"), indent=1, default=float)
    bad = [c for c in checks if not c["ok"]]
    print(f"{len(checks) - len(bad)} of {len(checks)} checks hold")
    for c in bad:
        print("FAILS:", c["what"], c["detail"])
    sys.exit(1 if bad else 0)
