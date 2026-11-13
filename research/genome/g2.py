"""G2: coherent-Ising-machine solvers (SimCIM, AHC; cim.py) against simulated annealing and parallel tempering
(solvers.py), all tuned with the same budget, by the optimal time to solution TTS* on one core, on haplotype phasing
(16-256 sites, 3 instances per size) and on Boev et al.'s assembly QUBO (one variable per overlap edge).
    python g2.py exact <out>            # exact optima -> <out>/exact.json
    python g2.py run <solver> <out>     # tune on separate instances, then evaluate -> <out>/<solver>.json
    python g2.py merge <out> <json>     # checks and fits -> the run's record
"""
import itertools
import json
import math
import os
import sys
import time

import numpy as np

import cim as C
import instances as I
import solvers as S

HERE = os.path.dirname(os.path.abspath(__file__))
HAP_SIZES = (16, 32, 64, 128, 256)  # short reads (2-6 sites): a narrow band, solved exactly by dynamic programming
LONG_SIZES = (64, 128, 256, 512)     # long reads (10-40 sites): the band is too wide for that
LONG_SPAN = (10, 40)
HAP_SEEDS = (0, 1, 2)                 # instance seed = size + offset; offset 0 is G0's instance
ASM_READS = (6, 8, 12, 16)
TUNE = {"hap": [(64, 1000), (64, 1001), (64, 1002)], "hapL": [(128, 2000), (128, 2001), (128, 2002)],
        "asm": [(8, 1008)]}
LENGTHS = {"sa": [10, 32, 100, 316, 1000, 3162, 10000], "pt": [10, 32, 100, 316, 1000, 3162, 10000],
           "simcim": [100, 316, 1000, 3162, 10000, 31623], "ahc": [100, 316, 1000, 3162, 10000, 31623]}
GRID = {
    "sa": [{"b0": a, "b1": b} for a, b in itertools.product((0.5, 1.0), (0.25, 0.5, 1.0, 2.0))],
    "pt": [{"replicas": r, "b1": b} for r, b in itertools.product((8, 16), (0.5, 1.0, 2.0, 4.0))],
    "simcim": [{"amp": a, "zeta": z, "sigma": s} for a, z, s in itertools.product((0.1, 0.3), (0.3, 1.0), (0.03, 0.15))],
    "ahc": [{"p": p, "eps": e, "beta": b} for p, e, b in itertools.product((0.0, 0.5), (0.3, 1.0), (0.1, 0.4))],
}
RUNS_TUNE, RUNS_EVAL = 20, 30
POINT_BUDGET = 60.0


# --- instances: each as an Ising model (for the CIM solvers) and a QUBO (for SA and PT), with a common energy ------

def hap_instance(n, seed, long=False):
    h = I.haplotype(n, seed=seed, span=LONG_SPAN) if long else I.haplotype(n, seed=seed)
    Q, c = I.ising_to_qubo(h.J)
    return {"kind": "hapL" if long else "hap", "n": n, "seed": seed, "J": h.J, "h": {}, "offset": 0.0, "Q": Q,
            "c": c, "nvar": n, "obj": h}


def asm_instance(reads, seed):
    g = I.read_fasta(os.path.join(HERE, "data", "phix174.fasta"))
    a = I.assembly(g, reads, seed=seed)
    Q, c, E = I.layout_qubo_edges(a)
    J, h, off = I.qubo_to_ising(Q, c)
    return {"kind": "asm", "n": reads, "seed": seed, "J": J, "h": h, "offset": off, "Q": Q, "c": c, "nvar": len(E),
            "obj": a}


def instances():
    for n in HAP_SIZES:
        for o in HAP_SEEDS:
            yield f"hap-{n}-{n + o}", hap_instance(n, n + o)
    for n in LONG_SIZES:
        for o in HAP_SEEDS:
            yield f"long-{n}-{n + o}", hap_instance(n, n + o, long=True)
    for r in ASM_READS:
        yield f"asm-{r}-{r}", asm_instance(r, r)


def tuning_instances():
    for n, s in TUNE["hap"]:
        yield hap_instance(n, s)
    for n, s in TUNE["hapL"]:
        yield hap_instance(n, s, long=True)
    for r, s in TUNE["asm"]:
        yield asm_instance(r, s)


def exact(inst):
    if inst["kind"] in ("hap", "hapL"):
        # 8 workers: on one, CP-SAT once stopped at a far-from-optimal solution after 300 s on an instance it proves
        # in 0.1 s with 8; the exact phase runs alone, before the solvers
        e = S.exact_haplotype(inst["obj"], 600.0, workers=8)
        out = {"status": e["status"], "seconds": e["seconds"], "energy": float(e["energy"]),
               "switch_errors": int(inst["obj"].switch_errors(e["spins"]))}
        if inst["kind"] == "hap":
            dp = S.exact_banded(inst["n"], inst["J"])
            out["dp"] = {"seconds": dp["seconds"], "energy": dp["energy"], "band": dp["band"]}
        return out
    e = S.exact_layout(inst["obj"], 60.0)
    return {"status": e["status"], "seconds": e["seconds"], "energy": -e["overlap"] / 100.0}


# --- running a solver ------------------------------------------------------------------------------------------------

def make(solver, inst, params):
    """A function (length, seed) -> energy in the instance's common units."""
    if solver in ("sa", "pt"):
        d, Sm = S.dense(inst["Q"], inst["nvar"])
        b0, b1 = S.default_betas(d, Sm)
        betas = (b0 * params.get("b0", 1.0), b1 * params["b1"])
        if solver == "sa":
            return lambda t, s: S.simulated_annealing(d, Sm, inst["c"], t, s, betas=betas)[0]
        return lambda t, s: S.parallel_tempering(d, Sm, inst["c"], t, s, replicas=params["replicas"], betas=betas)[0]
    m = C.Ising(inst["nvar"], inst["J"], inst["h"])
    f = m.simcim if solver == "simcim" else m.ahc
    return lambda t, s: f(t, s, **params)[0] + inst["offset"]


def tts(fn, target, lengths, runs):
    fn(lengths[0], 0)  # compile
    points = []
    for t in lengths:
        ok, secs = 0, []
        t0 = time.perf_counter()
        for r in range(runs):
            a = time.perf_counter()
            e = fn(t, 7919 * t + r)
            secs.append(time.perf_counter() - a)
            ok += e <= target + 1e-6
        p = ok / runs
        run = float(np.median(secs))
        v = run if p == 1 else (math.inf if p == 0 else run * math.log(0.01) / math.log(1 - p))
        points.append({"length": t, "p": p, "run_seconds": run, "tts": v})
        if p == 1 or time.perf_counter() - t0 > POINT_BUDGET:
            break
    best = min(points, key=lambda q: q["tts"])
    edge = best["length"] in (lengths[0], lengths[-1]) and math.isfinite(best["tts"])
    return {"tts_star": best["tts"], "at_length": best["length"], "at_grid_edge": edge, "points": points}


def geo(xs):
    xs = [x for x in xs]
    return math.inf if any(not math.isfinite(x) for x in xs) else math.exp(sum(math.log(x) for x in xs) / len(xs))


def run(solver, out):
    tune_set = list(tuning_instances())
    targets = [exact(t)["energy"] for t in tune_set]
    tuning = []
    for params in GRID[solver]:
        scores = {"hap": [], "hapL": [], "asm": []}
        for inst, tgt in zip(tune_set, targets):
            r = tts(make(solver, inst, params), tgt, LENGTHS[solver], RUNS_TUNE)
            scores[inst["kind"]].append(r["tts_star"])
        tuning.append({"params": params, **{k: geo(v) for k, v in scores.items()}})
        print(solver, "tune", params, {k: f"{tuning[-1][k]:.3g}" for k in scores}, flush=True)
    best = {k: min(tuning, key=lambda t: t[k])["params"] for k in ("hap", "hapL", "asm")}
    ex = json.load(open(os.path.join(out, "exact.json")))
    res = {}
    for name, inst in instances():
        r = tts(make(solver, inst, best[inst["kind"]]), ex[name]["energy"], LENGTHS[solver], RUNS_EVAL)
        res[name] = r
        print(solver, name, f"TTS* {r['tts_star']:.3g} at {r['at_length']}", flush=True)
    json.dump({"solver": solver, "tuning": tuning, "chosen": best, "results": res},
              open(os.path.join(out, f"{solver}.json"), "w"), indent=1, default=float)


def merge(out, dest):
    ex = json.load(open(os.path.join(out, "exact.json")))
    solvers = {s: json.load(open(os.path.join(out, f"{s}.json"))) for s in GRID}
    checks = []
    for name, e in ex.items():
        checks.append({"what": f"{name}: the exact optimum is proven", "ok": e["status"] == "optimal",
                       "detail": f"{e['status']}, {e['seconds']:.2f} s"})
        if "dp" in e:
            checks.append({"what": f"{name}: the banded dynamic programme gives the same optimum",
                           "ok": abs(e["dp"]["energy"] - e["energy"]) < 1e-6,
                           "detail": f"{e['dp']['energy']} in {e['dp']['seconds']:.3f} s, band {e['dp']['band']}"})
    # scaling: median TTS* per size over the largest three sizes of each family, fitted as n^q
    fits = {}
    for fam, prefix, sizes in (("short", "hap", (64, 128, 256)), ("long", "long", (128, 256, 512))):
        fits[fam] = {}
        for s, rec in list(solvers.items()) + [("exact", None)]:
            med = []
            for n in sizes:
                if rec is None:
                    vals = [ex[f"{prefix}-{n}-{n + o}"]["seconds"] for o in HAP_SEEDS]
                else:
                    vals = [rec["results"][f"{prefix}-{n}-{n + o}"]["tts_star"] for o in HAP_SEEDS]
                med.append(float(np.median(vals)))
            if all(math.isfinite(m) for m in med):
                q = np.polyfit(np.log(sizes), np.log(med), 1)[0]
                fits[fam][s] = {"sizes": list(sizes), "median": med, "factor_per_doubling": float(2 ** q)}
            else:
                fits[fam][s] = {"sizes": list(sizes), "median": med, "note": "not reached at some size"}
    json.dump({"exact": ex, "solvers": solvers, "fits": fits, "checks": checks}, open(dest, "w"), indent=1,
              default=float)
    bad = [c for c in checks if not c["ok"]]
    print(f"{len(checks) - len(bad)} of {len(checks)} checks hold")
    return not bad


if __name__ == "__main__":
    mode, out = sys.argv[1], sys.argv[-1] if sys.argv[1] != "merge" else sys.argv[2]
    if mode == "exact":
        ex = {}
        for name, inst in instances():
            ex[name] = exact(inst)
            print(name, ex[name], flush=True)
        json.dump(ex, open(os.path.join(out, "exact.json"), "w"), indent=1)
    elif mode == "run":
        run(sys.argv[2], out)
    else:
        sys.exit(0 if merge(sys.argv[2], sys.argv[3]) else 1)
