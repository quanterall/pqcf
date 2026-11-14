"""The quantum optics head-to-head's tables (README.md): every tool's answers against the closed forms (mpmath, 50
digits) and against each other, and the times. Writes the run's JSON, RESULTS.md, and a line per tool and task to
history.csv; exits 1 if a tool's answer is wrong by more than its tolerance.
    python compare.py out_dir results.json RESULTS.md history.csv meta.json tasks.json
"""
import json
import math
import os
import sys

import mpmath as mp

mp.mp.dps = 50
TOOLS = ["pqcf", "qutip-default", "qutip-tight"]
# an answer more than this far from the closed form is an error, not a precision
TOLERANCE = 1e-4


def coherent_weights(n, alpha):
    """|alpha>'s Poisson weights truncated at n photons and renormalised."""
    w = [mp.mpf(alpha) ** (2 * k) / mp.factorial(k) for k in range(n + 1)]
    s = sum(w)
    return [x / s for x in w]


def references(t):
    r = {}
    rb = t["rabi"]
    r["rabi"] = [mp.cos(rb["g"] * mp.mpf(x)) ** 2 for x in rb["times"]]
    rv = t["revival"]
    w = coherent_weights(rv["n"] - 1, rv["alpha"])
    r["revival"] = [sum(p * mp.cos(rv["g"] * mp.sqrt(k + 1) * mp.mpf(x)) ** 2 for k, p in enumerate(w))
                    for x in rv["times"]]
    cd = t["cavity_decay"]
    w = coherent_weights(cd["n"], cd["alpha"])
    n0 = sum(k * p for k, p in enumerate(w))
    r["cavity_decay"] = [n0 * mp.exp(-mp.mpf(cd["kappa"]) * mp.mpf(x)) for x in cd["times"]]
    fl = t["fluorescence"]
    g = mp.mpf(fl["gamma"])
    r["fluorescence"] = [(mp.mpf(om) ** 2 / 4) / (mp.mpf(dl) ** 2 + g ** 2 / 4 + mp.mpf(om) ** 2 / 2)
                         for om, dl in fl["cases"]]
    return r


def digits(err):
    return 99 if err == 0 else min(99, int(math.floor(-math.log10(float(err)))))


def main():
    out, res_path, md_path, hist_path, meta_path, tasks_path = sys.argv[1:7]
    tasks = json.load(open(tasks_path))
    meta = json.load(open(meta_path))
    runs = {k: json.load(open(os.path.join(out, k + ".json"))) for k in TOOLS if os.path.exists(os.path.join(out, k + ".json"))}
    refs = references(tasks)
    ok = True
    table = {}
    for name in ["rabi", "revival", "cavity_decay", "fluorescence"]:
        for k, run in runs.items():
            got = run["tasks"][name]
            err = max(abs(mp.mpf(v) - ref) for v, ref in zip(got["values"], refs[name]))
            table[(name, k)] = {"digits": digits(err), "error": float(err), "seconds": got["seconds"]}
            if err > TOLERANCE or len(got["values"]) != len(refs[name]):
                ok = False
    # the driven damped model has no closed form: QuTiP against PQCF at 128 bits
    sizes = {}
    for k, run in runs.items():
        for row in run["tasks"]["damped_jc"]:
            sizes.setdefault(row["size"], {})[k] = row
    jc = []
    for n in sorted(sizes):
        row = {"size": n, "dim": 2 * (n + 1)}
        base = sizes[n].get("pqcf")
        for k in runs:
            if k in sizes[n]:
                row[k] = {"seconds": sizes[n][k]["seconds"]}
                if base and k != "pqcf":
                    d = max(abs(mp.mpf(a) - mp.mpf(b)) for a, b in zip(sizes[n][k]["values"], base["values"]))
                    row[k]["diff"] = float(d)
                    if d > TOLERANCE:
                        ok = False
        jc.append(row)
    result = {"meta": meta, "versions": {k: r["version"] for k, r in runs.items()},
              "closed_form": [{"task": t, "tool": k, **v} for (t, k), v in table.items()], "damped_jc": jc}
    json.dump(result, open(res_path, "w"), indent=1)
    write_md(md_path, result, runs, table, jc)
    new = not os.path.exists(hist_path)
    with open(hist_path, "a") as f:
        if new:
            f.write("date,host,commit,tool,task,digits,seconds\n")
        for (t, k), v in table.items():
            f.write(f"{meta['date']},{meta['host']},{meta['commit']},{k},{t},{v['digits']},{v['seconds']:.6g}\n")
    sys.exit(0 if ok else 1)


def secs(s):
    return f"{s:.3g}"


def write_md(path, result, runs, table, jc):
    m = result["meta"]
    names = {k: runs[k]["tool"] for k in runs}
    L = ["# Quantum optics head-to-head: latest run", "",
         f"**Machine:** {m['host']}, {m['cpu']}, one core ({m['core']}). **Date:** {m['date']}. "
         f"**PQCF:** {m['commit']}. **QuTiP:** {result['versions'].get('qutip-default', '-')}.", "",
         "How the tasks are set and checked: [README.md](README.md). Every run: [`results/`](results), "
         "[history.csv](history.csv).", "",
         "## Against the closed forms", "",
         "Correct decimal digits (the worst point of each task) and seconds (the median of repeated runs).", "",
         "| task | " + " | ".join(names[k] for k in runs) + " |",
         "| --- | " + " | ".join("---" for _ in runs) + " |"]
    for t in ["rabi", "revival", "cavity_decay", "fluorescence"]:
        cells = [f"{table[(t, k)]['digits']} digits, {secs(table[(t, k)]['seconds'])} s" for k in runs]
        L.append(f"| {t} | " + " | ".join(cells) + " |")
    L += ["", "## Driven, damped Jaynes-Cummings model", "",
          "No closed form: QuTiP's answers against PQCF's at 128 bits (the largest difference of <a^dagger a> and "
          "P_e at t = 5), and the seconds of one run. Sizes are the photon cutoff N; the Hilbert space has 2(N + 1) "
          "states and the Liouvillian 4(N + 1)^2.", "",
          "| N | dimension | " + " | ".join(names[k] for k in runs) + " |",
          "| --- | --- | " + " | ".join("---" for _ in runs) + " |"]
    for row in jc:
        cells = []
        for k in runs:
            if k not in row:
                cells.append("-")
            elif "diff" in row[k]:
                cells.append(f"{secs(row[k]['seconds'])} s (differs by {row[k]['diff']:.1e})")
            else:
                cells.append(f"{secs(row[k]['seconds'])} s")
        L.append(f"| {row['size']} | {row['dim']} | " + " | ".join(cells) + " |")
    open(path, "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
