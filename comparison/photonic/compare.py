"""Checks the answers of the photonic head-to-head against each other and writes the tables (README.md).
    python compare.py <dir with inputs.json and <tool>.json> <record.json> <RESULTS.md> <history.csv> <meta.json>
Exits 1 if any answer disagrees.
"""
import csv
import json
import math
import os
import sys
from itertools import combinations_with_replacement

import numpy as np
from scipy.stats import chisquare
import thewalrus as tw

work, record_path, results_md, history_csv, meta_path = sys.argv[1:]
meta = json.load(open(meta_path))
inp = json.load(open(os.path.join(work, "inputs.json")))
TOOLS = ["pqcf", "thewalrus", "piquasso", "perceval", "deepquantum", "strawberryfields"]
NAMES = {"pqcf": "PQCF", "thewalrus": "The Walrus", "piquasso": "Piquasso", "perceval": "Perceval",
         "deepquantum": "DeepQuantum", "strawberryfields": "Strawberry Fields"}
res = {t: json.load(open(os.path.join(work, f"{t}.json"))) for t in TOOLS
       if os.path.exists(os.path.join(work, f"{t}.json"))}
checks = []  # (task, what, ok, detail)


def check(task, what, ok, detail):
    checks.append({"task": task, "what": what, "ok": bool(ok), "detail": detail})


def by_size(tool, task):
    return {r["size"]: r for r in res.get(tool, {}).get(task, [])}


# --- answers -----------------------------------------------------------------------------------------------------

# permanents and torontonians: every program against The Walrus, the reference both PQCF and Piquasso test against
for task, tol in (("permanent", 1e-9), ("torontonian", 1e-7)):
    ref = by_size("thewalrus", task)
    for t in TOOLS:
        if t == "thewalrus":
            continue
        for size, r in by_size(t, task).items():
            if size not in ref or "value" not in r:
                continue
            a, b = r["value"], ref[size]["value"]
            a = complex(*a) if isinstance(a, list) else a
            b = complex(*b) if isinstance(b, list) else b
            rel = abs(a - b) / abs(b)
            check(task, f"{NAMES[t]}, size {size}", rel <= tol, f"relative difference {rel:.1e} (limit {tol:.0e})")

# boson sampling: every program's 20,000 samples against the exact output distribution, chi-square
bc = inp["boson_check"]
U = np.array([[complex(*z) for z in row] for row in bc["unitary"]])
n, m = bc["photons"], bc["modes"]
pats = [tuple(np.bincount(c, minlength=m)) for c in combinations_with_replacement(range(m), n)]
exact = {}
for p in pats:
    rows = [i for i in range(m) for _ in range(p[i])]
    sub = U[np.ix_(rows, list(range(n)))]
    exact[p] = abs(tw.perm(sub)) ** 2 / math.prod(math.factorial(k) for k in p)
assert abs(sum(exact.values()) - 1) < 1e-12
for t in TOOLS:
    c = res.get(t, {}).get("boson_check")
    if not c:
        continue
    got = {tuple(x["pattern"]): x["count"] for x in c}
    outside = sum(v for k, v in got.items() if k not in exact)
    obs = [got.get(p, 0) for p in pats]
    stat, pval = chisquare(obs, [exact[p] * bc["shots"] for p in pats])
    check("boson_sampling", f"{NAMES[t]}, {n} photons in {m} modes, {bc['shots']} samples",
          pval > 1e-3 and outside == 0, f"chi-square p = {pval:.3f}, {len(pats) - 1} degrees of freedom")

# Fock: the same truncation must give the same state; per-mode truncations against each other
fk = {t: {r["cutoff"]: r for r in res[t].get("fock", [])} for t in res if res[t].get("fock")}
for t in fk:
    for cut, r in fk[t].items():
        if t == "pqcf" or r["truncation"] != "total" or cut not in fk.get("pqcf", {}):
            continue
        d = max(abs(a - b) for a, b in zip(r["mean_photons"], fk["pqcf"][cut]["mean_photons"]))
        check("fock", f"{NAMES[t]} against PQCF, total photons <= {cut}", d < 1e-9,
              f"largest difference of a mode's mean photon number {d:.1e}")
if "strawberryfields" in fk and "deepquantum" in fk:
    for cut in fk["strawberryfields"]:
        if cut in fk["deepquantum"]:
            d = max(abs(a - b) for a, b in zip(fk["strawberryfields"][cut]["mean_photons"],
                                               fk["deepquantum"][cut]["mean_photons"]))
            check("fock", f"DeepQuantum against Strawberry Fields, {cut} photons per mode", d < 1e-6,
                  f"largest difference of a mode's mean photon number {d:.1e}")
# the truncation error of each result: its distance from PQCF's largest cutoff
best = fk["pqcf"][max(fk["pqcf"])]["mean_photons"] if "pqcf" in fk else None

# --- tables ------------------------------------------------------------------------------------------------------


def secs(x):
    if x is None:
        return "—"
    return f"{x:.2e}" if x < 1e-3 else (f"{x:.4f}" if x < 1 else f"{x:.2f}")


def table(task, sizes, label, value):
    tools = [t for t in TOOLS if res.get(t, {}).get(task)]
    lines = [f"| {label} | " + " | ".join(NAMES[t] for t in tools) + " |",
             "| --- | " + " | ".join("---" for _ in tools) + " |"]
    for s in sizes:
        cells = []
        p = by_size("pqcf", task).get(s)
        for t in tools:
            r = by_size(t, task).get(s)
            if r is None:
                cells.append("not run (budget)")
            elif r.get("not_run"):
                cells.append(f"not run ({r['not_run']})")
            else:
                v = value(r)
                if t != "pqcf" and p and not p.get("not_run"):
                    ratio = value(r) / value(p) if task != "boson_sampling" else value(p) / value(r)
                    cells.append(f"{secs(v) if task != 'boson_sampling' else f'{v:,.0f}'} ({ratio:.2g}×)")
                else:
                    cells.append(secs(v) if task != "boson_sampling" else f"{v:,.0f}")
        lines.append(f"| {s} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


ver = ", ".join(f"{NAMES[t]} {res[t]['version']}" for t in TOOLS if t in res)
md = [f"# Photonic head-to-head: latest run",
      "",
      f"{meta['date']} on {meta['host']} ({meta['cpu']}), every program pinned to one core with one thread. PQCF "
      f"at commit {meta['commit']}. {ver}. Written by `run.sh`; the method is in [README.md](README.md).",
      "",
      f"**Answers: {sum(c['ok'] for c in checks)} of {len(checks)} checks agree.**",
      "",
      "## Permanent of an n × n block of a random unitary (seconds; in brackets, the time relative to PQCF's)",
      "",
      table("permanent", [c["n"] for c in inp["permanent"]], "n", lambda r: r["seconds"]),
      "",
      "## Torontonian of a k-mode Gaussian state, every mode clicking (seconds)",
      "",
      table("torontonian", [c["k"] for c in inp["torontonian"]], "k", lambda r: r["seconds"]),
      "",
      "## Boson sampling: n photons in n² modes, exact samples per second (in brackets, PQCF's rate relative to "
      "the program's)",
      "",
      table("boson_sampling", [c["photons"] for c in inp["boson_sampling"]], "photons",
            lambda r: r["samples_per_second"]),
      "",
      "## Fock simulation: one layer of a continuous-variable neural network on 10 modes (seconds)",
      "",
      "PQCF and Piquasso keep every state with at most c photons in total; Strawberry Fields and DeepQuantum keep at "
      "most c photons in each mode. *Truncation error*: the largest difference of a mode's mean photon number from "
      "PQCF's largest cutoff.",
      "",
      "| program | truncation | cutoff | states | seconds | truncation error |",
      "| --- | --- | --- | --- | --- | --- |"]
for t in TOOLS:
    for cut, r in sorted(fk.get(t, {}).items()):
        err = max(abs(a - b) for a, b in zip(r["mean_photons"], best)) if best else float("nan")
        md.append(f"| {NAMES[t]} | {r['truncation']} | {cut} | {r['dimension']:,} | {secs(r['seconds'])} | "
                  f"{err:.1e} |")
md += ["", "## Every check", "", "| task | what | agrees | detail |", "| --- | --- | --- | --- |"]
md += [f"| {c['task']} | {c['what']} | {'yes' if c['ok'] else '**NO**'} | {c['detail']} |" for c in checks]
open(results_md, "w").write("\n".join(md) + "\n")

json.dump({**meta, "versions": {t: res[t]["version"] for t in res}, "results": res, "checks": checks},
          open(record_path, "w"), indent=1)

# history: one row per program, task and size
new = not os.path.exists(history_csv)
with open(history_csv, "a", newline="") as f:
    w = csv.writer(f)
    if new:
        w.writerow(["date", "host", "commit", "program", "version", "task", "size", "seconds"])
    for t in res:
        for task in ("permanent", "torontonian", "boson_sampling"):
            for r in res[t].get(task, []):
                if r.get("seconds") is not None:
                    w.writerow([meta["date"], meta["host"], meta["commit"], t, res[t]["version"], task, r["size"],
                                f"{r['seconds']:.6g}"])
        for r in res[t].get("fock", []):
            w.writerow([meta["date"], meta["host"], meta["commit"], t, res[t]["version"],
                        f"fock_{r['truncation'].replace(' ', '_')}", r["cutoff"], f"{r['seconds']:.6g}"])

bad = [c for c in checks if not c["ok"]]
print(f"{len(checks) - len(bad)} of {len(checks)} checks agree")
for c in bad:
    print("DISAGREES:", c["task"], c["what"], c["detail"])
sys.exit(1 if bad else 0)
