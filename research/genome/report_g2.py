"""Tables of a G2 run (g2.py merge's JSON) as Markdown.
    python report_g2.py results/g2-<host>-<date>.json > RESULTS_G2.md"""
import json
import math
import os
import sys

path = sys.argv[1]
r = json.load(open(path))
NAMES = {"sa": "SA", "pt": "PT", "simcim": "SimCIM", "ahc": "AHC (CIM)"}
ORDER = ["sa", "pt", "simcim", "ahc"]


def s(x):
    if x is None or not math.isfinite(x):
        return "not reached"
    return f"{x:.2e}" if x < 1e-2 else (f"{x:.3f}" if x < 10 else f"{x:.1f}")


def cell(res):
    v = s(res["tts_star"])
    return v + (" *" if res.get("at_grid_edge") else "")


ok = sum(c["ok"] for c in r["checks"])
out = ["# G2: coherent-Ising-machine solvers against SA and PT, latest run", "",
       f"From [`results/{os.path.basename(path)}`](results/{os.path.basename(path)}), written by `g2_run.sh`. The "
       "method is in [README.md](README.md). TTS\\* in seconds on one core: the expected time to reach the exact "
       "optimum with 99% probability, at the best run length. \\* marks a best run length at the edge of the grid "
       "(an upper bound).", "", f"**Checks: {ok} of {len(r['checks'])} hold.**", "",
       "## Haplotype phasing from short reads (2–6 sites per read; 3 instances per size)", "",
       "The couplings form a band 5 sites wide, so a dynamic programme over the band solves it exactly in linear time "
       "(the *DP* column, one core).", "",
       "| sites | instance | exact: CP-SAT, 8 workers | exact: DP, 1 core | " + " | ".join(NAMES[x] for x in ORDER) + " |",
       "| --- | --- | --- | --- | " + " | ".join("---" for _ in ORDER) + " |"]
for name, e in r["exact"].items():
    if name.startswith("hap"):
        _, n, seed = name.split("-")
        out.append(f"| {n} | {seed} | {s(e['seconds'])} | {s(e['dp']['seconds'])} | " +
                   " | ".join(cell(r["solvers"][x]["results"][name]) for x in ORDER) + " |")
out += ["", "## Haplotype phasing from long reads (10–40 sites per read; 3 instances per size)", "",
        "The band is up to 39 sites wide: too wide for the dynamic programme above.", "",
        "| sites | instance | exact: CP-SAT, 8 workers | switch errors of the optimum | " +
        " | ".join(NAMES[x] for x in ORDER) + " |",
        "| --- | --- | --- | --- | " + " | ".join("---" for _ in ORDER) + " |"]
for name, e in r["exact"].items():
    if name.startswith("long"):
        _, n, seed = name.split("-")
        out.append(f"| {n} | {seed} | {s(e['seconds'])} | {e['switch_errors']} | " +
                   " | ".join(cell(r["solvers"][x]["results"][name]) for x in ORDER) + " |")
out += ["", "## How the time grows (median over the 3 instances per size)", "",
        "| family | method | median time at the three largest sizes | factor per doubling of the sites |",
        "| --- | --- | --- | --- |"]
for fam, label in (("short", "short reads, 64 / 128 / 256 sites"), ("long", "long reads, 128 / 256 / 512 sites")):
    for x in ["exact"] + ORDER:
        f = r["fits"][fam][x]
        nm = "exact (CP-SAT, 8 workers)" if x == "exact" else NAMES[x] + " TTS\\*"
        med = " / ".join(s(v) for v in f["median"])
        out.append(f"| {label} | {nm} | {med} | " + (f"{f['factor_per_doubling']:.1f}" if "factor_per_doubling" in f
                                                       else "—") + " |")
out += ["", "## φX174 assembly, Boev et al.'s QUBO with one variable per overlap edge", "",
        "| reads | native exact | " + " | ".join(NAMES[x] for x in ORDER) + " |",
        "| --- | --- | " + " | ".join("---" for _ in ORDER) + " |"]
for name, e in r["exact"].items():
    if name.startswith("asm"):
        out.append(f"| {name.split('-')[1]} | {s(e['seconds'])} | " +
                   " | ".join(cell(r["solvers"][x]["results"][name]) for x in ORDER) + " |")
out += ["", "## Tuning: the parameters chosen (best geometric-mean TTS\\* on separate tuning instances)", "",
        "| solver | short reads | long reads | assembly | grid |", "| --- | --- | --- | --- | --- |"]
for x in ORDER:
    c = r["solvers"][x]["chosen"]
    out.append(f"| {NAMES[x]} | {c['hap']} | {c['hapL']} | {c['asm']} | {len(r['solvers'][x]['tuning'])} settings |")
out += ["", "## Every check", "", "| check | holds | detail |", "| --- | --- | --- |"]
out += [f"| {c['what']} | {'yes' if c['ok'] else '**NO**'} | {c['detail']} |" for c in r["checks"]]
print("\n".join(out))
