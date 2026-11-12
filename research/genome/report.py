"""Tables of a G0 run (g0.py's JSON) as Markdown.
    python report.py results/g0-<host>-<date>.json > RESULTS.md"""
import json
import math
import os
import sys

path = sys.argv[1]
r = json.load(open(path))


def s(x):
    if x is None or (isinstance(x, float) and math.isinf(x)):
        return "not reached"
    return f"{x:.2e}" if x < 1e-2 else (f"{x:.3f}" if x < 10 else f"{x:.1f}")


def tts(block):
    t = block["tts_star"]
    return "not reached" if math.isinf(t) else f"{s(t)} ({block['at_sweeps']} sweeps)"


def qubo_cell(q):
    c = q["cpsat_on_qubo"]
    return ("optimum, " if c["reached_optimum"] else "**not reached**, ") + \
        (f"proven, {s(c['seconds'])} s" if c["status"] == "optimal" else f"{c['status']} at the {s(c['seconds'])} s limit")


name = os.path.basename(path)
ok = sum(c["ok"] for c in r["checks"])
out = [f"# G0: classical baselines on genomic instances, latest run",
       "",
       f"From [`results/{name}`](results/{name}), written by `run.sh` (one core, one thread). The method is in "
       "[README.md](README.md). Times in seconds. TTS\\* is the optimal time to solution: the shortest expected time to "
       "reach the exact optimum with 99% probability, over run lengths of 10 to 10,000 sweeps.",
       "",
       f"**Checks: {ok} of {len(r['checks'])} hold.**",
       "",
       "## Haplotype phasing (Ising, from simulated reads: coverage 10, 2–6 sites per read, 5% errors)",
       "",
       "| sites | Ising couplings | exact (CP-SAT) | switch errors of the optimum | SA TTS\\* | PT TTS\\* |",
       "| --- | --- | --- | --- | --- | --- |"]
for h in r["haplotype"]:
    out.append(f"| {h['size']} | {h['terms']} | {s(h['exact']['seconds'])} ({h['exact']['status']}) | "
               f"{h['exact']['switch_errors']} | {tts(h['sa'])} | {tts(h['pt'])} |")
out += ["",
        "## φX174 assembly: the layout of reads as a Hamiltonian path of largest overlap",
        "",
        "Reads of 100 bp, coverage 4, from a random segment of the genome. *Native*: CP-SAT's circuit constraint on "
        "the overlap graph. *Positions*: Boev et al.'s QUBO with n² variables. *Edges*: one variable per overlap edge.",
        "",
        "| reads | overlap edges | native exact | positions: variables | CP-SAT on it | SA TTS\\* | PT TTS\\* "
        "| edges: variables | CP-SAT on it | SA TTS\\* | PT TTS\\* |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
for a in r["assembly"]:
    p, e = a["qubo_positions"], a["qubo_edges"]
    out.append(f"| {a['size']} | {a['overlap_edges']} | {s(a['exact_native']['seconds'])} | {p['variables']} | "
               f"{qubo_cell(p)} | {tts(p['sa'])} | {tts(p['pt'])} | {e['variables']} | {qubo_cell(e)} | {tts(e['sa'])} | "
               f"{tts(e['pt'])} |")
out += ["",
        "## Dense modules in the STRING DNA-repair network (45 genes, links with score ≥ 0.9)",
        "",
        "| k | links in the densest module | the module | exact (CP-SAT) | CP-SAT on the QUBO | SA TTS\\* | PT TTS\\* |",
        "| --- | --- | --- | --- | --- | --- | --- |"]
for m in r["modules"]:
    ex = m["exact"]
    if m.get("open"):
        proven = ex["status"] == "optimal"
        out.append(f"| {m['size']} | {ex['edges']} of {m['size'] * (m['size'] - 1) // 2}"
                   f"{'' if proven else ' (best found, not proven)'} | {', '.join(ex['module'])} | "
                   f"{s(ex['seconds'])} ({ex['status']}; beyond the 60 s limit of the other rows, so not timed "
                   "against the heuristics) | — | — | — |")
        continue
    out.append(f"| {m['size']} | {ex['edges']} of {m['size'] * (m['size'] - 1) // 2} | {', '.join(ex['module'])} | "
               f"{s(ex['seconds'])} ({ex['status']}) | {qubo_cell(m)} | {tts(m['sa'])} | {tts(m['pt'])} |")
out += ["", "## Every check", "", "| check | holds | detail |", "| --- | --- | --- |"]
out += [f"| {c['what']} | {'yes' if c['ok'] else '**NO**'} | {c['detail']} |" for c in r["checks"]]
print("\n".join(out))
