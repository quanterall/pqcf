"""Tables of a G1 run (g1.py's JSON) as Markdown.
    python report_g1.py results/g1-<host>-<date>.json > RESULTS_G1.md"""
import json
import os
import sys

path = sys.argv[1]
r = json.load(open(path))
NAMES = {"gbs": "GBS (PQCF)", "oh": "Oh et al.", "uniform": "uniform"}
GRAPHS = {"planted": "Planted graph: 30 vertices, a dense 10-vertex subgraph of low-degree vertices (as Arrazola & "
                     "Bromley's Fig. 3)",
          "string": "STRING DNA-repair network: 45 genes, links with score ≥ 0.9"}


def pct(x):
    return f"{100 * x:.0f}%"


ok = sum(c["ok"] for c in r["checks"])
out = ["# G1: Gaussian boson sampling for dense modules and cliques, latest run", "",
       f"From [`results/{os.path.basename(path)}`](results/{os.path.basename(path)}), written by `g1_run.sh`. The "
       "method is in [README.md](README.md). Every sampler gives post-selected k-vertex samples. *After s samples*: "
       "over 10 disjoint blocks of 200 samples, the fraction of blocks whose best sample so far is optimal.", "",
       f"**Checks: {ok} of {len(r['checks'])} hold.**", ""]
for g in ("planted", "string"):
    rows = [x for x in r["results"] if x["graph"] == g]
    if not rows:
        continue
    out += [f"## {GRAPHS[g]}", "",
            "### Densest k-subgraph by random search: reached the optimum after s samples", "",
            "| k | optimum (links) | sampler | mean density of a sample | s = 1 | 10 | 50 | 200 |",
            "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for x in rows:
        cp = x["gbs"]["checkpoints"]
        idx = [cp.index(s) for s in (1, 10, 50, 200)]
        for lab in ("gbs", "oh", "uniform"):
            p = x[lab]["dks_p_optimum"]
            out.append(f"| {x['k']} | {x['densest']['edges']} of {x['densest']['of']} | {NAMES[lab]} | "
                       f"{x[lab]['mean_density_per_sample']:.2f} | " + " | ".join(pct(p[i]) for i in idx) + " |")
    out += ["", f"### Maximum clique (size {rows[0]['max_clique']}) by shrink-and-search: found after s samples", "",
            "| k | sampler | s = 1 | 10 | 50 | 200 |", "| --- | --- | --- | --- | --- | --- |"]
    for x in rows:
        cp = x["gbs"]["checkpoints"]
        idx = [cp.index(s) for s in (1, 10, 50, 200)]
        for lab in ("gbs", "oh", "uniform"):
            p = x[lab]["clique_p_maximum"]
            out.append(f"| {x['k']} | {NAMES[lab]} | " + " | ".join(pct(p[i]) for i in idx) + " |")
    if any("planted_found" in x for x in rows):
        x = next(x for x in rows if "planted_found" in x)
        out += ["", "**The planted subgraph itself, per sample:** " + ", ".join(
            f"{NAMES[l]} {100 * v:.2f}%" for l, v in x["planted_found"].items()) + "."]
    out += ["", "### What a sample costs", "",
            "| k | GBS: simulated shots per kept sample | GBS: PQCF CPU-seconds per kept sample | GBS: mean photon "
            "number, cap | Oh et al.: seconds per sample | uniform: seconds per sample |",
            "| --- | --- | --- | --- | --- | --- |"]
    for x in rows:
        c, pt = x["cost"], x["sampling_point"]
        out.append(f"| {x['k']} | {c['pqcf_gbs']['shots'] / c['pqcf_gbs']['kept']:,.0f} | "
                   f"{c['pqcf_gbs']['cpu_seconds_per_sample']:.3f} | {pt['mean_photons']:.2f}, {pt['cap']} | "
                   f"{c['oh']['seconds_per_sample']:.1e} | {c['uniform']['seconds_per_sample']:.1e} |")
    out.append("")
if all("equal_cost" in x for x in r["results"]):
    out += ["## At equal cost: samples and time to reach the densest k-subgraph", "",
            "GBS: the samples for a 50% chance, read from the curves above (so at their resolution: 1, 2, 5, 10, 20, "
            "50, 100, 200), and PQCF's CPU time to simulate them. Oh et al. and uniform: the median number of samples "
            "until the optimum first appears, over 20 runs of up to 2,000,000 samples, and its time on one core. "
            "*Break-even rate*: the post-selected samples per second a GBS device would need to match Oh et al.'s "
            "sampler on one core.", "",
            "| graph | k | GBS samples | Oh et al. samples | uniform samples | simulated GBS, CPU s | Oh et al., s | "
            "uniform, s | break-even rate |",
            "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for x in r["results"]:
        e = x["equal_cost"]
        g, o, u = e["gbs"], e["oh"], e["uniform"]
        f = lambda v, fmt: "not reached" if v is None else fmt.format(v)
        be = (f"{g['samples_for_50pct'] / o['median_seconds']:,.0f} per s"
              if g["samples_for_50pct"] and o["median_seconds"] else "—")
        gs = "> 200 (under 50% within the run's 200)" if g["samples_for_50pct"] is None else f"{g['samples_for_50pct']:,}"
        out.append(f"| {x['graph']} | {x['k']} | {gs} | "
                   f"{f(o['median_samples'], '{:,.0f}')} | {f(u['median_samples'], '{:,.0f}')} | "
                   f"{'—' if g['cpu_seconds_for_50pct'] is None else format(g['cpu_seconds_for_50pct'], '.3g')} | "
                   f"{f(o['median_seconds'], '{:.2g}')} | "
                   f"{f(u['median_seconds'], '{:.2g}')} | {be} |")
    out.append("")
mods = [x for x in r["results"] if x["graph"] == "string"]
if mods:
    out += ["## The densest modules of the DNA-repair network", "", "| k | links | genes |", "| --- | --- | --- |"]
    out += [f"| {x['k']} | {x['densest']['edges']} of {x['densest']['of']} | {', '.join(x['densest']['module'])} |"
            for x in mods]
    out.append("")
out += ["## Every check", "", "| check | holds | detail |", "| --- | --- | --- |"]
out += [f"| {c['what']} | {'yes' if c['ok'] else '**NO**'} | {c['detail']} |" for c in r["checks"]]
print("\n".join(out))
