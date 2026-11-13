"""The worked example's results as Markdown.
    python report_example.py example.json > EXAMPLE.md"""
import json
import math
import sys

r = json.load(open(sys.argv[1]))
ok = sum(c["ok"] for c in r["checks"])


def secs(x):
    return "never reached" if x is None or not math.isfinite(x) else (f"{x:.3f} s" if x < 10 else f"{x:.1f} s")


out = ["# Worked example: a virus, a bacterium and a human chromosome", "",
       "Results of `run.sh`, on public genomes from NCBI. What each step means is in [README.md](README.md).", "",
       f"**Checks: {ok} of {len(r['checks'])} hold.**", "",
       "## 1. Assembling genomes from reads", "",
       "Simulated error-free reads, laid out by the exact solver on the overlap graph, then joined.", "",
       "| genome | bases | reads | overlaps | layout (exact) | assembly = reference | QUBO for an annealer: n² encoding / edge encoding |",
       "| --- | --- | --- | --- | --- | --- | --- |"]
for a in r["assembly"]:
    out.append(f"| {a['genome']} | {a['bases']:,} | {a['reads']:,} × {a['read_length']} bp | {a['overlap_edges']:,} | "
               f"{a['layout_seconds']:.2f} s, {a['layout_status']} | {'yes' if a['assembly_equals_reference'] else '**no**'} | "
               f"{a['qubo_positions_variables']:,} / {a['qubo_edges_variables']:,} variables |")
p = r["phasing"]
sr, lr = p["short_reads"], p["long_reads"]
out += ["", "## 2. Phasing a simulated individual on human chromosome 21", "",
        f"{p['region']}: {p['heterozygous_sites']} heterozygous sites placed on the real sequence (about one per "
        "1,500 bases), in a random phase; reads from either haplotype.", "",
        "| reads | fragments joining sites | phase blocks | largest block (sites) | sites phased |",
        "| --- | --- | --- | --- | --- |",
        f"| {sr['kind']} | {sr['fragments_linking_sites']:,} | {sr['blocks']} | {sr['largest']} | {sr['sites_phased']} |",
        f"| {lr['kind']} | {lr['fragments_linking_sites']:,} | {lr['blocks']} | {lr['largest']} | {lr['sites_phased']} |",
        "", "The long-read phasing as an Ising problem:", "",
        "| method | time to the optimum | phasing |", "| --- | --- | --- |",
        f"| exact (CP-SAT, 8 workers) | {lr['exact']['seconds']:.2f} s ({lr['exact']['status']}) | "
        f"{lr['exact']['switch_errors']} switch errors, {lr['exact']['mec']} allele corrections |",
        f"| parallel tempering, G2's tuned settings | {secs(lr['pt']['tts_star'])} | — |",
        f"| CIM-derived AHC, G2's tuned settings | {secs(lr['ahc']['tts_star'])} | — |", ""]
m = r["modules"]
out += ["## 3. Dense modules in E. coli's DNA replication and repair network", "",
        f"STRING 12.5, {m['genes']} genes, {m['links']} links with score ≥ 0.9; GBS simulated by PQCF.", "",
        "| k | densest module | links | reached within 200 samples: GBS / Oh et al. / uniform | PQCF CPU-seconds per GBS sample |",
        "| --- | --- | --- | --- | --- |"]
for row in m["rows"]:
    pr = [f"{100 * row[l]['dks_p_optimum'][-1]:.0f}%" for l in ("gbs", "oh", "uniform")]
    out.append(f"| {row['k']} | {', '.join(row['densest']['module'])} | {row['densest']['edges']} of "
               f"{row['densest']['of']} | {' / '.join(pr)} | {row['pqcf_cpu_seconds_per_sample']:.3f} |")
e = r["encoding"]
out += ["", "## 4. The genomes as quantum states", "", f"From G4 (`{e['source']}`).", "",
        "| genome | bases | qubits | Schmidt values needed (error ≤ 1e-6) | a random sequence of the same length |",
        "| --- | --- | --- | --- | --- |",
        f"| SARS-CoV-2 | {e['SARS-CoV-2']['bases']:,} | {e['SARS-CoV-2']['qubits']} | largest bond "
        f"{e['SARS-CoV-2']['full'][-1]['max_bond']} | {e['SARS-CoV-2 random']['full'][-1]['max_bond']} |",
        f"| E. coli K-12 (middle cut) | {e['E. coli']['bases']:,} | {e['E. coli']['qubits']} | "
        f"{e['E. coli']['needed']['1e-06']:,} | (G4: the same) |",
        f"| human chromosome 21 (middle cut) | {e['human chr21']['bases']:,} | {e['human chr21']['qubits']} | "
        f"{e['human chr21']['needed']['1e-06']:,} | (every block of the chromosome) |",
        "", "## Every check", "", "| check | holds | detail |", "| --- | --- | --- |"]
out += [f"| {c['what']} | {'yes' if c['ok'] else '**NO**'} | {c['detail']} |" for c in r["checks"]]
print("\n".join(out))
