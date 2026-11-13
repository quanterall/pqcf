"""A worked example on public genomes: what the genome research (G0-G4) does with a virus, a bacterium and a human
chromosome. See README.md in this directory.
    python example.py prepare <out>    # assembly, phasing, encoding; the GBS device specs for PQCF
    python example.py analyse <out>    # after PQCF has sampled: the dense-module search -> <out>/example.json
"""
import json
import math
import os
import random
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

import cim as C          # noqa: E402
import g1                # noqa: E402
import g2                # noqa: E402
import gbs_graph as G    # noqa: E402
import instances as I    # noqa: E402
import solvers as S      # noqa: E402

DATA = os.path.join(ROOT, "data")
GBS_KS = (4, 6, 8)
GBS_TARGET = 2000


def fasta(acc):
    return I.read_fasta(os.path.join(DATA, f"{acc}.fasta"))


# --- 1. assembly from reads -----------------------------------------------------------------------------------------

def tiled_reads(seq, read_len, coverage, min_overlap, seed):
    """Error-free reads covering seq end to end: evenly spaced starts with a random jitter, shuffled."""
    rng = random.Random(seed)
    n = math.ceil(len(seq) * coverage / read_len)
    step = (len(seq) - read_len) / (n - 1)
    assert read_len - 5 * step / 3 >= min_overlap
    starts = sorted({0, len(seq) - read_len} |
                    {int(round(i * step + rng.uniform(-step / 3, step / 3))) for i in range(1, n - 1)})
    order = list(range(len(starts)))
    rng.shuffle(order)
    reads = [seq[starts[i]:starts[i] + read_len] for i in order]
    truth = sorted(range(len(reads)), key=lambda k: starts[order[k]])
    return reads, truth


def overlap_graph(reads, min_overlap=40, seed_len=30):
    """Exact suffix-prefix overlaps, through an index of read prefixes: u -> v of weight k when the last k bases of
    u are the first k of v (the longest such k)."""
    index = {}
    for v, r in enumerate(reads):
        index.setdefault(r[:seed_len], []).append(v)
    edges = {}
    for u, r in enumerate(reads):
        L = len(r)
        for j in range(1, L - min_overlap + 1):
            for v in index.get(r[j:j + seed_len], ()):
                if v != u and (u, v) not in edges and r[j:] == reads[v][:L - j]:
                    edges[(u, v)] = L - j
    return edges


def assemble(name, seq, read_len=300, coverage=6, seed=1):
    t0 = time.perf_counter()
    reads, truth = tiled_reads(seq, read_len, coverage, 40, seed)
    edges = overlap_graph(reads)
    t_graph = time.perf_counter() - t0
    asm = I.Assembly(reads, truth, edges, seq)
    ex = S.exact_layout(asm, 600.0)
    path = ex["path"]
    contig = reads[path[0]] + "".join(reads[v][edges[(u, v)]:] for u, v in zip(path, path[1:]))
    out = {"genome": name, "bases": len(seq), "reads": len(reads), "read_length": read_len, "coverage": coverage,
           "overlap_edges": len(edges), "overlap_seconds": t_graph, "layout_status": ex["status"],
           "layout_seconds": ex["seconds"], "layout_is_true_order": path == truth,
           "assembled_bases": len(contig), "assembly_equals_reference": contig == seq,
           "qubo_positions_variables": len(reads) ** 2, "qubo_edges_variables": len(edges)}
    print("assembly", name, {k: out[k] for k in ("reads", "overlap_edges", "layout_status", "layout_seconds",
                                                 "assembly_equals_reference")}, flush=True)
    return out


# --- 2. phasing a simulated individual on a human chromosome region --------------------------------------------------

def simulate_individual(region, snp_rate, seed):
    """Heterozygous sites at random positions of the real sequence, each with the reference base on one haplotype
    and another base on the other, in a random phase: h_i = +1 if haplotype 1 carries the reference."""
    rng = random.Random(seed)
    pos = sorted(p for p in range(len(region)) if rng.random() < snp_rate)
    phase = np.array([rng.choice((1, -1)) for _ in pos])
    return np.array(pos), phase


def reads_over_sites(pos, phase, length, n_reads, region_len, error, seed, insert=None):
    """Reads (or read pairs, with `insert`) from a random haplotype: the alleles they show at the sites they cover,
    each wrong with probability `error`."""
    rng = random.Random(seed)
    frags = []
    for _ in range(n_reads):
        s = rng.randint(0, region_len - (insert or length))
        spans = [(s, s + length)] if insert is None else [(s, s + length), (s + insert - length, s + insert)]
        side = rng.choice((1, -1))
        f = {}
        for a, b in spans:
            lo, hi = np.searchsorted(pos, a), np.searchsorted(pos, b)
            for i in range(lo, hi):
                f[int(i)] = (-1 if rng.random() < error else 1) * side * int(phase[i])
        if len(f) >= 2:
            frags.append(f)
    return frags


def blocks(n_sites, frags):
    """Phase blocks: sites joined by reads covering more than one of them (union-find)."""
    parent = list(range(n_sites))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for f in frags:
        s = sorted(f)
        for a, b in zip(s, s[1:]):
            parent[find(a)] = find(b)
    sizes = {}
    for i in range(n_sites):
        r = find(i)
        sizes[r] = sizes.get(r, 0) + 1
    big = [v for v in sizes.values() if v > 1]
    return {"blocks": len(big), "largest": max(sizes.values()), "sites_phased": sum(big)}


def phasing(seed=3):
    chrom = fasta("NC_000021.9")
    start, length = 30_000_000, 1_000_000
    region = chrom[start:start + length]
    assert "N" not in region
    pos, phase = simulate_individual(region, 1 / 1500, seed)
    n = len(pos)
    out = {"region": f"chr21:{start + 1:,}-{start + length:,} (GRCh38)", "heterozygous_sites": n}
    short = reads_over_sites(pos, phase, 150, int(30 * length / 300), length, 0.001, seed + 1, insert=500)
    long_ = reads_over_sites(pos, phase, 15000, int(10 * length / 15000), length, 0.01, seed + 2)
    out["short_reads"] = {"kind": "paired 2 x 150 bp, insert 500 bp, coverage 30, 0.1% errors",
                          "fragments_linking_sites": len(short), **blocks(n, short)}
    out["long_reads"] = {"kind": "15 kb, coverage 10, 1% errors", "fragments_linking_sites": len(long_),
                         **blocks(n, long_)}
    # the long-read Ising problem: exact, then G2's tuned heuristics
    hap = I.Haplotype(phase, long_)
    J = {}
    for f in long_:
        s = sorted(f)
        for a in range(len(s)):
            for b in range(a + 1, len(s)):
                i, j = s[a], s[b]
                J[(i, j)] = J.get((i, j), 0) + f[i] * f[j]
    hap.J = {k: v for k, v in J.items() if v}
    ex = S.exact_haplotype(hap, 600.0, workers=8)
    out["long_reads"]["exact"] = {"status": ex["status"], "seconds": ex["seconds"],
                                  "switch_errors": hap.switch_errors(ex["spins"]), "mec": hap.mec(ex["spins"]),
                                  "band": max(j - i for (i, j) in hap.J)}
    inst = {"kind": "hapL", "J": hap.J, "h": {}, "offset": 0.0, "nvar": n}
    inst["Q"], inst["c"] = I.ising_to_qubo(hap.J)
    tuned = {"pt": {"replicas": 8, "b1": 4.0}, "ahc": {"p": 0.0, "eps": 0.3, "beta": 0.4}}   # G2's long-read choice
    for s, params in tuned.items():
        r = g2.tts(g2.make(s, inst, params), float(ex["energy"]), g2.LENGTHS[s], 20)
        out["long_reads"][s] = {"params": params, "tts_star": r["tts_star"], "at_length": r["at_length"]}
    print("phasing", {k: out[k] for k in ("heterozygous_sites",)}, out["short_reads"], out["long_reads"], flush=True)
    return out


# --- 3. dense modules in E. coli's DNA replication and repair network, with PQCF's GBS --------------------------------

def ecoli_network():
    net = I.string_network(os.path.join(DATA, "string_ecoli_replication_repair.tsv"), 0.9)
    n = len(net.genes)
    A = np.zeros((n, n), int)
    for (i, j) in net.edges:
        A[i, j] = A[j, i] = 1
    return net, A


def gbs_specs(out, chunks):
    net, A = ecoli_network()
    jobs = []
    for k in GBS_KS:
        # at most 1e-6 of the probability above PQCF's photon cap (its truncation is step by step: G1)
        pt = G.sampling_point(A, k, tail=1e-6)
        e = G.encode(A, pt["mean_photons"])
        jobs.append({"name": f"ecoli-k{k}", "modes": len(A), "r": e["r"], "U": g1.cx(e["U"]), "k": k,
                     "target": GBS_TARGET, "max_shots": int(40 * GBS_TARGET / pt["p_k"]), "seed": 77 * k,
                     "max_photons": pt["cap"], "point": pt})
    for c in range(chunks):
        d = os.path.join(out, f"chunk-{c}")
        os.makedirs(d, exist_ok=True)
        json.dump([{**j, "target": -(-j["target"] // chunks), "max_shots": -(-j["max_shots"] // chunks),
                    "seed": j["seed"] + 100 * c} for j in jobs], open(os.path.join(d, "specs.json"), "w"))
    json.dump({j["name"]: j["point"] for j in jobs}, open(os.path.join(out, "gbs_plan.json"), "w"))


def gbs_analyse(out):
    gbs = g1.merge(out)
    net, A = ecoli_network()
    import networkx as nx
    omega = max(len(c) for c in nx.find_cliques(nx.from_numpy_array(A)))
    rng = np.random.default_rng(5)
    rows = []
    for k in GBS_KS:
        ex = S.exact_dks(net, k, 120.0)
        g = gbs[f"ecoli-k{k}"]
        oh, _ = G.oh_sampler(A, k, GBS_TARGET, rng)
        un = G.uniform_sampler(len(A), k, GBS_TARGET, rng)
        row = {"k": k, "densest": {"edges": ex["edges"], "of": k * (k - 1) // 2, "status": ex["status"],
                                   "module": [net.genes[i] for i in ex["chosen"]]},
               "gbs_samples": len(g["samples"]), "pqcf_cpu_seconds_per_sample": g["seconds"] / max(1, len(g["samples"])),
               "max_truncated": max(g["truncated"]), "max_clique": omega}
        for lab, samp in (("gbs", g["samples"][:GBS_TARGET]), ("oh", oh), ("uniform", un)):
            row[lab] = g1.search_stats(A, samp, ex["edges"], omega, rng)
        rows.append(row)
        print("gbs", k, ex["edges"], row["densest"]["module"],
              {lab: row[lab]["dks_p_optimum"][-1] for lab in ("gbs", "oh", "uniform")}, flush=True)
    return {"genes": len(net.genes), "links": len(net.edges), "rows": rows}


# --- 4. encoding: from G4's latest run --------------------------------------------------------------------------------

def encoding():
    res = sorted(f for f in os.listdir(os.path.join(ROOT, "results")) if f.startswith("g4-"))
    r = json.load(open(os.path.join(ROOT, "results", res[-1])))
    small = {x["sequence"]: x for x in r["small"]}
    large = {x["sequence"]: x for x in r["large"]}
    return {"source": res[-1], "SARS-CoV-2": small["SARS-CoV-2"], "SARS-CoV-2 random": small["random, 29903 bases"],
            "E. coli": large["E. coli K-12"], "human chr21": large["human chromosome 21 (N removed)"]}


def main():
    mode, out = sys.argv[1], sys.argv[2]
    if mode == "prepare":
        rec = {"assembly": [assemble("SARS-CoV-2 (NC_045512.2), whole genome", fasta("NC_045512.2")),
                            assemble("E. coli K-12 (NC_000913.3), 340,001-390,000 (the lac operon region)",
                                     fasta("NC_000913.3")[340_000:390_000])],
               "phasing": phasing(), "encoding": encoding()}
        json.dump(rec, open(os.path.join(out, "prepared.json"), "w"), indent=1, default=float)
        gbs_specs(out, int(sys.argv[3]))
    else:
        rec = json.load(open(os.path.join(out, "prepared.json")))
        rec["modules"] = gbs_analyse(out)
        checks = []
        for a in rec["assembly"]:
            checks.append({"what": f"{a['genome']}: the assembly equals the reference", "ok": a["assembly_equals_reference"]
                           and a["layout_status"] == "optimal", "detail": f"{a['assembled_bases']:,} bases"})
        lr = rec["phasing"]["long_reads"]
        checks.append({"what": "chr21 region, long reads: the exact phasing is proven",
                       "ok": lr["exact"]["status"] == "optimal", "detail": f"{lr['exact']['seconds']:.1f} s"})
        for row in rec["modules"]["rows"]:
            checks.append({"what": f"E. coli network, k = {row['k']}: densest module proven, PQCF's truncation negligible",
                           "ok": row["densest"]["status"] == "optimal" and row["max_truncated"] < 1e-4,
                           "detail": f"{row['densest']['status']}, truncated {row['max_truncated']:.1e}"})
        rec["checks"] = checks
        json.dump(rec, open(os.path.join(out, "example.json"), "w"), indent=1, default=float)
        bad = [c for c in checks if not c["ok"]]
        print(f"{len(checks) - len(bad)} of {len(checks)} checks hold")
        sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
