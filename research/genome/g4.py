"""G4: how far a genome's quantum encoding (mps.py) compresses: the bond dimension a matrix product state needs, for
phiX174 and SARS-CoV-2 in full and for E. coli K-12 and human chromosome 21 at their middle cut, against a random
sequence of the same length (no structure) and a periodic one (all structure). Reproduces Creevey et al.'s phiX174
result (bond dimension 98 for an error below 1e-6) first.
    python g4.py <out>     # -> <out>/g4.json
"""
import json
import math
import os
import random
import sys
import time

import numpy as np

import instances as I
import mps as M

HERE = os.path.dirname(os.path.abspath(__file__))
EPS = (1e-2, 1e-4, 1e-6)


def controls(L, seed):
    rng = random.Random(seed)
    rand = "".join(rng.choice("ACGT") for _ in range(L))
    unit = "".join(rng.choice("ACGT") for _ in range(100))
    periodic = (unit * (L // 100 + 1))[:L]
    return {"random": rand, "periodic (a 100-base unit repeated)": periodic}


def full(name, seq):
    n = M.qubits(len(seq))
    psi = M.state(seq)
    rows = []
    for eps in EPS:
        t0 = time.perf_counter()
        _, bonds, fid = M.tt_svd(psi, n, eps)
        rows.append({"eps": eps, "max_bond": max(bonds), "bonds": bonds, "infidelity": max(0.0, 1 - fid),
                     "parameters": M.parameters(bonds), "seconds": time.perf_counter() - t0})
    print(name, len(seq), n, [(r["eps"], r["max_bond"], r["parameters"]) for r in rows], flush=True)
    return {"sequence": name, "bases": len(seq), "qubits": n, "amplitudes": 1 << n, "full": rows}


def middle(name, seq):
    n = M.qubits(len(seq))
    p = n - 2
    cut = n // 2                          # position qubits on the left of the cut: the most balanced cut
    t0 = time.perf_counter()
    w = M.cut_spectrum(seq, cut)
    secs = time.perf_counter() - t0
    w = w / w.sum()
    nz = w[w > 1e-15]
    left, right = cut, n - cut
    out = {"sequence": name, "bases": len(seq), "qubits": n, "cut_after_qubits": cut,
           "max_possible": 1 << min(left, right), "seconds": secs,
           "entropy_bits": float(-(nz * np.log2(nz)).sum()), "max_entropy_bits": min(left, right),
           "rank": int((w > 1e-12).sum()),
           "needed": {str(e): M.bond_for_error(w, e) for e in EPS}}
    print(name, len(seq), n, out["needed"], "of", out["max_possible"], f"{secs:.1f}s", flush=True)
    return out


def main(out):
    rec = {"checks": [], "small": [], "large": []}
    g = {"phiX174": M.clean(I.read_fasta(os.path.join(HERE, "data", "phix174.fasta"))),
         "SARS-CoV-2": M.clean(I.read_fasta(os.path.join(HERE, "data", "NC_045512.2.fasta")))}
    for name, seq in g.items():
        rec["small"].append(full(name, seq))
        for cname, cseq in controls(len(seq), 1).items():
            rec["small"].append(full(f"{cname}, {len(seq)} bases", cseq))
    phix = rec["small"][0]["full"][-1]
    rec["checks"].append({"what": "phiX174: bond dimension 98 for an error below 1e-6, as Creevey et al. report",
                          "ok": phix["max_bond"] == 98 and phix["infidelity"] < 1e-6,
                          "detail": f"max bond {phix['max_bond']}, infidelity {phix['infidelity']:.1e}"})
    # the middle-cut spectrum against the full SVD, on SARS-CoV-2
    s = g["SARS-CoV-2"]
    n = M.qubits(len(s))
    w = M.cut_spectrum(s, n // 2, dtype=np.float64)
    S = np.linalg.svd(M.state(s).reshape(1 << (n // 2), -1), compute_uv=False)
    k = min(len(w), len(S))
    rec["checks"].append({"what": "the block-overlap spectrum equals the SVD of the state (SARS-CoV-2, middle cut)",
                          "ok": bool(np.allclose(w[:k], S[:k] ** 2, atol=1e-9)),
                          "detail": f"largest difference {np.abs(w[:k] - S[:k] ** 2).max():.1e}"})
    big = {"E. coli K-12": "NC_000913.3", "human chromosome 21 (N removed)": "NC_000021.9"}
    for name, acc in big.items():
        seq = M.clean(I.read_fasta(os.path.join(HERE, "data", f"{acc}.fasta")))
        rec["large"].append(middle(name, seq))
        if acc == "NC_000913.3":
            for cname, cseq in controls(len(seq), 2).items():
                rec["large"].append(middle(f"{cname}, {len(seq):,} bases", cseq))
    json.dump(rec, open(os.path.join(out, "g4.json"), "w"), indent=1)
    bad = [c for c in rec["checks"] if not c["ok"]]
    print(f"{len(rec['checks']) - len(bad)} of {len(rec['checks'])} checks hold")
    for c in bad:
        print("FAILS:", c["what"], c["detail"])
    return not bad


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1]) else 1)
