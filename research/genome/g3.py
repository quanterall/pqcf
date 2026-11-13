"""G3: the GBS graph kernel against classical features, for telling which organism a window of DNA comes from
(SARS-CoV-2, E. coli K-12, human chromosome 21). See gbs_kernel.py and README.md.
    python g3.py specs <out>       # PQCF's check: two window graphs as GBS devices -> <out>/chunk-0/specs.json
    python g3.py analyse <out>     # the check against PQCF's samples, then the classification -> <out>/g3.json
"""
import json
import math
import os
import sys
import time

import numpy as np
from scipy.stats import chisquare
from sklearn.model_selection import GridSearchCV, RepeatedStratifiedKFold, StratifiedKFold, cross_val_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

import gbs_graph as G
import gbs_kernel as K
import instances as I

HERE = os.path.dirname(os.path.abspath(__file__))
GENOMES = {"SARS-CoV-2": "NC_045512.2", "E. coli K-12": "NC_000913.3", "human chr21": "NC_000021.9"}
LENGTHS = (100, 250, 500)
PER_CLASS = 150
CHECK_N = 4
CHECK_TARGET = 4000


def genomes():
    return {name: I.read_fasta(os.path.join(HERE, "data", f"{acc}.fasta")) for name, acc in GENOMES.items()}


def check_graphs():
    g = genomes()
    return {"check-SARS-CoV-2": K.graph(K.windows(g["SARS-CoV-2"], 250, 1, 7)[0]),
            "check-human": K.graph(K.windows(g["human chr21"], 250, 1, 7)[0])}


def specs(out):
    jobs = []
    for name, A in check_graphs().items():
        pt = G.sampling_point(A, CHECK_N)
        e = G.encode(A, pt["mean_photons"])
        jobs.append({"name": name, "modes": 16, "r": e["r"],
                     "U": [[[float(z.real), float(z.imag)] for z in row] for row in e["U"]], "k": CHECK_N,
                     "target": CHECK_TARGET, "max_shots": int(30 * CHECK_TARGET / pt["p_k"]), "seed": 31,
                     "max_photons": pt["cap"], "collisions": True, "point": pt})
    d = os.path.join(out, "chunk-0")
    os.makedirs(d, exist_ok=True)
    json.dump(jobs, open(os.path.join(d, "specs.json"), "w"))
    print(len(jobs), "check specs")


def check_against_pqcf(out):
    gbs = {g["name"]: g for g in json.load(open(os.path.join(out, "chunk-0", "gbs.json")))}
    checks = []
    for name, A in check_graphs().items():
        pats, probs = K.exact_pattern_probabilities(A, CHECK_N)
        index = {tuple(np.bincount(p, minlength=16)): i for i, p in enumerate(pats)}
        obs = np.zeros(len(pats))
        for s in gbs[name]["samples"]:
            obs[index[tuple(s)]] += 1
        n = obs.sum()
        exp = probs * n
        big = exp >= 5
        o, e = np.append(obs[big], obs[~big].sum()), np.append(exp[big], exp[~big].sum())
        pv = chisquare(o, e).pvalue
        checks.append({"what": f"{name}: PQCF's simulated GBS device gives the exact pattern probabilities, "
                               f"{CHECK_N} photons", "ok": bool(pv > 1e-3),
                       "detail": f"chi-square p = {pv:.3f}, {len(o) - 1} degrees of freedom, {int(n)} samples"})
    return checks


def classify(X, y, seed=0):
    """Accuracy of an RBF support-vector machine, its C and gamma chosen by an inner 3-fold search, over 5-fold
    cross-validation repeated 3 times."""
    keep = X.std(axis=0) > 1e-12
    X = X[:, keep]
    inner = GridSearchCV(make_pipeline(StandardScaler(), SVC(kernel="rbf")),
                         {"svc__C": [1, 10, 100], "svc__gamma": ["scale", 0.01, 0.1]},
                         cv=StratifiedKFold(3, shuffle=True, random_state=seed))
    scores = cross_val_score(inner, X, y, cv=RepeatedStratifiedKFold(n_splits=5, n_repeats=3, random_state=seed))
    return {"mean": float(scores.mean()), "sd": float(scores.std()), "features": int(keep.sum())}


def analyse(out):
    checks = check_against_pqcf(out)
    g = genomes()
    of = K.OrbitFeatures(16, (4, 6, 8))
    results = []
    for L in LENGTHS:
        # balanced classes: as many windows per class as the smallest genome (SARS-CoV-2, 29,903 bases) holds
        per = min(PER_CLASS, min(len(seq) // L - 1 for seq in g.values()))
        wins, labels = [], []
        for c, (name, seq) in enumerate(g.items()):
            ws = K.windows(seq, L, per, 1000 + L)
            checks.append({"what": f"{name}, {L} bp: {per} windows without N", "ok": len(ws) == per,
                           "detail": f"{len(ws)}"})
            wins += ws
            labels += [c] * len(ws)
        y = np.array(labels)
        graphs = [K.graph(w) for w in wins]
        feats, secs = {}, {}
        for fname, fn in (("GBS orbits (N = 4, 6, 8)", lambda i: of(graphs[i])),
                          ("closed walks Tr(A^k), k = 2..8", lambda i: K.walk_features(graphs[i])),
                          ("spectrum of A", lambda i: K.spectrum_features(graphs[i])),
                          ("three-letter words (no graph)", lambda i: K.kmer_features(wins[i]))):
            t0 = time.perf_counter()
            feats[fname] = np.array([fn(i) for i in range(len(wins))])
            secs[fname] = (time.perf_counter() - t0) / len(wins)
        row = {"length": L, "windows": len(wins), "per_class": per, "classes": list(g), "methods": {}}
        for fname, X in feats.items():
            r = classify(X, y)
            r["seconds_per_window"] = secs[fname]
            row["methods"][fname] = r
            print(L, fname, f"{r['mean']:.3f} ± {r['sd']:.3f}", flush=True)
        results.append(row)
    rec = {"results": results, "checks": checks, "feature_names": of.names()}
    json.dump(rec, open(os.path.join(out, "g3.json"), "w"), indent=1)
    bad = [c for c in checks if not c["ok"]]
    print(f"{len(checks) - len(bad)} of {len(checks)} checks hold")
    for c in bad:
        print("FAILS:", c["what"], c["detail"])
    return not bad


if __name__ == "__main__":
    if sys.argv[1] == "specs":
        specs(sys.argv[2])
    else:
        sys.exit(0 if analyse(sys.argv[2]) else 1)
