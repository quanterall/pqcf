"""Tables of a G3 run as Markdown.
    python report_g3.py results/g3-<host>-<date>.json > RESULTS_G3.md"""
import json
import os
import sys

path = sys.argv[1]
r = json.load(open(path))
ok = sum(c["ok"] for c in r["checks"])
methods = list(r["results"][0]["methods"])
out = ["# G3: the GBS graph kernel against classical features, latest run", "",
       f"From [`results/{os.path.basename(path)}`](results/{os.path.basename(path)}), written by `g3_run.sh`. The "
       "method is in [README.md](README.md). Task: which of " + ", ".join(r["results"][0]["classes"]) +
       " a window of DNA comes from, with the same number of windows per class (up to 150; fewer for long windows, "
       "as the SARS-CoV-2 genome holds only 29,903 bases). Accuracy of an RBF support-vector machine (C and gamma "
       "chosen by an inner cross-validation) over 5-fold cross-validation repeated 3 times; chance is 33%.", "",
       f"**Checks: {ok} of {len(r['checks'])} hold.**", "",
       "| window | windows per class | " + " | ".join(methods) + " |",
       "| --- | --- | " + " | ".join("---" for _ in methods) + " |"]
for row in r["results"]:
    out.append(f"| {row['length']} bp | {row.get('per_class', '')} | " + " | ".join(
        f"{100 * row['methods'][m]['mean']:.1f}% ± {100 * row['methods'][m]['sd']:.1f}" for m in methods) + " |")
out += ["", "## Features and their cost (one core, per window)", "", "| features | number used | seconds per window |",
        "| --- | --- | --- |"]
row = r["results"][-1]
out += [f"| {m} | {row['methods'][m]['features']} | {row['methods'][m]['seconds_per_window']:.2e} |" for m in methods]
out += ["", "## Every check", "", "| check | holds | detail |", "| --- | --- | --- |"]
out += [f"| {c['what']} | {'yes' if c['ok'] else '**NO**'} | {c['detail']} |" for c in r["checks"]]
print("\n".join(out))
