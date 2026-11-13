"""The variant benchmark's results as Markdown.
    python3 report.py results/<host>-<date>.json > RESULTS.md"""
import json
import os
import sys

path = sys.argv[1]
r = json.load(open(path))
t = r["timing"]
ok = sum(c["ok"] for c in r["checks"])
s8, s1, un = t["sliced8"], t["sliced1"], t["unsliced"]
out = [f"# Finding the variants of GIAB HG002, {r['chrom']}: latest run", "",
       f"From [`results/{os.path.basename(path)}`](results/{os.path.basename(path)}), written by `run.sh` on "
       f"{r['host']}, {r['date']}. The method is in [README.md](README.md).", "",
       f"**Checks: {ok} of {len(r['checks'])} hold.**", "",
       "## Time", "",
       "| step | tasks | workers | seconds |", "| --- | --- | --- | --- |",
       f"| alignment of the reads (minimap2) | {s8['chunks']} chunks | 1 at a time, many threads | {s8['align_seconds']:.1f} |",
       f"| variant calling, sliced | {s8['slices']} slices | 8 | {s8['call_seconds']:.1f} |",
       f"| variant calling, sliced | {s1['slices']} slices | 1 | {s1['call_seconds']:.1f} |",
       f"| variant calling, unsliced | 1 | 1 | {un['call_seconds']:.1f} |",
       f"| gathering the slices | — | — | {s8['gather_seconds']:.2f} |", "",
       f"Calling speed-up with 8 workers: **{s1['call_seconds'] / s8['call_seconds']:.1f}×** over 1 worker, "
       f"**{un['call_seconds'] / s8['call_seconds']:.1f}×** over the unsliced call. Slowest slice: "
       f"{max(s8['call_task_seconds']):.1f} s, median {sorted(s8['call_task_seconds'])[len(s8['call_task_seconds']) // 2]:.1f} s.",
       "", "## Accuracy against the GIAB v4.2.1 truth set (confident regions)", "",
       "| variants | found (true positives) | false positives | missed | precision | recall | F1 |",
       "| --- | --- | --- | --- | --- | --- | --- |"]
for k, name in (("snps", "SNPs"), ("indels", "indels")):
    a = r["accuracy"][k]
    out.append(f"| {name} | {a['tp_baseline']:,.0f} | {a['fp']:,.0f} | {a['fn']:,.0f} | {100 * a['precision']:.2f}% | "
               f"{100 * a['recall']:.2f}% | {100 * a['f1']:.2f}% |")
out += ["", "## Every check", "", "| check | holds | detail |", "| --- | --- | --- |"]
out += [f"| {c['what']} | {'yes' if c['ok'] else '**NO**'} | {c['detail']} |" for c in r["checks"]]
print("\n".join(out))
