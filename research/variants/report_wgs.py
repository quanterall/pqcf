"""The whole-genome benchmark's results as Markdown.
    python3 report_wgs.py results/wgs-<host>-<date>.json > RESULTS_WGS.md"""
import json
import os
import statistics
import sys

path = sys.argv[1]
r = json.load(open(path))
t = r["timing"]
w, one, un = t["wgs"], t["chr20_sliced1"], t["chr20_unsliced"]
ok = sum(c["ok"] for c in r["checks"])
hm = lambda s: f"{s / 3600:.1f} h" if s >= 3600 else f"{s / 60:.1f} min"
cts = w["call_task_seconds"]
out = ["# Finding every variant of GIAB HG002's whole genome: latest run", "",
       f"From [`results/{os.path.basename(path)}`](results/{os.path.basename(path)}), written by `run_wgs.sh` on "
       f"{r['host']}, {r['date']}: one machine, {r['threads']} threads for alignment, {r['slots']} calling workers. "
       "The method is in [README.md](README.md).", "", f"**Checks: {ok} of {len(r['checks'])} hold.**", "",
       "## Time", "", "| step | tasks | time |", "| --- | --- | --- |",
       f"| alignment (minimap2) | {w['chunks']} chunks of 6 million read pairs | {hm(w['align_seconds'])} |",
       f"| variant calling (bcftools) | {w['slices']:,} slices of 1 Mb, {r['slots']} workers | {hm(w['call_seconds'])} |",
       f"| gathering | — | {w['gather_seconds']:.1f} s |",
       f"| **total** | | **{hm(w['align_seconds'] + w['call_seconds'] + w['gather_seconds'])}** |", "",
       (f"Alignment per chunk: median {statistics.median(w['align_task_seconds']):.0f} s. " if w['align_task_seconds']
        else f"Alignment: {w.get('align_note', '')} ") + f"Calling per slice: median "
       f"{statistics.median(cts):.1f} s, slowest {max(cts):.0f} s. On chromosome 20, calling with 1 worker took "
       f"{one['call_seconds']:.0f} s and one unsliced call {un['call_seconds']:.0f} s.", "",
       "## Accuracy against the GIAB v4.2.1 truth set (chr1–22, confident regions)", "",
       "| variants | found (true positives) | false positives | missed | precision | recall | F1 |",
       "| --- | --- | --- | --- | --- | --- | --- |"]
for k, name in (("snps", "SNPs"), ("indels", "indels")):
    a = r["accuracy"][k]
    out.append(f"| {name} | {a['tp_baseline']:,.0f} | {a['fp']:,.0f} | {a['fn']:,.0f} | {100 * a['precision']:.2f}% | "
               f"{100 * a['recall']:.2f}% | {100 * a['f1']:.2f}% |")
out += ["", f"**Variant records:** {r['records']:,} in all, by chromosome: " +
        ", ".join(f"{c} {n:,}" for c, n in r["records_per_chromosome"].items()) + ".",
        "", "## Every check", "", "| check | holds | detail |", "| --- | --- | --- |"]
out += [f"| {c['what']} | {'yes' if c['ok'] else '**NO**'} | {c['detail']} |" for c in r["checks"]]
print("\n".join(out))
