"""The whole-genome benchmark's checks and accuracy (run_wgs.sh). On chromosome 20, with the whole run's
alignments: the whole run's calls (SLOTS workers) against 1 worker, byte for byte, and against one unsliced call,
identical at every site with depth <= 255 (bcftools subsamples deeper sites from a run-wide random stream; see
README.md). Then RTG vcfeval against the GIAB HG002 v4.2.1 truth set (chr1-22, its confident regions), SNPs and
indels apart.
    python3 evaluate_wgs.py <work> <ref> <out.json> <host> <date> <slots> <threads>"""
import bisect
import gzip
import json
import subprocess
import sys

work, ref, out, host, date, slots, threads = sys.argv[1:]
D = "/v/data"
TRUTH = f"{D}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
BED = f"{D}/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"


def sh(cmd, capture=False):
    r = subprocess.run(cmd, shell=True, check=True, capture_output=capture, text=True)
    return r.stdout if capture else None


def body_gz(path):
    with gzip.open(path, "rt") as f:
        return [l for l in f if not l.startswith("#")]


def depth(l):
    return int(next(x for x in l.split("\t")[7].split(";") if x.startswith("DP=")).split("=")[1])


bed20 = [(int(s), int(e)) for c, s, e, *_ in (l.split("\t") for l in open(BED)) if c == "chr20"]
starts20 = [s for s, e in bed20]


def confident20(pos):
    i = bisect.bisect_right(starts20, pos - 1) - 1
    return i >= 0 and bed20[i][0] < pos <= bed20[i][1]


def vcfeval(kind):
    base = f"{work}/eval-{kind}"
    t, c = f"{work}/truth-{kind}.vcf.gz", f"{work}/calls-{kind}.vcf.gz"
    flag = {"snps": "-v snps", "indels": "-v indels"}[kind]
    sh(f"bcftools view {flag} -Oz -o {t} {TRUTH} && tabix -f -p vcf {t}")
    sh(f"bcftools view {flag} -Oz -o {c} {work}/wgs.vcf.gz && tabix -f -p vcf {c}")
    sh(f"rm -rf {base} && rtg vcfeval --threads {threads} -b {t} -c {c} -e {BED} -t {D}/{ref}.sdf -o {base} > /dev/null 2>&1")
    rows = [l.split() for l in open(f"{base}/summary.txt") if l.strip()]
    none = next(r for r in rows if r[0] == "None")
    return dict(zip(["tp_baseline", "tp_call", "fp", "fn", "precision", "recall", "f1"], [float(x) for x in none[1:]]))


timing = {n: json.load(open(f"{work}/{n}.vcf.gz.timing.json")) for n in ("wgs", "chr20_sliced1", "chr20_unsliced")}
whole20 = [l for l in sh(f"bcftools view -H -r chr20 {work}/wgs.vcf.gz", capture=True).splitlines(keepends=True)]
one20, un20 = body_gz(f"{work}/chr20_sliced1.vcf.gz"), body_gz(f"{work}/chr20_unsliced.vcf.gz")
diff = set(one20) ^ set(un20)
checks = [
    {"what": f"chromosome 20: the whole run ({slots} workers) and 1 worker give byte-identical calls",
     "ok": whole20 == one20, "detail": f"{len(whole20):,} and {len(one20):,} records"},
    {"what": "chromosome 20: sliced and unsliced calls are identical at every site with depth ≤ 255",
     "ok": [l for l in one20 if depth(l) <= 255] == [l for l in un20 if depth(l) <= 255],
     "detail": f"{sum(depth(l) <= 255 for l in one20):,} records"},
    {"what": "chromosome 20: every sliced/unsliced difference is at a site deeper than 255 reads",
     "ok": all(depth(l) > 255 for l in diff),
     "detail": f"{len(diff):,} differing records; smallest depth among them {min((depth(l) for l in diff), default=0)}"},
    {"what": "chromosome 20: no sliced/unsliced difference lies in the truth set's confident regions",
     "ok": not any(confident20(int(l.split(chr(9))[1])) for l in diff),
     "detail": f"{sum(confident20(int(l.split(chr(9))[1])) for l in diff)} in confident regions"},
]
counts = sh(f"bcftools index -s {work}/wgs.vcf.gz", capture=True)
per_chrom = {c: int(n) for c, _, n in (l.split("\t") for l in counts.splitlines())}
acc = {k: vcfeval(k) for k in ("snps", "indels")}
rec = {"host": host, "date": date, "slots": int(slots), "threads": int(threads), "records": sum(per_chrom.values()),
       "records_per_chromosome": per_chrom, "timing": timing, "accuracy": acc, "checks": checks}
json.dump(rec, open(out, "w"), indent=1)
for c in checks:
    print(("ok   " if c["ok"] else "FAIL ") + c["what"] + ": " + c["detail"])
print({k: (v["precision"], v["recall"], v["f1"]) for k, v in acc.items()})
sys.exit(0 if all(c["ok"] for c in checks) else 1)
