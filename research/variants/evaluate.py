"""The benchmark's checks and accuracy (run.sh): the sliced calls against the unsliced ones and across worker
counts, then RTG vcfeval against the GIAB HG002 truth set in its confident regions, for SNPs and indels apart.
    python3 evaluate.py <work> <chrom> <ref> <out.json> <host> <date>"""
import gzip
import hashlib
import json
import os
import subprocess
import sys

work, chrom, ref, out, host, date = sys.argv[1:]
D = "/v/data"


def body(path):
    """The VCF's records (header lines differ by command and date)."""
    with gzip.open(path, "rt") as f:
        return [l for l in f if not l.startswith("#")]


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True)


def vcfeval(calls, kind):
    base = f"{work}/eval-{kind}"
    t, c = f"{work}/truth-{kind}.vcf.gz", f"{work}/calls-{kind}.vcf.gz"
    flag = {"snps": "-v snps", "indels": "-v indels"}[kind]
    sh(f"bcftools view {flag} -Oz -o {t} {D}/truth.{chrom}.vcf.gz && tabix -f -p vcf {t}")
    sh(f"bcftools view {flag} -Oz -o {c} {calls} && tabix -f -p vcf {c}")
    sh(f"rm -rf {base}")
    sh(f"rtg vcfeval -b {t} -c {c} -e {D}/truth.{chrom}.bed -t {D}/{ref}.sdf -o {base} > /dev/null 2>&1")
    rows = [l.split() for l in open(f"{base}/summary.txt") if l.strip()]
    head, none = rows[0], next(r for r in rows if r[0] == "None")
    return dict(zip(["threshold", "tp_baseline", "tp_call", "fp", "fn", "precision", "recall", "f1"],
                    [none[0]] + [float(x) for x in none[1:]]))


timing = {n: json.load(open(f"{work}/{n}.vcf.gz.timing.json")) for n in ("sliced8", "sliced1", "unsliced")}
b8, b1, bu = (body(f"{work}/{n}.vcf.gz") for n in ("sliced8", "sliced1", "unsliced"))
checks = [
    {"what": "8 workers and 1 worker give byte-identical calls", "ok": b8 == b1,
     "detail": f"{len(b8):,} and {len(b1):,} records, sha256 {hashlib.sha256(''.join(b8).encode()).hexdigest()[:16]}"},
]
# Sliced against unsliced. bcftools' error model keeps a random 255 of the bases at a deeper site, from one random
# stream for the whole run, so at such sites the result depends on where the run started: there the two can differ.
# Everywhere else they must be identical, and no difference may touch the truth set's confident regions.
def depth(l):
    return int(next(x for x in l.split("\t")[7].split(";") if x.startswith("DP=")).split("=")[1])


def confident(pos, bed=[(int(s), int(e)) for c, s, e, *_ in (l.split("\t") for l in open(f"{D}/truth.{chrom}.bed"))]):
    import bisect
    i = bisect.bisect_right([s for s, e in bed], pos - 1) - 1
    return i >= 0 and bed[i][0] < pos <= bed[i][1]


diff = set(b8) ^ set(bu)
shallow8, shallowu = [l for l in b8 if depth(l) <= 255], [l for l in bu if depth(l) <= 255]
checks += [
    {"what": "sliced and unsliced calls are identical at every site with depth ≤ 255", "ok": shallow8 == shallowu,
     "detail": f"{len(shallow8):,} records"},
    {"what": "every sliced/unsliced difference is at a site deeper than 255 reads (bcftools' random subsampling)",
     "ok": all(depth(l) > 255 for l in diff),
     "detail": f"{len(diff):,} differing records of {len(b8) + len(bu):,}; smallest depth among them "
               f"{min((depth(l) for l in diff), default=0)}"},
    {"what": "no sliced/unsliced difference lies in the truth set's confident regions",
     "ok": not any(confident(int(l.split("\t")[1])) for l in diff),
     "detail": f"{sum(confident(int(l.split(chr(9))[1])) for l in diff)} in confident regions"},
]
acc = {k: vcfeval(f"{work}/sliced8.vcf.gz", k) for k in ("snps", "indels")}
rec = {"host": host, "date": date, "chrom": chrom, "records": len(b8), "timing": timing, "accuracy": acc,
       "checks": checks}
json.dump(rec, open(out, "w"), indent=1)
for c in checks:
    print(("ok   " if c["ok"] else "FAIL ") + c["what"] + ": " + c["detail"])
print({k: (v["precision"], v["recall"], v["f1"]) for k, v in acc.items()})
sys.exit(0 if all(c["ok"] for c in checks) else 1)
