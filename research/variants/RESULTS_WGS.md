# Finding every variant of GIAB HG002's whole genome: latest run

From [`results/wgs-dgx-spark-2026-11-13T2237Z.json`](results/wgs-dgx-spark-2026-11-13T2237Z.json), written by `run_wgs.sh` on dgx-spark, 2026-11-13T2237Z: one machine, 10 threads for alignment, 8 calling workers. The method is in [README.md](README.md).

**Checks: 4 of 4 hold.**

## Time

| step | tasks | time |
| --- | --- | --- |
| alignment (minimap2) | 70 chunks of 6 million read pairs | 3.0 h |
| variant calling (bcftools) | 3,102 slices of 1 Mb, 8 workers | 25.5 min |
| gathering | — | 18.4 s |
| **total** | | **3.4 h** |

Alignment: alignment ran in a first attempt that stopped later at a calling slice (bcftools norm on an IUPAC code in GRCh38, chr3:16902879; fixed with norm -c w); its 70 chunk alignments were kept and the run resumed at calling (SKIP_ALIGN). Alignment time from the work directory's file times: 19:33:28 to 22:31:21 UTC. Calling per slice: median 3.9 s, slowest 71 s. On chromosome 20, calling with 1 worker took 266 s and one unsliced call 248 s.

## Accuracy against the GIAB v4.2.1 truth set (chr1–22, confident regions)

| variants | found (true positives) | false positives | missed | precision | recall | F1 |
| --- | --- | --- | --- | --- | --- | --- |
| SNPs | 3,334,821 | 29,131 | 31,553 | 99.13% | 99.06% | 99.10% |
| indels | 507,436 | 13,423 | 17,961 | 97.42% | 96.58% | 97.00% |

**Variant records:** 4,935,855 in all, by chromosome: chr1 377,473, chr2 377,060, chr3 330,119, chr4 343,549, chr5 297,248, chr6 306,104, chr7 279,034, chr8 251,941, chr9 208,336, chr10 252,137, chr11 237,372, chr12 231,342, chr13 187,808, chr14 154,624, chr15 143,903, chr16 145,136, chr17 143,168, chr18 137,880, chr19 115,607, chr20 126,428, chr21 79,970, chr22 76,693, chrX 120,396, chrY 12,527.

## Every check

| check | holds | detail |
| --- | --- | --- |
| chromosome 20: the whole run (8 workers) and 1 worker give byte-identical calls | yes | 126,428 and 126,428 records |
| chromosome 20: sliced and unsliced calls are identical at every site with depth ≤ 255 | yes | 125,011 records |
| chromosome 20: every sliced/unsliced difference is at a site deeper than 255 reads | yes | 1,659 differing records; smallest depth among them 262 |
| chromosome 20: no sliced/unsliced difference lies in the truth set's confident regions | yes | 0 in confident regions |
