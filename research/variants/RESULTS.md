# Finding the variants of GIAB HG002, chr20: latest run

From [`results/dgx-spark-2026-11-13T1611Z.json`](results/dgx-spark-2026-11-13T1611Z.json), written by `run.sh` on dgx-spark, 2026-11-13T1611Z. The method is in [README.md](README.md).

**Checks: 4 of 4 hold.**

## Time

| step | tasks | workers | seconds |
| --- | --- | --- | --- |
| alignment of the reads (minimap2) | 8 chunks | 1 at a time, many threads | 301.8 |
| variant calling, sliced | 65 slices | 8 | 29.2 |
| variant calling, sliced | 65 slices | 1 | 218.2 |
| variant calling, unsliced | 1 | 1 | 204.8 |
| gathering the slices | — | — | 0.46 |

Calling speed-up with 8 workers: **7.5×** over 1 worker, **7.0×** over the unsliced call. Slowest slice: 7.7 s, median 3.2 s.

## Accuracy against the GIAB v4.2.1 truth set (confident regions)

| variants | found (true positives) | false positives | missed | precision | recall | F1 |
| --- | --- | --- | --- | --- | --- | --- |
| SNPs | 70,851 | 614 | 511 | 99.14% | 99.28% | 99.21% |
| indels | 10,847 | 308 | 409 | 97.24% | 96.37% | 96.80% |

## Every check

| check | holds | detail |
| --- | --- | --- |
| 8 workers and 1 worker give byte-identical calls | yes | 126,066 and 126,066 records, sha256 5fe94058edef1b07 |
| sliced and unsliced calls are identical at every site with depth ≤ 255 | yes | 124,789 records |
| every sliced/unsliced difference is at a site deeper than 255 reads (bcftools' random subsampling) | yes | 1,402 differing records of 252,132; smallest depth among them 261 |
| no sliced/unsliced difference lies in the truth set's confident regions | yes | 0 in confident regions |
