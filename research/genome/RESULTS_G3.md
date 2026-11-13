# G3: the GBS graph kernel against classical features, latest run

From [`results/g3-dgx-spark-2026-11-13T1130Z.json`](results/g3-dgx-spark-2026-11-13T1130Z.json), written by `g3_run.sh`. The method is in [README.md](README.md). Task: which of SARS-CoV-2, E. coli K-12, human chr21 a window of DNA comes from, with the same number of windows per class (up to 150; fewer for long windows, as the SARS-CoV-2 genome holds only 29,903 bases). Accuracy of an RBF support-vector machine (C and gamma chosen by an inner cross-validation) over 5-fold cross-validation repeated 3 times; chance is 33%.

**Checks: 11 of 11 hold.**

| window | windows per class | GBS orbits (N = 4, 6, 8) | closed walks Tr(A^k), k = 2..8 | spectrum of A | three-letter words (no graph) |
| --- | --- | --- | --- | --- | --- |
| 100 bp | 150 | 47.4% ± 4.2 | 46.5% ± 4.4 | 51.0% ± 4.0 | 83.5% ± 4.0 |
| 250 bp | 118 | 56.1% ± 5.8 | 52.1% ± 4.9 | 60.7% ± 6.6 | 91.1% ± 2.6 |
| 500 bp | 58 | 60.4% ± 8.8 | 58.2% ± 3.8 | 69.9% ± 8.3 | 96.0% ± 2.7 |

## Features and their cost (one core, per window)

| features | number used | seconds per window |
| --- | --- | --- |
| GBS orbits (N = 4, 6, 8) | 25 | 7.22e-02 |
| closed walks Tr(A^k), k = 2..8 | 7 | 1.47e-05 |
| spectrum of A | 16 | 1.07e-05 |
| three-letter words (no graph) | 64 | 6.87e-05 |

## Every check

| check | holds | detail |
| --- | --- | --- |
| check-SARS-CoV-2: PQCF's simulated GBS device gives the exact pattern probabilities, 4 photons | yes | chi-square p = 0.482, 215 degrees of freedom, 4000 samples |
| check-human: PQCF's simulated GBS device gives the exact pattern probabilities, 4 photons | yes | chi-square p = 0.119, 177 degrees of freedom, 4000 samples |
| SARS-CoV-2, 100 bp: 150 windows without N | yes | 150 |
| E. coli K-12, 100 bp: 150 windows without N | yes | 150 |
| human chr21, 100 bp: 150 windows without N | yes | 150 |
| SARS-CoV-2, 250 bp: 118 windows without N | yes | 118 |
| E. coli K-12, 250 bp: 118 windows without N | yes | 118 |
| human chr21, 250 bp: 118 windows without N | yes | 118 |
| SARS-CoV-2, 500 bp: 58 windows without N | yes | 58 |
| E. coli K-12, 500 bp: 58 windows without N | yes | 58 |
| human chr21, 500 bp: 58 windows without N | yes | 58 |
