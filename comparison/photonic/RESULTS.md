# Photonic head-to-head: latest run

2026-11-12T1658Z on dgx-spark (Cortex-A725+Cortex-X925), every program pinned to one core with one thread. PQCF at commit ef79196. PQCF 1.4.0, The Walrus 0.22.0, Piquasso 8.0.2, Perceval 1.3.1, DeepQuantum 4.5.0, Strawberry Fields 0.23.0. Written by `run.sh`; the method is in [README.md](README.md).

**Answers: 51 of 51 checks agree.**

*Not run*: a smaller size took more than 120 s or could not run, so the program's larger sizes were skipped.

## Permanent of an n × n block of a random unitary (seconds; in brackets, the time relative to PQCF's)

| n | PQCF | The Walrus | Piquasso | Perceval | DeepQuantum |
| --- | --- | --- | --- | --- | --- |
| 12 | 6.00e-05 | 7.78e-05 (1.3×) | 2.31e-04 (3.8×) | 3.20e-05 (0.53×) | 0.0032 (53×) |
| 16 | 0.0014 | 0.0014 (1×) | 0.0014 (0.99×) | 7.14e-04 (0.5×) | 0.0597 (42×) |
| 20 | 0.0291 | 0.0262 (0.9×) | 0.0216 (0.74×) | 0.0166 (0.57×) | 2.86 (98×) |
| 24 | 0.5889 | 0.4736 (0.8×) | 0.3916 (0.66×) | 0.3375 (0.57×) | 58.58 (99×) |
| 26 | 2.56 | 2.01 (0.79×) | 1.80 (0.7×) | 1.48 (0.58×) | 260.10 (1e+02×) |
| 28 | 11.07 | 8.62 (0.78×) | 7.74 (0.7×) | 6.44 (0.58×) | not run |

## Torontonian of a k-mode Gaussian state, every mode clicking (seconds)

| k | PQCF | The Walrus | Piquasso | DeepQuantum |
| --- | --- | --- | --- | --- |
| 8 | 5.70e-05 | 1.88e-04 (3.3×) | 2.39e-05 (0.42×) | 0.0015 (27×) |
| 12 | 8.65e-04 | 0.0037 (4.3×) | 5.29e-04 (0.61×) | 0.0166 (19×) |
| 16 | 0.0191 | 0.0726 (3.8×) | 0.0095 (0.5×) | 0.8684 (45×) |
| 20 | 0.3955 | 1.40 (3.6×) | 0.1707 (0.43×) | 24.96 (63×) |
| 22 | 1.78 | 6.18 (3.5×) | 0.7137 (0.4×) | 114.27 (64×) |
| 24 | 8.00 | 27.05 (3.4×) | 2.97 (0.37×) | 546.11 (68×) |

## Boson sampling: n photons in n² modes, exact samples per second (in brackets, PQCF's rate relative to the program's)

| photons | PQCF | Piquasso | Perceval | DeepQuantum |
| --- | --- | --- | --- | --- |
| 4 | 1,324,328 | 4,624 (2.9e+02×) | 19,973 (66×) | 33,686 (39×) |
| 8 | 121,855 | 334 (3.6e+02×) | 1,997 (61×) | not run (1.06e+10 output patterns) |
| 12 | 5,700 | 121 (47×) | 199 (29×) | not run |
| 16 | 264 | 36.9 (7.2×) | 19.7 (13×) | not run |
| 20 | 12.8 | 3.16 (4×) | not run (RuntimeError: Fock state size is capped from 0 to 256 modes) | not run |
| 24 | 0.643 | 0.144 (4.4×) | not run | not run |

## Fock simulation: one layer of a continuous-variable neural network on 10 modes (seconds)

PQCF and Piquasso keep every state with at most c photons in total; Strawberry Fields and DeepQuantum keep at most c photons in each mode. *Truncation error*: the largest difference of a mode's mean photon number from PQCF's largest cutoff.

| program | truncation | cutoff | states | seconds | truncation error |
| --- | --- | --- | --- | --- | --- |
| PQCF | total | 4 | 1,001 | 0.0053 | 1.9e-04 |
| PQCF | total | 6 | 8,008 | 0.0435 | 8.4e-06 |
| PQCF | total | 8 | 43,758 | 0.2333 | 2.9e-07 |
| PQCF | total | 10 | 184,756 | 1.10 | 0.0e+00 |
| Piquasso | total | 4 | 1,001 | 0.0061 | 1.9e-04 |
| Piquasso | total | 6 | 8,008 | 0.0187 | 8.4e-06 |
| Piquasso | total | 8 | 43,758 | 0.0799 | 2.9e-07 |
| Piquasso | total | 10 | 184,756 | 0.3481 | 5.2e-16 |
| DeepQuantum | per mode | 3 | 59,049 | 0.0179 | 8.4e-04 |
| DeepQuantum | per mode | 4 | 1,048,576 | 0.1456 | 4.9e-04 |
| Strawberry Fields | per mode | 3 | 59,049 | 0.8192 | 8.4e-04 |
| Strawberry Fields | per mode | 4 | 1,048,576 | 10.90 | 4.9e-04 |

## Every check

| task | what | agrees | detail |
| --- | --- | --- | --- |
| permanent | PQCF, size 12 | yes | relative difference 2.7e-14 (limit 1e-08) |
| permanent | PQCF, size 16 | yes | relative difference 1.0e-13 (limit 1e-08) |
| permanent | PQCF, size 20 | yes | relative difference 5.5e-12 (limit 1e-08) |
| permanent | PQCF, size 24 | yes | relative difference 2.5e-10 (limit 1e-08) |
| permanent | PQCF, size 26 | yes | relative difference 3.8e-10 (limit 1e-08) |
| permanent | PQCF, size 28 | yes | relative difference 4.6e-09 (limit 1e-08) |
| permanent | Piquasso, size 12 | yes | relative difference 3.5e-15 (limit 1e-08) |
| permanent | Piquasso, size 16 | yes | relative difference 9.2e-14 (limit 1e-08) |
| permanent | Piquasso, size 20 | yes | relative difference 6.3e-12 (limit 1e-08) |
| permanent | Piquasso, size 24 | yes | relative difference 2.5e-10 (limit 1e-08) |
| permanent | Piquasso, size 26 | yes | relative difference 3.9e-10 (limit 1e-08) |
| permanent | Piquasso, size 28 | yes | relative difference 4.8e-09 (limit 1e-08) |
| permanent | Perceval, size 12 | yes | relative difference 5.5e-14 (limit 1e-08) |
| permanent | Perceval, size 16 | yes | relative difference 3.8e-13 (limit 1e-08) |
| permanent | Perceval, size 20 | yes | relative difference 2.6e-11 (limit 1e-08) |
| permanent | Perceval, size 24 | yes | relative difference 2.8e-10 (limit 1e-08) |
| permanent | Perceval, size 26 | yes | relative difference 8.1e-10 (limit 1e-08) |
| permanent | Perceval, size 28 | yes | relative difference 7.4e-09 (limit 1e-08) |
| permanent | DeepQuantum, size 12 | yes | relative difference 5.2e-13 (limit 1e-08) |
| permanent | DeepQuantum, size 16 | yes | relative difference 1.1e-12 (limit 1e-08) |
| permanent | DeepQuantum, size 20 | yes | relative difference 2.5e-10 (limit 1e-08) |
| permanent | DeepQuantum, size 24 | yes | relative difference 3.3e-10 (limit 1e-08) |
| permanent | DeepQuantum, size 26 | yes | relative difference 2.4e-10 (limit 1e-08) |
| torontonian | PQCF, size 8 | yes | relative difference 6.0e-15 (limit 1e-07) |
| torontonian | PQCF, size 12 | yes | relative difference 2.4e-15 (limit 1e-07) |
| torontonian | PQCF, size 16 | yes | relative difference 1.6e-13 (limit 1e-07) |
| torontonian | PQCF, size 20 | yes | relative difference 5.4e-14 (limit 1e-07) |
| torontonian | PQCF, size 22 | yes | relative difference 5.0e-13 (limit 1e-07) |
| torontonian | PQCF, size 24 | yes | relative difference 6.2e-13 (limit 1e-07) |
| torontonian | Piquasso, size 8 | yes | relative difference 3.0e-16 (limit 1e-07) |
| torontonian | Piquasso, size 12 | yes | relative difference 1.6e-13 (limit 1e-07) |
| torontonian | Piquasso, size 16 | yes | relative difference 1.6e-12 (limit 1e-07) |
| torontonian | Piquasso, size 20 | yes | relative difference 1.6e-11 (limit 1e-07) |
| torontonian | Piquasso, size 22 | yes | relative difference 4.2e-12 (limit 1e-07) |
| torontonian | Piquasso, size 24 | yes | relative difference 3.6e-11 (limit 1e-07) |
| torontonian | DeepQuantum, size 8 | yes | relative difference 9.7e-15 (limit 1e-07) |
| torontonian | DeepQuantum, size 12 | yes | relative difference 5.0e-14 (limit 1e-07) |
| torontonian | DeepQuantum, size 16 | yes | relative difference 3.1e-14 (limit 1e-07) |
| torontonian | DeepQuantum, size 20 | yes | relative difference 9.6e-13 (limit 1e-07) |
| torontonian | DeepQuantum, size 22 | yes | relative difference 9.8e-13 (limit 1e-07) |
| torontonian | DeepQuantum, size 24 | yes | relative difference 1.4e-11 (limit 1e-07) |
| boson_sampling | PQCF, 3 photons in 6 modes, 20000 samples | yes | chi-square p = 0.610, 55 degrees of freedom |
| boson_sampling | Piquasso, 3 photons in 6 modes, 20000 samples | yes | chi-square p = 0.052, 55 degrees of freedom |
| boson_sampling | Perceval, 3 photons in 6 modes, 20000 samples | yes | chi-square p = 0.175, 55 degrees of freedom |
| boson_sampling | DeepQuantum, 3 photons in 6 modes, 20000 samples | yes | chi-square p = 0.011, 55 degrees of freedom |
| fock | Piquasso against PQCF, total photons <= 4 | yes | largest difference of a mode's mean photon number 3.8e-17 |
| fock | Piquasso against PQCF, total photons <= 6 | yes | largest difference of a mode's mean photon number 7.3e-17 |
| fock | Piquasso against PQCF, total photons <= 8 | yes | largest difference of a mode's mean photon number 2.9e-16 |
| fock | Piquasso against PQCF, total photons <= 10 | yes | largest difference of a mode's mean photon number 5.2e-16 |
| fock | DeepQuantum against Strawberry Fields, 3 photons per mode | yes | largest difference of a mode's mean photon number 1.5e-07 |
| fock | DeepQuantum against Strawberry Fields, 4 photons per mode | yes | largest difference of a mode's mean photon number 3.0e-07 |
