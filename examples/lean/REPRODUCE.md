# Reproducing the Lean 4 certificates of PQCF's published results

This page takes you from nothing to Lean 4 checking PQCF's published results, on your machine. Part 1 covers the
symbolic engine: the polynomial gcds (task T7) and factorisations (task T8) of the published comparison with FORM
(`comparison/form`). What the certificates are and what they prove: [docs/lean-certificates.md](../../docs/lean-certificates.md).

## What you need

- **Docker**, about 15 GB of disk: PQCF's development image (2.5 GB) and lean-check's image (about 12 GB: Lean
  4.36.0, the Lean REPL, Mathlib's prebuilt cache, and lean-check's certificate library).
- **git**, and a network connection for the first builds (the certificate checks themselves run without network).
- A Linux machine, x86_64 or aarch64, with about 16 GB of memory (a Lean worker holds Mathlib).

## Steps

```sh
git clone https://github.com/quanterall/pqcf.git
git clone https://github.com/quanterall/erlang-lean4.git
cd pqcf
docker build -t pqcf-dev .                                   # a few minutes
(cd ../erlang-lean4 && docker build -t lean-check docker)    # 15-30 minutes the first time
LEAN_CHECK_SRC=../erlang-lean4 DEMO=certify_published sh examples/lean/run.sh
```

The script starts a lean-check node and a PQCF node on an internal Docker network (no outside access) and joins
them by Erlang distribution. The PQCF node regenerates the published problems from their fixed seeds
(`compare:gcd_problems("nontrivial")`, `compare:factor_problems("nontrivial")`: T. Ueda's polybench defaults, 5
variables, basic polynomials of 22-30 terms of degree 22-30), computes each result with PQCF's engine, writes it
as a Lean statement with its proof (`pq_lean`), and sends it to the lean-check node, which checks it with Lean.

## What you should see

One line per result, then a summary:

```
T7 gcd 1: 575 input terms, Lean source 177902 bytes -> valid in 66.4 s, axioms [<<"Classical.choice">>,
                                                                                <<"Quot.sound">>,
                                                                                <<"propext">>]
...
T8 factor 5: 623 input terms, Lean source 94151 bytes -> valid in 29.0 s, axioms [...]

summary: 15 of 15 valid
```

Our run (25 Nov 2026, aarch64, one Lean worker) is in [results/part1-2026-11-25.txt](results/part1-2026-11-25.txt):
the ten gcds (572-751 input terms) in 58-99 s each, the five factorisations (572-722 terms) in 27-39 s each, 16
minutes in all, every one valid from the standard axioms.

`valid` means Lean's kernel accepted the certificate: for a gcd G of A and B, that A = G A' and B = G B' with the
cofactors PQCF computed; for a factorisation, that the product of the factors (with their multiplicities and the
constant) is the input polynomial. The axioms listed are Lean's standard ones (`propext`, `Classical.choice`,
`Quot.sound`) and nothing else.

## How to compare with the published results

The problems are the published ones (fixed seeds); every line must say `valid`, and the summary `15 of 15 valid`.
The times depend on your machine.

## Tested

From fresh clones of both public repositories on 25 Nov 2026 (PQCF dd472fc, erlang-lean4 137489b), on an aarch64
machine (Cortex-X925 and Cortex-A725, 20 cores): both images built from the clones (with Docker's layer cache: the
Mathlib download was not repeated), then the script: 15 of 15 valid, the same Lean sources byte for byte as our run,
28-104 s per certificate, 16.5 minutes in all.
