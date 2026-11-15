# Reproducing the three-loop QCD β-function

This page takes you from nothing to the results in [README.md](README.md): the β-function of QCD to three loops from
the gluon and ghost self-energies in Landau gauge, with FORM, Forcer and PQCF, and the check that it is the known
one (Tarasov, Vladimirov, Zharkov 1980; Larin and Vermaseren, [arXiv:hep-ph/9302208](https://arxiv.org/abs/hep-ph/9302208)).

## What you need

- **Docker**, and about 8 GB of disk for three images: PQCF's development image (2.5 GB), the FORM image (FORM 5 and
  FLINT from their sources) and the Forcer image (FORM plus Forcer), the last two sharing most of their layers.
- **git** and a network connection: the images, and FORM's colour package (`color.h`), fetched at a fixed commit and
  checked by its SHA-256.
- Any Linux machine, x86_64 or aarch64, with 8 cores for the time below (the three-loop gluon self-energy is the long
  step).

## Steps

**Reproduce the published run**, at the commit it was tested at (see "Tested" below):

```sh
git clone https://github.com/quanterall/pqcf.git
cd pqcf
git checkout eef0e6f
sh comparison/beta/run.sh 0-7
```

**Run the latest** version: the same commands without `git checkout eef0e6f`. The answers must still be
the published ones; what changed since is in [CHANGELOG.md](../../CHANGELOG.md).

`0-7` are the cores (a Docker cpuset). [run.sh](run.sh) builds the three images and then, in a temporary directory:

1. **FORM + Forcer** compute the bare gluon and ghost self-energies at one, two and three loops ([se.frm](se.frm):
   FORM's own diagram generator, QCD's Feynman rules, the colour factors by FORM's colour package, Dirac traces in d
   dimensions, then Forcer);
2. **PQCF** reduces every one- and two-loop integral itself (FORM writes the integrands, `beta_pqcf:reduce/3`
   reduces them to master integrals with exact coefficients), and **Forcer** gives the values of PQCF's masters;
3. **PQCF** renormalises both and extracts the β-function ([renorm.erl](renorm.erl), [beta_run.erl](beta_run.erl)):
   Z_3 and Z̃_3 order by order, Z_g² = 1/(Z̃_3² Z_3), β from its single poles, the higher poles checked against the
   lower orders.

It rewrites the files of [results/](results) and checks the three β coefficients line by line against the known ones.
On our machine it takes about 15 minutes on 8 cores, most of it the three-loop gluon.

## What you should see

FORM's progress, one line per self-energy (`FORM + Forcer: glu, 3 loops`, ...), PQCF's reductions, then the β
coefficients as maps (colour factor → coefficient), and at the end:

```
beta functions: as known
```

The script exits with 1 if the β-functions differ from the known ones or if either renormalisation is not finite.

## How to compare with the published results

- **The results files** must be unchanged: `git status comparison/beta/results` shows nothing. The three-loop result
  ([renormalisation_form_forcer_3loops.txt](results/renormalisation_form_forcer_3loops.txt)) has
  `beta_2 = ... + 2857/54*cA^3`, and for SU(3), 2857/2 − 5033 n_f/18 + 325 n_f²/54; the two-loop result computed with
  PQCF's own reductions ([renormalisation_pqcf_2loops.txt](results/renormalisation_pqcf_2loops.txt)) has β₀ and β₁.
- Every coefficient is an exact rational; there is no tolerance.

## What this does not show

At three loops PQCF does not reduce the integrals itself: Forcer does (README.md, "Where PQCF stands"). The
three-loop β-function here checks the pipeline (diagrams, Feynman rules, colour, renormalisation by PQCF) against the
known result; PQCF's own reductions are checked at one and two loops.

## Tested

From a clean clone on 15 November 2026 (commit 92e2aeb with this script, as committed in eef0e6f), aarch64 (Cortex-X925
and Cortex-A725), cores 9-16: `beta functions: as known`, every results file unchanged, 870 seconds.
