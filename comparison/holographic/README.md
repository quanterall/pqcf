# Holographic QCD in PQCF: how to reproduce

Seven computations of decision record [SY2](../../docs/decisions/SY2-symbolic-foundations.md), each from the model's
action or equations to a published number. Each script writes its result file in this directory; the committed files
are the results of the commit named at their top.

## What you need

Docker, and this repository. Everything runs in PQCF's development image, built from the repository's
[Dockerfile](../../Dockerfile):

```sh
git clone https://github.com/quanterall/pqcf.git && cd pqcf
git checkout d77ead6
docker build -t pqcf-dev .
```

`git checkout d77ead6` **reproduces the published runs**: it is the commit of this page, with the code they were
tested at (see below). To **run the latest** version, leave it out; the results must still be the committed ones,
and what changed since is in [CHANGELOG.md](../../CHANGELOG.md).

Then each script runs from the repository root with the command at the top of its source file, of this form
(`<script>` one of the names in the table below):

```sh
docker run --rm -u $(id -u):$(id -g) -e PQCF_COMMIT=$(git rev-parse --short HEAD) -v "$PWD":/w -w /w pqcf-dev \
  sh -c 'rebar3 compile >/dev/null && erlc -o /tmp comparison/holographic/<script>.erl &&
  erl -noshell -pa _build/default/lib/*/ebin -pa /tmp -eval "<script>:main(\"comparison/holographic\"), halt()."'
```

It rewrites the result file. `git diff` then shows whether anything changed: only the commit line, the date and the
times should. The first run also compiles PQCF, which takes about a minute.

## The computations

| script | result file | model and paper | checked against | time (one core unless stated) |
| --- | --- | --- | --- | --- |
| `walls` | [RESULTS.md](RESULTS.md) | soft wall ([hep-ph/0602229](https://arxiv.org/abs/hep-ph/0602229)), hard wall ([hep-ph/0501128](https://arxiv.org/abs/hep-ph/0501128)): vector mesons from the action | the exact spectra, 30 digits | 12 s |
| `hardwall` | [HARDWALL.md](HARDWALL.md) | hard wall, the paper's two fits: masses, decay constants, f_pi, m_pi | the paper's Table II | 8 s |
| `sakai` | [SAKAI.md](SAKAI.md) | Sakai-Sugimoto ([hep-th/0412141](https://arxiv.org/abs/hep-th/0412141)): vector and axial mesons | the paper's eigenvalues | 16 s |
| `brane` | [BRANE.md](BRANE.md) | the AdS_5 black brane ([hep-th/0002230](https://arxiv.org/abs/hep-th/0002230)): free energy by holographic renormalisation | F = -pi^2 N^2 T^4/8, exactly | 1 s |
| `weyl` | [WEYL.md](WEYL.md) | the Weyl anomaly ([hep-th/9806087](https://arxiv.org/abs/hep-th/9806087)) on S^4 and S^2 x S^2 | a = c = N^2/4, exactly | 4 s |
| `gubser` | [GUBSER.md](GUBSER.md) | Gubser-Nellore ([0804.0434](https://arxiv.org/abs/0804.0434)): speed of sound against temperature | the paper's parameters, limits and figures | 4 min |
| `ihqcd` | [IHQCD.md](IHQCD.md) | improved holographic QCD ([0707.1349](https://arxiv.org/abs/0707.1349)): scalar and tensor glueballs | the paper's mass ratios | 40 min on 16 cores; quick check 4 min |

**Tested from a fresh clone** of this repository on 15 November 2026 (commit 73c0246, a new image built without cache,
the commands above unchanged): every script reproduced its committed result file, differing only in the commit line,
the date and the times.

## Improved holographic QCD: the quick check

The full run of `ihqcd` computes about 30 spectra (several starting values, cutoffs, a fit of b0, the first 12
states) and takes about 40 minutes on 16 cores. The quick check computes the two ratios the paper fits, R00 and R20,
in about 4 minutes, and compares the masses with the committed raw results
([ihqcd-raw.etf](ihqcd-raw.etf), Erlang's external term format):

```sh
docker run --rm -u $(id -u):$(id -g) -v "$PWD":/w -w /w pqcf-dev \
  sh -c 'rebar3 compile >/dev/null && erlc -o /tmp comparison/holographic/ihqcd.erl &&
  erl -noshell -pa _build/default/lib/*/ebin -pa /tmp -eval "ihqcd:quick(\"comparison/holographic\"), halt()."'
```

It prints:

```
m^2: 0++ 1.8131136152, 6.3763556293; 2++ 3.5214364768
committed (ihqcd-raw.etf): 0++ 1.8131136152, 6.3763556293; 2++ 3.5214364768
R00 = 1.8753 (paper 1.866), R20 = 1.3936 (paper 1.393)
PASS: equal to the committed results to 1e-7
```

`ihqcd:report_only("comparison/holographic")` rewrites IHQCD.md from the raw results without computing, in a second.

## Independent checks in double precision

Three of the computations are also solved by a separate program: SciPy's DOP853, in the paper's own coordinates,
sharing no code or method with PQCF. Each script reads PQCF's committed result file and exits 1 if they differ by
more than its tolerance. They run in the oracle image ([test/oracle](../../test/oracle/Dockerfile): SciPy 1.17.1),
in seconds:

```sh
docker build -t pqcf-oracle test/oracle
for s in sakai hardwall gubser; do
  docker run --rm -v "$PWD":/w -w /w pqcf-oracle python comparison/holographic/scipy_$s.py
done
```

| script | checks | tolerance | agreement |
| --- | --- | --- | --- |
| [scipy_sakai.py](scipy_sakai.py) | the first four Sakai-Sugimoto eigenvalues (SAKAI.md) | relative 1e-9 | 3e-11 or better |
| [scipy_hardwall.py](scipy_hardwall.py) | m_a1 and the chiral-limit f_pi of both models (HARDWALL.md) | relative 1e-7, the published rounding | m_a1 2e-8, f_pi 1e-9 |
| [scipy_gubser.py](scipy_gubser.py) | c_s^2 for phi_H =< 5 (GUBSER.md) | 2e-6 | 5e-7, the table's rounding |

The hard-wall check found an error in PQCF on 15 November 2026: its f_pi was 1e-6 low, because a Frobenius series
whose terms have gaps was summed too few terms. Fixed (2c971f2), with a test; HARDWALL.md says what changed.

## What these results are, and are not

- **Precision.** Soft wall, hard wall and Sakai-Sugimoto are computed to 25-30 digits; brane and Weyl are exact. Gubser-Nellore
  and improved holographic QCD are integrated at 12 digits and compared with papers that give 2-4 digits: they are
  numerical results at that precision.
- **Open questions** are stated in each file: the hard wall's f_pi (the paper keeps the quark mass) and one a1 mass;
  the last 0.6 % of the improved holographic QCD ratios.
- **Independent cross-checks.** The SciPy solutions quoted in SAKAI.md, HARDWALL.md and GUBSER.md were run
  separately and are not yet scripts in this repository. The checks against SymPy and mpmath of the engine's series,
  arithmetic and tensors are: `test/oracle/` writes the fixtures that `rebar3 eunit` compares against.
- **History.** IHQCD.md keeps the first run (a = 5/16, 6-15 % low) as a record of the error, and the decision record
  explains it.
