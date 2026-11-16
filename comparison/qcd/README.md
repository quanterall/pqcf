# QCD's renormalisation group: PQCF against mpmath and CRunDec

`pq_qcd` holds the beta function and the quark mass anomalous dimension of QCD to five loops, exactly (rational
numbers and zeta values, polynomials in the number of flavours n_f), and runs alpha_s and MS-bar quark masses at
arbitrary precision; it decouples alpha_s at a heavy-quark threshold (MS-bar mass) to three loops. Conventions are
RunDec's ([arXiv:1703.03751](https://arxiv.org/abs/1703.03751)): a = alpha_s/pi,
da/d ln mu^2 = -sum beta_i a^(i+2), d ln m/d ln mu^2 = -sum gamma_i a^(i+1).

Sources: beta_3 [hep-ph/9701390](https://arxiv.org/abs/hep-ph/9701390), beta_4
[arXiv:1606.08659](https://arxiv.org/abs/1606.08659) and [arXiv:1701.01404](https://arxiv.org/abs/1701.01404);
gamma_3 [hep-ph/9703278](https://arxiv.org/abs/hep-ph/9703278) and [hep-ph/9703284](https://arxiv.org/abs/hep-ph/9703284),
gamma_4 [arXiv:1402.6611](https://arxiv.org/abs/1402.6611); decoupling [hep-ph/9706430](https://arxiv.org/abs/hep-ph/9706430).

## The checks (`test/pq_qcd_tests.erl`)

| check | against | agreement |
| --- | --- | --- |
| beta_0 .. beta_3 | van Ritbergen, Vermaseren and Larin's formula for any gauge group, with SU(3)'s invariants computed by `pq_lie` | exact, term by term |
| beta_0 .. beta_4, gamma_0 .. gamma_4, n_f = 3..6 | a separate transcription in mpmath (`test/oracle/gt_qcd.py`) | 38 digits |
| alpha_s(mu) and m(mu), 10 cases, 1 to 5 loops, n_f = 3..6 | mpmath's ODE solver at 45 digits | 28 digits |
| the same | CRunDec 3.1.2, `AsmMSrunexact` | 1e-8 |
| alpha_s decoupling down and up, 8 cases | mpmath; CRunDec `DecAsDownMS`, `DecAsUpMS` | 28 digits; 1e-14 |
| one loop | the closed form | 55 digits |

## About CRunDec's numbers

CRunDec 3.1.2 (Herren and Steinhauser; MIT licence; its authors' repository at commit 31e0c7f, built in
[crundec/](crundec/Dockerfile)) is the reference program of the field. Two of its functions are not meant to be exact,
which matters when it is compared at more than double precision:

- `AlphasExact` and `AsmMSrunexact` integrate with adaptive Runge-Kutta steps at a tolerance of 1e-10 per step: their
  results differ from the exact solution by about 1e-10 to 1e-9 (e.g. one-loop alpha_s(10 GeV) from 0.118 at
  91.1876 GeV: 0.17308363604581 against the closed form's 0.17308363622084).
- `mMS2mMS` evaluates RunDec's expanded formula for the mass, not the exact integral of gamma/beta; at five loops it
  differs from the exact running by about 1e-5 relative (2.85763 against 2.85765 from m_b(m_b) = 4.18 to M_Z).

## Reproduce

```sh
docker build -t pqcf-dev . && docker build -t pqcf-oracle test/oracle && docker build -t pqcf-crundec comparison/qcd/crundec
sh test/oracle/gt_qcd.sh            # a few seconds; rewrites test/fixtures/gt_qcd.json
git diff --stat test/fixtures       # no change expected
docker run --rm pqcf-dev rebar3 eunit --module=pq_qcd_tests
```
