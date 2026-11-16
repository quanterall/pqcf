# One-loop integrals: PQCF against LoopTools and mpmath

`pq_oneloop` computes the one-loop one- and two-point functions (A0, B0, B1, B00, B11 and dB0/dp^2) at arbitrary
precision for real kinematics, in LoopTools' conventions (finite parts, Delta = 0, scale mu^2).

The references are independent of each other and of PQCF:

- **LoopTools 2.16** (T. Hahn and M. Perez-Victoria, [arXiv:hep-ph/9807565](https://arxiv.org/abs/hep-ph/9807565)),
  built in `looptools/` from the release tarball, checked by its SHA-256; double precision;
- **mpmath**: the Feynman-parameter integrals at 40 digits, split at the zeros of the denominator.

`test/oracle/gt_oneloop.sh` draws 60 configurations (p^2, m1^2, m2^2), exact fractions fixed by a seed, over every
region: space-like, below the pseudo-threshold, between the thresholds, above threshold, one or both masses zero,
p^2 = 0. It checks that the two references agree (1e-11 relative) and writes `test/fixtures/gt_oneloop.json`.
`test/pq_oneloop_tests.erl` then requires PQCF to equal mpmath to 35 digits (B0, B1, B11) and LoopTools to 1e-11 (all
five functions), and to refuse the points where a function is not finite (dB0 at threshold or on shell with a
massless line, the tensor coefficients at p^2 = 0).

## Reproduce

```sh
git clone https://github.com/quanterall/pqcf && cd pqcf
docker build -t pqcf-dev .
docker build -t pqcf-oracle test/oracle                     # mpmath
docker build -t pqcf-looptools comparison/oneloop/looptools  # LoopTools 2.16, about a minute
sh test/oracle/gt_oneloop.sh                                 # regenerates the fixture; prints "60 configurations; ..."
git diff --stat test/fixtures/gt_oneloop.json                # no change expected
docker run --rm pqcf-dev rebar3 eunit --module=pq_oneloop_tests
```

## A finding on the way

LoopTools' default three-point function (version "a", from the FF library) returns a complex value,
C0(1, 1, 1; 1, 1, 1) = 6.87 + 34.2 i, at a point below every threshold, where C0 is real; it reports
"lambda(p1,p2,p3) < 0, unphysical configuration". Its alternative version "b" (Denner's code, `LTVERSION=48`)
gives the real value -0.6712531057480048. The three- and four-point functions in PQCF (next part of F3) are compared
with version "b", and the points where the two versions differ are recorded.
