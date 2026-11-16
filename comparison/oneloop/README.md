# One-loop integrals: PQCF against LoopTools, QCDLoop and mpmath

`pq_oneloop` computes the one-loop one- and two-point functions (A0, B0, B1, B00, B11 and dB0/dp^2), the
three-point function C0 and the four-point function D0 at arbitrary precision for real kinematics, in LoopTools'
conventions (finite parts, Delta = 0, scale mu^2).

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

## The three-point function and the two external programs

In our builds on this aarch64 machine (an x86 build was not available to compare), both external programs are wrong
in some configurations where the Kallen function lambda(p1^2, p2^2, p3^2) is negative; the cause (a limit of the
methods or the platform) is not established:

- LoopTools 2.16: C0 returns wrong values when lambda < 0, e.g. C0(1, 1, 1; 1, 1, 1) = 6.87 + 34.2 i (FF reports
  "lambda(p1,p2,p3) < 0, unphysical configuration"; the value is the real -0.671253105748004686), and the same call
  from Fortran does not return. `LTVERSION` (the alternative versions) changes nothing. Its two-point functions are
  right (they agree with mpmath to 1e-11) and so is C0 for lambda > 0.
- QCDLoop 2.1.0 (comparison/oneloop/qcdloop): right below threshold and for lambda > 0, but wrong above threshold
  with lambda < 0, e.g. C0(9/2, 7, 20; 1, 1/4, 2) = 0.42341 - 0.01664 i, where two independent mpmath computations
  (contour deformation at two values of kappa) and PQCF give 0.408381951022 - 0.152701852704 i.

An earlier version of this page (and of the CHANGELOG) said LoopTools' default C0 was wrong and its Denner version
right; that was a misreading of a debug line: no version of our LoopTools build returns the right value there.

The C0 reference is therefore mpmath (`test/oracle/gt_c0.py`): contour deformation for every configuration, and the
direct integral at 30 digits where D > 0 on the whole triangle. The values of QCDLoop and LoopTools are recorded in
`test/fixtures/gt_c0.json` with a flag for each disagreement.

In the 40 configurations of `test/oracle/gt_c0.sh` (drawn by a seed over space-like, below- and above-threshold
kinematics, both signs of lambda, massless lines):

| reference | configurations | PQCF agrees |
| --- | --- | --- |
| mpmath, direct integral (D > 0) | 9 | to 28 digits |
| mpmath, contour deformation, stable under kappa = 0.05 / 0.1 | 28 | to 1e-12 (largest difference 1.3e-14) |
| two of deformation (kappa = 0.1), QCDLoop, LoopTools agreeing to 1e-12 | 11 of the other 12 | to 1e-12 |
| the deformation alone (both programs wrong, deformation not stable) | 1 | to 1e-6 (in fact 1e-16) |

QCDLoop disagrees with PQCF in 6 configurations (and by 7e-9 in one more), LoopTools in 7, always with
lambda < 0; the deformation is unstable (its two values differ by more than 1e-12, up to 0.14) in 12, mostly with
massless lines, where the corners of the triangle are singular. Where the deformation is unstable and the two
programs agree, they agree with PQCF to about 1e-15. (`gt_c0.sh`'s own summary counts 12 disagreements for each
program: its flags compare with the deformation, which is itself off in some of those configurations.)

## The four-point function

`test/oracle/gt_d0.sh` writes `test/fixtures/gt_d0.json`: three Euclidean configurations by mpmath's direct 3D
integral (20 digits) and four above threshold by LoopTools and QCDLoop. PQCF equals mpmath to 16 digits and LoopTools
to about 1e-15 in all seven; QCDLoop returns `nan` for two of the four above threshold (D0(2, 3, 1, 4, 20, -5;
1, 2, 1, 3) and D0(1, 2, 3, 4, 30, -10; 2, 1, 2, 1)) and agrees with the others in the other two. D0 is slow above
threshold (45 s to 5 minutes at 18 digits on 20 cores); the test suite checks three of the four.
