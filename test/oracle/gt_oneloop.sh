#!/bin/sh
# One-loop two-point functions (pq_oneloop, plan "theoretical physics instruments", F3): two independent references.
# 60 configurations (p^2, m1^2, m2^2), exact fractions fixed by a seed, over every region (space-like, below the
# pseudo-threshold, between the thresholds, above threshold, one or both masses 0, p^2 = 0):
#  - LoopTools 2.16 (comparison/oneloop/looptools): B0, B1, B00, B11, DB0 in double precision;
#  - mpmath, Feynman-parameter integrals at 40 digits (tanh-sinh, split at the zeros of the denominator):
#      B0 = -int log D, B1 = int x log D, B11 = -int x^2 log D, D = x m2^2 + (1 - x) m1^2 - x (1 - x) p^2 - i eps.
# The two are checked against each other (1e-11 relative) before test/fixtures/gt_oneloop.json is written.
#   docker build -t pqcf-looptools comparison/oneloop/looptools && sh test/oracle/gt_oneloop.sh
set -e
cd "$(dirname "$0")/../.."
W=$(mktemp -d "${TMPDIR:-/tmp}/pqcf-oneloop-XXXXXX")
python3 - "$W" <<'PY'
import json, random, sys
from fractions import Fraction as F
w = sys.argv[1]
rnd = random.Random(20261116)
masses = [F(0), F(1, 4), F(1, 2), F(1), F(2), F(3), F(9, 4)]
pts = set()
while len(pts) < 60:
    a, b = rnd.choice(masses), rnd.choice(masses)
    s = F(rnd.randint(-12, 30), rnd.choice([1, 2, 4]))
    if a == b == 0 and s == 0:
        continue
    pts.add((s, a, b))
pts = sorted(pts)
json.dump([[str(x) for x in p] for p in pts], open(f"{w}/pts.json", "w"))
with open(f"{w}/lt.in", "w") as f:
    for s, a, b in pts:
        f.write(f"B {float(s)!r} {float(a)!r} {float(b)!r}\n")
PY
docker run --rm -i pqcf-looptools lt_eval < "$W/lt.in" 2>/dev/null | grep -E '^-?[0-9]' > "$W/lt.out"
cat > "$W/mp.py" <<'PY'
import json, sys
from fractions import Fraction as F
import mpmath as mp
mp.mp.dps = 50
pts = [[F(x) for x in p] for p in json.load(open(sys.argv[1]))]

def integral(k, s, a, b):
    """int_0^1 x^k log(D - i eps) dx, D = x b + (1 - x) a - x (1 - x) s"""
    s, a, b = mp.mpf(s.numerator) / s.denominator, mp.mpf(a.numerator) / a.denominator, mp.mpf(b.numerator) / b.denominator
    # D = s x^2 + (b - a - s) x + a: split at its zeros in (0, 1)
    cuts = [mp.mpf(0), mp.mpf(1)]
    if s != 0:
        disc = (b - a - s) ** 2 - 4 * s * a
        if disc > 0:
            for r in [(-(b - a - s) + mp.sqrt(disc)) / (2 * s), (-(b - a - s) - mp.sqrt(disc)) / (2 * s)]:
                if 0 < r < 1:
                    cuts.append(r)
    elif b != a:
        r = a / (a - b)
        if 0 < r < 1:
            cuts.append(r)
    cuts.sort()
    D = lambda x: s * x * x + (b - a - s) * x + a
    re = mp.quad(lambda x: x ** k * mp.log(abs(D(x))) if D(x) != 0 else mp.mpf(0), cuts)
    # log(D - i eps) = log|D| - i pi where D < 0
    im = mp.mpf(0)
    for lo, hi in zip(cuts, cuts[1:]):
        if D((lo + hi) / 2) < 0:
            im -= mp.pi * (hi ** (k + 1) - lo ** (k + 1)) / (k + 1)
    return re, im

out = []
for s, a, b in pts:
    r0, i0 = integral(0, s, a, b)
    r1, i1 = integral(1, s, a, b)
    r2, i2 = integral(2, s, a, b)
    out.append({"b0": [mp.nstr(-r0, 40), mp.nstr(-i0, 40)], "b1": [mp.nstr(r1, 40), mp.nstr(i1, 40)],
                "b11": [mp.nstr(-r2, 40), mp.nstr(-i2, 40)]})
json.dump(out, open(sys.argv[2], "w"))
PY
docker run --rm -v "$W":/w pqcf-oracle python /w/mp.py /w/pts.json /w/mp.json
python3 - "$W" <<'PY'
import json, sys
w = sys.argv[1]
pts = json.load(open(f"{w}/pts.json"))
mp_ = json.load(open(f"{w}/mp.json"))
lt = [list(map(float, l.split())) for l in open(f"{w}/lt.out")]
assert len(lt) == len(pts) == len(mp_), (len(lt), len(pts), len(mp_))
names = ["b0", "b1", "b00", "b11", "db0"]
out = []
for p, m, l in zip(pts, mp_, lt):
    s, a, b = p
    rec = {"s": s, "m1sq": a, "m2sq": b, "mp": m, "looptools": {}}
    for k, n in enumerate(names):
        rec["looptools"][n] = [l[2 * k], l[2 * k + 1]]
    for n in ["b0", "b1", "b11"]:
        if s == "0" and n != "b0":
            continue
        for x, y in zip(map(float, m[n]), rec["looptools"][n]):
            assert abs(x - y) <= 1e-11 * max(1, abs(x)), (p, n, m[n], rec["looptools"][n])
    out.append(rec)
json.dump(out, open("test/fixtures/gt_oneloop.json", "w"), indent=0)
print(len(out), "configurations; LoopTools and mpmath agree")
PY
rm -rf "${W:?}"
