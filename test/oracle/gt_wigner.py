"""Wigner 3j, 6j and 9j symbols and Clebsch-Gordan coefficients by SymPy (sympy.physics.wigner), the reference for
pq_lie (plan "theoretical physics instruments", F1). Every value v is written as its sign and its exact square v^2
(a rational for 3j, 6j and Clebsch-Gordan; 9j as a 40-digit decimal, since a 9j symbol is a sum of square roots).
The angular momenta are given doubled (2j, 2m). Writes test/fixtures/gt_wigner.json.
    docker run --rm -v "$PWD":/src -w /src pqcf-oracle python test/oracle/gt_wigner.py
"""
import itertools
import json
import random

from sympy import Rational, expand, sign
from sympy.physics.wigner import clebsch_gordan, wigner_3j, wigner_6j, wigner_9j

H = lambda x: Rational(x, 2)  # noqa: E731


def sq(v):
    """the sign and the exact square of a signed square root of a rational"""
    s2 = expand(v * v)
    if not s2.is_Rational:
        raise ValueError(f"not a square root of a rational: {v}")
    return [int(sign(v)), str(s2)]


def ok3(j1, j2, j3):
    return abs(j1 - j2) <= j3 <= j1 + j2 and (j1 + j2 + j3) % 2 == 0


out = {"3j": [], "6j": [], "9j": [], "cg": []}
# every 3j with 2j up to 6 (j up to 3), all m (zero ones included: selection rules)
for j1, j2, j3 in itertools.product(range(7), repeat=3):
    if not ok3(j1, j2, j3):
        continue
    for m1 in range(-j1, j1 + 1, 2):
        for m2 in range(-j2, j2 + 1, 2):
            m3 = -m1 - m2
            if abs(m3) > j3:
                continue
            out["3j"].append([j1, j2, j3, m1, m2, m3] + sq(wigner_3j(H(j1), H(j2), H(j3), H(m1), H(m2), H(m3))))
rnd = random.Random(20261116)
# a random sample of 3j and Clebsch-Gordan with 2j up to 20
while len(out["cg"]) < 300:
    j1, j2 = rnd.randint(0, 20), rnd.randint(0, 20)
    j3 = rnd.randrange(abs(j1 - j2), j1 + j2 + 1, 2)
    m1, m2 = rnd.randrange(-j1, j1 + 1, 2), rnd.randrange(-j2, j2 + 1, 2)
    if abs(m1 + m2) > j3:
        continue
    out["3j"].append([j1, j2, j3, m1, m2, -m1 - m2] + sq(wigner_3j(H(j1), H(j2), H(j3), H(m1), H(m2), H(-m1 - m2))))
    out["cg"].append([j1, m1, j2, m2, j3, m1 + m2] + sq(clebsch_gordan(H(j1), H(j2), H(j3), H(m1), H(m2), H(m1 + m2))))
# 6j: all with 2j up to 6, and a random sample up to 16
for js in itertools.product(range(7), repeat=6):
    j1, j2, j3, j4, j5, j6 = js
    if ok3(j1, j2, j3) and ok3(j1, j5, j6) and ok3(j4, j2, j6) and ok3(j4, j5, j3):
        out["6j"].append(list(js) + sq(wigner_6j(*[H(x) for x in js])))
while len(out["6j"]) < 2500:
    js = [rnd.randint(0, 16) for _ in range(6)]
    j1, j2, j3, j4, j5, j6 = js
    if ok3(j1, j2, j3) and ok3(j1, j5, j6) and ok3(j4, j2, j6) and ok3(j4, j5, j3):
        out["6j"].append(js + sq(wigner_6j(*[H(x) for x in js])))
# 9j: a random sample with 2j up to 6, as decimals
while len(out["9j"]) < 150:
    js = [rnd.randint(0, 6) for _ in range(9)]
    j1, j2, j3, j4, j5, j6, j7, j8, j9 = js
    if ok3(j1, j2, j3) and ok3(j4, j5, j6) and ok3(j7, j8, j9) and ok3(j1, j4, j7) and ok3(j2, j5, j8) and ok3(j3, j6, j9):
        v = wigner_9j(*[H(x) for x in js], prec=60)
        out["9j"].append(js + [str(v.evalf(45))])
json.dump(out, open("test/fixtures/gt_wigner.json", "w"))
print({k: len(v) for k, v in out.items()})
