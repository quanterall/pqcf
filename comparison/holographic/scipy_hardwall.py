"""An independent check of HARDWALL.md in double precision: SciPy's DOP853, in the paper's coordinate z.

  docker build -t pqcf-oracle test/oracle
  docker run --rm -v "$PWD":/w -w /w pqcf-oracle python comparison/holographic/scipy_hardwall.py

Erlich, Katz, Son and Stephanov (hep-ph/0501128), with z_m = 1 and v(z) = mu z + s3 z^3 (mu = m_q z_m,
s3 = sigma z_m^3), g5^2 = 4 pi^2:
- a_1 (eq. 16): psi'' - psi'/z + (q^2 - g5^2 v^2/z^2) psi = 0, psi ~ z^2 at 0, psi'(1) = 0; m_a1^2 = q^2, the first.
- f_pi (eq. 20, q^2 = 0): A'' - A'/z - g5^2 (v/z)^2 A = 0, A'(1) = 0, and f_pi^2 = -A'(z)/(g5^2 z A(0)) at z -> 0.
  In (A, A'/z) the chiral limit (mu = 0) is regular at z = 0, so it is integrated from z = 1 down to 0 itself.
PQCF solves the same problems in u = z/z_m in Schrödinger form, by Frobenius series and Taylor-series integration
(hardwall.erl); nothing is shared.

Tolerance: m_a1 and the chiral-limit f_pi equal PQCF's (HARDWALL.md, 8 digits) to a relative 1e-7, the rounding of
the published values. Exits 1 otherwise. With m_q kept, f_pi depends on the cutoff eps (the O(m_q^2) contact term):
printed for eps = 1e-2 ... 1e-4, for the comparison with the paper; not checked.
"""
import math
import re
import sys

from scipy.integrate import solve_ivp
from scipy.optimize import brentq

G5SQ = 4 * math.pi ** 2
TOL = 1e-7
MODELS = {"A": (323.0, 2.29, 327.0), "B": (346.0, 2.3, 308.0)}   # 1/z_m, m_q, sigma^(1/3), in MeV


def params(name):
    zinv, mq, sig = MODELS[name]
    return zinv, mq / zinv, (sig / zinv) ** 3


def axial_end(q2, mu, s3, eps=1e-6):
    def f(z, y):
        w = mu + s3 * z * z                        # v/z
        return [y[1], y[1] / z - (q2 - G5SQ * w * w) * y[0]]
    s = solve_ivp(f, (eps, 1.0), [eps * eps, 2 * eps], method="DOP853", rtol=1e-13, atol=1e-30)
    return s.y[1, -1]                              # psi'(1)


def m_a1(name):
    zinv, mu, s3 = params(name)
    q2 = brentq(axial_end, 10.0, 22.0, args=(mu, s3), xtol=1e-14, rtol=1e-15)
    return zinv * math.sqrt(q2)


def f_pi(name, chiral=True, eps=1e-5):
    zinv, mu, s3 = params(name)
    if chiral:
        # in (A, B = A'/z) the chiral-limit equation is regular at 0: A' = z B, B' = g5^2 s3^2 z^3 A; integrated from
        # z = 1 to 0, f_pi^2 = -B(0) / (g5^2 A(0))
        def f(z, y):
            return [z * y[1], G5SQ * s3 * s3 * z ** 3 * y[0]]
        s = solve_ivp(f, (1.0, 0.0), [1.0, 0.0], method="DOP853", rtol=1e-14, atol=1e-30)
        return zinv * math.sqrt(-s.y[1, -1] / (s.y[0, -1] * G5SQ))

    # with m_q, B' = g5^2 (v/z)^2 A / z diverges at 0 (the logarithm of the contact term): stopped at eps
    def f(z, y):
        w = mu + s3 * z * z
        return [z * y[1], G5SQ * w * w * y[0] / z]
    s = solve_ivp(f, (1.0, eps), [1.0, 0.0], method="DOP853", rtol=1e-13, atol=1e-30)
    return zinv * math.sqrt(-s.y[1, -1] / (s.y[0, -1] * G5SQ))


def pqcf(name, label):
    text = open("comparison/holographic/HARDWALL.md").read()
    m = re.search(r"^\| %s \| %s \| `([^`]+)`" % (name, re.escape(label)), text, re.M)
    return float(m.group(1))


def main():
    ok = True
    print(f"{'':10} {'observable':22} {'PQCF':>12} {'SciPy':>14} {'relative':>9}")
    for name in MODELS:
        for label, value in (("m_a1", m_a1(name)), ("f_pi (chiral limit)", f_pi(name))):
            ref = pqcf(name, label)
            rel = abs(value - ref) / ref
            ok &= rel < TOL
            print(f"model {name:4} {label:22} {ref:12.7g} {value:14.9g} {rel:9.1e}")
    for name in MODELS:
        vals = ", ".join(f"{f_pi(name, chiral=False, eps=e):.2f}" for e in (1e-2, 1e-3, 1e-4))
        print(f"model {name}: f_pi with m_q kept, eps = 1e-2, 1e-3, 1e-4 (in units of z_m): {vals} MeV")
    print(f"agree to {TOL:g}" if ok else "DIFFER")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
