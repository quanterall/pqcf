"""An independent check of GUBSER.md in double precision: SciPy's DOP853 on the same equations.

  docker build -t pqcf-oracle test/oracle
  docker run --rm -v "$PWD":/w -w /w pqcf-oracle python comparison/holographic/scipy_gubser.py

Gubser and Nellore (arXiv:0804.0434): the master equation (34) for G = dA/dphi, from the horizon series (36) at
phi_H - 1e-8 down to phi = 1e-6, with the integrals of (45) and (51) carried along; c_s^2 = d log T / d log s by
central differences in phi_H (h = 1e-3); V = -12 cosh(gamma phi) + b phi^2, gamma = 0.606, b = 2.06. The same
steps as gubser.erl, in floating point instead of PQCF's arbitrary-precision Taylor-series method (pq_jet).

Tolerance: c_s^2 equal to PQCF's (GUBSER.md, 6 decimals) to 2e-6 for phi_H =< 5; they agree to 5e-7, the table's
rounding. Beyond that a double-precision solution depends on its tolerance (the 0/0 at the horizon): its values at
phi_H = 6 and 7 are printed for two tolerances and not checked. Exits 1 otherwise.
"""
import math
import re
import sys

from scipy.integrate import solve_ivp

GM, B = 0.606, 2.06
DM4 = -2 + math.sqrt(4 + 2 * (B - 6 * GM * GM))       # Delta - 4, eq. (63)
TOL = 2e-6


def pot(phi):
    v = -12 * math.cosh(GM * phi) + B * phi * phi
    v1 = -12 * GM * math.sinh(GM * phi) + 2 * B * phi
    v2 = -12 * GM * GM * math.cosh(GM * phi) + 2 * B
    return v, v1, v2


def rhs(phi, y):
    g, gp, _, _ = y
    v, v1, v2 = pot(phi)
    r = v / (3 * v1)
    rp = 1 / 3 - v * v2 / (3 * v1 * v1)
    p = g + r
    a = gp / g + 1 / (6 * g) - 4 * g - gp / p
    num = a * gp / p + gp * gp / (g * g) + gp / (6 * g * g) + 4 * gp - gp * (gp + rp) / (p * p)
    return [gp, g * p / r * num, g - 1 / (DM4 * phi), 1 / (6 * g)]


def thermo(phi_h, rtol):
    v, v1, v2 = pot(phi_h)
    g0 = -v / (3 * v1)
    g1 = (v * v2 / (v1 * v1) - 1) / 6
    d = 1e-8
    s = solve_ivp(rhs, (phi_h - d, 1e-6), [g0 - d * g1, g1, 0.0, 0.0], method="DOP853", rtol=rtol, atol=1e-14)
    i1, i2 = s.y[2, -1], s.y[3, -1]
    j1 = -i1 + d * (g0 - 1 / (DM4 * phi_h))
    j2 = -i2 + d / (6 * g0)
    log_s = 3 / DM4 * math.log(phi_h) + 3 * j1
    log_t = 1 / DM4 * math.log(phi_h) + math.log(v / pot(0)[0]) + j1 + j2
    return log_s, log_t


def cs2(phi_h, rtol=1e-12, h=1e-3):
    s1, t1 = thermo(phi_h - h, rtol)
    s2, t2 = thermo(phi_h + h, rtol)
    return (t2 - t1) / (s2 - s1)


def pqcf():
    text = open("comparison/holographic/GUBSER.md").read()
    return {float(m.group(1)): float(m.group(2))
            for m in re.finditer(r"^\| ([0-9.]+) \| [0-9.]+ \| ([0-9.]+) \|", text, re.M)}


def main():
    ok = True
    print(f"{'phi_H':>6} {'PQCF':>10} {'SciPy':>12} {'difference':>11}")
    ref = pqcf()
    for phi_h, c_ref in sorted(ref.items()):
        if phi_h > 5:
            continue
        c = cs2(phi_h)
        ok &= abs(c - c_ref) < TOL
        print(f"{phi_h:6.2f} {c_ref:10.6f} {c:12.8f} {c - c_ref:11.1e}")
    for phi_h in (6.0, 7.0):
        if phi_h in ref:
            print(f"{phi_h:6.2f} {ref[phi_h]:10.6f}  SciPy, rtol 1e-10 and 1e-12: {cs2(phi_h, 1e-10):.6f}, "
                  f"{cs2(phi_h, 1e-12):.6f} (not checked)")
    print(f"agree to {TOL:g} for phi_H =< 5" if ok else "DIFFER")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
