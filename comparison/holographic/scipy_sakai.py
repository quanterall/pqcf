"""An independent check of SAKAI.md in double precision: SciPy's DOP853, in the original coordinate z.

  docker build -t pqcf-oracle test/oracle
  docker run --rm -v "$PWD":/w -w /w pqcf-oracle python comparison/holographic/scipy_sakai.py

Sakai and Sugimoto (hep-th/0412141), eq. (3.28): -K^(1/3) d_z (K d_z psi) = lambda psi, K = 1 + z^2. The vectors are
even in z (psi'(0) = 0), the axial vectors odd (psi(0) = 0). At large z the solutions go as 1 and 1/z; the
normalisable one decays, so at a cutoff Z we impose psi + z psi' = 0, the condition 1/z satisfies, and extrapolate the
eigenvalue in Z: Richardson over Z = 200 ... 12800, in the powers Z^(-5/3), Z^(-7/3), ... of the cutoff's error. It
shares nothing with PQCF's method (Frobenius series in u = K^(1/3), Taylor-series integration, a Wronskian match).

Tolerance: the first four eigenvalues equal PQCF's (SAKAI.md) to a relative 1e-9; they agree to 3e-11 or better (the
extrapolation's accuracy, which falls with n). Exits 1 otherwise.
"""
import re
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

TOL = 1e-9


def mismatch(lam, parity, Z):
    # y = (psi, K psi'): psi' = y1 / K, (K psi')' = -lambda psi / K^(1/3)
    def f(z, y):
        k = 1 + z * z
        return [y[1] / k, -lam * y[0] / k ** (1 / 3)]
    y0 = [1.0, 0.0] if parity == "even" else [0.0, 1.0]
    s = solve_ivp(f, (0, Z), y0, method="DOP853", rtol=1e-13, atol=1e-15)
    psi, kdpsi = s.y[0, -1], s.y[1, -1]
    return psi + Z * kdpsi / (1 + Z * Z)


def eigenvalue(guess, parity, Z):
    a, b = guess * 0.9, guess * 1.1
    return brentq(mismatch, a, b, args=(parity, Z), xtol=1e-15, rtol=1e-15)


def extrapolated(guess, parity):
    zs = [200.0 * 2 ** k for k in range(7)]
    lams = [eigenvalue(guess, parity, z) for z in zs]
    # Richardson: the decaying solution is 1/z times a series in z^(-2/3), so the cutoff's error is a series in
    # Z^(-5/3), Z^(-7/3), Z^(-3), ...; each pass removes one term (Z doubles from one value to the next)
    table = [lams]
    for k in range(1, len(zs)):
        prev, r = table[-1], 2 ** (1 + 2 * k / 3)
        table.append([(r * prev[i + 1] - prev[i]) / (r - 1) for i in range(len(prev) - 1)])
    return table[-1][0], lams[-1]


def pqcf_values():
    text = open("comparison/holographic/SAKAI.md").read()
    return [(m.group(1), float(m.group(2))) for m in
            re.finditer(r"^\| \d+ \| (vector|axial) \| (?:even|odd) \| `([^`]+)`", text, re.M)]


def main():
    ok = True
    print(f"{'n':>2} {'meson':7} {'PQCF':>20} {'SciPy':>20} {'relative':>10}")
    for n, (kind, ref) in enumerate(pqcf_values()[:4], 1):
        lam, _ = extrapolated(ref, "even" if kind == "vector" else "odd")
        rel = abs(lam - ref) / ref
        ok &= rel < TOL
        print(f"{n:>2} {kind:7} {ref:20.12f} {lam:20.12f} {rel:10.1e}")
    print(f"agree to {TOL:g}" if ok else "DIFFER")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
