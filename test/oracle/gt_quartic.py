"""The quartic oscillator H = p^2 + x^4: its first four eigenvalues by mpmath, independently of PQCF (plan "theoretical
physics instruments", QM2). On the half-line [0, B]: even states with psi'(0) = 0, odd ones with psi(0) = 0, and
psi(B) = 0 at B = 7 (the eigenfunctions fall like exp(-x^3/3): 1e-50 there). Shooting with mpmath's Taylor-series ODE
solver (odefun) and findroot, at 45 digits. Writes test/fixtures/gt_quartic.json.
    docker run --rm -v "$PWD":/src -w /src pqcf-oracle python test/oracle/gt_quartic.py
"""
import json

import mpmath as mp

mp.mp.dps = 45
B = 7


def end_value(E, even):
    y0 = [mp.mpf(1), mp.mpf(0)] if even else [mp.mpf(0), mp.mpf(1)]
    f = mp.odefun(lambda x, y: [y[1], (x ** 4 - E) * y[0]], 0, y0)
    return f(B)[0]


out = []
for n, guess in enumerate(["1.0604", "3.7997", "7.4557", "11.6447"]):
    even = n % 2 == 0
    # psi(B) is of size e^(B^3/3) away from an eigenvalue: bracket the root and judge convergence on E, not on psi(B)
    E = mp.findroot(lambda e: end_value(e, even), (mp.mpf(guess) - mp.mpf("0.001"), mp.mpf(guess) + mp.mpf("0.001")),
                    solver="illinois", tol=mp.mpf(10) ** -84, verify=False)
    out.append({"n": n, "even": even, "E": mp.nstr(E, 40)})
    print(n, mp.nstr(E, 40), flush=True)
json.dump({"B": B, "states": out}, open("test/fixtures/gt_quartic.json", "w"), indent=1)
