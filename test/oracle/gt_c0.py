"""The one-loop three-point function C0 by mpmath, independently of PQCF (plan "theoretical physics instruments", F3),
for test/fixtures/gt_c0.json. Every configuration by contour deformation (Soper, hep-ph/9910292; Binoth and Heinrich,
hep-ph/0004013): the triangle mapped to the unit square (x = s, y = s u), z_k = x_k - i kappa x_k (1 - x_k) dD/dx_k
with its exact complex Jacobian; kept only if kappa = 0.05 and 0.1 agree to 1e-12. Where D > 0 on the whole triangle
(no threshold inside) also the direct real integral at 30 digits, which must agree with the deformed one.
The configurations are read from argv[1] (JSON, exact fractions as strings); the results go to argv[2].
    docker run --rm -v "$PWD":/w -w /w pqcf-oracle python test/oracle/gt_c0.py in.json out.json
"""
import json
import sys
from fractions import Fraction

import mpmath as mp


def coeffs(p1, p2, p3, m1, m2, m3):
    return p1, p2, p3 - p1 - p2, m2 - m1 - p1, m3 - m2 + p1 - p3, m1


def deformed(args, kappa):
    a, b, c, d, e, f = coeffs(*args)
    D = lambda x, y: a*x*x + b*y*y + c*x*y + d*x + e*y + f

    def integrand(s, u):
        # F(s, u) = D(s, s u), its gradient and Hessian
        gs = 2*a*s + 2*b*s*u*u + 2*c*s*u + d + e*u
        gu = 2*b*s*s*u + c*s*s + e*s
        gss = 2*a + 2*b*u*u + 2*c*u
        gsu = 4*b*s*u + 2*c*s + e
        guu = 2*b*s*s
        zs = s - 1j*kappa*s*(1 - s)*gs
        zu = u - 1j*kappa*u*(1 - u)*gu
        j11 = 1 - 1j*kappa*((1 - 2*s)*gs + s*(1 - s)*gss)
        j12 = -1j*kappa*s*(1 - s)*gsu
        j21 = -1j*kappa*u*(1 - u)*gsu
        j22 = 1 - 1j*kappa*((1 - 2*u)*gu + u*(1 - u)*guu)
        return -zs*(j11*j22 - j12*j21)/D(zs, zs*zu)
    return mp.quad(integrand, [0, 1], [0, 1])


def positive(args):
    a, b, c, d, e, f = coeffs(*args)
    n = 60
    vals = [a*x*x + b*y*y + c*x*y + d*x + e*y + f
            for i in range(n + 1) for j in range(i + 1) for x, y in [(mp.mpf(i)/n, mp.mpf(j)/n)]]
    return min(vals) > mp.mpf('0.02')


def direct(args):
    a, b, c, d, e, f = coeffs(*args)
    return -mp.quad(lambda x: mp.quad(lambda y: 1/(a*x*x + b*y*y + c*x*y + d*x + e*y + f), [0, x]), [0, 1])


out = []
for cfg in json.load(open(sys.argv[1])):
    args = [mp.mpf(Fraction(v).numerator)/Fraction(v).denominator for v in cfg]
    mp.mp.dps = 24
    v1, v2 = deformed(args, mp.mpf('0.05')), deformed(args, mp.mpf('0.1'))
    stable = abs(v1 - v2) <= mp.mpf(10)**-12*max(1, abs(v2))
    rec = {"args": cfg, "deformed": [mp.nstr(v2.real, 20), mp.nstr(v2.imag, 20)], "stable": bool(stable)}
    if positive(args):
        mp.mp.dps = 34
        v = direct(args)
        assert abs(v - v2) <= mp.mpf(10)**-12*max(1, abs(v)), (cfg, v, v2)
        rec["direct"] = mp.nstr(v, 30)
    out.append(rec)
    print(cfg, rec["deformed"], "stable" if stable else "UNSTABLE", "direct" if "direct" in rec else "", flush=True)
json.dump(out, open(sys.argv[2], "w"), indent=0)
