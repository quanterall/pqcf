"""QCD's renormalisation group by mpmath, independently of PQCF (plan "theoretical physics instruments", QCD1).
The beta and gamma_m coefficients transcribed again from the papers (Baikov, Chetyrkin and Kuehn, arXiv:1606.08659 and
arXiv:1402.6611; van Ritbergen, Vermaseren and Larin, hep-ph/9701390; Chetyrkin, hep-ph/9703278) as Fractions with
zeta values; alpha_s and m by mpmath's Taylor-series ODE solver at 45 digits (da/dt = -beta(a), dln m/dt = -gamma(a),
t = ln mu^2); the decoupling at an MS-bar threshold to three loops (Chetyrkin, Kniehl and Steinhauser, hep-ph/9706430)
summed directly. Reads the cases from argv[1], writes argv[2].
    docker run --rm -v "$PWD":/w -w /w pqcf-oracle python test/oracle/gt_qcd.py in.json out.json
"""
import json
import sys
from fractions import Fraction

import mpmath as mp

mp.mp.dps = 45


def F(a, b=1):
    return mp.mpf(a)/b

z = {k: mp.zeta(k) for k in range(3, 8)}


def beta(nf):
    n = F(nf)
    return [
        (F(11) - F(2, 3)*n)/4,
        (F(102) - F(38, 3)*n)/16,
        (F(2857, 2) - F(5033, 18)*n + F(325, 54)*n**2)/64,
        (F(149753, 6) + 3564*z[3] + (-F(1078361, 162) - F(6508, 27)*z[3])*n + (F(50065, 162) + F(6472, 81)*z[3])*n**2
         + F(1093, 729)*n**3)/256,
        (F(8157455, 16) + F(621885, 2)*z[3] - F(88209, 2)*z[4] - 288090*z[5]
         + (-F(336460813, 1944) - F(4811164, 81)*z[3] + F(33935, 6)*z[4] + F(1358995, 27)*z[5])*n
         + (F(25960913, 1944) + F(698531, 81)*z[3] - F(10526, 9)*z[4] - F(381760, 81)*z[5])*n**2
         + (-F(630559, 5832) - F(48722, 243)*z[3] + F(1618, 27)*z[4] + F(460, 9)*z[5])*n**3
         + (F(1205, 2916) - F(152, 81)*z[3])*n**4)/1024]


def gamma(nf):
    n = F(nf)
    return [
        mp.mpf(1),
        (F(202, 3) - F(20, 9)*n)/16,
        (1249 + (-F(2216, 27) - F(160, 3)*z[3])*n - F(140, 81)*n**2)/64,
        (F(4603055, 162) + F(135680, 27)*z[3] - 8800*z[5]
         + (-F(91723, 27) - F(34192, 9)*z[3] + 880*z[4] + F(18400, 9)*z[5])*n
         + (F(5242, 243) + F(800, 9)*z[3] - F(160, 3)*z[4])*n**2 + (-F(332, 243) + F(64, 27)*z[3])*n**3)/256,
        (F(99512327, 162) + F(46402466, 243)*z[3] + 96800*z[3]**2 - F(698126, 9)*z[4] - F(231757160, 243)*z[5]
         + 242000*z[6] + 412720*z[7]
         + n*(-F(150736283, 1458) - F(12538016, 81)*z[3] - F(75680, 9)*z[3]**2 + F(2038742, 27)*z[4]
              + F(49876180, 243)*z[5] - F(638000, 9)*z[6] - F(1820000, 27)*z[7])
         + n**2*(F(1320742, 729) + F(2010824, 243)*z[3] + F(46400, 27)*z[3]**2 - F(166300, 27)*z[4]
                 - F(264040, 81)*z[5] + F(92000, 27)*z[6])
         + n**3*(F(91865, 1458) + F(12848, 81)*z[3] + F(448, 9)*z[4] - F(5120, 27)*z[5])
         + n**4*(-F(260, 243) - F(320, 243)*z[3] + F(64, 27)*z[4]))/1024]


def mpf(x):
    x = Fraction(x)
    return mp.mpf(x.numerator)/x.denominator


def run(as0, mu0, mu1, nf, loops, m0):
    B = beta(nf)[:loops]
    G = gamma(nf)[:loops]
    t0, t1 = mp.log(mpf(mu0)**2), mp.log(mpf(mu1)**2)
    # odefun integrates forward only: downward in mu, in -t
    s = 1 if t1 >= t0 else -1
    f = mp.odefun(lambda t, y: [-s*sum(b*y[0]**(i + 2) for i, b in enumerate(B)),
                                -s*sum(g*y[0]**(i + 1) for i, g in enumerate(G))], s*t0, [mpf(as0)/mp.pi, mp.mpf(0)])
    a1, lnm = f(s*t1)
    return a1*mp.pi, mpf(m0)*mp.exp(lnm)


def dec(direction, als, mh, mu, nl, loops):
    a = mpf(als)/mp.pi
    L = mp.log(mpf(mu)**2/mpf(mh)**2)
    c1 = -L/6
    c2 = mp.mpf(11)/72 - 11*L/24 + L**2/36
    c3 = (mp.mpf(564731)/124416 - 955*L/576 + 53*L**2/576 - L**3/216
          + nl*(-mp.mpf(2633)/31104 + 67*L/576 - L**2/36) - mp.mpf(82043)/27648*z[3])
    cs = [c1, c2, c3] if direction == "down" else [-c1, 2*c1**2 - c2, -5*c1**3 + 5*c1*c2 - c3]
    return mpf(als)*(1 + sum(c*a**(k + 1) for k, c in enumerate(cs[:loops - 1])))


cases = json.load(open(sys.argv[1]))
out = {"coefficients": [], "running": [], "decoupling": []}
for nf in range(3, 7):
    out["coefficients"].append({"nf": nf, "beta": [mp.nstr(b, 40) for b in beta(nf)],
                                "gamma_m": [mp.nstr(g, 40) for g in gamma(nf)]})
for c in cases["running"]:
    a1, m1 = run(c["as0"], c["mu0"], c["mu1"], c["nf"], c["loops"], c["m0"])
    out["running"].append(dict(c, alphas=mp.nstr(a1, 30), mass=mp.nstr(m1, 30)))
    print("running", c, mp.nstr(a1, 20), mp.nstr(m1, 20), flush=True)
for c in cases["decoupling"]:
    v = dec(c["direction"], c["as"], c["mh"], c["mu"], c["nl"], c["loops"])
    out["decoupling"].append(dict(c, value=mp.nstr(v, 30)))
json.dump(out, open(sys.argv[2], "w"), indent=0)
