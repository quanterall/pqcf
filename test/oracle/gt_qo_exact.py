"""Exact steady states and spectra of small master equations, by SymPy (quantum optics plan, Q3): the reference for
pq_qo_exact. The models are dense matrices in the Fock basis (square roots kept exact), the Lindblad right-hand side
written as matrix products, rho solved exactly; independent of pq_qop, of the superoperator and of the rescaled basis
f_n = sqrt(n!) |n> that pq_qo_exact uses. Writes test/fixtures/gt_qo_exact.json: every value as [Re, Im] strings.
    docker run --rm -v "$PWD":/src -w /src pqcf-oracle python test/oracle/gt_qo_exact.py
"""
import json

import sympy as sp

I = sp.I


def destroy(n):
    return sp.Matrix(n + 1, n + 1, lambda i, j: sp.sqrt(j) if i == j - 1 else 0)


def kron(*ms):
    out = ms[0]
    for m in ms[1:]:
        out = sp.kronecker_product(out, m)
    return out


def steady(h, ls):
    """rho with -i[H, rho] + sum (L rho L^+ - 1/2 {L^+ L, rho}) = 0 and Tr rho = 1, exactly"""
    d = h.shape[0]
    xs = sp.symbols(f"r0:{d * d}")
    rho = sp.Matrix(d, d, xs)
    rhs = -I * (h * rho - rho * h)
    for l in ls:
        k = l.H * l
        rhs += l * rho * l.H - (k * rho + rho * k) / 2
    eqs = list(rhs)[1:] + [rho.trace() - 1]
    sol = sp.solve(eqs, xs, dict=True)[0]
    return rho.subs(sol)


def liouvillian(h, ls):
    """the superoperator on vec(rho) (row by row), from the matrices"""
    d = h.shape[0]
    cols = []
    for k in range(d * d):
        e = sp.zeros(d, d)
        e[k // d, k % d] = 1
        out = -I * (h * e - e * h)
        for l in ls:
            kk = l.H * l
            out += l * e * l.H - (kk * e + e * kk) / 2
        cols.append(list(out))
    return sp.Matrix(cols).T


def spectrum(h, ls, a, b, rho, ws):
    """Re Tr[A (-(L + i w))^-1 (B rho - <B> rho)] at each w"""
    d = h.shape[0]
    lv = liouvillian(h, ls)
    x = b * rho - (b * rho).trace() * rho
    xv = sp.Matrix(list(x))
    out = []
    for w in ws:
        y = (lv + I * w * sp.eye(d * d)).LUsolve(-xv)
        ym = sp.Matrix(d, d, list(y))
        out.append(num((a * ym).trace())[0])
    return out


def num(z):
    """[Re, Im] as exact rationals; anything else is an error (the steady states of these models are rational)"""
    z = sp.radsimp(sp.expand(sp.simplify(z)))
    re, im = sp.simplify(sp.re(z)), sp.simplify(sp.im(z))
    if not (re.is_Rational and im.is_Rational):
        raise ValueError(f"not rational: {z}")
    return [str(re), str(im)]


def two_level(omega, delta, gamma):
    sm = sp.Matrix([[0, 1], [0, 0]])          # |1><2|, level 1 the ground state (index 0)
    h = omega / 2 * (sm + sm.T) - delta * sm.T * sm
    return sm, h, [sp.sqrt(gamma) * sm]


def main():
    cases = []
    r = sp.Rational
    # resonance fluorescence: steady state and the incoherent spectrum
    for om, dl, g in [(1, 0, 1), (3, 0, 1), (2, 1, 1), (r(1, 2), 2, 1), (20, 0, 1)]:
        om, dl, g = r(om), r(dl), r(g)
        sm, h, ls = two_level(om, dl, g)
        rho = steady(h, ls)
        ws = [r(1, 3), 1, r(5, 2), 4, 20]
        spec = spectrum(h, ls, sm.T, sm, rho, ws)
        cases.append({"model": "two_level", "params": {"omega": str(om), "delta": str(dl), "gamma": str(g)},
                      "pe": num(rho[1, 1]), "sm": num((sm * rho).trace()),
                      "spectrum": [[str(w), s] for w, s in zip(ws, spec)]})
    # a driven Kerr cavity with loss, truncated at N photons
    for n, e, chi, kappa, det in [(3, 1, 1, 1, 0), (4, r(1, 2), r(1, 3), 2, 1)]:
        e, chi, kappa, det = r(e), r(chi), r(kappa), r(det)
        a = destroy(n)
        ad = a.T
        h = det * ad * a + chi * ad * ad * a * a + e * (a + ad)
        rho = steady(h, [sp.sqrt(kappa) * a])
        cases.append({"model": "kerr", "params": {"n": n, "e": str(e), "chi": str(chi), "kappa": str(kappa),
                                                  "delta": str(det)},
                      "n": num((ad * a * rho).trace()), "a": num((a * rho).trace()),
                      "n2": num((ad * ad * a * a * rho).trace())})
    # the Jaynes-Cummings model with cavity drive and both losses, the cavity truncated at 2 photons
    n = 2
    a = kron(destroy(n), sp.eye(2))
    sm = kron(sp.eye(n + 1), sp.Matrix([[0, 1], [0, 0]]))
    g, e, kappa, gamma = r(1), r(1, 2), r(1, 2), r(1, 4)
    h = g * (a.T * sm + a * sm.T) + e * (a + a.T)
    rho = steady(h, [sp.sqrt(kappa) * a, sp.sqrt(gamma) * sm])
    cases.append({"model": "jc", "params": {"n": n, "g": str(g), "e": str(e), "kappa": str(kappa), "gamma": str(gamma)},
                  "n": num((a.T * a * rho).trace()), "pe": num((sm.T * sm * rho).trace()), "a": num((a * rho).trace())})
    json.dump({"cases": cases}, open("test/fixtures/gt_qo_exact.json", "w"), indent=1)


if __name__ == "__main__":
    main()
