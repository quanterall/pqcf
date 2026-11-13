"""G2: coherent-Ising-machine (CIM) solvers, simulated on one core with Numba, for Ising problems
E(s) = -sum_{i<j} J_ij s_i s_j - sum_i h_i s_i, s in {-1, +1}^n.

- SimCIM (Tiunov, Ulanov, Lvovsky, arXiv:1901.08927, Eqs. 5-6): x_i <- clip(x_i + dx_i, -1, 1) with
  dx = v(t) x + zeta J x + noise, the pump-loss v(t) ramped by a tanh law, and momentum 0.9.
- AHC, the CIM with amplitude-heterogeneity correction (Leleu, Yamamoto, McMahon, arXiv:1810.12565, Eqs. 1-2):
  dx_i/dt = (p - 1) x_i - x_i^3 + eps e_i sum_j J_ij x_j,  de_i/dt = -beta (x_i^2 - a) e_i, with a constant target
  amplitude a (the paper also modulates a; not done here), integrated by Euler steps.

Local fields h are handled with one extra spin s_0 coupled to s_i by h_i; the solution is read with s_0 = +1.
Both report the lowest energy among the sign patterns visited (checked every `every` iterations).
Couplings are stored sparse (CSR): haplotype instances have about 10 couplings per spin.
"""
import numpy as np
from numba import njit


def with_ancilla(n, J, h):
    """(J, h) on n spins -> symmetric CSR couplings on n + 1 spins, spin n the ancilla."""
    rows = [[] for _ in range(n + 1)]
    for (i, j), w in J.items():
        rows[i].append((j, w))
        rows[j].append((i, w))
    for i, w in (h or {}).items():
        if w:
            rows[i].append((n, w))
            rows[n].append((i, w))
    ptr = np.zeros(n + 2, np.int64)
    for i, r in enumerate(rows):
        ptr[i + 1] = ptr[i] + len(r)
    idx = np.array([j for r in rows for j, _ in r], np.int64)
    val = np.array([w for r in rows for _, w in r], np.float64)
    return ptr, idx, val


def spectral_norm(ptr, idx, val, iters=200):
    n = len(ptr) - 1
    x = np.random.default_rng(0).standard_normal(n)
    lam = 1.0
    for _ in range(iters):
        y = _matvec(ptr, idx, val, x)
        lam = np.linalg.norm(y)
        x = y / max(lam, 1e-300)
    return lam


@njit(cache=True)
def _matvec(ptr, idx, val, x):
    n = len(ptr) - 1
    y = np.zeros(n)
    for i in range(n):
        s = 0.0
        for k in range(ptr[i], ptr[i + 1]):
            s += val[k] * x[idx[k]]
        y[i] = s
    return y


@njit(cache=True)
def _energy_signs(ptr, idx, val, x):
    # E = -1/2 sum_ij J_ij s_i s_j over the symmetric CSR, with the ancilla's sign folded in (s_i * s_anc)
    n = len(ptr) - 1
    anc = 1.0 if x[n - 1] >= 0 else -1.0
    e = 0.0
    for i in range(n):
        si = (1.0 if x[i] >= 0 else -1.0) * anc
        if i == n - 1:
            si = 1.0
        for k in range(ptr[i], ptr[i + 1]):
            j = idx[k]
            sj = (1.0 if x[j] >= 0 else -1.0) * anc
            if j == n - 1:
                sj = 1.0
            e -= 0.5 * val[k] * si * sj
    return e


@njit(cache=True)
def _simcim(ptr, idx, val, iters, amp, zeta, sigma, every, seed):
    np.random.seed(seed)
    n = len(ptr) - 1
    x = np.zeros(n)
    m = np.zeros(n)
    best = 1e300
    bestx = np.zeros(n)
    for t in range(iters):
        v = amp * np.tanh(3.0 * (2.0 * t / max(iters - 1, 1) - 1.0))
        f = _matvec(ptr, idx, val, x)
        for i in range(n):
            dx = v * x[i] + zeta * f[i] + sigma * np.random.standard_normal()
            m[i] = 0.9 * m[i] + 0.1 * dx
            xi = x[i] + m[i]
            x[i] = 1.0 if xi > 1.0 else (-1.0 if xi < -1.0 else xi)
        if (t + 1) % every == 0 or t == iters - 1:
            e = _energy_signs(ptr, idx, val, x)
            if e < best:
                best = e
                bestx[:] = x
    return best, bestx


@njit(cache=True)
def _ahc(ptr, idx, val, iters, p, eps, beta, a, dt, every, seed):
    np.random.seed(seed)
    n = len(ptr) - 1
    x = 0.01 * np.random.standard_normal(n)
    e = np.ones(n)
    best = 1e300
    bestx = np.zeros(n)
    for t in range(iters):
        f = _matvec(ptr, idx, val, x)
        for i in range(n):
            xi = x[i]
            x[i] = xi + dt * ((p - 1.0) * xi - xi * xi * xi + eps * e[i] * f[i])
            e[i] = e[i] + dt * (-beta * (xi * xi - a) * e[i])
            if e[i] > 1e6:
                e[i] = 1e6
            elif e[i] < 1e-6:
                e[i] = 1e-6
        if (t + 1) % every == 0 or t == iters - 1:
            en = _energy_signs(ptr, idx, val, x)
            if en < best:
                best = en
                bestx[:] = x
    return best, bestx


def spins(x):
    """The solution's spins on the original n variables, the ancilla's sign folded in."""
    s = np.where(x >= 0, 1, -1)
    return s[:-1] * s[-1]


class Ising:
    def __init__(self, n, J, h=None):
        self.n = n
        self.ptr, self.idx, self.val = with_ancilla(n, J, h)
        self.norm = spectral_norm(self.ptr, self.idx, self.val)

    def simcim(self, iters, seed, amp=0.2, zeta=0.5, sigma=0.05, every=10):
        e, x = _simcim(self.ptr, self.idx, self.val, iters, amp, zeta / self.norm, sigma, every, seed)
        return e, spins(x)

    def ahc(self, iters, seed, p=0.5, eps=0.5, beta=0.2, a=1.0, dt=0.1, every=10):
        e, x = _ahc(self.ptr, self.idx, self.val, iters, p, eps / self.norm, beta, a, dt, every, seed)
        return e, spins(x)
