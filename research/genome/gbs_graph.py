"""G1: Gaussian boson sampling (GBS) for dense subgraphs and cliques, and the classical samplers it is compared with.

- encode(): a graph's adjacency matrix A as a GBS device, by the Takagi decomposition A = U diag(d) U^T: squeezing
  tanh r_i = c d_i into an interferometer U, so the output's GBS matrix is B = c A and a collision-free photon
  pattern S has probability proportional to |haf(A_S)|^2 (Arrazola & Bromley, arXiv:1803.10730). c sets the mean
  photon number. The device itself is simulated by PQCF (g1_pqcf.erl).
- oh_sampler(): the quantum-inspired classical sampler of Oh, Fefferman, Jiang, Quesada (arXiv:2302.00536): for
  a nonnegative A, draw k/2 edges with probability proportional to their weight and keep the draw if they are
  disjoint, so a k-vertex set S appears with probability proportional to haf(A_S) (not its square).
- uniform k-subsets.
- random search for the densest k-subgraph, and shrink-and-search for the maximum clique (Banchi et al.,
  arXiv:1902.00462), the same procedures for every sampler.
"""
import math
from itertools import combinations

import numpy as np


# --- the GBS encoding -------------------------------------------------------------------------------------------------

def takagi_real_symmetric(A):
    """A = U diag(d) U^T with d >= 0 and U unitary, for a real symmetric A: from A = V diag(e) V^T, d = |e| and the
    columns of V with a negative eigenvalue multiplied by i (i^2 = -1 restores the sign)."""
    e, V = np.linalg.eigh(A)
    U = V.astype(complex) * np.where(e < 0, 1j, 1.0)
    return np.abs(e), U


def mean_photons(d, c):
    t = c * d
    return float(np.sum(t ** 2 / (1 - t ** 2)))


def total_photon_distribution(d, c, K):
    """P(total photons = 0..K) of the encoded state: the interferometer conserves photon number, so it is the
    convolution of the single-mode squeezed vacua, P(2m) = C(2m, m) / 4^m tanh^2m r / cosh r."""
    p = np.zeros(K + 1)
    p[0] = 1.0
    for t in c * d:
        s = np.zeros(K + 1)
        for m in range(K // 2 + 1):
            s[2 * m] = math.comb(2 * m, m) / 4 ** m * t ** (2 * m) * math.sqrt(1 - t * t)
        p = np.convolve(p, s)[:K + 1]
    return p


def sampling_point(A, k, tail=1e-5):
    """The mean photon number and photon cap to simulate the device at: PQCF's cap truncates step by step, which is
    exact only when little probability lies above it, so among caps k+2..k+6 and mean photon numbers whose mass
    above the cap is below `tail`, the one with the largest P(total = k)."""
    d, _ = takagi_real_symmetric(np.asarray(A, float))
    best = None
    for mean in np.geomspace(0.05, 2 * k, 120):
        c = encode_scale(d, mean)
        for K in (k + 2, k + 4, k + 6):
            p = total_photon_distribution(d, c, K + 40)
            above = max(0.0, 1.0 - p[:K + 1].sum())
            if above < tail and (best is None or p[k] > best["p_k"]):
                best = {"mean_photons": float(mean), "cap": K, "p_k": float(p[k]), "mass_above_cap": float(above)}
    return best


def encode_scale(d, mean):
    lo, hi = 0.0, (1 - 1e-12) / d.max()
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if mean_photons(d, mid) < mean else (lo, mid)
    return (lo + hi) / 2


def encode(A, mean):
    """Squeezing parameters and interferometer for B = c A with the given mean photon number."""
    d, U = takagi_real_symmetric(np.asarray(A, float))
    c = encode_scale(d, mean)
    r = np.arctanh(c * d)
    assert np.allclose(U @ np.diag(c * d) @ U.T, c * np.asarray(A, float), atol=1e-9)
    return {"c": c, "r": r.tolist(), "U": U, "mean_photons": mean_photons(d, c)}


# --- hafnians and the exact distributions (small graphs, for checks) ---------------------------------------------------

def hafnian(M):
    """By expansion along the first row; fine for the 4-10 vertex subgraphs of the checks."""
    n = len(M)
    if n == 0:
        return 1.0
    if n % 2:
        return 0.0
    total = 0.0
    rest = list(range(1, n))
    for j in rest:
        if M[0][j]:
            keep = [x for x in rest if x != j]
            total += M[0][j] * hafnian([[M[a][b] for b in keep] for a in keep])
    return total


def exact_distribution(A, k, power):
    """P(S) proportional to haf(A_S)^power over the k-subsets S: power 2 for GBS post-selected on k collision-free
    photons, 1 for Oh et al.'s sampler."""
    n = len(A)
    subsets = list(combinations(range(n), k))
    w = np.array([hafnian([[A[i][j] for j in S] for i in S]) ** power for S in subsets])
    return subsets, w / w.sum()


# --- the classical samplers -------------------------------------------------------------------------------------------

def oh_sampler(A, k, n_samples, rng, max_tries=10_000_000):
    """k/2 edges drawn with probability proportional to A_ij, kept when pairwise disjoint (a matching): the vertex
    set appears with probability proportional to haf(A_S). Returns the samples and the number of draws."""
    A = np.asarray(A, float)
    iu = np.triu_indices(len(A), 1)
    w = A[iu]
    edges = np.stack(iu, axis=1)[w > 0]
    p = w[w > 0] / w[w > 0].sum()
    m = k // 2
    out, draws = [], 0
    while len(out) < n_samples and draws < max_tries:
        batch = rng.choice(len(edges), size=(4096, m), p=p)
        draws += len(batch)
        for row in batch:
            verts = edges[row].ravel()
            if len(set(verts.tolist())) == k:
                out.append(tuple(sorted(verts.tolist())))
                if len(out) == n_samples:
                    break
    return out, draws


def uniform_sampler(n, k, n_samples, rng):
    return [tuple(sorted(rng.choice(n, size=k, replace=False).tolist())) for _ in range(n_samples)]


# --- the searches -----------------------------------------------------------------------------------------------------

def edges_in(A, S):
    return int(sum(A[i][j] > 0 for a, i in enumerate(S) for j in S[a + 1:]))


def random_search_curve(A, samples, checkpoints):
    """The best number of edges among the first s samples, for each checkpoint s."""
    best, out, it = -1, [], iter(samples)
    done = 0
    for s in checkpoints:
        while done < s:
            S = next(it)
            best = max(best, edges_in(A, S))
            done += 1
        out.append(best)
    return out


def shrink_and_search(A, S, rng, steps=None):
    """Banchi et al.'s post-processing: remove the vertex of lowest degree within S until S is a clique (ties at
    random), then grow it greedily by vertices adjacent to all of it (ties at random)."""
    A = np.asarray(A) > 0
    S = list(S)
    while True:
        sub = A[np.ix_(S, S)]
        deg = sub.sum(axis=1)
        if deg.min() == len(S) - 1:
            break
        low = np.flatnonzero(deg == deg.min())
        S.pop(int(rng.choice(low)))
    while True:
        cand = [v for v in range(len(A)) if v not in S and A[v, S].all()]
        if not cand:
            return tuple(sorted(S))
        S.append(int(rng.choice(cand)))


# --- instances --------------------------------------------------------------------------------------------------------

def planted_graph(n=30, k=10, p_in=0.95, p_out=0.4, seed=0):
    """A random graph with a planted dense k-subgraph whose vertices have lower degree than the rest, the kind of
    graph Arrazola & Bromley use (their Fig. 3): degree-based greedy methods miss it. Planted vertices connect
    inside with p_in and to the rest only sparsely."""
    rng = np.random.default_rng(seed)
    A = np.zeros((n, n), int)
    planted = set(range(k))
    for i in range(n):
        for j in range(i + 1, n):
            if i in planted and j in planted:
                p = p_in
            elif i in planted or j in planted:
                p = 0.05
            else:
                p = p_out
            A[i, j] = A[j, i] = int(rng.random() < p)
    perm = rng.permutation(n)
    A = A[np.ix_(perm, perm)]
    planted_at = sorted(int(np.flatnonzero(perm == v)[0]) for v in planted)
    return A, planted_at
