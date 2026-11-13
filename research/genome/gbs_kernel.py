"""G3: the Gaussian-boson-sampling graph kernel (Schuld, Brádler, Israel, Su, Gupt, arXiv:1905.12646) on graphs made
from DNA, against classical features, for telling which organism a stretch of DNA comes from.

- A window of DNA as a graph: the 16 two-letter words are the nodes, and the weight of the edge {xy, yz} is the
  number of times the three-letter word xyz (or its reverse zyx) occurs, divided by the window's length; no
  self-loops.
- The GBS feature vector: for a device programmed with the graph's adjacency matrix A, the probabilities of the
  photon-number "orbits" (the patterns up to a permutation of the modes, e.g. [2, 1, 1]) given N photons in total,
  N = 4, 6, 8. Given N, p(n) is proportional to |haf(A_n)|^2 / prod n_i!, independently of the device's squeezing
  scale, so the features are computed here exactly, summing over every pattern with N photons (`orbit_features`).
  PQCF's simulated device is used to check them (g3.py).
- Classical features: the closed-walk counts Tr(A^k), k = 2..8 (the classical analogue of the matching counts GBS
  measures), the sorted spectrum of A, and, without any graph, the 64 three-letter-word frequencies (the standard
  alignment-free genome classifier).
"""
import itertools
import math
import random

import numpy as np
from numba import njit

BASES = "ACGT"
WORDS2 = ["".join(p) for p in itertools.product(BASES, repeat=2)]
WORDS3 = ["".join(p) for p in itertools.product(BASES, repeat=3)]


def windows(genome, length, count, seed):
    """`count` non-overlapping windows of `length` bases at random positions, skipping any with N."""
    rng = random.Random(seed)
    starts = list(range(0, len(genome) - length, length))
    rng.shuffle(starts)
    out = []
    for s in starts:
        w = genome[s:s + length]
        if set(w) <= set(BASES):
            out.append(w)
            if len(out) == count:
                break
    return out


def graph(window):
    A = np.zeros((16, 16))
    for i in range(len(window) - 2):
        x, y, z = window[i:i + 3]
        u, v = WORDS2.index(x + y), WORDS2.index(y + z)
        if u != v:
            A[u, v] += 1
            A[v, u] += 1
    return A / len(window)


def kmer_features(window, k=3):
    words = WORDS3 if k == 3 else ["".join(p) for p in itertools.product(BASES, repeat=k)]
    idx = {w: i for i, w in enumerate(words)}
    f = np.zeros(len(words))
    for i in range(len(window) - k + 1):
        f[idx[window[i:i + k]]] += 1
    return f / max(1, f.sum())


def walk_features(A, kmax=8):
    out, P = [], np.eye(len(A))
    for k in range(1, kmax + 1):
        P = P @ A
        if k >= 2:
            out.append(np.trace(P))
    return np.array(out)


def spectrum_features(A):
    return np.sort(np.linalg.eigvalsh(A))


# --- exact GBS orbit probabilities ------------------------------------------------------------------------------------

def _partitions_of_patterns(m, N):
    """Every multiset of N modes out of m (as index lists with repeats), its orbit (sorted multiplicities), and
    1 / prod n_i!."""
    pats, orbit_ids, weights, orbits = [], [], [], {}
    for combo in itertools.combinations_with_replacement(range(m), N):
        counts = {}
        for c in combo:
            counts[c] = counts.get(c, 0) + 1
        key = tuple(sorted(counts.values(), reverse=True))
        orbits.setdefault(key, len(orbits))
        pats.append(combo)
        orbit_ids.append(orbits[key])
        weights.append(1.0 / math.prod(math.factorial(v) for v in counts.values()))
    return np.array(pats, np.int64), np.array(orbit_ids, np.int64), np.array(weights), orbits


@njit(cache=True)
def _haf(M, mask, n, memo, seen):
    if mask == 0:
        return 1.0
    if seen[mask]:
        return memo[mask]
    i = 0
    while not (mask >> i) & 1:
        i += 1
    rest = mask & ~(1 << i)
    total = 0.0
    j = i + 1
    while j < n:
        if (rest >> j) & 1 and M[i, j] != 0.0:
            total += M[i, j] * _haf(M, rest & ~(1 << j), n, memo, seen)
        j += 1
    seen[mask] = True
    memo[mask] = total
    return total


@njit(cache=True)
def _orbit_sums(A, pats, orbit_ids, weights, n_orbits):
    N = pats.shape[1]
    out = np.zeros(n_orbits)
    M = np.zeros((N, N))
    memo = np.zeros(1 << N)
    seen = np.zeros(1 << N, np.bool_)
    for p in range(pats.shape[0]):
        for a in range(N):
            for b in range(N):
                M[a, b] = A[pats[p, a], pats[p, b]]
        seen[:] = False
        h = _haf(M, (1 << N) - 1, N, memo, seen)
        out[orbit_ids[p]] += weights[p] * h * h
    return out


class OrbitFeatures:
    """The GBS feature vector for m-mode graphs: orbit probabilities given N photons, for each N in `totals`."""

    def __init__(self, m=16, totals=(4, 6, 8)):
        self.tables = [(N,) + _partitions_of_patterns(m, N) for N in totals]

    def names(self):
        return [f"N={N}:{list(k)}" for N, _, _, _, orbits in self.tables for k in orbits]

    def __call__(self, A):
        feats = []
        for N, pats, ids, w, orbits in self.tables:
            s = _orbit_sums(np.asarray(A, float), pats, ids, w, len(orbits))
            feats.append(s / s.sum() if s.sum() > 0 else s)
        return np.concatenate(feats)


def exact_pattern_probabilities(A, N):
    """p(pattern | N photons) for every pattern, for the check against PQCF's sampled device."""
    pats, ids, w, orbits = _partitions_of_patterns(len(A), N)
    probs = np.zeros(len(pats))
    for p, idx in enumerate(pats):
        M = np.asarray(A, float)[np.ix_(idx, idx)]
        h = _orbit_sums(np.asarray(A, float), pats[p:p + 1], np.zeros(1, np.int64), w[p:p + 1], 1)[0]
        probs[p] = h
    return pats, probs / probs.sum()
