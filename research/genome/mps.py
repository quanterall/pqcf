"""G4: genomes as matrix product states (MPS), as in Creevey, Hassan, McCafferty, Hollenberg, Strelchuk
(arXiv:2508.06184): |psi> = L^(-1/2) sum_i |i>|b_i>, the position i on ceil(log2 L) qubits (most significant
first) followed by the base b_i on 2 qubits (A = 00, T = 01, G = 10, C = 11). An MPS's bond dimension at a cut is
the number of Schmidt values kept there; a state with no structure needs 2^min(left, right) of them.

- tt_svd(): the MPS of a whole state by sequential SVDs, each truncated so that the discarded weight is at most
  eps / (n - 1); the fidelity of the result is then computed exactly by contracting it. For genomes up to ~2^20
  amplitudes (phiX174, SARS-CoV-2).
- cut_spectrum(): the Schmidt values at one cut of a large genome's state without building the state. Splitting the
  position bits into a prefix p (the first `cut` bits) and a suffix, the state is the matrix M[p, (s, b)] =
  L^(-1/2) [base at p 2^k + s is b], and M M^T [p, q] = (the number of positions where blocks p and q of 2^k
  bases agree) / L: the squared Schmidt values are the eigenvalues of that block-overlap matrix.
"""
import math

import numpy as np

CODE = {"A": 0, "T": 1, "G": 2, "C": 3}


def clean(seq):
    """The sequence without unknown bases (N): the encoding has one of four bases at every position."""
    return "".join(c for c in seq if c in CODE)


def qubits(L):
    return max(1, math.ceil(math.log2(L))) + 2


def state(seq):
    L = len(seq)
    p = qubits(L) - 2
    psi = np.zeros(1 << (p + 2))
    idx = (np.arange(L) << 2) | np.array([CODE[c] for c in seq])
    psi[idx] = 1.0 / math.sqrt(L)
    return psi


def tt_svd(psi, n, eps):
    """Left-canonical MPS of psi (2^n amplitudes) with discarded weight <= eps / (n - 1) at each cut; returns the
    cores, the bond dimensions and the exact fidelity |<psi|mps>|^2."""
    cores, bonds = [], []
    rest = psi.reshape(1, -1)
    chi = 1
    for k in range(n - 1):
        rest = rest.reshape(chi * 2, -1)
        U, S, Vt = np.linalg.svd(rest, full_matrices=False)
        w = S ** 2
        tail = np.cumsum(w[::-1])[::-1]           # tail[j] = weight discarded if we keep j values
        keep = len(S)
        budget = eps / (n - 1)
        for j in range(1, len(S) + 1):
            if j == len(S) or tail[j] <= budget:
                keep = j
                break
        cores.append(U[:, :keep].reshape(chi, 2, keep))
        rest = S[:keep, None] * Vt[:keep]
        chi = keep
        bonds.append(keep)
    cores.append(rest.reshape(chi, 2, 1))
    approx = contract(cores)
    fid = float((psi @ approx) ** 2 / (approx @ approx))
    return cores, bonds, fid


def contract(cores):
    v = cores[0].reshape(2, -1)
    for c in cores[1:]:
        v = (v @ c.reshape(c.shape[0], -1)).reshape(-1, c.shape[2])
    return v.reshape(-1)


def parameters(bonds):
    dims = [1] + list(bonds) + [1]
    return sum(2 * dims[i] * dims[i + 1] for i in range(len(dims) - 1))


def cut_spectrum(seq, cut, dtype=np.float32):
    """Squared Schmidt values (descending) at the cut after the first `cut` position qubits, from the block-overlap
    matrix; the state's position register has p = qubits(L) - 2 bits."""
    L = len(seq)
    p = qubits(L) - 2
    k = p - cut                                   # bits within a block
    blocks, size = 1 << cut, 1 << k
    codes = np.full(blocks * size, -1, np.int8)
    codes[:L] = [CODE[c] for c in seq]
    codes = codes.reshape(blocks, size)
    gram = np.zeros((blocks, blocks), dtype)
    for b in range(4):
        X = (codes == b).astype(dtype)
        gram += X @ X.T
    ev = np.linalg.eigvalsh(gram.astype(np.float64) / L)[::-1]
    return np.clip(ev, 0, None)


def bond_for_error(w, eps):
    """Schmidt values needed at a cut for discarded weight <= eps (w: squared values, descending, summing to 1)."""
    tail = np.cumsum(w[::-1])[::-1]
    for j in range(1, len(w) + 1):
        if j == len(w) or tail[j] <= eps:
            return j
    return len(w)
