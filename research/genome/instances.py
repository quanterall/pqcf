"""Problem instances from genomic data (G0 of docs/research/genome/README.md). Every instance is written two ways:
its native formulation, which a classical exact solver works on directly, and the binary quadratic form (QUBO /
Ising) that annealers, coherent Ising machines and QAOA need. Both are kept so the cost of the translation can be
measured, as Maurizio & Mazzola ask (arXiv:2507.04111).

A QUBO here is (Q, c): minimise x^T Q x + c over x in {0,1}^n, Q upper triangular as a dict {(i, j): q}, i <= j.
"""
import csv
import random
from dataclasses import dataclass, field

import numpy as np

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def read_fasta(path):
    with open(path) as f:
        return "".join(line.strip() for line in f if not line.startswith(">")).upper()


# --- haplotype phasing ---------------------------------------------------------------------------------------------

@dataclass
class Haplotype:
    """Two haplotypes over n heterozygous sites, h and -h (h in {+1, -1}^n), and the fragments sequencing reads give:
    each fragment a {site: allele} map of the sites it covers, alleles +1 / -1, read from one haplotype with
    errors. The Ising model is the one HapCUT's Max-Cut relaxation optimises (Bansal & Bafna 2008): for sites i, j
    on the same fragment, a_i a_j says whether they are in phase, and E(s) = -sum_ij W_ij s_i s_j with
    W_ij = sum over fragments of a_i a_j."""
    truth: np.ndarray
    fragments: list
    J: dict = field(default_factory=dict)

    def energy(self, s):
        return -sum(w * s[i] * s[j] for (i, j), w in self.J.items())

    def mec(self, s):
        # minimum error correction: per fragment, the alleles to change to fit h or -h, whichever is fewer
        total = 0
        for f in self.fragments:
            wrong = sum(1 for i, a in f.items() if a != s[i])
            total += min(wrong, len(f) - wrong)
        return total

    def switch_errors(self, s):
        # adjacent sites phased differently from the truth (up to the global flip)
        rel = s * self.truth
        return int(np.sum(rel[1:] != rel[:-1]))


def haplotype(n_sites, coverage=10, span=(2, 6), error=0.05, seed=0):
    rng = random.Random(seed)
    truth = np.array([rng.choice((1, -1)) for _ in range(n_sites)])
    fragments = []
    n_frag = int(coverage * n_sites / (sum(span) / 2))
    for _ in range(n_frag):
        length = rng.randint(*span)
        start = rng.randint(0, n_sites - length)
        side = rng.choice((1, -1))
        fragments.append({i: (-1 if rng.random() < error else 1) * side * int(truth[i])
                          for i in range(start, start + length)})
    J = {}
    for f in fragments:
        sites = sorted(f)
        for a in range(len(sites)):
            for b in range(a + 1, len(sites)):
                i, j = sites[a], sites[b]
                J[(i, j)] = J.get((i, j), 0) + f[i] * f[j]
    return Haplotype(truth, fragments, {k: v for k, v in J.items() if v})


def ising_to_qubo(J, h=None, n=None):
    """E(s) = -sum J_ij s_i s_j - sum h_i s_i with s = 2x - 1, as a QUBO."""
    h = h or {}
    Q, c = {}, 0.0
    for (i, j), w in J.items():
        # -w (2xi - 1)(2xj - 1) = -4w xi xj + 2w xi + 2w xj - w
        Q[(i, j)] = Q.get((i, j), 0) - 4 * w
        Q[(i, i)] = Q.get((i, i), 0) + 2 * w
        Q[(j, j)] = Q.get((j, j), 0) + 2 * w
        c -= w
    for i, v in h.items():
        Q[(i, i)] = Q.get((i, i), 0) - 2 * v
        c += v
    return Q, c


# --- genome assembly (overlap-layout-consensus) ---------------------------------------------------------------------

@dataclass
class Assembly:
    """Reads of a genome segment, their order along it (the truth), and the overlap graph: an edge u -> v of weight
    the length of the longest suffix of u equal (up to max_mismatch) to a prefix of v, at least min_overlap. The
    layout step is the Hamiltonian path of largest total overlap (Boev et al., arXiv:2004.06719)."""
    reads: list
    order: list
    edges: dict
    segment: str

    @property
    def n(self):
        return len(self.reads)


def overlap(a, b, min_overlap, max_mismatch):
    for k in range(min(len(a), len(b)) - 1, min_overlap - 1, -1):
        if sum(x != y for x, y in zip(a[-k:], b[:k])) <= max_mismatch:
            return k
    return 0


def assembly(genome, n_reads, read_len=100, coverage=4, error=0.0, min_overlap=20, max_mismatch=0, seed=0,
             shuffle=True):
    rng = random.Random(seed)
    seg_len = int(n_reads * read_len / coverage)
    start = rng.randint(0, len(genome) - seg_len - 1)
    segment = genome[start:start + seg_len]
    # read starts evenly spaced with a random jitter of up to a third of the spacing, the first at 0 and the last at
    # the segment's end: consecutive reads overlap by at least read_len - 5/3 spacing, so the reads cover the segment
    step = (seg_len - read_len) / (n_reads - 1)
    assert read_len - 5 * step / 3 >= min_overlap, "coverage too low for min_overlap"
    starts = [0] + [int(round(i * step + rng.uniform(-step / 3, step / 3))) for i in range(1, n_reads - 1)] \
        + [seg_len - read_len]
    starts = sorted(starts)
    reads = []
    for s in starts:
        r = list(segment[s:s + read_len])
        for i in range(len(r)):
            if rng.random() < error:
                r[i] = rng.choice([b for b in "ACGT" if b != r[i]])
        reads.append("".join(r))
    idx = list(range(n_reads))
    if shuffle:
        rng.shuffle(idx)
    reads = [reads[i] for i in idx]
    order = sorted(range(n_reads), key=lambda k: starts[idx[k]])
    edges = {}
    for u in range(n_reads):
        for v in range(n_reads):
            if u != v:
                k = overlap(reads[u], reads[v], min_overlap, max_mismatch)
                if k:
                    edges[(u, v)] = k
    return Assembly(reads, order, edges, segment)


def layout_qubo_positions(asm, penalty=None, reward=1.0):
    """Boev et al. Eq. (1) with an overlap reward: x[v, p] = read v at position p of the path; N^2 variables.
    Minimise A (sum_v (1 - sum_p x_vp)^2 + sum_p (1 - sum_v x_vp)^2 + sum_{(u,v) not an edge} sum_p x_up x_v,p+1)
    - reward sum_{(u,v) edge} w_uv/read_len sum_p x_up x_v,p+1. With weights at most 1, a violated constraint saves
    at most 2 in reward, so A = 3 keeps every optimum a Hamiltonian path."""
    n = asm.n
    A = penalty if penalty is not None else 3.0 * reward
    L = max(len(r) for r in asm.reads)
    var = lambda v, p: v * n + p
    Q, c = {}, 0.0

    def add(i, j, w):
        if i > j:
            i, j = j, i
        Q[(i, j)] = Q.get((i, j), 0) + w

    # (1 - sum x)^2 = 1 - 2 sum x + sum x^2 + 2 sum_{a<b} x_a x_b, with x^2 = x
    for group in ([[var(v, p) for p in range(n)] for v in range(n)] + [[var(v, p) for v in range(n)] for p in range(n)]):
        c += A
        for a in range(len(group)):
            add(group[a], group[a], -A)
            for b in range(a + 1, len(group)):
                add(group[a], group[b], 2 * A)
    for u in range(n):
        for v in range(n):
            if u == v:
                continue
            w = -reward * asm.edges[(u, v)] / L if (u, v) in asm.edges else A
            for p in range(n - 1):
                add(var(u, p), var(v, p + 1), w)
    return Q, c


def layout_qubo_edges(asm, penalty=None, reward=1.0):
    """Boev et al. Eq. (2) with an overlap reward: y[e] = edge e in the path; one variable per edge. Each vertex has
    at most one outgoing and one incoming path edge, and the path has n - 1 edges: A (sum_u (sum_out y)^2 - sum_out
    y) etc. would only forbid two; Boev's (1 - sum)^2 form needs a start and an end vertex. Here: at most one out
    and one in, by the pair penalty A y_e y_f for edges sharing a tail or a head, and exactly n - 1 edges by
    A (n - 1 - sum y)^2. Without acyclicity a set of edges can also form cycles; the overlap graph of reads from a
    linear segment has none when overlaps are unique, which the solver results are checked for."""
    n = asm.n
    A = penalty if penalty is not None else 3.0 * reward
    L = max(len(r) for r in asm.reads)
    E = sorted(asm.edges)
    Q, c = {}, 0.0

    def add(i, j, w):
        if i > j:
            i, j = j, i
        Q[(i, j)] = Q.get((i, j), 0) + w

    for a in range(len(E)):
        add(a, a, -reward * asm.edges[E[a]] / L)
        for b in range(a + 1, len(E)):
            if E[a][0] == E[b][0] or E[a][1] == E[b][1]:
                add(a, b, A)
    # A (m - sum y)^2 = A m^2 - 2 A m sum y + A sum y + 2 A sum_{a<b} y_a y_b
    m = n - 1
    c += A * m * m
    for a in range(len(E)):
        add(a, a, A - 2 * A * m)
        for b in range(a + 1, len(E)):
            add(a, b, 2 * A)
    return Q, c, E


# --- dense modules in a gene network ---------------------------------------------------------------------------------

@dataclass
class Network:
    genes: list
    edges: dict  # {(i, j): score}, i < j

    def density_of(self, chosen):
        s = set(chosen)
        return sum(1 for (i, j) in self.edges if i in s and j in s)


def string_network(path, threshold=0.9):
    genes, edges = [], {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if float(row["score"]) < threshold:
                continue
            a, b = row["preferredName_A"], row["preferredName_B"]
            for g in (a, b):
                if g not in genes:
                    genes.append(g)
            i, j = sorted((genes.index(a), genes.index(b)))
            edges[(i, j)] = max(edges.get((i, j), 0), float(row["score"]))
    return Network(genes, edges)


def dks_qubo(net, k, penalty=None):
    """Densest k-subgraph: maximise the edges among k chosen genes. Minimise -sum_{ij edge} x_i x_j +
    P (k - sum x)^2; adding a vertex gains at most its degree, so P = max degree + 1 keeps the optimum at size k."""
    n = len(net.genes)
    deg = [0] * n
    for (i, j) in net.edges:
        deg[i] += 1
        deg[j] += 1
    P = penalty if penalty is not None else max(deg) + 1
    Q, c = {}, P * k * k
    for (i, j) in net.edges:
        Q[(i, j)] = Q.get((i, j), 0) - 1
    for i in range(n):
        Q[(i, i)] = Q.get((i, i), 0) + P - 2 * P * k
        for j in range(i + 1, n):
            Q[(i, j)] = Q.get((i, j), 0) + 2 * P
    return Q, c


def qubo_value(Q, c, x):
    return c + sum(w * x[i] * x[j] for (i, j), w in Q.items())


def qubo_to_ising(Q, c):
    """x = (1 + s)/2: x^T Q x + c = -sum_{i<j} J_ij s_i s_j - sum_i h_i s_i + offset."""
    J, h, off = {}, {}, c
    for (i, j), w in Q.items():
        if i == j:
            h[i] = h.get(i, 0) - w / 2
            off += w / 2
        else:
            J[(i, j)] = J.get((i, j), 0) - w / 4
            h[i] = h.get(i, 0) - w / 4
            h[j] = h.get(j, 0) - w / 4
            off += w / 4
    return J, h, off
