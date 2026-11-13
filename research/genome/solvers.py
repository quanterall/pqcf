"""Solvers for the G0 instances: exact (OR-Tools CP-SAT) on the native formulation and on the QUBO, and the
classical heuristics every photonic solver will be measured against, simulated annealing and parallel tempering,
single-threaded with Numba."""
import time

import numpy as np
from numba import njit
from ortools.sat.python import cp_model

SCALE = 1000  # CP-SAT needs integer objectives: weights are scaled and rounded


def _solve(model, time_limit, workers=1):
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = time_limit
    s.parameters.num_workers = workers
    t0 = time.perf_counter()
    status = s.Solve(model)
    return s, status, time.perf_counter() - t0


def _status(status):
    return {cp_model.OPTIMAL: "optimal", cp_model.FEASIBLE: "feasible"}.get(status, "none")


# --- exact, native ------------------------------------------------------------------------------------------------

def exact_haplotype(hap, time_limit=60.0, workers=1):
    """max sum W_ij s_i s_j, s_i = 1 - 2 b_i, s_0 = +1; z_ij = b_i xor b_j, s_i s_j = 1 - 2 z_ij."""
    n = len(hap.truth)
    m = cp_model.CpModel()
    b = [m.NewBoolVar(f"b{i}") for i in range(n)]
    m.Add(b[0] == 0)
    terms = []
    for (i, j), w in hap.J.items():
        z = m.NewBoolVar(f"z{i}_{j}")
        m.AddBoolXOr([b[i], b[j], z.Not()])
        terms.append(2 * int(w) * z)
    m.Minimize(sum(terms))
    s, status, t = _solve(m, time_limit, workers)
    spins = np.array([1 - 2 * s.Value(x) for x in b]) if status in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None
    return {"status": _status(status), "seconds": t, "spins": spins,
            "energy": hap.energy(spins) if spins is not None else None}


def exact_layout(asm, time_limit=60.0):
    """The Hamiltonian path of largest total overlap, natively: CP-SAT's circuit constraint over the reads and a
    dummy node n (its arcs start and end the path)."""
    n = asm.n
    m = cp_model.CpModel()
    arcs, lits = [], {}
    for (u, v), w in asm.edges.items():
        lit = m.NewBoolVar(f"e{u}_{v}")
        lits[(u, v)] = lit
        arcs.append((u, v, lit))
    for v in range(n):
        arcs.append((n, v, m.NewBoolVar(f"s{v}")))
        arcs.append((v, n, m.NewBoolVar(f"t{v}")))
    m.AddCircuit(arcs)
    m.Maximize(sum(w * lits[e] for e, w in asm.edges.items()))
    s, status, t = _solve(m, time_limit)
    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return {"status": _status(status), "seconds": t, "path": None, "overlap": None}
    nxt = {u: v for (u, v), lit in lits.items() if s.Value(lit)}
    start = (set(range(n)) - set(nxt.values())).pop()
    path = [start]
    while path[-1] in nxt:
        path.append(nxt[path[-1]])
    return {"status": _status(status), "seconds": t, "path": path,
            "overlap": sum(asm.edges[(path[i], path[i + 1])] for i in range(n - 1))}


def exact_dks(net, k, time_limit=60.0):
    n = len(net.genes)
    m = cp_model.CpModel()
    x = [m.NewBoolVar(f"x{i}") for i in range(n)]
    m.Add(sum(x) == k)
    ys, inc = [], [[] for _ in range(n)]
    for (i, j) in net.edges:
        y = m.NewBoolVar(f"y{i}_{j}")
        m.AddImplication(y, x[i])
        m.AddImplication(y, x[j])
        ys.append(y)
        inc[i].append(y)
        inc[j].append(y)
    # a chosen gene has at most k - 1 links inside the module, an unchosen one none: tightens the bound
    for i in range(n):
        m.Add(sum(inc[i]) <= (k - 1) * x[i])
    m.Add(sum(ys) <= k * (k - 1) // 2)
    m.Maximize(sum(ys))
    s, status, t = _solve(m, time_limit)
    chosen = [i for i in range(n) if s.Value(x[i])] if status in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None
    return {"status": _status(status), "seconds": t, "chosen": chosen,
            "edges": net.density_of(chosen) if chosen is not None else None}


# --- exact, on the QUBO (to check a QUBO's penalties: its optimum must be the native one) ---------------------------

def exact_qubo(Q, c, n, time_limit=60.0):
    m = cp_model.CpModel()
    x = [m.NewBoolVar(f"x{i}") for i in range(n)]
    terms = []
    for (i, j), w in Q.items():
        wi = int(round(w * SCALE))
        if not wi:
            continue
        if i == j:
            terms.append(wi * x[i])
        else:
            y = m.NewBoolVar(f"y{i}_{j}")
            m.AddBoolAnd([x[i], x[j]]).OnlyEnforceIf(y)
            m.AddBoolOr([x[i].Not(), x[j].Not()]).OnlyEnforceIf(y.Not())
            terms.append(wi * y)
    m.Minimize(sum(terms))
    s, status, t = _solve(m, time_limit)
    xs = np.array([s.Value(v) for v in x]) if status in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None
    return {"status": _status(status), "seconds": t, "x": xs,
            "value": c + sum(w * xs[i] * xs[j] for (i, j), w in Q.items()) if xs is not None else None}


# --- heuristics on a QUBO -------------------------------------------------------------------------------------------

def dense(Q, n):
    """The QUBO as a diagonal d and a symmetric off-diagonal matrix S: x^T Q x = d.x + x^T S x / ... with
    E(x) = sum_i d_i x_i + sum_{i<j} q_ij x_i x_j, so flipping x_i changes E by (1 - 2 x_i)(d_i + sum_j S_ij x_j)
    where S_ij = S_ji = q_ij."""
    d = np.zeros(n)
    S = np.zeros((n, n))
    for (i, j), w in Q.items():
        if i == j:
            d[i] += w
        else:
            S[i, j] += w
            S[j, i] += w
    return d, S


def default_betas(d, S):
    # dwave-neal's default range: the hottest beta accepts the largest single-flip change half the time, the coldest
    # accepts the smallest non-zero change with probability 1%
    big = np.max(np.abs(d) + np.abs(S).sum(axis=1))
    nz = np.abs(np.concatenate([d, S[np.triu_indices_from(S, 1)]]))
    small = nz[nz > 1e-12].min()
    return np.log(2) / big, np.log(100) / small


@njit(cache=True)
def _sa(d, S, sweeps, b0, b1, seed):
    np.random.seed(seed)
    n = d.shape[0]
    x = np.random.randint(0, 2, n).astype(np.float64)
    f = S @ x
    e = d @ x + 0.5 * x @ (S @ x)
    best, bestx = e, x.copy()
    ratio = (b1 / b0) ** (1.0 / max(sweeps - 1, 1))
    beta = b0
    for _ in range(sweeps):
        for i in range(n):
            de = (1.0 - 2.0 * x[i]) * (d[i] + f[i])
            if de <= 0.0 or np.random.random() < np.exp(-beta * de):
                s = 1.0 - 2.0 * x[i]
                x[i] = 1.0 - x[i]
                for j in range(n):
                    f[j] += s * S[j, i]
                e += de
                if e < best:
                    best = e
                    bestx[:] = x
        beta *= ratio
    return best, bestx


@njit(cache=True)
def _pt(d, S, sweeps, betas, seed):
    np.random.seed(seed)
    n, R = d.shape[0], betas.shape[0]
    X = np.random.randint(0, 2, (R, n)).astype(np.float64)
    F = np.empty((R, n))
    E = np.empty(R)
    for r in range(R):
        F[r] = S @ X[r]
        E[r] = d @ X[r] + 0.5 * X[r] @ F[r]
    best = E.min()
    bestx = X[np.argmin(E)].copy()
    for _ in range(sweeps):
        for r in range(R):
            for i in range(n):
                de = (1.0 - 2.0 * X[r, i]) * (d[i] + F[r, i])
                if de <= 0.0 or np.random.random() < np.exp(-betas[r] * de):
                    s = 1.0 - 2.0 * X[r, i]
                    X[r, i] = 1.0 - X[r, i]
                    for j in range(n):
                        F[r, j] += s * S[j, i]
                    E[r] += de
                    if E[r] < best:
                        best = E[r]
                        bestx[:] = X[r]
        for r in range(R - 1):
            a = (betas[r + 1] - betas[r]) * (E[r + 1] - E[r])
            if a >= 0.0 or np.random.random() < np.exp(a):
                tmp = X[r].copy(); X[r] = X[r + 1]; X[r + 1] = tmp
                tmp = F[r].copy(); F[r] = F[r + 1]; F[r + 1] = tmp
                te = E[r]; E[r] = E[r + 1]; E[r + 1] = te
    return best, bestx


def simulated_annealing(d, S, c, sweeps, seed, betas=None):
    b0, b1 = betas or default_betas(d, S)
    best, x = _sa(d, S, sweeps, b0, b1, seed)
    return best + c, x


def parallel_tempering(d, S, c, sweeps, seed, replicas=16, betas=None):
    b0, b1 = betas or default_betas(d, S)
    ladder = np.geomspace(b0, b1, replicas)
    best, x = _pt(d, S, sweeps, ladder, seed)
    return best + c, x


# --- exact by dynamic programming over a band (haplotype phasing from reads spanning at most w + 1 sites) ---------

def exact_banded(n, J, h=None):
    """min E(s) = -sum J_ij s_i s_j - sum h_i s_i when every coupling has |i - j| <= w: a Viterbi pass over the
    last w spins, O(n 2^w w). Short reads couple only sites within one read, so w is the longest read span - 1;
    the standard phasing tools exploit the same structure (WhatsHap's dynamic programme is exponential only in the
    read coverage)."""
    t0 = time.perf_counter()
    w = max((j - i for (i, j) in J), default=0)
    nb = [[] for _ in range(n)]          # couplings of site i to earlier sites, as (offset, weight)
    for (i, j), v in J.items():
        nb[j].append((j - i, v))
    h = h or {}
    S = 1 << w                            # state: the last w spins, bit k = spin at position i-1-k (1 = -1)
    cost = np.full(S, np.inf)
    cost[0] = 0.0 if w else 0.0
    # positions before 0 do not exist: start from site 0 with every state over "virtual" earlier sites free of
    # couplings, so all states cost 0 and only real couplings count
    cost[:] = 0.0
    back = []
    for i in range(n):
        new = np.full(S, np.inf)
        arg = np.zeros(S, np.int64)
        for st in range(S):
            if not np.isfinite(cost[st]):
                continue
            for b in (0, 1):
                si = -1.0 if b else 1.0
                e = -h.get(i, 0.0) * si
                for off, v in nb[i]:
                    sj = -1.0 if (st >> (off - 1)) & 1 else 1.0
                    e -= v * si * sj
                ns = ((st << 1) | b) & (S - 1) if w else 0
                c = cost[st] + e
                if c < new[ns]:
                    new[ns], arg[ns] = c, st * 2 + b
        back.append(arg)
        cost = new
    st = int(np.argmin(cost))
    best = float(cost[st])
    spins = np.zeros(n, int)
    for i in range(n - 1, -1, -1):
        a = back[i][st]
        spins[i] = -1 if a & 1 else 1
        st = a >> 1
    return {"energy": best, "spins": spins, "seconds": time.perf_counter() - t0, "band": w}
