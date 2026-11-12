"""The peers' side of the photonic head-to-head (README.md): every task of inputs.json that a program supports,
through its public functions, timed the same way as photonic_pqcf.erl; writes <tool>.json.
    python peers.py thewalrus|piquasso|perceval|deepquantum|strawberryfields <dir>
Each program runs in its own process, pinned to one core with one thread.
"""
import json
import math
import os
import sys
import time

import numpy as np

BUDGET = 120.0  # a size slower than this is the last one run for that task


def timed(f):
    # the median of repeated runs: once if it takes a second or more, else up to 7 times within about two seconds
    t0 = time.perf_counter()
    v = f()
    t1 = time.perf_counter() - t0
    ts = [t1]
    while len(ts) < 7 and t1 < 1.0 and t1 * (len(ts) + 1) < 2.0:
        t0 = time.perf_counter()
        f()
        ts.append(time.perf_counter() - t0)
    return sorted(ts)[(len(ts) + 1) // 2 - 1], v


def once(f):
    t0 = time.perf_counter()
    v = f()
    return time.perf_counter() - t0, v


def sized(cases, key, f):
    # one untimed call on the smallest case first: just-in-time compilation (Numba, PyTorch) is not timed
    f(cases[0])
    out = []
    for c in cases:
        r = f(c)
        r["size"] = c[key]
        out.append(r)
        if r.get("not_run") or r["seconds"] > BUDGET:
            break
    return out


def cmat(rows):
    return np.array([[complex(re, im) for re, im in row] for row in rows])


def occupation(photons, modes):
    return [1] * photons + [0] * (modes - photons)


def counts(patterns):
    c = {}
    for p in patterns:
        p = tuple(int(x) for x in p)
        c[p] = c.get(p, 0) + 1
    return [{"pattern": list(p), "count": n} for p, n in sorted(c.items())]


def mean_photons(probs, modes):
    # probs: an array of shape (d,) * modes, the probability of every per-mode occupation
    d = probs.shape[0]
    n = np.arange(d)
    out = []
    for i in range(modes):
        marg = probs.sum(axis=tuple(j for j in range(modes) if j != i))
        out.append(float((n * marg).sum()))
    return out, float(probs.sum())


def xpxp_real(O):
    # O = I - Q^-1 in the (a, a^dagger) ordering, as the quadratures x = (a + a^dagger)/sqrt 2, p = (a - a^dagger)/(i
    # sqrt 2), interleaved per mode; a change of basis within each mode leaves every determinant of the torontonian
    k = O.shape[0] // 2
    W = np.zeros((2 * k, 2 * k), dtype=complex)
    for j in range(k):
        W[2 * j, j], W[2 * j, k + j] = 1 / math.sqrt(2), 1 / math.sqrt(2)
        W[2 * j + 1, j], W[2 * j + 1, k + j] = -1j / math.sqrt(2), 1j / math.sqrt(2)
    R = W @ O @ np.linalg.inv(W)
    assert np.abs(R.imag).max() < 1e-10, np.abs(R.imag).max()
    return np.ascontiguousarray(R.real)


# --- The Walrus ---------------------------------------------------------------------------------------------------

def thewalrus(inp):
    import thewalrus as tw

    def perm(c):
        A = cmat(c["matrix"])
        best = None
        for method in ("glynn", "ryser", "bbfg"):
            t, v = timed(lambda: tw.perm(A, method=method))
            if best is None or t < best[0]:
                best = (t, v, method)
        return {"seconds": best[0], "value": [best[1].real, best[1].imag], "method": best[2]}

    def tor(c):
        O = cmat(c["matrix"])
        t, v = timed(lambda: tw.tor(O))
        return {"seconds": t, "value": complex(v).real}

    return {"version": tw.__version__, "permanent": sized(inp["permanent"], "n", perm),
            "torontonian": sized(inp["torontonian"], "k", tor)}


# --- Piquasso -----------------------------------------------------------------------------------------------------

def piquasso(inp):
    import piquasso as pq
    from piquasso._math.permanent import permanent
    from piquasso._math.torontonian import torontonian

    def perm(c):
        A = cmat(c["matrix"])
        ones = np.ones(len(A), dtype=np.int32)
        t, v = timed(lambda: permanent(A, ones, ones))
        return {"seconds": t, "value": [complex(v).real, complex(v).imag]}

    def tor(c):
        R = xpxp_real(cmat(c["matrix"]))
        t, v = timed(lambda: torontonian(R))
        return {"seconds": t, "value": float(v)}

    def program(occ, U):
        with pq.Program() as p:
            pq.Q() | pq.StateVector(occ)
            pq.Q() | pq.Interferometer(U)
            pq.Q() | pq.ParticleNumberMeasurement()
        return p

    def sample(c, shots, seed=1):
        U = cmat(c["unitary"])
        sim = pq.SamplingSimulator(d=c["modes"], config=pq.Config(seed_sequence=seed))
        p = program(occupation(c["photons"], c["modes"]), U)
        return lambda: sim.execute(p, shots=shots).samples

    def boson(c):
        t, _ = once(sample(c, c["shots"]))
        return {"seconds": t, "shots": c["shots"], "samples_per_second": c["shots"] / t}

    def fock(f):
        out = []
        for cut in f["total_cutoffs"]:
            with pq.Program() as p:
                pq.Q() | pq.Vacuum()
                for g in f["gates"]:
                    name, *a = g
                    if name == "phaseshift":
                        pq.Q(a[0]) | pq.Phaseshifter(phi=a[1])
                    elif name == "beamsplitter":
                        pq.Q(a[0], a[1]) | pq.Beamsplitter(theta=a[2], phi=a[3])
                    elif name == "squeezing":
                        pq.Q(a[0]) | pq.Squeezing(r=a[1], phi=a[2])
                    elif name == "displacement":
                        pq.Q(a[0]) | pq.Displacement(r=a[1], phi=a[2])
                    elif name == "kerr":
                        pq.Q(a[0]) | pq.Kerr(xi=a[1])
            sim = pq.PureFockSimulator(d=f["modes"], config=pq.Config(cutoff=cut + 1))
            if cut == f["total_cutoffs"][0]:
                sim.execute(p)
            t, st = timed(lambda: sim.execute(p).state)
            pm = st.fock_probabilities_map
            means = [float(sum(occ[i] * pr for occ, pr in pm.items())) for i in range(f["modes"])]
            out.append({"cutoff": cut, "truncation": "total", "dimension": math.comb(cut + f["modes"], f["modes"]),
                        "seconds": t, "mean_photons": means, "norm": float(st.norm)})
            if t > BUDGET:
                break
        return out

    bc = inp["boson_check"]
    return {"version": pq.__version__, "permanent": sized(inp["permanent"], "n", perm),
            "torontonian": sized(inp["torontonian"], "k", tor),
            "boson_sampling": sized(inp["boson_sampling"], "photons", boson),
            "boson_check": counts(sample(bc, bc["shots"])()), "fock": fock(inp["fock"])}


# --- Perceval -----------------------------------------------------------------------------------------------------

def perceval(inp):
    import exqalibur as xq
    import perceval as pcvl
    from perceval.algorithm import Sampler

    def perm(c):
        A = np.ascontiguousarray(cmat(c["matrix"]))
        t, v = timed(lambda: xq.permanent_cx(A, "glynn"))
        return {"seconds": t, "value": [v.real, v.imag], "method": "exqalibur glynn"}

    def sampler(c):
        proc = pcvl.Processor("CliffordClifford2017", pcvl.Unitary(pcvl.Matrix(cmat(c["unitary"]))))
        proc.with_input(pcvl.BasicState(occupation(c["photons"], c["modes"])))
        return Sampler(proc)

    def boson(c):
        s = sampler(c)
        t, _ = once(lambda: s.samples(c["shots"]))
        return {"seconds": t, "shots": c["shots"], "samples_per_second": c["shots"] / t}

    bc = inp["boson_check"]
    check = [list(b) for b in sampler(bc).samples(bc["shots"])["results"]]
    return {"version": pcvl.__version__, "permanent": sized(inp["permanent"], "n", perm),
            "boson_sampling": sized(inp["boson_sampling"], "photons", boson), "boson_check": counts(check)}


# --- DeepQuantum --------------------------------------------------------------------------------------------------

DQ_MAX_OUTPUTS = 200_000  # DeepQuantum's boson sampler first computes the probability of every output pattern


def deepquantum(inp):
    import torch
    import deepquantum as dq
    torch.set_num_threads(1)

    def perm(c):
        A = torch.tensor(cmat(c["matrix"]))
        t, v = timed(lambda: dq.permanent(A))
        v = complex(v)
        return {"seconds": t, "value": [v.real, v.imag]}

    def tor(c):
        O = torch.tensor(cmat(c["matrix"]))
        t, v = timed(lambda: dq.torontonian(O))
        return {"seconds": t, "value": complex(v).real}

    def circuit(c):
        cir = dq.QumodeCircuit(nmode=c["modes"], init_state=occupation(c["photons"], c["modes"]), basis=True)
        cir.any(torch.tensor(cmat(c["unitary"])), wires=list(range(c["modes"])))
        return cir

    def sample(c, shots):
        cir = circuit(c)

        def run():
            cir()
            return cir.measure(shots=shots)
        return run

    def boson(c):
        outputs = math.comb(c["modes"] + c["photons"] - 1, c["photons"])
        if outputs > DQ_MAX_OUTPUTS:
            return {"not_run": f"{outputs:.3g} output patterns", "seconds": None}
        t, _ = once(sample(c, c["shots"]))
        return {"seconds": t, "shots": c["shots"], "samples_per_second": c["shots"] / t}

    def check():
        bc = inp["boson_check"]
        res = sample(bc, bc["shots"])()
        pats = []
        for state, n in res.items():
            occ = [int(x) for x in (state.state.tolist() if hasattr(state, "state") else state)]
            pats.extend([occ] * int(n))
        return counts(pats)

    def fock(f):
        out = []
        for d in f["per_mode_cutoffs"]:
            def run():
                cir = dq.QumodeCircuit(nmode=f["modes"], init_state="vac", cutoff=d, backend="fock", basis=False)
                for g in f["gates"]:
                    name, *a = g
                    if name == "phaseshift":
                        cir.ps(a[0], inputs=a[1])
                    elif name == "beamsplitter":
                        cir.bs([a[0], a[1]], inputs=[a[2], a[3]])
                    elif name == "squeezing":
                        cir.s(a[0], r=a[1], theta=a[2])
                    elif name == "displacement":
                        cir.d(a[0], r=a[1], theta=a[2])
                    elif name == "kerr":
                        cir.k(a[0], inputs=a[1])
                return cir()
            if d == f["per_mode_cutoffs"][0]:
                run()
            t, psi = timed(run)
            probs = (psi.abs() ** 2).reshape((d,) * f["modes"]).numpy()
            means, norm = mean_photons(probs, f["modes"])
            out.append({"cutoff": d, "truncation": "per mode", "dimension": d ** f["modes"], "seconds": t,
                        "mean_photons": means, "norm": norm})
            if t > BUDGET:
                break
        return out

    return {"version": dq.__version__, "permanent": sized(inp["permanent"], "n", perm),
            "torontonian": sized(inp["torontonian"], "k", tor),
            "boson_sampling": sized(inp["boson_sampling"], "photons", boson), "boson_check": check(),
            "fock": fock(inp["fock"])}


# --- Strawberry Fields (its own image: NumPy 1.x) -----------------------------------------------------------------

def strawberryfields(inp):
    import strawberryfields as sf
    from strawberryfields import ops

    def fock(f):
        out = []
        for d in f["per_mode_cutoffs"]:
            prog = sf.Program(f["modes"])
            with prog.context as q:
                for g in f["gates"]:
                    name, *a = g
                    if name == "phaseshift":
                        ops.Rgate(a[1]) | q[a[0]]
                    elif name == "beamsplitter":
                        ops.BSgate(a[2], a[3]) | (q[a[0]], q[a[1]])
                    elif name == "squeezing":
                        ops.Sgate(a[1], a[2]) | q[a[0]]
                    elif name == "displacement":
                        ops.Dgate(a[1], a[2]) | q[a[0]]
                    elif name == "kerr":
                        ops.Kgate(a[1]) | q[a[0]]

            def run():
                return sf.Engine("fock", backend_options={"cutoff_dim": d}).run(prog).state
            if d == f["per_mode_cutoffs"][0]:
                run()
            t, st = timed(run)
            means = [float(st.mean_photon(i)[0]) for i in range(f["modes"])]
            out.append({"cutoff": d, "truncation": "per mode", "dimension": d ** f["modes"], "seconds": t,
                        "mean_photons": means, "norm": float(st.trace())})
            if t > BUDGET:
                break
        return out

    return {"version": sf.__version__, "fock": fock(inp["fock"])}


if __name__ == "__main__":
    tool, out = sys.argv[1], sys.argv[2]
    inp = json.load(open(os.path.join(out, "inputs.json")))
    rec = {"thewalrus": thewalrus, "piquasso": piquasso, "perceval": perceval, "deepquantum": deepquantum,
           "strawberryfields": strawberryfields}[tool](inp)
    rec["tool"] = tool
    json.dump(rec, open(os.path.join(out, f"{tool}.json"), "w"))
    print(tool, "done")
