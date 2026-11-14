"""QuTiP's side of the quantum optics head-to-head (README.md): every task of tasks.json, timed, at QuTiP's default
tolerances and at tight ones. Writes qutip-default.json and qutip-tight.json in the format of qo_pqcf.erl.
    python qutip_side.py tasks.json out_dir
"""
import json
import math
import sys
import time
import warnings

warnings.filterwarnings("ignore", message="matplotlib not found")
import numpy as np  # noqa: E402
import qutip as qt  # noqa: E402

# the tightest tolerances QuTiP's default integrator (scipy's zvode, Adams) accepts without failing on these problems
TIGHT = {"atol": 1e-14, "rtol": 1e-12, "nsteps": 10**7}
BUDGET = 120.0


def timed(f):
    """The median of repeated runs, as in comparison/photonic: once if it takes a second or more, else up to 7 times
    within about two seconds."""
    f()  # untimed warm-up
    ts, v = [], None
    while True:
        t0 = time.perf_counter()
        v = f()
        ts.append(time.perf_counter() - t0)
        if len(ts) >= 7 or ts[0] >= 1.0 or ts[0] * (len(ts) + 1) >= 2.0:
            break
    return sorted(ts)[(len(ts) - 1) // 2], v


def coherent_amps(n, alpha):
    """|alpha> truncated at n photons and renormalised, as pq_lindblad:coherent/3."""
    a = np.array([alpha ** k / math.sqrt(math.factorial(k)) for k in range(n + 1)], dtype=complex)
    return a / np.linalg.norm(a)


def jc_ops(n, g):
    a = qt.tensor(qt.destroy(n + 1), qt.qeye(2))
    # the atom's level 1 (ground) is basis 0, level 2 (excited) basis 1: sigma- = |1><2|
    sm = qt.tensor(qt.qeye(n + 1), qt.basis(2, 0) * qt.basis(2, 1).dag())
    return a, sm, g * (a.dag() * sm + a * sm.dag())


def rabi(t, opts):
    a, sm, h = jc_ops(t["n"], t["g"])
    psi0 = qt.tensor(qt.basis(t["n"] + 1, 0), qt.basis(2, 1))
    return lambda: list(qt.sesolve(h, psi0, [0] + t["times"], e_ops=[sm.dag() * sm], options=opts).expect[0][1:])


def revival(t, opts):
    n = t["n"]
    a, sm, h = jc_ops(n, t["g"])
    field = qt.Qobj(np.concatenate([coherent_amps(n - 1, t["alpha"]), [0]]))
    psi0 = qt.tensor(field, qt.basis(2, 1))
    times = t["times"]
    return lambda: list(qt.sesolve(h, psi0, times, e_ops=[sm.dag() * sm], options=opts).expect[0])


def cavity_decay(t, opts):
    n = t["n"]
    a = qt.destroy(n + 1)
    psi0 = qt.Qobj(coherent_amps(n, t["alpha"]))
    c = [math.sqrt(t["kappa"]) * a]
    return lambda: list(qt.mesolve(0 * a.dag() * a, psi0, [0] + t["times"], c_ops=c, e_ops=[a.dag() * a],
                                   options=opts).expect[0][1:])


def fluorescence(t, _opts):
    sm = qt.basis(2, 0) * qt.basis(2, 1).dag()
    pe = sm.dag() * sm

    def run():
        out = []
        for om, dl in t["cases"]:
            h = 0.5 * om * (sm + sm.dag()) - dl * pe
            rho = qt.steadystate(h, [math.sqrt(t["gamma"]) * sm])
            out.append(qt.expect(pe, rho))
        return out
    return run


def damped_jc(t, n, opts):
    a, sm, h = jc_ops(n, t["g"])
    h = h + t["drive"] * (a + a.dag())
    c = [math.sqrt(t["kappa"]) * a, math.sqrt(t["gamma"]) * sm]
    psi0 = qt.tensor(qt.basis(n + 1, 0), qt.basis(2, 0))
    return lambda: [float(x[-1]) for x in qt.mesolve(h, psi0, [0, t["t"]], c_ops=c, e_ops=[a.dag() * a, sm.dag() * sm],
                                                     options=opts).expect]


def side(tasks, opts):
    out = {}
    for name, f in [("rabi", rabi), ("revival", revival), ("cavity_decay", cavity_decay),
                    ("fluorescence", fluorescence)]:
        s, v = timed(f(tasks[name], opts))
        out[name] = {"values": [float(np.real(x)) for x in v], "seconds": s}
    rows = []
    t = tasks["damped_jc"]
    for n in t["sizes"]:
        s, v = timed(damped_jc(t, n, opts))
        rows.append({"size": n, "values": v, "seconds": s})
        if s > BUDGET:
            break
    out["damped_jc"] = rows
    return out


def main():
    tasks = json.load(open(sys.argv[1]))
    for label, opts in [("default", {}), ("tight", TIGHT)]:
        res = {"tool": f"QuTiP ({label})", "version": qt.__version__, "tasks": side(tasks, opts)}
        with open(f"{sys.argv[2]}/qutip-{label}.json", "w") as f:
            json.dump(res, f)


if __name__ == "__main__":
    main()
