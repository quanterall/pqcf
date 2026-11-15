"""QuTiP's PIQS (Shammah et al., PRA 98, 063815 (2018), arXiv:1805.05129) on three Dicke-basis problems, for PQCF's
pq_dicke tests (quantum optics plan, Q5). Writes test/fixtures/piqs_dicke.json.
    docker run --rm -v "$PWD":/w -w /w pqcf-qutip python comparison/dicke/piqs_check.py
"""
import json
import warnings

warnings.filterwarnings("ignore")
import numpy as np  # noqa: E402
import qutip as qt  # noqa: E402
from qutip import piqs  # noqa: E402

OPTS = {"atol": 1e-12, "rtol": 1e-12, "nsteps": 10**7}


def run(N, rates, H, rho0, times):
    system = piqs.Dicke(N, **rates)
    if H is not None:
        system.hamiltonian = H
    L = system.liouvillian()
    jz = piqs.jspin(N, "z")
    jp, jm = piqs.jspin(N, "+"), piqs.jspin(N, "-")
    res = qt.mesolve(L, rho0, [0.0] + list(times), e_ops=[jz, jp * jm], options=OPTS)
    return [float(np.real(x)) for x in res.expect[0][1:]], [float(np.real(x)) for x in res.expect[1][1:]]


out = {"qutip": qt.__version__}
# A: N = 6, every process, no Hamiltonian, all atoms excited
rA = dict(emission=0.5, dephasing=0.3, pumping=0.2, collective_emission=1.0, collective_dephasing=0.4,
          collective_pumping=0.1)
tA = [0.25, 0.5, 1.0, 2.0]
jzA, jjA = run(6, rA, None, piqs.excited(6), tA)
out["A"] = {"N": 6, "rates": rA, "t": tA, "jz": jzA, "jpjm": jjA}
# B: N = 5, H = omega Jz + omega_x Jx, local and collective emission, from the ground state
rB = dict(emission=0.25, collective_emission=0.5)
HB = 1.0 * piqs.jspin(5, "z") + 0.75 * piqs.jspin(5, "x")
tB = [0.5, 1.0, 2.0]
jzB, jjB = run(5, rB, HB, piqs.ground(5), tB)
out["B"] = {"N": 5, "rates": rB, "omega": 1.0, "omega_x": 0.75, "t": tB, "jz": jzB, "jpjm": jjB}
# C: the superradiant burst, N = 100, collective emission only, all excited: J+J- on a grid around the peak
tC = [round(0.002 * k, 3) for k in range(1, 61)]
jzC, jjC = run(100, dict(collective_emission=1.0), None, piqs.excited(100), tC)
out["C"] = {"N": 100, "rates": {"collective_emission": 1.0}, "t": tC, "jz": jzC, "jpjm": jjC}
json.dump(out, open("test/fixtures/piqs_dicke.json", "w"), indent=1)
k = int(np.argmax(jjC))
print(f"A jz {jzA}\nB jz {jzB}\nC peak of J+J- at t = {tC[k]}: {jjC[k]:.6g} (N^2/4 = 2500; ln N / N = {np.log(100) / 100:.4f})")
