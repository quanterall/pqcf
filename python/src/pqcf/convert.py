"""Converters between pqcf.Circuit and Perceval and Piquasso, in the conventions measured by
test/oracle/converter_conventions.py (the same as pq_convert on the Erlang side):

- Perceval: our B(theta, phi) is BS.Rx(theta = 2 theta, phi_bl = pi/2 - phi, phi_br = phi - pi/2) on adjacent modes
  (m2 = m1 + 1); PS(phi) is our phase shift. Other components, and beamsplitters on other mode pairs, are carried
  exactly as unitaries (their computed matrices). Perceval circuits are passive.
- Piquasso: Beamsplitter, Phaseshifter, Squeezing, Displacement and Kerr are ours one for one (Kerr's xi is our
  kappa); Interferometer is our unitary; its covariance is twice ours.

Both libraries are imported only when a converter is called. Parameters must be numbers (bind symbols first).
"""

from __future__ import annotations

import math

from .circuit import Circuit
from .linear import gate_matrix, _num


def _np():
    import numpy as np
    return np


# Perceval

def from_perceval(circuit, input_state=None) -> Circuit:
    """a pqcf.Circuit from a perceval Circuit (and optionally a BasicState input, as a Fock pattern)"""
    import perceval as pcvl  # noqa: F401
    np = _np()
    n = circuit.m
    occ = [int(x) for x in input_state] if input_state is not None else None
    out = Circuit(n, fock=occ)
    for modes, comp in circuit:
        modes = list(modes)
        u = np.array(comp.compute_unitary())
        if len(modes) == 1:
            out.phaseshift(modes[0], float(np.angle(u[0, 0])))
        else:
            out.unitary(modes, u)
    return out


def to_perceval(circuit):
    """a perceval Circuit from a passive pqcf.Circuit (measurements dropped; the input is not part of it)"""
    import perceval as pcvl
    np = _np()
    doc = circuit.to_json() if isinstance(circuit, Circuit) else circuit
    n = doc["n_modes"]
    pc = pcvl.Circuit(n)
    for ins in doc["instructions"]:
        op = ins["op"]
        if op.startswith("measure_"):
            continue
        if op == "beamsplitter" and ins["m2"] == ins["m1"] + 1:
            t, p = _num(ins["theta"]), _num(ins["phi"])
            pc.add(ins["m1"], pcvl.BS.Rx(theta=2 * t, phi_bl=math.pi / 2 - p, phi_br=p - math.pi / 2))
        elif op == "phaseshift":
            pc.add(ins["mode"], pcvl.PS(_num(ins["phi"])))
        elif op in ("beamsplitter", "unitary"):
            # exactly, as a unitary on the contiguous span of its modes
            ms = [ins["m1"], ins["m2"]] if op == "beamsplitter" else ins["modes"]
            lo, hi = min(ms), max(ms)
            g = gate_matrix(ins, n)[lo:hi + 1, lo:hi + 1]
            pc.add(lo, pcvl.Unitary(pcvl.Matrix(np.array(g))))
        else:
            raise ValueError(f"{op} has no Perceval equivalent (Perceval circuits are passive)")
    return pc


# Piquasso

_PIQ = {"beamsplitter": "Beamsplitter", "phaseshift": "Phaseshifter", "squeezing": "Squeezing",
        "displacement": "Displacement", "kerr": "Kerr"}


def to_piquasso(circuit):
    """a piquasso Program from a pqcf.Circuit (vacuum, Fock or Gaussian input; the supported gates and
    measurements); returns (program, d)"""
    import piquasso as pq
    np = _np()
    doc = circuit.to_json() if isinstance(circuit, Circuit) else circuit
    n = doc["n_modes"]
    inp = doc.get("input", {"kind": "vacuum"})
    with pq.Program() as prog:
        if inp["kind"] == "fock":
            state = pq.StateVector(inp["occupation"]) if hasattr(pq, "StateVector") else pq.NumberState(inp["occupation"])
            pq.Q(*range(n)) | state * 1.0
        elif inp["kind"] == "gaussian":
            mean = np.array(inp["mean"])
            cov = np.array(inp["cov"]).reshape(2 * n, 2 * n) * 2  # Piquasso's covariance is twice ours
            # xxpp -> xpxp
            perm = [k // 2 + (k % 2) * n for k in range(2 * n)]
            pq.Q(*range(n)) | pq.Mean(mean[perm]) | pq.Covariance(cov[np.ix_(perm, perm)])
        for ins in doc["instructions"]:
            op = ins["op"]
            if op == "beamsplitter":
                pq.Q(ins["m1"], ins["m2"]) | pq.Beamsplitter(theta=_num(ins["theta"]), phi=_num(ins["phi"]))
            elif op == "phaseshift":
                pq.Q(ins["mode"]) | pq.Phaseshifter(phi=_num(ins["phi"]))
            elif op in ("squeezing", "displacement"):
                gate = getattr(pq, _PIQ[op])
                pq.Q(ins["mode"]) | gate(r=_num(ins["r"]), phi=_num(ins["phi"]))
            elif op == "kerr":
                pq.Q(ins["mode"]) | pq.Kerr(xi=_num(ins["kappa"]))
            elif op == "unitary":
                k = len(ins["modes"])
                u = np.array([complex(z["re"], z["im"]) for z in ins["matrix"]]).reshape(k, k)
                pq.Q(*ins["modes"]) | pq.Interferometer(u)
            elif op == "measure_pnr":
                pq.Q(*ins["modes"]) | pq.ParticleNumberMeasurement()
            elif op == "measure_threshold":
                pq.Q(*ins["modes"]) | pq.ThresholdMeasurement()
            elif op == "measure_homodyne":
                pq.Q(ins["mode"]) | pq.HomodyneMeasurement(phi=_num(ins["phi"]))
            else:
                raise ValueError(f"{op} has no Piquasso equivalent here")
    return prog, n


def from_piquasso(program, d: int) -> Circuit:
    """a pqcf.Circuit from a piquasso Program on d modes (its gates and measurements; a preparation of a number
    state becomes the Fock input)"""
    import piquasso as pq  # noqa: F401
    np = _np()
    out = Circuit(d)
    for ins in program.instructions:
        name = type(ins).__name__
        modes = list(ins.modes)
        p = ins.params
        if name == "Beamsplitter":
            out.beamsplitter(modes[0], modes[1], float(p["theta"]), float(p["phi"]))
        elif name == "Phaseshifter":
            out.phaseshift(modes[0], float(p["phi"]))
        elif name in ("Squeezing", "Displacement"):
            getattr(out, name.lower())(modes[0], float(p["r"]), float(p["phi"]))
        elif name == "Kerr":
            out.kerr(modes[0], float(p["xi"]))
        elif name == "Interferometer":
            out.unitary(modes, np.array(p["matrix"]))
        elif name == "ParticleNumberMeasurement":
            out.measure_pnr(modes)
        elif name == "ThresholdMeasurement":
            out.measure_threshold(modes)
        elif name == "HomodyneMeasurement":
            out.measure_homodyne(modes[0], float(p.get("phi", 0.0)))
        elif name in ("StateVector", "NumberState"):
            occ = [int(x) for x in np.array(p.get("occupation_numbers", p.get("occupation_number", ())))]
            if len(occ) != d:
                raise ValueError(f"a number-state preparation on all {d} modes, got {occ}")
            out.input = {"kind": "fock", "occupation": occ}
        else:
            raise ValueError(f"{name} has no pqcf equivalent here")
    return out
