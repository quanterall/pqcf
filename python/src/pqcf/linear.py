"""The single-photon unitary of a passive circuit in PQCF's conventions (FS 3.3), with numpy: the reference the
converters are tested against. B(theta, phi) = [[cos t, -e^(-i phi) sin t], [e^(i phi) sin t, cos t]] on (m1, m2)
(column j the image of mode j), R(phi) = e^(i phi); the circuit's unitary is G_K ... G_1."""

from __future__ import annotations

import numpy as np

from .circuit import Circuit

PASSIVE = ("beamsplitter", "phaseshift", "unitary")
MEASURE = ("measure_pnr", "measure_threshold", "measure_homodyne")


def _num(p) -> float:
    if isinstance(p, (int, float)) and not isinstance(p, bool):
        return float(p)
    raise ValueError(f"a numeric parameter is needed here, not {p!r} (bind the symbols first)")


def gate_matrix(ins: dict, n: int) -> np.ndarray:
    g = np.eye(n, dtype=complex)
    op = ins["op"]
    if op == "beamsplitter":
        a, b = ins["m1"], ins["m2"]
        t, p = _num(ins["theta"]), _num(ins["phi"])
        c, s = np.cos(t), np.sin(t)
        g[a, a], g[a, b] = c, -np.exp(-1j * p) * s
        g[b, a], g[b, b] = np.exp(1j * p) * s, c
    elif op == "phaseshift":
        g[ins["mode"], ins["mode"]] = np.exp(1j * _num(ins["phi"]))
    elif op == "unitary":
        ms = ins["modes"]
        k = len(ms)
        u = np.array([complex(z["re"], z["im"]) for z in ins["matrix"]]).reshape(k, k)
        for i, mi in enumerate(ms):
            for j, mj in enumerate(ms):
                g[mi, mj] = u[i, j]
    else:
        raise ValueError(f"{op} is not passive")
    return g


def unitary(circuit) -> np.ndarray:
    """U = G_K ... G_1 of the passive instructions (measurements ignored)"""
    doc = circuit.to_json() if isinstance(circuit, Circuit) else circuit
    n = doc["n_modes"]
    u = np.eye(n, dtype=complex)
    for ins in doc["instructions"]:
        if ins["op"] in MEASURE:
            continue
        u = gate_matrix(ins, n) @ u
    return u
