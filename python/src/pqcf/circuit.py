"""Circuits in PQCF's JSON interchange format (FS Appendix B, version 1): a small builder whose `to_json()` is the
document the service reads. Parameters are numbers or expressions in the circuit's symbols, written as text in
PQCF's expression syntax ("pi/4", "2*theta", "sqrt(eta)")."""

from __future__ import annotations

import copy
import json
from typing import Any, Iterable, Sequence

Param = "float | int | str"


def _param(p):
    if isinstance(p, bool) or not isinstance(p, (int, float, str)):
        raise TypeError(f"a parameter is a number or an expression string, not {p!r}")
    return p


def _modes(ms: Iterable[int]) -> list[int]:
    out = [int(m) for m in ms]
    if any(m < 0 for m in out):
        raise ValueError("modes are non-negative")
    return out


def _complex(z) -> dict:
    z = complex(z)
    return {"re": z.real, "im": z.imag}


class Circuit:
    """A circuit of `n_modes` modes. The input is the vacuum, a Fock pattern (`fock=[...]`) or a Gaussian state
    (`gaussian=(mean, cov)`, xxpp order). Every gate method returns the circuit, so calls chain."""

    def __init__(self, n_modes: int, *, fock: Sequence[int] | None = None,
                 gaussian: tuple[Sequence[float], Sequence[Sequence[float]]] | None = None,
                 hbar: float | None = None):
        if n_modes < 1:
            raise ValueError("at least one mode")
        if fock is not None and gaussian is not None:
            raise ValueError("one input: fock or gaussian")
        self.n_modes = int(n_modes)
        if fock is not None:
            if len(fock) != self.n_modes:
                raise ValueError("one occupation per mode")
            self.input: dict = {"kind": "fock", "occupation": [int(x) for x in fock]}
        elif gaussian is not None:
            mean, cov = gaussian
            self.input = {"kind": "gaussian", "mean": [float(x) for x in mean],
                          "cov": [float(x) for row in cov for x in row],
                          "hbar": float(hbar if hbar is not None else 2)}
        else:
            self.input = {"kind": "vacuum"}
        self.hbar = hbar
        self.instructions: list[dict] = []
        self.symbols: dict[str, str] = {}
        self.params: dict[str, Any] = {}

    # declarations
    def symbol(self, name: str, domain: str = "real") -> "Circuit":
        """declares a symbol's domain: real, positive, nonnegative or complex"""
        if domain not in ("real", "positive", "nonnegative", "complex"):
            raise ValueError(f"domain {domain!r}")
        self.symbols[name] = domain
        return self

    def param(self, name: str, default) -> "Circuit":
        """a default value (a number or an expression) for a symbol"""
        self.params[name] = _param(default)
        return self

    def _add(self, d: dict) -> "Circuit":
        self.instructions.append(d)
        return self

    # gates (FS 3.3 conventions)
    def beamsplitter(self, m1: int, m2: int, theta, phi=0) -> "Circuit":
        return self._add({"op": "beamsplitter", "m1": int(m1), "m2": int(m2), "theta": _param(theta), "phi": _param(phi)})

    def phaseshift(self, mode: int, phi) -> "Circuit":
        return self._add({"op": "phaseshift", "mode": int(mode), "phi": _param(phi)})

    def displacement(self, mode: int, r, phi=0) -> "Circuit":
        return self._add({"op": "displacement", "mode": int(mode), "r": _param(r), "phi": _param(phi)})

    def squeezing(self, mode: int, r, phi=0) -> "Circuit":
        return self._add({"op": "squeezing", "mode": int(mode), "r": _param(r), "phi": _param(phi)})

    def kerr(self, mode: int, kappa) -> "Circuit":
        return self._add({"op": "kerr", "mode": int(mode), "kappa": _param(kappa)})

    def loss(self, mode: int, eta) -> "Circuit":
        return self._add({"op": "loss", "mode": int(mode), "eta": _param(eta)})

    def thermal_loss(self, mode: int, eta, nth) -> "Circuit":
        return self._add({"op": "thermal_loss", "mode": int(mode), "eta": _param(eta), "nth": _param(nth)})

    def unitary(self, modes: Sequence[int], matrix) -> "Circuit":
        """a numeric passive unitary on `modes` (rows of complex numbers, or a numpy array)"""
        rows = [[_complex(z) for z in row] for row in matrix]
        k = len(rows)
        if k != len(modes) or any(len(r) != k for r in rows):
            raise ValueError("a k x k matrix for k modes")
        return self._add({"op": "unitary", "modes": _modes(modes), "matrix": [z for r in rows for z in r]})

    # measurements
    def measure_pnr(self, modes: Sequence[int]) -> "Circuit":
        return self._add({"op": "measure_pnr", "modes": _modes(modes)})

    def measure_threshold(self, modes: Sequence[int]) -> "Circuit":
        return self._add({"op": "measure_threshold", "modes": _modes(modes)})

    def measure_homodyne(self, mode: int, phi=0) -> "Circuit":
        return self._add({"op": "measure_homodyne", "mode": int(mode), "phi": _param(phi)})

    # the document
    def to_json(self) -> dict:
        doc = {"pq": 1, "n_modes": self.n_modes, "input": copy.deepcopy(self.input),
               "symbols": dict(self.symbols), "params": dict(self.params),
               "instructions": copy.deepcopy(self.instructions)}
        if self.hbar is not None:
            doc["hbar"] = self.hbar
        return doc

    def dumps(self) -> str:
        return json.dumps(self.to_json())

    @classmethod
    def from_json(cls, doc: dict | str) -> "Circuit":
        if isinstance(doc, str):
            doc = json.loads(doc)
        if doc.get("pq") != 1:
            raise ValueError("a PQCF circuit document, version 1")
        c = cls(int(doc["n_modes"]), hbar=doc.get("hbar"))
        c.input = copy.deepcopy(doc.get("input", {"kind": "vacuum"}))
        c.symbols = dict(doc.get("symbols", {}))
        c.params = dict(doc.get("params", {}))
        c.instructions = copy.deepcopy(doc.get("instructions", []))
        return c

    def __repr__(self) -> str:
        return f"Circuit({self.n_modes} modes, {len(self.instructions)} instructions, input {self.input['kind']})"


# observables of pq_gradient
def probability(modes: Sequence[int], pattern: Sequence[int]) -> dict:
    return {"probability": {"modes": _modes(modes), "pattern": [int(x) for x in pattern]}}


def mean_photon_number(modes: Sequence[int]) -> dict:
    return {"mean_photon_number": {"modes": _modes(modes)}}


def weighted_sum(terms: Iterable[tuple[Any, dict]]) -> dict:
    """terms (weight, observable); a weight is a number or an expression in the circuit's symbols"""
    return {"weighted_sum": [[_param(w), o] for w, o in terms]}
