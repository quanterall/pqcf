"""The client of PQCF's JSON service (pq_service): run, sample, exact probabilities, evolve and gradients from
Python. Two transports:

- stdio (the default): the client starts `erl ... -s pq_service stdio` as a subprocess and talks to it over its
  standard input and output, one JSON document per line. Nothing listens on the network.
- HTTP: `connect(url="http://127.0.0.1:PORT")` to a service started with `-s pq_service http`.

Only the Python standard library is needed; numpy is optional (pqcf.convert, pqcf.linear).
"""

from __future__ import annotations

import glob
import itertools
import json
import os
import shlex
import subprocess
import threading
import urllib.request
from dataclasses import dataclass, field
from typing import Any, Sequence

from .circuit import Circuit


class PqcfError(Exception):
    """An error returned by the service: `category` is PQCF's error category (FS 10.2), `reason` the reason as
    JSON, `detail` its printed form."""

    def __init__(self, category: str, detail: str, reason: Any = None):
        super().__init__(f"{category}: {detail}")
        self.category = category
        self.detail = detail
        self.reason = reason


@dataclass
class Result:
    """The samples of a run (FS 4.5): one row per accepted shot, columns in `measured` order."""
    kind: str
    measured: list[int]
    samples: list[list]
    metadata: dict = field(default_factory=dict)


def _home() -> str:
    return os.environ.get("PQCF_HOME", os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")))


def default_command(transport: str = "stdio", home: str | None = None) -> list[str]:
    """`erl` with the compiled PQCF on its code path (PQCF_HOME, else the checkout this package sits in);
    PQCF_COMMAND overrides it entirely"""
    if os.environ.get("PQCF_COMMAND") and transport == "stdio":
        return shlex.split(os.environ["PQCF_COMMAND"])
    home = home or _home()
    paths = sorted(glob.glob(os.path.join(home, "_build", "default", "lib", "*", "ebin")))
    if not paths:
        raise PqcfError("unavailable", f"no compiled PQCF under {home} (run rebar3 compile, or set PQCF_HOME)")
    cmd = [os.environ.get("PQCF_ERL", "erl"), "-noshell"]
    for p in paths:
        cmd += ["-pa", p]
    return cmd + ["-s", "pq_service", transport]


class _Stdio:
    def __init__(self, command: Sequence[str], timeout: float):
        self.proc = subprocess.Popen(list(command), stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                     stderr=subprocess.DEVNULL if os.environ.get("PQCF_QUIET") else None)
        self.lock = threading.Lock()
        line = self.proc.stdout.readline()
        if not line:
            raise PqcfError("unavailable", f"the service did not start: {' '.join(command[:2])} ...")
        if json.loads(line).get("ready") != 1:
            raise PqcfError("unavailable", f"unexpected first line from the service: {line[:200]!r}")

    def call(self, req: dict) -> dict:
        data = (json.dumps(req) + "\n").encode()
        with self.lock:
            if self.proc.poll() is not None:
                raise PqcfError("unavailable", "the service has stopped")
            self.proc.stdin.write(data)
            self.proc.stdin.flush()
            line = self.proc.stdout.readline()
        if not line:
            raise PqcfError("unavailable", "the service closed its output")
        return json.loads(line)

    def close(self):
        if self.proc.poll() is None:
            self.proc.stdin.close()
            try:
                self.proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                self.proc.kill()
                self.proc.wait()
        self.proc.stdout.close()


class _Http:
    def __init__(self, url: str, timeout: float):
        self.url = url.rstrip("/") + "/api"
        self.timeout = timeout

    def call(self, req: dict) -> dict:
        r = urllib.request.Request(self.url, data=json.dumps(req).encode(), method="POST",
                                   headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(r, timeout=self.timeout) as resp:
            return json.loads(resp.read())

    def close(self):
        pass


def _circuit(c) -> dict:
    if isinstance(c, Circuit):
        return c.to_json()
    if isinstance(c, dict):
        return c
    if isinstance(c, str):
        return json.loads(c)
    raise TypeError("a pqcf.Circuit, a circuit document (dict) or its JSON text")


def _options(opts: dict) -> dict:
    out = {}
    for k, v in opts.items():
        if v is None:
            continue
        if k == "gram":
            v = [[[complex(z).real, complex(z).imag] for z in row] for row in v]
        elif k == "dark_counts" and isinstance(v, dict):
            v = {str(m): x for m, x in v.items()}
        out[k] = v
    return out


class Client:
    """One connection to the service. Use `pqcf.connect()`; close it, or use it as a context manager."""

    def __init__(self, transport):
        self._t = transport
        self._ids = itertools.count(1)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def close(self):
        self._t.close()

    def request(self, op: str, **fields) -> Any:
        """a raw request: the `ok` value, or PqcfError"""
        rid = next(self._ids)
        resp = self._t.call({"id": rid, "op": op, **fields})
        if "error" in resp:
            e = resp["error"]
            raise PqcfError(e.get("category", "unknown"), e.get("detail", ""), e.get("reason"))
        if resp.get("id") != rid:
            raise PqcfError("pq_bug", f"response id {resp.get('id')!r} for request {rid}")
        return resp["ok"]

    def version(self) -> dict:
        return self.request("version")

    def validate(self, circuit) -> dict:
        return self.request("validate", circuit=_circuit(circuit))

    def backend_for(self, circuit, **opts) -> str:
        return self.request("backend_for", circuit=_circuit(circuit), options=_options(opts))

    def run(self, circuit, **opts) -> Result:
        r = self.request("run", circuit=_circuit(circuit), options=_options(opts))
        return Result(r["kind"], r["measured"], r["samples"], r["metadata"])

    def sample(self, circuit, shots: int, **opts) -> Result:
        r = self.request("sample", circuit=_circuit(circuit), options=_options({"shots": shots, **opts}))
        return Result(r["kind"], r["measured"], r["samples"], r["metadata"])

    def probabilities(self, circuit, **opts) -> dict[tuple, float]:
        """exact probabilities of the measured modes' outcomes, keyed by pattern tuples"""
        r = self.request("probabilities", circuit=_circuit(circuit), options=_options(opts))
        return {tuple(p): v for p, v in r["probabilities"]}

    def probabilities_with_metadata(self, circuit, **opts) -> tuple[dict[tuple, float], dict]:
        r = self.request("probabilities", circuit=_circuit(circuit), options=_options(opts))
        return {tuple(p): v for p, v in r["probabilities"]}, r["metadata"]

    def evolve(self, circuit, **opts) -> dict:
        """the output state: {"kind": "gaussian", "mean", "cov"} | {"kind": "fock_state", "amps": [[re, im], ...]}
        | {"kind": "fock_dm", "rho": ...}"""
        return self.request("evolve", circuit=_circuit(circuit), options=_options(opts))["state"]

    def gradient(self, circuit, observable: dict, binding: dict, **opts) -> dict[str, float]:
        """d observable / d symbol at `binding`; mode = "numeric" (default) | "symbolic" | "adjoint" |
        "parameter_shift"; observables from pqcf.probability, mean_photon_number, weighted_sum"""
        r = self.request("gradient", circuit=_circuit(circuit), observable=observable, binding=binding,
                         options=_options(opts))
        return r["gradient"]


def connect(url: str | None = None, *, command: Sequence[str] | None = None, home: str | None = None,
            timeout: float = 600.0) -> Client:
    """a client: over HTTP when `url` is given, else a stdio subprocess (`command`, or `erl` with PQCF_HOME's build)"""
    if url is not None:
        return Client(_Http(url, timeout))
    return Client(_Stdio(command or default_command("stdio", home), timeout))
