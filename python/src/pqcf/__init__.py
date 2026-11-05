"""pqcf: the Python client of PQCF, the Photonic Quantum Computing Framework (github.com/quanterall/pqcf).

    import pqcf
    c = pqcf.Circuit(2, fock=[1, 1]).beamsplitter(0, 1, "pi/4").measure_pnr([0, 1])
    with pqcf.connect() as q:
        q.probabilities(c)          # {(0, 2): 0.5, (2, 0): 0.5, ...}

The simulation runs in PQCF (Erlang/OTP and Rust); this package builds circuits in its JSON format, talks to its
JSON service (pq_service) and converts to and from Perceval and Piquasso.
"""

from .circuit import Circuit, mean_photon_number, probability, weighted_sum
from .client import Client, PqcfError, Result, connect, default_command

__version__ = "1.2.0"
__all__ = ["Circuit", "Client", "PqcfError", "Result", "connect", "default_command", "probability",
           "mean_photon_number", "weighted_sum", "__version__"]
