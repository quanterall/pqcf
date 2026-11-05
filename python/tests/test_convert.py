"""The converters against Perceval and Piquasso themselves: every converted circuit's unitary (or Gaussian moments)
equals what the library computes. Skipped without numpy, perceval or piquasso; they are in the oracle image:
    docker run --rm -v "$PWD":/src -w /src pqcf-oracle python -m unittest discover -s python/tests"""

import math
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

try:
    import numpy as np
    import perceval as pcvl
except ImportError:  # pragma: no cover
    pcvl = None
try:
    import piquasso as pq
except ImportError:  # pragma: no cover
    pq = None

from pqcf import Circuit  # noqa: E402


def random_unitary(k, seed):
    rng = np.random.default_rng(seed)
    z = rng.normal(size=(k, k)) + 1j * rng.normal(size=(k, k))
    q, r = np.linalg.qr(z)
    return q * (np.diag(r) / abs(np.diag(r)))


@unittest.skipIf(pcvl is None, "perceval not installed")
class Perceval(unittest.TestCase):
    def test_from_perceval(self):
        from pqcf.convert import from_perceval
        from pqcf.linear import unitary
        c = (pcvl.Circuit(4) // (0, pcvl.BS.Rx(theta=0.7, phi_tl=0.1, phi_br=0.3)) // (1, pcvl.PS(0.4))
             // (2, pcvl.BS.H(theta=1.1)) // (1, pcvl.Unitary(pcvl.Matrix(random_unitary(3, 1)))))
        q = from_perceval(c, pcvl.BasicState([1, 0, 1, 0]))
        self.assertTrue(np.allclose(unitary(q), np.array(c.compute_unitary()), atol=1e-12))
        self.assertEqual(q.input, {"kind": "fock", "occupation": [1, 0, 1, 0]})

    def test_to_perceval(self):
        from pqcf.convert import to_perceval
        from pqcf.linear import unitary
        q = (Circuit(4).beamsplitter(0, 1, 0.45, 1.2).phaseshift(2, 0.3).beamsplitter(3, 1, 0.8, -0.4)
             .unitary([1, 2, 3], random_unitary(3, 2)).beamsplitter(2, 3, 0.2, 0.0).measure_pnr([0, 1, 2, 3]))
        c = to_perceval(q)
        self.assertTrue(np.allclose(np.array(c.compute_unitary()), unitary(q), atol=1e-12))

    def test_round_trip(self):
        from pqcf.convert import from_perceval, to_perceval
        from pqcf.linear import unitary
        q = Circuit(3).beamsplitter(0, 1, 0.3, 0.7).beamsplitter(1, 2, 1.0, 0.0).phaseshift(0, 2.0)
        self.assertTrue(np.allclose(unitary(from_perceval(to_perceval(q))), unitary(q), atol=1e-12))

    def test_refusals(self):
        from pqcf.convert import to_perceval
        with self.assertRaises(ValueError):
            to_perceval(Circuit(1).squeezing(0, 0.3))
        with self.assertRaises(ValueError):
            to_perceval(Circuit(1).phaseshift(0, "a"))


@unittest.skipIf(pq is None, "piquasso not installed")
class Piquasso(unittest.TestCase):
    def test_round_trip(self):
        from pqcf.convert import from_piquasso, to_piquasso
        q = (Circuit(2).squeezing(0, 0.3, 0.2).displacement(1, 0.5, 0.4).beamsplitter(0, 1, 0.6, 0.1)
             .phaseshift(1, 0.7).kerr(0, 0.2).measure_pnr([0, 1]))
        prog, d = to_piquasso(q)
        self.assertEqual(from_piquasso(prog, d).to_json()["instructions"], q.to_json()["instructions"])

    def test_gaussian_moments(self):
        """a Gaussian circuit converted to Piquasso: its mean and covariance are ours (FS 3.3, hbar = 2) with
        Piquasso's covariance twice ours"""
        from pqcf.convert import to_piquasso
        r, p, a, ap, t, tp = 0.3, 0.2, 0.5, 0.4, 0.6, 0.1
        q = Circuit(2).squeezing(0, r, p).displacement(1, a, ap).beamsplitter(0, 1, t, tp)
        prog, d = to_piquasso(q)
        st = pq.GaussianSimulator(d=d, config=pq.Config(hbar=2)).execute(prog).state
        # ours: squeezing S on mode 0, displacement of mode 1 by alpha, then the beamsplitter's symplectic map
        ch, sh = math.cosh(r), math.sinh(r)
        s0 = np.array([[ch - math.cos(p) * sh, -math.sin(p) * sh], [-math.sin(p) * sh, ch + math.cos(p) * sh]])
        cov = np.eye(4)  # xpxp: x0 p0 x1 p1
        cov[0:2, 0:2] = s0 @ s0.T
        mean = np.array([0, 0, 2 * a * math.cos(ap), 2 * a * math.sin(ap)])
        b = np.array([[math.cos(t), -np.exp(-1j * tp) * math.sin(t)], [np.exp(1j * tp) * math.sin(t), math.cos(t)]])
        # the passive map on (x, p) pairs: [[Re B, -Im B], [Im B, Re B]] in xxpp, reordered to xpxp
        big = np.block([[b.real, -b.imag], [b.imag, b.real]])
        perm = [0, 2, 1, 3]
        m = big[np.ix_(perm, perm)]
        self.assertTrue(np.allclose(st.xpxp_mean_vector, m @ mean, atol=1e-12))
        self.assertTrue(np.allclose(st.xpxp_covariance_matrix, 2 * (m @ cov @ m.T), atol=1e-12))


if __name__ == "__main__":
    unittest.main()
