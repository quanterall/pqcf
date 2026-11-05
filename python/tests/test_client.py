"""The Python client against a running service (stdio and HTTP): closed forms, the golden file GT-P01 (thewalrus
0.22.0, test/fixtures/gt_p01.json), gradients, the noise options, errors, and the same samples over both
transports. Standard library only:  python3 -m unittest discover -s python/tests"""

import json
import math
import os
import subprocess
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "python", "src"))
os.environ.setdefault("PQCF_HOME", ROOT)
os.environ.setdefault("PQCF_QUIET", "1")

import pqcf  # noqa: E402
from pqcf import Circuit, PqcfError  # noqa: E402


def hom():
    return Circuit(2, fock=[1, 1]).beamsplitter(0, 1, "pi/4", 0).measure_pnr([0, 1])


class Stdio(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.q = pqcf.connect()

    @classmethod
    def tearDownClass(cls):
        cls.q.close()

    def test_version(self):
        v = self.q.version()
        self.assertEqual(v["service"], 1)
        self.assertTrue(v["pq"].startswith("1."))

    def test_hom(self):
        p = self.q.probabilities(hom())
        self.assertAlmostEqual(p[(2, 0)], 0.5, places=14)
        self.assertAlmostEqual(p[(0, 2)], 0.5, places=14)
        self.assertLess(p.get((1, 1), 0.0), 1e-15)

    def test_samples_reproducible(self):
        a = self.q.sample(hom(), shots=200, seed=7)
        b = self.q.sample(hom(), shots=200, seed=7)
        self.assertEqual(a.samples, b.samples)
        self.assertEqual(a.measured, [0, 1])
        self.assertEqual(a.kind, "pnr")
        self.assertNotIn([1, 1], a.samples)
        self.assertEqual(a.metadata["seed"], 7)

    def test_gt_p01(self):
        """every boson-sampling probability of the golden file small enough to enumerate, through the client"""
        with open(os.path.join(ROOT, "test", "fixtures", "gt_p01.json")) as f:
            cases = json.load(f)["probabilities"]
        done = 0
        for c in cases:
            m = len(c["unitary"])
            n = sum(c["input"])
            if math.comb(n + m - 1, n) > 20000:
                continue
            u = [[complex(re, im) for re, im in row] for row in c["unitary"]]
            circ = Circuit(m, fock=c["input"]).unitary(list(range(m)), u).measure_pnr(list(range(m)))
            p = self.q.probabilities(circ, backend="boson")
            got, want = p.get(tuple(c["output"]), 0.0), c["probability"]
            self.assertLessEqual(abs(got - want), 1e-12 + 1e-9 * abs(want), (m, c["input"], c["output"]))
            done += 1
        self.assertGreater(done, 0)

    def test_gradients(self):
        c = Circuit(2, fock=[1, 0]).beamsplitter(0, 1, "t", 0)
        obs = pqcf.probability([0], [1])
        want = -math.sin(2 * 0.3)
        for mode in ("numeric", "symbolic", "adjoint", "parameter_shift"):
            g = self.q.gradient(c, obs, {"t": 0.3}, mode=mode)
            self.assertAlmostEqual(g["t"], want, places=6 if mode == "numeric" else 12, msg=mode)
        w = pqcf.weighted_sum([("t*t", obs), (2, pqcf.mean_photon_number([1]))])
        g = self.q.gradient(c, w, {"t": 0.3}, mode="adjoint")
        # d/dt [t^2 cos^2 t + 2 sin^2 t]
        self.assertAlmostEqual(g["t"], 2 * 0.3 * math.cos(0.3) ** 2 - 0.09 * math.sin(0.6) + 2 * math.sin(0.6), places=12)

    def test_symbols_and_bindings(self):
        c = Circuit(2, fock=[1, 0]).symbol("t").param("t", 0.0).beamsplitter(0, 1, "t").measure_pnr([0])
        self.assertAlmostEqual(self.q.probabilities(c)[(1,)], 1.0, places=14)
        p = self.q.probabilities(c, bindings={"t": math.pi / 3})
        self.assertAlmostEqual(p[(1,)], 0.25, places=14)
        self.assertEqual(self.q.validate(Circuit(1).phaseshift(0, "a"))["free_symbols"], ["a"])

    def test_gaussian_state(self):
        r = 0.4
        s = self.q.evolve(Circuit(1).squeezing(0, r, 0))
        self.assertEqual(s["kind"], "gaussian")
        self.assertAlmostEqual(s["cov"][0][0], math.exp(-2 * r), places=13)
        self.assertAlmostEqual(s["cov"][1][1], math.exp(2 * r), places=13)
        self.assertEqual(self.q.backend_for(Circuit(1).squeezing(0, r).measure_homodyne(0)), "gaussian")

    def test_noise_options(self):
        p = self.q.probabilities(hom(), visibility=0.6)
        self.assertAlmostEqual(p[(1, 1)], 0.2, places=14)
        g = [[1, 0.6], [0.6, 1]]
        self.assertAlmostEqual(self.q.probabilities(hom(), gram=g)[(1, 1)], (1 - 0.36) / 2, places=14)
        thr = Circuit(1, fock=[0]).measure_threshold([0])
        self.assertAlmostEqual(self.q.probabilities(thr, dark_counts=0.2, cutoff=1)[(1,)], 1 - math.exp(-0.2), places=14)
        self.assertAlmostEqual(self.q.probabilities(thr, dark_counts={0: 0.2}, cutoff=1)[(1,)], 1 - math.exp(-0.2), places=14)
        lossy = Circuit(1, fock=[2]).loss(0, 0.5).measure_pnr([0])
        self.assertEqual(self.q.probabilities(lossy, dm_storage="blocks"), self.q.probabilities(lossy))

    def test_postselection(self):
        c = Circuit(2, fock=[1, 1]).beamsplitter(0, 1, 0.3).measure_pnr([0, 1])
        r = self.q.sample(c, shots=500, seed=3, postselect={"op": "eq", "mode": 0, "value": 1})
        self.assertTrue(all(row[0] == 1 for row in r.samples))
        self.assertLess(r.metadata["shots_accepted"], 500)

    def test_errors(self):
        with self.assertRaises(PqcfError) as e:
            self.q.probabilities(Circuit(1, fock=[1]).phaseshift(0, "q").measure_pnr([0]))
        self.assertEqual(e.exception.category, "unbound_symbols")
        with self.assertRaises(PqcfError) as e:
            self.q.probabilities(hom(), not_an_option=1)
        self.assertEqual(e.exception.category, "bad_term")
        with self.assertRaises(PqcfError) as e:
            self.q.probabilities(Circuit(1, fock=[1]).loss(0, 1.5).measure_pnr([0]))
        self.assertEqual(e.exception.category, "domain_violation")
        with self.assertRaises(PqcfError) as e:
            self.q.request("no_such_op")
        self.assertEqual(e.exception.category, "bad_term")
        # the connection still works after errors
        self.assertAlmostEqual(self.q.probabilities(hom())[(2, 0)], 0.5, places=14)


class Http(unittest.TestCase):
    """the HTTP transport, on a free port on 127.0.0.1: the same answers and the same samples as stdio"""

    @classmethod
    def setUpClass(cls):
        env = dict(os.environ, PQ_HTTP_PORT="0")
        cls.proc = subprocess.Popen(pqcf.default_command("http"), stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, env=env)
        port = json.loads(cls.proc.stdout.readline())["listening"]
        cls.q = pqcf.connect(f"http://127.0.0.1:{port}")

    @classmethod
    def tearDownClass(cls):
        cls.proc.kill()
        cls.proc.wait()
        cls.proc.stdout.close()

    def test_same_as_stdio(self):
        with pqcf.connect() as s:
            self.assertEqual(self.q.probabilities(hom()), s.probabilities(hom()))
            a = self.q.sample(hom(), shots=300, seed=11)
            b = s.sample(hom(), shots=300, seed=11)
            self.assertEqual(a.samples, b.samples)

    def test_errors_over_http(self):
        with self.assertRaises(PqcfError) as e:
            self.q.probabilities(hom(), backend="nope")
        self.assertEqual(e.exception.category, "bad_term")

    def test_only_post_api(self):
        import urllib.error
        import urllib.request
        url = self.q._t.url
        with self.assertRaises(urllib.error.HTTPError) as e:
            urllib.request.urlopen(url, timeout=10)
        self.assertEqual(e.exception.code, 405)
        with self.assertRaises(urllib.error.HTTPError) as e:
            urllib.request.urlopen(url.replace("/api", "/other"), timeout=10)
        self.assertEqual(e.exception.code, 404)


class Builder(unittest.TestCase):
    def test_round_trip(self):
        c = hom().symbol("a", "positive").param("a", "1/2").kerr(0, "a")
        d = Circuit.from_json(c.dumps())
        self.assertEqual(d.to_json(), c.to_json())

    def test_shape_checks(self):
        with self.assertRaises(ValueError):
            Circuit(2, fock=[1])
        with self.assertRaises(ValueError):
            Circuit(2).unitary([0, 1], [[1, 0]])
        with self.assertRaises(TypeError):
            Circuit(1).phaseshift(0, None)


if __name__ == "__main__":
    unittest.main()
