"""S7 fixtures (FS 8.5, AN 2.1): the draws derived from the stream and the stream order of the boson samplers,
from an implementation that shares no code with the kernels (the ChaCha8 of chacha8_ref.py). The kernels must
reproduce these bit for bit (normals to 1e-15: they go through the math library). Writes test/fixtures/gt_s7.json.

    python3 test/oracle/gt_s7.py
"""
import cmath
import itertools
import json
import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from chacha8_ref import rng_u64  # noqa: E402


class Stream:
    def __init__(self, seed, stream):
        self.seed, self.stream, self.pos, self.buf = seed, stream, 0, []

    def u64(self):
        if self.pos >= len(self.buf):
            self.buf = rng_u64(self.seed, self.stream, len(self.buf) + 64)
        self.pos += 1
        return self.buf[self.pos - 1]

    def uniform(self):
        return (self.u64() >> 11) * 2.0 ** -53

    def below(self, s):                      # Lemire, 64-bit
        m = self.u64() * s
        low = m & (2 ** 64 - 1)
        if low < s:
            t = (2 ** 64 - s) % s
            while low < t:
                m = self.u64() * s
                low = m & (2 ** 64 - 1)
        return m >> 64

    def normal_pair(self):                   # Box-Muller, cosine first, u1 = 0 rejected
        u1 = self.uniform()
        while u1 == 0.0:
            u1 = self.uniform()
        u2 = self.uniform()
        r = math.sqrt(-2.0 * math.log(u1))
        a = 2.0 * math.pi * u2
        return r * math.cos(a), r * math.sin(a)

    def shuffle(self, v):                    # Fisher-Yates, i = len-1 .. 1
        for i in range(len(v) - 1, 0, -1):
            j = self.below(i + 1)
            v[i], v[j] = v[j], v[i]
        return v


def perm(a):
    n = len(a)
    return sum(math.prod(a[i][s[i]] for i in range(n)) for s in itertools.permutations(range(n))) if n else 1


def categorical(w, u):
    total = math.fsum(w)
    cum = 0.0
    for j, x in enumerate(w):
        cum += x
        if cum > u * total:
            return j
    return max(j for j, x in enumerate(w) if x > 0)


def patterns(m, n):
    if m == 1:
        return [[n]]
    return [[k] + rest for k in range(n + 1) for rest in patterns(m - 1, n - k)]


def brute_force(U, inp, shots, seed, start):
    m, n = len(U), sum(inp)
    pats = patterns(m, n)
    def prob(t):
        rows = [j for j in range(m) for _ in range(t[j])]
        cols = [i for i in range(m) for _ in range(inp[i])]
        p = perm([[U[r][c] for c in cols] for r in rows])
        return abs(p) ** 2 / (math.prod(map(math.factorial, inp)) * math.prod(map(math.factorial, t)))
    probs = [prob(t) for t in pats]
    return [pats[categorical(probs, Stream(seed, start + i).uniform())] for i in range(shots)]


def clifford_clifford(U, inp, shots, seed, start):
    m = len(U)
    out = []
    for i in range(shots):
        s = Stream(seed, start + i)
        cols = s.shuffle([c for c in range(m) for _ in range(inp[c])])
        ys = []
        for k in range(1, len(cols) + 1):
            p = [perm([[U[y][cols[c]] for c in range(k) if c != l] for y in ys]) for l in range(k)]
            w = [abs(sum(U[j][cols[l]] * p[l] for l in range(k))) ** 2 for j in range(m)]
            ys.append(categorical(w, s.uniform()))
        out.append([ys.count(j) for j in range(m)])
    return out


# a fixed 3-mode unitary: the discrete Fourier transform times a phase on mode 1
w3 = cmath.exp(2j * math.pi / 3)
U = [[w3 ** (j * k) / math.sqrt(3) * (cmath.exp(0.4j) if j == 1 else 1) for k in range(3)] for j in range(3)]

fx = {"test": "S7", "generator": "test/oracle/gt_s7.py",
      "draws": [], "unitary": [[[z.real, z.imag] for z in row] for row in U], "samplers": []}
for seed, stream in [(0, 0), (42, 7), (2 ** 64 - 1, 3)]:
    st = Stream(seed, stream)
    rec = {"seed": seed, "stream": stream,
           "uniform": [st.uniform() for _ in range(16)],
           "below": [[s, [st.below(s) for _ in range(8)]] for s in (1, 2, 3, 7, 10, 1000, 2 ** 32 + 1, 2 ** 63 + 5, 2 ** 64 - 1)],
           "normal": [list(st.normal_pair()) for _ in range(8)],
           "shuffle": st.shuffle(list(range(10)))}
    fx["draws"].append(rec)
for inp in ([1, 1, 1], [2, 1, 0], [0, 3, 0]):
    fx["samplers"].append({"input": inp, "seed": 11, "stream_start": 4, "shots": 40,
                           "brute_force": brute_force(U, inp, 40, 11, 4),
                           "clifford_clifford": clifford_clifford(U, inp, 40, 11, 4)})
path = pathlib.Path(__file__).resolve().parents[1] / "fixtures" / "gt_s7.json"
path.write_text(json.dumps(fx, indent=1) + "\n")
print("wrote", path)
