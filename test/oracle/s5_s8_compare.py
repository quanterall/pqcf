"""Spikes S5 and S8: the kernels' matrices against closed forms in mpmath at 80 digits.

  displacement <m|D(alpha)|n> = sqrt(n!/m!) alpha^(m-n) e^(-|alpha|^2/2) L_n^(m-n)(|alpha|^2)       (m >= n)
                              = sqrt(m!/n!) (-conj alpha)^(n-m) e^(-|alpha|^2/2) L_m^(n-m)(|alpha|^2) (m < n)
  squeezing    S = exp(-t a^dag^2 / 2) sech(r)^(n + 1/2) exp(conj(t) a^2 / 2), t = e^(i phi) tanh r, so
               <m|S|n> = sqrt(m! n! sech r) sum_k sech(r)^k / k! (-t/2)^i / i! (conj(t)/2)^j / j!,
               i = (m - k)/2, j = (n - k)/2 non-negative integers
  beamsplitter the sector matrix exp(theta (e^(i phi) a b^dag - e^(-i phi) a^dag b)) by mpmath.expm

    docker run --rm -v "$PWD":/w -w /w pqcf-oracle python test/oracle/s5_s8_compare.py s5_s8.json
"""
import json
import sys

import mpmath as mp

mp.mp.dps = 80
data = json.load(open(sys.argv[1]))


def cm(rows, d):
    return [[complex(*rows[i * d + j]) for j in range(d)] for i in range(d)]


def worst(got, ref):
    return max(abs(complex(got[i][j]) - complex(ref[i][j])) for i in range(len(ref)) for j in range(len(ref)))


def lag(n, k, x):
    # the generalised Laguerre polynomial as its finite sum (exact at this precision)
    return mp.fsum((-1) ** i * mp.binomial(n + k, n - i) * x ** i / mp.factorial(i) for i in range(n + 1))


CACHE = {}


def disp_ref(alpha, c):
    key = ("d", c, str(alpha))
    if key not in CACHE:
        CACHE[key] = _disp_ref(alpha, c)
    return CACHE[key]


def _disp_ref(alpha, c):
    x = abs(alpha) ** 2
    e = mp.e ** (-x / 2)
    R = [[0] * (c + 1) for _ in range(c + 1)]
    for m in range(c + 1):
        for n in range(c + 1):
            if m >= n:
                R[m][n] = mp.sqrt(mp.factorial(n) / mp.factorial(m)) * alpha ** (m - n) * e * lag(n, m - n, x)
            else:
                R[m][n] = mp.sqrt(mp.factorial(m) / mp.factorial(n)) * (-mp.conj(alpha)) ** (n - m) * e * lag(m, n - m, x)
    return R


def sq_ref(r, phi, c):
    key = ("s", c, str(r), str(phi))
    if key not in CACHE:
        CACHE[key] = _sq_ref(r, phi, c)
    return CACHE[key]


def _sq_ref(r, phi, c):
    t = mp.e ** (1j * phi) * mp.tanh(r)
    sech = 1 / mp.cosh(r)
    R = [[0] * (c + 1) for _ in range(c + 1)]
    f = [mp.factorial(k) for k in range(c + 1)]
    for m in range(c + 1):
        for n in range(c + 1):
            if (m - n) % 2:
                continue
            s = 0
            for k in range(m % 2, min(m, n) + 1, 2):
                i, j = (m - k) // 2, (n - k) // 2
                s += sech ** k / f[k] * (-t / 2) ** i / f[i] * (mp.conj(t) / 2) ** j / f[j]
            R[m][n] = mp.sqrt(f[m] * f[n] * sech) * s
    return R


def bs_ref(n, theta, phi):
    d = n + 1
    G = mp.zeros(d, d)
    for j in range(d):           # basis j = n_a, n_b = N - j
        if j >= 1:               # a b^dag: |j> -> sqrt(j (N - j + 1)) |j - 1>
            G[j - 1, j] += mp.e ** (1j * phi) * mp.sqrt(j * (n - j + 1))
        if j < n:                # a^dag b: |j> -> sqrt((j + 1)(N - j)) |j + 1>
            G[j + 1, j] -= mp.e ** (-1j * phi) * mp.sqrt((j + 1) * (n - j))
    E = mp.expm(theta * G)
    return [[E[i, k] for k in range(d)] for i in range(d)]


bs_worst = 0.0
s5 = {}
for item in data:
    if item["kind"] == "bs":
        d = item["n"] + 1
        bs_worst = max(bs_worst, worst(cm(item["m"], d), bs_ref(item["n"], mp.mpf(item["theta"]), mp.mpf(item["phi"]))))
    else:
        c = item["c"]
        if item["kind"] in ("disp", "disp_lag"):
            alpha = mp.mpf(item["a"]) * mp.e ** (1j * mp.mpf(item["phi"]))
            err = worst(cm(item["m"], c + 1), disp_ref(alpha, c))
            s5[(item["kind"], c, item["a"])] = err
        else:
            err = worst(cm(item["m"], c + 1), sq_ref(mp.mpf(item["r"]), mp.mpf(item["phi"]), c))
            s5[(item["kind"], c, item["r"])] = err
        print(item["kind"], c, item.get("a", item.get("r")), f"{err:.1e}", flush=True)

print(f"S8: largest entry error of the beamsplitter sectors N = 0..20, three angle pairs: {bs_worst:.1e}")
json.dump({"s8_worst": bs_worst, "s5": [[k[0], k[1], k[2], v] for k, v in s5.items()]},
          open("s5_s8_errors.json", "w"), indent=1)
