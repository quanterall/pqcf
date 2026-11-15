"""The superradiant burst exactly, independently of PIQS and PQCF: N atoms, collective emission only, from all excited.
The populations of the symmetric ladder j = N/2 obey dp_m/dt = -gamma a_m p_m + gamma a_(m+1) p_(m+1),
a_m = (j + m)(j - m + 1); p(t) = expm(A t) p(0), by mpmath at 40 digits. <J+ J-> = sum_m a_m p_m.
Writes test/fixtures/superradiance_mp.json.
    docker run --rm -v "$PWD":/w -w /w pqcf-oracle python comparison/dicke/superradiance_mp.py
"""
import json

import mpmath as mp

mp.mp.dps = 40
N, gamma = 100, 1
j = mp.mpf(N) / 2
ms = [j - k for k in range(N + 1)]                  # m = j .. -j
a = [(j + m) * (j - m + 1) for m in ms]
A = mp.zeros(N + 1, N + 1)
for k in range(N + 1):
    A[k, k] = -gamma * a[k]
    if k > 0:
        A[k, k - 1] = gamma * a[k - 1]              # from m + 1 (index k - 1) down to m
p0 = mp.zeros(N + 1, 1)
p0[0] = 1
ts = ["0.01", "0.03", "0.046", "0.05", "0.08", "0.12"]
out = {"N": N, "t": ts, "jpjm": [], "jz": []}
for t in ts:
    p = mp.expm(A * mp.mpf(t)) * p0
    out["jpjm"].append(mp.nstr(sum(a[k] * p[k] for k in range(N + 1)), 35))
    out["jz"].append(mp.nstr(sum(ms[k] * p[k] for k in range(N + 1)), 35))
json.dump(out, open("test/fixtures/superradiance_mp.json", "w"), indent=1)
print(out["jpjm"][:3])
