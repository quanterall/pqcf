"""Spike S2: both crates' results against numpy (LAPACK), relative errors, condition numbers, timings."""
import json
import numpy as np

mats = {m["n"]: m for m in json.load(open("matrices.json"))}
res = json.load(open("results.json"))


def c(a):
    a = np.array(a)
    return a[..., 0] + 1j * a[..., 1]


rows = []
for r in res:
    n = r["n"]
    H, A, G = (c(mats[n][k]).reshape(n, n) for k in "HAG")
    ev = np.linalg.eigvalsh(H)
    L = np.linalg.cholesky(A)
    _, ld = np.linalg.slogdet(G)
    X = np.linalg.inv(G)
    kappa = np.linalg.cond(G)
    for name in ("nalgebra", "faer"):
        o = r[name]
        e_eig = np.max(np.abs(np.array(o["eig"]) - ev)) / np.max(np.abs(ev))
        e_res = o["eig_residual"] / np.max(np.abs(ev))
        e_chol = np.max(np.abs(c(o["chol"]).reshape(n, n) - L)) / np.max(np.abs(L))
        e_ld = abs(o["logdet"] - ld) / max(1.0, abs(ld))
        e_inv = np.max(np.abs(c(o["inv"]).reshape(n, n) - X)) / np.max(np.abs(X))
        rows.append(dict(n=n, crate=name, eig=e_eig, eig_res=e_res, chol=e_chol, logdet=e_ld, inv=e_inv,
                         inv_over_kappa=e_inv / kappa, kappa=kappa, repeat=o["repeat_bitwise"],
                         **{"seconds_" + k: v for k, v in o["seconds"].items()}))

hdr = ("n", "crate", "eig", "eig_res", "chol", "logdet", "inv", "inv/kappa", "kappa", "repeat",
       "t_eig_us", "t_chol_us", "t_lu_us", "t_inv_us")
print(" ".join(f"{h:>9}" for h in hdr))
for x in rows:
    vals = [x["n"], x["crate"], f"{x['eig']:.1e}", f"{x['eig_res']:.1e}", f"{x['chol']:.1e}", f"{x['logdet']:.1e}",
            f"{x['inv']:.1e}", f"{x['inv_over_kappa']:.1e}", f"{x['kappa']:.0f}", x["repeat"],
            f"{1e6 * x['seconds_eig']:.1f}", f"{1e6 * x['seconds_chol']:.1f}", f"{1e6 * x['seconds_lu']:.1f}",
            f"{1e6 * x['seconds_inverse']:.1f}"]
    print(" ".join(f"{str(v):>9}" for v in vals))
json.dump(rows, open("compare.json", "w"), indent=1, default=float)
