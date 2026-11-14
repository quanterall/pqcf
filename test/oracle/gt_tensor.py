"""GT-TEN (decision record SY2, phase 4): canonical forms of Riemann monomials by SymPy 1.14's Butler-Portugal
implementation (sympy.tensor.tensor, canon_bp), for pq_tensor's tests: random products of 2 and 3 Riemann tensors,
fully contracted or with two free indices. Writes test/fixtures/gt_tensor.json.

    docker run --rm -u $(id -u):$(id -g) -v "$PWD":/w -w /w --entrypoint python pqcf-oracle test/oracle/gt_tensor.py
"""
import json
import pathlib
import random

from sympy.tensor.tensor import TensorHead, TensorIndexType, TensorSymmetry, tensor_indices

random.seed(20261114)
L = TensorIndexType("L", dummy_name="L")
R = TensorHead("R", [L] * 4, TensorSymmetry.riemann())
names = [f"i{k}" for k in range(1, 13)] + ["a", "b"]
idx = dict(zip(names, tensor_indices(" ".join(names), L)))

cases = []
for n in range(160):
    k = random.choice([2, 2, 3])
    slots = list(range(4 * k))
    random.shuffle(slots)
    labels = [None] * (4 * k)
    free = random.choice([0, 0, 2])
    for s, f in zip(slots[:free], ["a", "b"]):
        labels[s] = f
    rest = slots[free:]
    for j in range(0, len(rest), 2):
        labels[rest[j]] = labels[rest[j + 1]] = f"i{j // 2 + 1}"
    factors = [labels[4 * t:4 * t + 4] for t in range(k)]
    seen = set()
    expr = 1
    for fac in factors:
        args = []
        for name in fac:
            if name in ("a", "b") or name not in seen:
                args.append(idx[name])
                seen.add(name)
            else:
                args.append(-idx[name])
        expr = expr * R(*args)
    canon = expr.canon_bp()
    cases.append({"factors": factors, "sympy": str(canon)})

out = pathlib.Path(__file__).resolve().parents[1] / "fixtures" / "gt_tensor.json"
import sympy
out.write_text(json.dumps({"sympy": sympy.__version__, "cases": cases}, indent=1) + "\n")
print(f"{len(cases)} monomials -> {out}; zeros: {sum(c['sympy'] == '0' for c in cases)}; "
      f"classes: {len(set(c['sympy'].lstrip('-') for c in cases))}")
