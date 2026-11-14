"""GT-QOP (quantum optics plan, Q1): normal-ordered forms of random words of two bosonic modes by SymPy 1.14
(sympy.physics.quantum, normal_ordered_form, the modes independent), for pq_qop's tests. Writes
test/fixtures/gt_qop.json.

    docker run --rm -u $(id -u):$(id -g) -v "$PWD":/w -w /w --entrypoint python pqcf-oracle test/oracle/gt_qop.py
"""
import json
import pathlib
import random

import sympy
from sympy import Mul, Pow, Add, Integer
from sympy.physics.quantum import Dagger
from sympy.physics.quantum.boson import BosonOp
from sympy.physics.quantum.operatorordering import normal_ordered_form

random.seed(20261115)
a, b = BosonOp("a"), BosonOp("b")
letters = {"a": a, "ad": Dagger(a), "b": b, "bd": Dagger(b)}


def powers(term):
    """{(name, dagger): power} and the numeric coefficient of a normal-ordered term"""
    coef, rest = term.as_coeff_Mul()
    pw = {"ad": 0, "a": 0, "bd": 0, "b": 0}
    for f in Mul.make_args(rest):
        base, e = f.as_base_exp()
        if base == Integer(1):
            continue
        key = {Dagger(a): "ad", a: "a", Dagger(b): "bd", b: "b"}[base]
        pw[key] += int(e)
    return int(coef), pw


cases = []
for n in range(60):
    word = [random.choice(list(letters)) for _ in range(random.randint(2, 7))]
    expr = Mul(*[letters[w] for w in word])
    no = sympy.expand(normal_ordered_form(expr, independent=True))
    terms = []
    for t in Add.make_args(no):
        c, pw = powers(t)
        terms.append([c, pw["ad"], pw["a"], pw["bd"], pw["b"]])
    cases.append({"word": word, "terms": sorted(terms)})

out = pathlib.Path(__file__).resolve().parents[1] / "fixtures" / "gt_qop.json"
out.write_text(json.dumps({"sympy": sympy.__version__, "cases": cases}, indent=1) + "\n")
print(f"{len(cases)} words -> {out}")
