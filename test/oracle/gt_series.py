"""GT-SER (decision record SY2, phase 2): series expansions by SymPy 1.14, their coefficients evaluated by mpmath to
30 digits, for pq_series's tests. The parameters are set to fixed rationals before expanding; pq_series keeps them
symbolic and its coefficients are evaluated at the same values. Writes test/fixtures/gt_series.json.

    docker run --rm -u $(id -u):$(id -g) -v "$PWD":/w -w /w --entrypoint python pqcf-oracle test/oracle/gt_series.py
"""
import json
import pathlib

import sympy as sp

x = sp.Symbol("x", positive=True)
a = sp.Symbol("a")
b = sp.Symbol("b")
PARAMS = {a: sp.Rational(3, 7), b: sp.Rational(-5, 4)}
R = sp.Rational

# (name, Erlang term of pq_series, SymPy expression, terms below x^n, expansion point)
CASES = [
    ("geometric", "{pow, {add, [1, {mul, [-1, {sym, x}]}]}, -1}", 1 / (1 - x), 6, 0),
    ("exp_sin", "{exp, {sin, {sym, x}}}", sp.exp(sp.sin(x)), 8, 0),
    ("cos_squared", "{pow, {cos, {sym, x}}, 2}", sp.cos(x) ** 2, 8, 0),
    ("sqrt", "{sqrt, {add, [1, {sym, x}]}}", sp.sqrt(1 + x), 6, 0),
    ("cube_root", "{pow, {add, [1, {mul, [2, {sym, x}]}]}, {rat, 1, 3}}", (1 + 2 * x) ** R(1, 3), 6, 0),
    ("sec", "{pow, {cos, {sym, x}}, -1}", 1 / sp.cos(x), 9, 0),
    ("bernoulli", "{mul, [{sym, x}, {pow, {add, [{exp, {sym, x}}, -1]}, -1}]}", x / (sp.exp(x) - 1), 9, 0),
    ("log_cos", "{log, {cos, {sym, x}}}", sp.log(sp.cos(x)), 9, 0),
    ("log_sum", "{log, {add, [{sym, x}, {pow, {sym, x}, 2}]}}", sp.log(x + x ** 2), 6, 0),
    ("log_sin", "{log, {sin, {sym, x}}}", sp.log(sp.sin(x)), 7, 0),
    ("x_log_x_over", "{mul, [{pow, {sym, x}, 2}, {log, {sym, x}}, {pow, {add, [1, {sym, x}]}, -1}]}",
     x ** 2 * sp.log(x) / (1 + x), 6, 0),
    ("log_shift", "{log, {add, [2, {sym, x}]}}", sp.log(2 + x), 5, 0),
    ("exp_shift", "{exp, {add, [1, {sym, x}]}}", sp.exp(1 + x), 5, 0),
    ("sin_shift", "{sin, {add, [{rat, 1, 3}, {sym, x}]}}", sp.sin(R(1, 3) + x), 5, 0),
    ("param_geometric", "{pow, {add, [1, {mul, [-1, {sym, a}, {sym, x}]}]}, -1}", 1 / (1 - a * x), 6, 0),
    ("param_power", "{pow, {add, [1, {sym, x}]}, {sym, a}}", (1 + x) ** a, 5, 0),
    ("param_exp", "{exp, {mul, [{sym, a}, {sym, x}, {pow, {add, [1, {mul, [{sym, b}, {sym, x}]}]}, -1}]}}",
     sp.exp(a * x / (1 + b * x)), 5, 0),
    ("besselj0", "{besselj, 0, {sym, x}}", sp.besselj(0, x), 10, 0),
    ("besselj3", "{besselj, 3, {sym, x}}", sp.besselj(3, x), 10, 0),
    ("besselj_m2", "{besselj, -2, {sym, x}}", sp.besselj(-2, x), 8, 0),
    ("besselj_half", "{besselj, {rat, 1, 2}, {sym, x}}", sp.besselj(R(1, 2), x), 6, 0),
    ("besselj_third", "{besselj, {rat, 1, 3}, {mul, [2, {sym, x}]}}", sp.besselj(R(1, 3), 2 * x), 5, 0),
    ("bessely0", "{bessely, 0, {sym, x}}", sp.bessely(0, x), 6, 0),
    ("bessely1", "{bessely, 1, {sym, x}}", sp.bessely(1, x), 5, 0),
    ("bessely2", "{bessely, 2, {sym, x}}", sp.bessely(2, x), 5, 0),
    ("laguerre", "{laguerre, 3, {sym, a}, {sym, x}}", sp.assoc_laguerre(3, a, x), 4, 0),
    ("hyp2f1", "{hyp2f1, {sym, a}, {sym, b}, {rat, 5, 2}, {sym, x}}", sp.hyper([a, b], [R(5, 2)], x), 5, 0),
    ("gamma_pole", "{gamma, {sym, x}}", sp.gamma(x), 3, 0),
    ("gamma_one", "{gamma, {add, [1, {sym, x}]}}", sp.gamma(1 + x), 4, 0),
    ("gamma_pole2", "{gamma, {add, [-2, {sym, x}]}}", sp.gamma(x - 2), 2, 0),
    ("at_point", "{pow, {sym, x}, -1}", 1 / x, 4, 2),
    ("at_infinity", "{mul, [{sym, x}, {pow, {add, [{sym, x}, 1]}, -1}]}", x / (x + 1), 4, "infinity"),
]


def terms(expr, n, at):
    if at == "infinity":
        t = sp.Symbol("t", positive=True)
        e = sp.series(expr.subs(x, 1 / t), t, 0, n).removeO()
        e = sp.expand(e)
        out = []
        for term in sp.Add.make_args(e):
            c, rest = term.as_independent(t, as_Add=False)
            pw, lg = split(rest, t)
            out.append((-pw, lg, c * (-1) ** lg))      # t^e log(t)^m = x^-e (-log x)^m
        return out
    if at != 0:
        e = sp.series(expr.subs(x, x + at), x, 0, n).removeO()
    else:
        e = sp.series(expr, x, 0, n).removeO()
    e = sp.expand(e)
    out = []
    for term in sp.Add.make_args(e):
        c, rest = term.as_independent(x, as_Add=False)
        pw, lg = split(rest, x)
        out.append((pw, lg, c))
    return out


def split(rest, v):
    pw, lg = sp.Integer(0), 0
    for f in sp.Mul.make_args(rest):
        base, e = f.as_base_exp()
        if base == v:
            pw += e
        elif base == sp.log(v):
            lg += int(e)
        elif f != 1:
            raise ValueError(f"unexpected factor {f}")
    return pw, lg


cases = []
for name, erl, expr, n, at in CASES:
    ex = expr.subs(PARAMS)
    ts = {}
    for pw, lg, c in terms(ex, n, at):
        key = (sp.Rational(pw), lg)
        ts[key] = ts.get(key, 0) + c
    rows = []
    for (pw, lg), c in sorted(ts.items()):
        v = sp.N(c, 30)
        if v != 0:
            rows.append([int(pw.p), int(pw.q), lg, str(sp.re(v)), str(sp.im(v))])
    cases.append({"name": name, "erl": erl, "n": n, "at": at if at in (0, "infinity") else str(at),
                  "params": {str(k): [int(v.p), int(v.q)] for k, v in PARAMS.items()}, "terms": rows})

out = pathlib.Path(__file__).resolve().parents[1] / "fixtures" / "gt_series.json"
out.write_text(json.dumps({"sympy": sp.__version__, "cases": cases}, indent=1) + "\n")
print(f"{len(cases)} cases -> {out}")
