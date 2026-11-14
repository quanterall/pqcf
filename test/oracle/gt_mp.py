"""GT-MP (decision record SY2, phase 2): values of constants, elementary and special functions, and zeros of Bessel
functions by mpmath 1.3 at 60 digits, for pq_mp's tests (compared to 48 digits). Writes test/fixtures/gt_mp.json.

    docker run --rm -u $(id -u):$(id -g) -v "$PWD":/w -w /w --entrypoint python pqcf-oracle test/oracle/gt_mp.py
"""
import json
import pathlib

import mpmath as m

m.mp.dps = 60
F = m.mpf


def q(a, b):
    return F(a) / b


# (Erlang expression for pq_mp:evaluate/3, mpmath value)
VALUES = [
    ("{const, pi}", m.pi),
    ("{const, euler_gamma}", m.euler),
    ("{log, 2}", m.ln2),
    ("{exp, {rat, -37, 3}}", m.exp(q(-37, 3))),
    ("{exp, 50}", m.exp(50)),
    ("{log, {rat, 1, 7}}", m.log(q(1, 7))),
    ("{log, 123456789}", m.log(123456789)),
    ("{sin, 100}", m.sin(100)),
    ("{cos, {rat, 22, 7}}", m.cos(q(22, 7))),
    ("{sqrt, 2}", m.sqrt(2)),
    ("{pow, 3, {rat, 1, 3}}", m.cbrt(3)),
    ("{pow, {rat, 5, 2}, {const, pi}}", q(5, 2) ** m.pi),
    ("{zeta, 2}", m.zeta(2)),
    ("{zeta, 3}", m.zeta(3)),
    ("{zeta, 7}", m.zeta(7)),
    ("{gamma, {rat, 1, 3}}", m.gamma(q(1, 3))),
    ("{gamma, {rat, -5, 2}}", m.gamma(q(-5, 2))),
    ("{gamma, 30}", m.gamma(30)),
    ("{gamma, {rat, 201, 2}}", m.gamma(q(201, 2))),
    ("{polygamma, 0, {rat, 1, 3}}", m.psi(0, q(1, 3))),
    ("{polygamma, 1, 1}", m.psi(1, 1)),
    ("{polygamma, 2, {rat, 5, 2}}", m.psi(2, q(5, 2))),
    ("{besselj, 0, {rat, 1, 2}}", m.besselj(0, q(1, 2))),
    ("{besselj, 1, 10}", m.besselj(1, 10)),
    ("{besselj, 3, {rat, 307, 10}}", m.besselj(3, q(307, 10))),
    ("{besselj, -2, 7}", m.besselj(-2, 7)),
    ("{besselj, {rat, 1, 3}, {rat, 5, 4}}", m.besselj(q(1, 3), q(5, 4))),
    ("{besselj, {rat, 5, 2}, 20}", m.besselj(q(5, 2), 20)),
    ("{bessely, 0, {rat, 3, 10}}", m.bessely(0, q(3, 10))),
    ("{bessely, 1, 12}", m.bessely(1, 12)),
    ("{bessely, 3, {rat, 5, 2}}", m.bessely(3, q(5, 2))),
    ("{laguerre, 5, {rat, 1, 2}, {rat, 33, 10}}", m.laguerre(5, q(1, 2), q(33, 10))),
    ("{laguerre, 8, 0, 20}", m.laguerre(8, 0, 20)),
    ("{hyp2f1, {rat, 1, 3}, {rat, 2, 3}, {rat, 3, 2}, {rat, 1, 2}}", m.hyp2f1(q(1, 3), q(2, 3), q(3, 2), q(1, 2))),
    ("{hyp2f1, 1, 1, 2, {rat, -1, 2}}", m.hyp2f1(1, 1, 2, q(-1, 2))),
    # coefficients pq_series produces for Y_0 and Gamma(x): (2/pi)(gamma - log 2), (pi^2/6 + gamma^2)/2
    ("{mul, [2, {pow, {const, pi}, -1}, {add, [{const, euler_gamma}, {mul, [-1, {log, 2}]}]}]}",
     2 / m.pi * (m.euler - m.ln2)),
    ("{mul, [{rat, 1, 2}, {add, [{zeta, 2}, {pow, {const, euler_gamma}, 2}]}]}", (m.zeta(2) + m.euler ** 2) / 2),
]

# (nu as an Erlang term, k, mpmath zero)
ZEROS = [("0", k, m.besseljzero(0, k)) for k in (1, 2, 3, 10)] + \
        [("1", k, m.besseljzero(1, k)) for k in (1, 5)] + \
        [("{rat, 1, 3}", 2, m.besseljzero(q(1, 3), 2)), ("4", 3, m.besseljzero(4, 3))]

out = pathlib.Path(__file__).resolve().parents[1] / "fixtures" / "gt_mp.json"
out.write_text(json.dumps({"mpmath": m.__version__, "digits": 60,
                           "values": [{"erl": e, "value": m.nstr(v, 60, min_fixed=1, max_fixed=0)} for e, v in VALUES],
                           "zeros": [{"nu": n, "k": k, "value": m.nstr(v, 60, min_fixed=1, max_fixed=0)} for n, k, v in ZEROS]},
                          indent=1) + "\n")
print(f"{len(VALUES)} values, {len(ZEROS)} zeros -> {out}")
