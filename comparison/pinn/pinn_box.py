"""A physics-informed neural network (PINN) for Feynman integrals, against classical solvers (decision record R1).

The one-loop massless box's masters M = (box, s-bubble, t-bubble) obey dM/dt = A(d, s, t) M, A derived exactly by
PQCF (box_de.erl: derivatives by external-momentum operators, reduced to the masters by IBP). Here:

  1. the equation is checked against the integrals themselves: the box by direct numerical integration over its
     Feynman parameters, the bubbles in closed form, an ODE solver (DOP853) carrying the box from t0 to other t;
  2. a PINN learns M(t) on t in [-10, -0.1] (s = -1) from the equation and the value at t0 = -1 alone, against the
     ODE solution: accuracy and time;
  3. a PINN over (t, d) at once, a surrogate for a range of dimensions, its boundary values given at a few d only.

The network and training follow TechUPGRADE's PINN-PONN reproduction (TU-TR-2026-002, atech-up/sim-manager/ponn/
pinn.py; reused with permission): a tanh input layer, residual blocks of two tanh layers, AdamW inside a look-ahead
(k = 5, alpha = 0.5); here followed by L-BFGS, as PINN practice suggests for a final polish. Run in TechUPGRADE's
phnn-trainer image (JAX, CUDA), read-only:

  docker run --rm --gpus all -e XLA_PYTHON_CLIENT_PREALLOCATE=false -e XLA_PYTHON_CLIENT_MEM_FRACTION=0.05 \\
      -v $PWD:/w -w /w --entrypoint python phnn-trainer:0.5.0 pinn_box.py box_de.txt
"""
import math
import sys
import time

import numpy as np
from scipy import integrate, special

import jax
import jax.numpy as jnp
import optax

jax.config.update("jax_enable_x64", True)

S = -1.0
T0, T_LO, T_HI = -1.0, -10.0, -0.1


def load_de(path):
    """A(d, s, t) as a function returning a 3 x 3 array (numpy or jax), from box_de.erl's file"""
    entries = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        i, j, expr = line.split(" ", 2)
        entries[(int(i), int(j))] = compile(expr.strip(), "A", "eval")

    def A(d, s, t, xp=np):
        env = {"d": d, "s": s, "t": t}
        return [[xp.asarray(eval(entries[(i, j)], {}, env)) * xp.ones_like(t) for j in range(3)] for i in range(3)]
    return A


# --- the integrals themselves (Feynman parameters; Euclidean region s, t < 0) ---

def bubble(d, q2):
    return special.gamma(2 - d / 2) * (-q2) ** (d / 2 - 2) * special.beta(d / 2 - 1, d / 2 - 1)


def _tanh_sinh(n):
    """tanh-sinh nodes on [0, 1] with their complements 1 - x (exact near 1) and weights: endpoint singularities
    (integrable powers) cost little"""
    h = 6.0 / n
    tau = (np.arange(n + 1) - n / 2) * h
    q = np.pi / 2 * np.sinh(tau)
    x = 1 / (1 + np.exp(-2 * q))
    xc = 1 / (1 + np.exp(2 * q))
    w = h * np.pi / 2 * np.cosh(tau) / (2 * np.cosh(q) ** 2)
    return x, xc, w


def box(d, s, t, n=256):
    """Gamma(4 - d/2) int over the simplex of F^(d/2 - 4), F = -s x1 x3 - t x2 x4; x1, x3 = a u, a (1 - u),
    x2, x4 = (1 - a) v, (1 - a)(1 - v), Jacobian a (1 - a); a tensor tanh-sinh rule, its error estimated by the
    rule with half the nodes"""
    def rule(n):
        x, xc, w = _tanh_sinh(n)
        a, ac, wa = x[:, None, None], xc[:, None, None], w[:, None, None]
        uu = (x * xc)[None, :, None]
        vv = (x * xc)[None, None, :]
        F = -s * a * a * uu - t * ac * ac * vv
        return float((wa * w[None, :, None] * w[None, None, :] * a * ac * F ** (d / 2 - 4)).sum())
    fine, coarse = rule(n), rule(n // 2)
    g = special.gamma(4 - d / 2)
    return g * fine, g * abs(fine - coarse)


def masters(d, s, t):
    return np.array([box(d, s, t)[0], bubble(d, s), bubble(d, t)])


# --- the classical solution: x = ln(-t), dM/dx = t A M ---

def ode(A, d, m0, xs):
    def rhs(x, m):
        t = -math.exp(x)
        a = np.array([[float(e) for e in row] for row in A(d, S, t)])
        return t * a @ m
    x0 = math.log(-T0)
    out = np.empty((len(xs), 3))
    lo, hi = xs[xs <= x0][::-1], xs[xs > x0]
    for part, sign in ((lo, -1), (hi, 1)):
        if len(part):
            sol = integrate.solve_ivp(rhs, (x0, part[-1]), m0, method="DOP853", rtol=1e-13, atol=0, dense_output=True)
            mask = np.isin(xs, part)
            out[mask] = sol.sol(xs[mask]).T
    return out


# --- the PINN (after TU-TR-2026-002) ---

def init(key, n_in, width=32, blocks=3, n_out=3):
    def dense(k, a, b):
        lim = math.sqrt(6 / (a + b))
        return (jax.random.uniform(k, (a, b), minval=-lim, maxval=lim), jnp.zeros(b))
    ks = jax.random.split(key, 2 + 2 * blocks)
    return {"inp": dense(ks[0], n_in, width),
            "blocks": [(dense(ks[1 + 2 * i], width, width), dense(ks[2 + 2 * i], width, width)) for i in range(blocks)],
            "out": dense(ks[-1], width, n_out)}


def net(th, z):
    h = jnp.tanh(z @ th["inp"][0] + th["inp"][1])
    for (W1, b1), (W2, b2) in th["blocks"]:
        h = h + jnp.tanh(jnp.tanh(h @ W1 + b1) @ W2 + b2)
    return h @ th["out"][0] + th["out"][1]


def train(loss, th, steps_adam, steps_lbfgs, lr=2e-3):
    adam = optax.lookahead(optax.adamw(optax.cosine_decay_schedule(lr, steps_adam), weight_decay=1e-6), sync_period=5,
                           slow_step_size=0.5)
    params = optax.LookaheadParams.init_synced(th)
    state = adam.init(params)
    vg = jax.jit(jax.value_and_grad(loss))

    @jax.jit
    def step(params, state):
        # the gradient at the fast weights; the look-ahead keeps the slow ones
        val, g = vg(params.fast)
        upd, state = adam.update(g, state, params)
        return optax.apply_updates(params, upd), state, val

    hist = []
    for i in range(steps_adam):
        params, state, val = step(params, state)
        if i % 1000 == 0:
            hist.append(float(val))
    th = params.slow
    lb = optax.lbfgs()
    st = lb.init(th)
    vg2 = optax.value_and_grad_from_state(loss)

    @jax.jit
    def lstep(th, st):
        val, g = vg2(th, state=st)
        upd, st = lb.update(g, st, th, value=val, grad=g, value_fn=loss)
        return optax.apply_updates(th, upd), st, val

    for i in range(steps_lbfgs):
        th, st, val = lstep(th, st)
    hist.append(float(loss(th)))
    return th, hist


def experiment_1(A, d, m0, xs_test, ode_test):
    """the PINN on t in [T_LO, T_HI] at one d, the value at T0 given"""
    x0, xlo, xhi = math.log(-T0), math.log(-T_LO), math.log(-T_HI)
    xc = jnp.linspace(min(xlo, xhi), max(xlo, xhi), 256)
    zc = (2 * (xc - xc[0]) / (xc[-1] - xc[0]) - 1)[:, None]
    scale = jnp.asarray(np.abs(m0))
    m0j = jnp.asarray(m0)

    def M(th, x, z):
        # the boundary value exactly: M(x0) = m0
        return m0j * (1 + (x - x0)[:, None] * net(th, z))

    def loss(th):
        Mx = M(th, xc, zc)
        # each output depends on its own point only: one forward-mode pass with a tangent of ones is dM/dx
        dM = jax.jvp(lambda xx: M(th, xx, (2 * (xx - xc[0]) / (xc[-1] - xc[0]) - 1)[:, None]), (xc,), (jnp.ones_like(xc),))[1]
        t = -jnp.exp(xc)
        a = A(d, S, t, xp=jnp)
        rhs = jnp.stack([t * sum(a[i][j] * Mx[:, j] for j in range(3)) for i in range(3)], -1)
        return jnp.mean(((dM - rhs) / scale) ** 2)

    t_start = time.time()
    th, hist = train(loss, init(jax.random.PRNGKey(42), 1), 20000, 3000)
    jax.block_until_ready(th)
    t_train = time.time() - t_start
    xt = jnp.asarray(xs_test)
    pred = np.asarray(M(th, xt, (2 * (xt - xc[0]) / (xc[-1] - xc[0]) - 1)[:, None]))
    rel = np.abs(pred - ode_test) / np.abs(ode_test)
    return t_train, hist, rel.max(0)


def experiment_2(A, ds_train, ds_test, m0_of_d, ode_of_d):
    """one PINN over (t, d): the equation everywhere, the value at T0 only at the training dimensions"""
    x0, xlo, xhi = math.log(-T0), math.log(-T_LO), math.log(-T_HI)
    dlo, dhi = min(ds_train), max(ds_train)
    xa = np.linspace(min(xlo, xhi), max(xlo, xhi), 64)
    da = np.linspace(dlo, dhi, 32)
    X, D = np.meshgrid(xa, da)
    xc, dc = jnp.asarray(X.ravel()), jnp.asarray(D.ravel())
    nz = lambda x, d: jnp.stack([2 * (x - xa[0]) / (xa[-1] - xa[0]) - 1, 2 * (d - dlo) / (dhi - dlo) - 1], -1)
    b0 = jnp.asarray(np.array([m0_of_d(d) for d in ds_train]))
    scale = jnp.abs(b0).mean(0)

    def M(th, x, d):
        return scale * net(th, nz(x, d))

    def loss(th):
        Mx = M(th, xc, dc)
        dM = jax.jvp(lambda xx: M(th, xx, dc), (xc,), (jnp.ones_like(xc),))[1]
        t = -jnp.exp(xc)
        a = A(dc, S, t, xp=jnp)
        rhs = jnp.stack([t * sum(a[i][j] * Mx[:, j] for j in range(3)) for i in range(3)], -1)
        res = jnp.mean(((dM - rhs) / scale) ** 2)
        bc = jnp.mean(((M(th, jnp.full(len(ds_train), x0), jnp.asarray(ds_train)) - b0) / scale) ** 2)
        return res + 10 * bc

    t_start = time.time()
    th, hist = train(loss, init(jax.random.PRNGKey(7), 2, width=48), 30000, 3000)
    t_train = time.time() - t_start
    xs = np.linspace(min(xlo, xhi), max(xlo, xhi), 41)
    errs = {}
    for d in ds_test:
        pred = np.asarray(M(th, jnp.asarray(xs), jnp.full(len(xs), d)))
        ref = ode_of_d(d, xs)
        errs[d] = (np.abs(pred - ref) / np.abs(ref)).max(0)
    return t_train, hist, errs


def main(path):
    A = load_de(path)
    d = 5.3
    print(f"devices: {jax.devices()}")
    # 1. the equation against the integrals
    t1 = time.time()
    m0 = masters(d, S, T0)
    t_m0 = time.time() - t1
    print(f"masters at t0 = {T0}, d = {d}: {m0}  (box by 3-d integration in {t_m0:.1f} s)")
    xs_check = np.log([-T_LO, 3.0, 0.3, -T_HI])
    t1 = time.time()
    ode_check = ode(A, d, m0, np.sort(xs_check))
    t_ode = time.time() - t1
    for x, row in zip(np.sort(xs_check), ode_check):
        t = -math.exp(x)
        direct = masters(d, S, t)
        print(f"  t = {t:7.3f}: ODE box {row[0]: .12e}, direct {direct[0]: .12e}, rel. diff {abs(row[0] / direct[0] - 1):.1e};"
              f" t-bubble rel. diff {abs(row[2] / direct[2] - 1):.1e}")
    print(f"  ODE solve (DOP853, rtol 1e-13): {t_ode * 1000:.1f} ms")
    # 2. the PINN at one d
    xs_test = np.linspace(math.log(-T_LO), math.log(-T_HI), 201)[::-1]
    xs_test = np.sort(xs_test)
    t1 = time.time()
    ode_test = ode(A, d, m0, xs_test)
    t_ode2 = time.time() - t1
    t_train, hist, rel = experiment_1(A, d, m0, xs_test, ode_test)
    print(f"PINN, d = {d}: trained in {t_train:.1f} s (loss {hist[0]:.2e} -> {hist[-1]:.2e}); max rel. error vs ODE on "
          f"201 points: box {rel[0]:.1e}, s-bubble {rel[1]:.1e}, t-bubble {rel[2]:.1e}; the ODE: {t_ode2 * 1000:.1f} ms")
    # 3. the surrogate over (t, d)
    ds_train = [4.6, 4.9, 5.2, 5.5, 5.8]
    ds_test = [4.75, 5.05, 5.3, 5.65]
    t1 = time.time()
    cache = {dd: masters(dd, S, T0) for dd in ds_train + ds_test}
    t_bc = time.time() - t1
    ode_of_d = lambda dd, xs: ode(A, dd, cache[dd], xs)
    t_train2, hist2, errs = experiment_2(A, ds_train, ds_test, lambda dd: cache[dd], ode_of_d)
    print(f"PINN over (t, d), d in [4.6, 5.8], boundary values at {ds_train} ({t_bc:.0f} s of integration): trained in "
          f"{t_train2:.1f} s (loss {hist2[0]:.2e} -> {hist2[-1]:.2e})")
    for dd, e in errs.items():
        kind = "trained" if dd in ds_train else "unseen"
        print(f"  d = {dd} ({kind}): max rel. error box {e[0]:.1e}, s-bubble {e[1]:.1e}, t-bubble {e[2]:.1e}")


if __name__ == "__main__":
    main(sys.argv[1])
