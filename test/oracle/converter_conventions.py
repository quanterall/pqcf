"""The conventions behind pq_convert (FS 3.3: converters adapt to our convention), measured against the libraries.
Asserts, so a library version that changes them fails here. Ours: B(theta, phi) = [[cos t, -e^(-i phi) sin t],
[e^(i phi) sin t, cos t]] (column j = image of mode j), R(phi) = e^(i phi), S(z) and D(alpha) of FS 3.3, Kerr
e^(i kappa n^2), vacuum covariance (hbar/2) I.

    docker run --rm pqcf-oracle python test/oracle/converter_conventions.py
"""
import numpy as np
import perceval as pcvl
import piquasso as pq


def ours_bs(t, p):
    return np.array([[np.cos(t), -np.exp(-1j * p) * np.sin(t)], [np.exp(1j * p) * np.sin(t), np.cos(t)]])


# Perceval BS.Rx and PS
t, tl, bl, tr, br = 0.7, 0.1, 0.2, 0.3, 0.4
U = np.array(pcvl.BS.Rx(theta=t, phi_tl=tl, phi_bl=bl, phi_tr=tr, phi_br=br).compute_unitary())
c, s = np.cos(t / 2), np.sin(t / 2)
assert np.allclose(U, [[np.exp(1j * (tl + tr)) * c, 1j * np.exp(1j * (tr + bl)) * s],
                       [1j * np.exp(1j * (tl + br)) * s, np.exp(1j * (bl + br)) * c]])
th, ph = 0.45, 1.2
assert np.allclose(np.array(pcvl.BS.Rx(theta=2 * th, phi_bl=np.pi / 2 - ph, phi_br=ph - np.pi / 2).compute_unitary()),
                   ours_bs(th, ph))
assert np.allclose(np.array(pcvl.PS(0.5).compute_unitary()), [[np.exp(0.5j)]])
# orientation: a circuit's unitary is the product in application order (later components on the left)
circ = pcvl.Circuit(2) // pcvl.PS(0.3) // pcvl.BS.Rx(theta=0.9)
assert np.allclose(np.array(circ.compute_unitary()),
                   np.array(pcvl.BS.Rx(theta=0.9).compute_unitary()) @ np.diag([np.exp(0.3j), 1]))

# Piquasso gates
cn, cf = pq.NumpyConnector(), pq.Config()
assert np.allclose(np.array(pq.Beamsplitter(theta=th, phi=ph)._get_passive_block(cn, cf)), ours_bs(th, ph))
assert np.allclose(np.array(pq.Phaseshifter(phi=0.5)._get_passive_block(cn, cf)), [[np.exp(0.5j)]])
r, p = 0.4, 0.7
with pq.Program() as prog:
    pq.Q(0) | pq.Squeezing(r=r, phi=p)
    pq.Q(0) | pq.Displacement(r=0.5, phi=p)
st = pq.GaussianSimulator(d=1, config=pq.Config(hbar=2)).execute(prog).state
ch, sh = np.cosh(r), np.sinh(r)
S = np.array([[ch - np.cos(p) * sh, -np.sin(p) * sh], [-np.sin(p) * sh, ch + np.cos(p) * sh]])
assert np.allclose(st.xpxp_covariance_matrix, 2 * (S @ S.T)), "Piquasso's covariance is twice ours"
assert np.allclose(st.xpxp_mean_vector, S @ [0, 0] + [2 * 0.5 * np.cos(p), 2 * 0.5 * np.sin(p)])
with pq.Program() as prog:
    pq.Q(0) | pq.NumberState([2]) * 1.0 if hasattr(pq, "NumberState") else pq.StateVector([2]) * 1.0
    pq.Q(0) | pq.Kerr(xi=0.3)
v = np.array(pq.PureFockSimulator(d=1, config=pq.Config(cutoff=4)).execute(prog).state.state_vector)
assert np.isclose(v[2], np.exp(1j * 0.3 * 4))
print("converter conventions hold: perceval", pcvl.__version__, "piquasso", pq.__version__)
