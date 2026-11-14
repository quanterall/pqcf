# The AdS_5 black brane's free energy by holographic renormalisation

Produced by `brane.erl` (decision record SY2, phase 5), PQCF 5eb5a60. Each step below is computed by PQCF; the expressions are its output.

## Set-up

- **Metric** (Euclidean): ds^2 = (L^2/z^2)(f dtau^2 + dx^2 + dy^2 + dw^2 + dz^2/f), with f = 1 - z^4/z_h^4.
- **Action:** I = -(1/16 pi G) int sqrt(g) (R - 2 Lambda) - (1/8 pi G) int sqrt(gamma) K
  + (1/8 pi G) int sqrt(gamma) 3/L, with the boundary at z = eps (de Haro, Skenderis, Solodukhin,
  [hep-th/0002230](https://arxiv.org/abs/hep-th/0002230)).
- **Free energy:** per unit volume, F = T I = I / beta.

## Steps

1. **Curvature** (`pq_geometry`, 0.15 s): R = `-20*l^-2`. The Kretschmann scalar is `40*l^-4 + 72*l^-4*z^8*zh^-8`.
2. **Volume element and bulk term:** sqrt(g) = `l^5*z^-5`. The bulk term integrated from eps to
   z_h is `1/8*pi^-1*eps^-4*g5^-1*l^3 - 1/8*pi^-1*g5^-1*l^3*zh^-4`.
3. **Extrinsic curvature:** K = `4*sqrt(1 - z^4*zh^-4)*l^-1 + 2*(1 - z^4*zh^-4)^(-1/2)*l^-1*z^4*zh^-4`.
4. **Behaviour as eps -> 0, term by term:**
   - bulk: `1/8*pi^-1*eps^-4*g5^-1*l^3 - 1/8*pi^-1*g5^-1*l^3*zh^-4`;
   - Gibbons-Hawking: `-1/2*pi^-1*eps^-4*g5^-1*l^3 + 1/4*pi^-1*g5^-1*l^3*zh^-4`;
   - counterterm: `3/8*pi^-1*eps^-4*g5^-1*l^3 - 3/16*pi^-1*g5^-1*l^3*zh^-4`.

   The eps^-4 divergences cancel.
5. **Finite limit:** F = `-1/16*pi^-1*g5^-1*l^3*zh^-4`.
6. **In field-theory units:** with T = 1/(pi z_h) and L^3/G = 2 N^2/pi, F = `-1/8*pi^2*n^2*t^4`.

## Result

The expected value is -pi^2 N^2 T^4 / 8, and PQCF's `equiv` answers: `{ok,true}`.
