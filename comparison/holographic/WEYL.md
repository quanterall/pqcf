# The holographic Weyl anomaly: a = c = N^2/4

Produced by `weyl.erl` (decision record SY2, phase 5), PQCF 0c4a889. Every value below is PQCF's output.

## Method

Henningson and Skenderis, [hep-th/9806087](https://arxiv.org/abs/hep-th/9806087); de Haro, Skenderis and
Solodukhin, [hep-th/0002230](https://arxiv.org/abs/hep-th/0002230).

1. **Bulk metric:** L^2 [drho^2/(4 rho^2) + g(rho)_ij dx^i dx^j / rho], with
   g(rho) = g0 (1 + rho m + rho^2 tau) for an Einstein boundary g0.
2. **Coefficients:** the 5-dimensional Einstein equations R_ab + (4/L^2) G_ab = 0 are solved order by order.
   The boundary components at rho^0 fix m, so g_(2) = m g0; the radial component at rho^0 then fixes tau,
   and with it Tr g_(4).
3. **Volume coefficient:** a_4 is the rho^2 coefficient of sqrt(det g(rho) / det g0) = (1 + rho m + rho^2 tau)^2.
4. **Anomaly:** the log divergence of the on-shell action under a constant rescaling of g0 gives
   <T> = -(L^3/(2 pi G)) a_4, compared with (c W^2 - a E_4)/(16 pi^2).
5. **Separating a from c:** S^4 has W = 0, so it gives a. S^2 x S^2 with equal radii has W^2 /= 0, so it then
   gives c. Both are in stereographic coordinates.
6. **Units:** L^3/G = 2 N^2/pi.

## Results

| boundary (radius r) | m | tau | <T> | E_4 | W^2 |
| --- | --- | --- | --- | --- | --- |
| S^4 | `-1/2*r^-2` | `1/16*r^-4` | `-3/16*pi^-1*g5^-1*l^3*r^-4` | `24*r^-4` | `0` |
| S^2 x S^2 | `-1/6*r^-2` | `1/144*r^-4` | `-1/48*pi^-1*g5^-1*l^3*r^-4` | `8*r^-4` | `16/3*r^-4` |

- a = `1/4*n^2`
- c = `1/4*n^2`

Against N^2/4, PQCF's `equiv` answers `{{ok,true},{ok,true}}`.

Times: S^4 2.6 s, S^2 x S^2 0.7 s.
