//! Unitary decomposition (AN 9, FR-C10): the rectangular mesh of Clements et al. (2016) with the beamsplitter of
//! FS 3.3 as the block, B(theta, phi) = [[cos t, -e^(-i phi) sin t], [e^(i phi) sin t, cos t]] on adjacent modes.
//!
//! Entries are nulled diagonal by diagonal, alternately by a block on two columns from the right (U <- U B^dagger)
//! and on two rows from the left (U <- B U), leaving a diagonal D:  L_k .. L_1 U R_1^dagger .. R_n^dagger = D.
//! The left blocks commute through D, B(theta, phi) D = D B(theta, phi + a - b) for D = diag(e^(ia), e^(ib)) on the
//! block's modes, so U = D L'_1 .. L'_k R_n .. R_1: in circuit order R_1 .. R_n, L'_k .. L'_1, then the phases.
//! k(k-1)/2 beamsplitters and k phases.
#![forbid(unsafe_code)] // NFR-04: kernels are safe Rust; unsafe only in the Rustler glue (lib.rs)
use crate::num::C64;

pub type Mesh = (Vec<(usize, f64, f64)>, Vec<f64>);

fn arg(z: C64) -> f64 {
    z.im.atan2(z.re)
}

fn bs(theta: f64, phi: f64) -> [C64; 4] {
    let (c, s) = (theta.cos(), theta.sin());
    [C64::new(c, 0.0), C64::new(-phi.cos() * s, phi.sin() * s), C64::new(phi.cos() * s, phi.sin() * s), C64::new(c, 0.0)]
}

/// rows (r, r+1) of u <- B u
fn left(u: &mut [C64], k: usize, r: usize, b: &[C64; 4]) {
    for c in 0..k {
        let (x, y) = (u[r * k + c], u[(r + 1) * k + c]);
        u[r * k + c] = b[0].mul(x).add(b[1].mul(y));
        u[(r + 1) * k + c] = b[2].mul(x).add(b[3].mul(y));
    }
}

/// columns (c, c+1) of u <- u B^dagger
fn right_dagger(u: &mut [C64], k: usize, c: usize, b: &[C64; 4]) {
    let conj = |z: C64| C64::new(z.re, -z.im);
    // B^dagger = [[conj b0, conj b2], [conj b1, conj b3]]
    let bd = [conj(b[0]), conj(b[2]), conj(b[1]), conj(b[3])];
    for r in 0..k {
        let (x, y) = (u[r * k + c], u[r * k + c + 1]);
        u[r * k + c] = x.mul(bd[0]).add(y.mul(bd[2]));
        u[r * k + c + 1] = x.mul(bd[1]).add(y.mul(bd[3]));
    }
}

/// The mesh of a k x k unitary (row-major): beamsplitters (m, theta, phi) on (m, m+1) in application order, then
/// the phases.
pub fn decompose(u0: &[C64], k: usize) -> Mesh {
    let mut u = u0.to_vec();
    let mut rights: Vec<(usize, f64, f64)> = Vec::new();
    let mut lefts: Vec<(usize, f64, f64)> = Vec::new();
    for i in 0..k.saturating_sub(1) {
        for j in 0..=i {
            if i % 2 == 0 {
                // null u[k-1-j][i-j] with columns (i-j, i-j+1): tan t e^(i phi) = u_rc / u_r,c+1
                let (r, c) = (k - 1 - j, i - j);
                let (a, b) = (u[r * k + c], u[r * k + c + 1]);
                let theta = a.abs().atan2(b.abs());
                let phi = arg(a) - arg(b);
                right_dagger(&mut u, k, c, &bs(theta, phi));
                rights.push((c, theta, phi));
            } else {
                // null u[k-1-i+j][j] with rows (k-2-i+j, k-1-i+j): tan t e^(i phi) = -u_(r+1),c / u_r,c
                let (r, c) = (k - 2 - i + j, j);
                let (a, b) = (u[r * k + c], u[(r + 1) * k + c]);
                let theta = b.abs().atan2(a.abs());
                let phi = std::f64::consts::PI + arg(b) - arg(a);
                left(&mut u, k, r, &bs(theta, phi));
                lefts.push((r, theta, phi));
            }
        }
    }
    let mut phases: Vec<f64> = (0..k).map(|j| arg(u[j * k + j])).collect();
    // U = L_1^dag .. L_k^dag D R_n .. R_1; push D to the far left: L^dag = B(-t, phi), B(-t, phi) D = D B(-t, phi + a - b)
    let mut moved: Vec<(usize, f64, f64)> = Vec::new();
    for &(r, t, p) in lefts.iter().rev() {
        moved.push((r, -t, p + phases[r] - phases[r + 1]));
    }
    let mut mesh: Vec<(usize, f64, f64)> = rights;
    mesh.extend(moved);
    for p in phases.iter_mut() {
        *p = p.rem_euclid(2.0 * std::f64::consts::PI);
    }
    (mesh, phases)
}

/// The unitary of a mesh (for tests and checks): the beamsplitters in order, then the phases.
pub fn reconstruct(mesh: &Mesh, k: usize) -> Vec<C64> {
    let mut u: Vec<C64> = (0..k * k).map(|i| if i / k == i % k { C64::ONE } else { C64::ZERO }).collect();
    for &(m, t, p) in &mesh.0 {
        left(&mut u, k, m, &bs(t, p));
    }
    for (j, &a) in mesh.1.iter().enumerate() {
        let z = C64::new(a.cos(), a.sin());
        for c in 0..k {
            u[j * k + c] = z.mul(u[j * k + c]);
        }
    }
    u
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Stream;
    use nalgebra::{Complex, DMatrix};

    fn haar(k: usize, seed: u64) -> Vec<C64> {
        let mut s = Stream::new(seed, 0);
        let z = DMatrix::from_fn(k, k, |_, _| {
            let (a, b) = s.normal_pair();
            Complex::new(a, b)
        });
        let qr = z.qr();
        let (q, r) = (qr.q(), qr.r());
        let mut u = vec![C64::ZERO; k * k];
        for j in 0..k {
            let ph = r[(j, j)] / r[(j, j)].norm();
            for i in 0..k {
                let x = q[(i, j)] * ph;
                u[i * k + j] = C64::new(x.re, x.im);
            }
        }
        u
    }

    fn err(a: &[C64], b: &[C64]) -> f64 {
        a.iter().zip(b).map(|(x, y)| (x.re - y.re).abs() + (x.im - y.im).abs()).fold(0.0, f64::max)
    }

    #[test]
    fn random_unitaries_up_to_32() {
        for k in 1..=32 {
            let u = haar(k, k as u64);
            let m = decompose(&u, k);
            assert_eq!(m.0.len(), k * (k - 1) / 2);
            let e = err(&reconstruct(&m, k), &u);
            assert!(e <= 1e-12 * k as f64, "k = {k}: {e:e}");
        }
    }

    #[test]
    fn identity_permutation_fourier() {
        let k = 6;
        let id: Vec<C64> = (0..k * k).map(|i| if i / k == i % k { C64::ONE } else { C64::ZERO }).collect();
        let perm: Vec<C64> = (0..k * k).map(|i| if (i / k + 1) % k == i % k { C64::ONE } else { C64::ZERO }).collect();
        let w = 2.0 * std::f64::consts::PI / k as f64;
        let four: Vec<C64> = (0..k * k)
            .map(|i| {
                let a = w * ((i / k) * (i % k)) as f64;
                C64::new(a.cos(), a.sin()).scale(1.0 / (k as f64).sqrt())
            })
            .collect();
        for u in [id, perm, four] {
            let m = decompose(&u, k);
            assert!(err(&reconstruct(&m, k), &u) <= 1e-12 * k as f64);
        }
    }
}
