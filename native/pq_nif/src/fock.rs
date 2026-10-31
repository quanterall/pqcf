//! Single- and two-mode Fock-space matrices (AN 7.2, 7.4), the subjects of spikes S5 and S8.
//!
//! Beamsplitter (AN 7.2): B(theta, phi) conserves N = n_a + n_b. In the sector N, with basis index j = n_a
//! (n_b = N - j), B(theta, 0) = exp(theta G), G = a b^dagger - a^dagger b, real antisymmetric tridiagonal:
//! G[j-1][j] = sqrt(j (N - j + 1)), G[j+1][j] = -sqrt((j + 1)(N - j)). With D = diag(i^j), T = D^dagger (-i G) D
//! is real symmetric tridiagonal with off-diagonal sqrt(j (N - j + 1)) and eigenvalues N, N-2, ..., -N, so
//! U_N(theta, 0) = D V diag(exp(i theta lambda)) V^T D^dagger with the eigenvalues snapped to the exact set, and
//! B(theta, phi) = R_b(phi) B(theta, 0) R_b(-phi): U_N(theta, phi)[j][k] = exp(i phi (k - j)) U_N(theta, 0)[j][k].
//!
//! Displacement and squeezing (AN 7.4, amended by spike S5, docs/decisions/S5-S8-fock-matrices.md): the column
//! recurrences of AN 7.4 lose accuracy exponentially with the cutoff. Displacement uses the Laguerre form along
//! each diagonal; squeezing the Bargmann recurrence where it is accurate, else the exponential of the generator in
//! a padded space, computed per parity by a tridiagonal eigendecomposition.
#![forbid(unsafe_code)] // NFR-04: kernels are safe Rust; unsafe only in the Rustler glue (lib.rs)
use crate::num::C64;
use nalgebra::{DMatrix, SymmetricEigen};

/// Largest deviation of a computed eigenvalue from the exact set before it counts as an internal error.
pub const SNAP_TOL: f64 = 1e-10;

#[derive(Debug, PartialEq)]
pub enum Error {
    /// pq_bug: an eigenvalue of T_N further than SNAP_TOL from the exact set
    EigenSnap { n: usize, deviation: f64 },
}

/// The (N+1) x (N+1) matrix of B(theta, phi) in the sector N, row-major, index j = n_a.
pub fn beamsplitter_sector(n: usize, theta: f64, phi: f64) -> Result<Vec<C64>, Error> {
    let d = n + 1;
    if n == 0 {
        return Ok(vec![C64::ONE]);
    }
    let mut t = DMatrix::<f64>::zeros(d, d);
    for j in 1..d {
        let s = ((j * (n - j + 1)) as f64).sqrt();
        t[(j - 1, j)] = s;
        t[(j, j - 1)] = s;
    }
    let e = SymmetricEigen::new(t);
    // order by eigenvalue, descending: N, N-2, ..., -N
    let mut idx: Vec<usize> = (0..d).collect();
    idx.sort_by(|&a, &b| e.eigenvalues[b].partial_cmp(&e.eigenvalues[a]).unwrap());
    let mut worst: f64 = 0.0;
    let exact: Vec<f64> = (0..d).map(|k| n as f64 - 2.0 * k as f64).collect();
    for (k, &i) in idx.iter().enumerate() {
        worst = worst.max((e.eigenvalues[i] - exact[k]).abs());
    }
    if worst > SNAP_TOL {
        return Err(Error::EigenSnap { n, deviation: worst });
    }
    // U0 = V diag(exp(i theta lambda)) V^T in the T basis
    let phases: Vec<C64> = exact.iter().map(|&l| C64::new((theta * l).cos(), (theta * l).sin())).collect();
    let mut u = vec![C64::ZERO; d * d];
    for a in 0..d {
        for b in 0..d {
            let mut acc = C64::ZERO;
            for (k, &i) in idx.iter().enumerate() {
                acc = acc.add(phases[k].scale(e.eigenvectors[(a, i)] * e.eigenvectors[(b, i)]));
            }
            // back from T to -iG: D = diag(i^j), entry times i^a * conj(i^b) = i^(a - b); then the phase phi
            let rot = i_pow(a as i64 - b as i64);
            let ph = phi * (b as f64 - a as f64);
            u[a * d + b] = acc.mul(rot).mul(C64::new(ph.cos(), ph.sin()));
        }
    }
    Ok(u)
}

/// i^k, exactly
fn i_pow(k: i64) -> C64 {
    match k.rem_euclid(4) {
        0 => C64::new(1.0, 0.0),
        1 => C64::new(0.0, 1.0),
        2 => C64::new(-1.0, 0.0),
        _ => C64::new(0.0, -1.0),
    }
}

/// ln k!, by summation (deterministic, no lgamma)
fn ln_fact(k: usize) -> f64 {
    (2..=k).map(|x| (x as f64).ln()).sum()
}

/// <m|D(alpha)|n> through the Laguerre form, stable: for each k = m - n >= 0 the diagonal
/// f_n = D[n+k][n] = sqrt(n!/(n+k)!) alpha^k e^(-x/2) L_n^(k)(x), x = |alpha|^2, follows the three-term recurrence
///   (n+1) f_(n+1) = (2n+1+k-x) sqrt((n+1)/(n+k+1)) f_n - (n+k) sqrt((n+1) n / ((n+k+1)(n+k))) f_(n-1);
/// the upper triangle by D[n][m] = (-1)^(m-n) conj(D[m][n]).
pub fn displacement_matrix(alpha: C64, c: usize) -> Vec<C64> {
    let d = c + 1;
    let mut m = vec![C64::ZERO; d * d];
    let x = alpha.re * alpha.re + alpha.im * alpha.im;
    let r = x.sqrt();
    let arg = alpha.im.atan2(alpha.re);
    for k in 0..d {
        // f_0 = alpha^k e^(-x/2) / sqrt(k!)
        let f0 = if r == 0.0 {
            if k == 0 { C64::ONE } else { C64::ZERO }
        } else {
            let mag = (k as f64 * r.ln() - x / 2.0 - 0.5 * ln_fact(k)).exp();
            C64::new(mag * (k as f64 * arg).cos(), mag * (k as f64 * arg).sin())
        };
        let kf = k as f64;
        let mut prev = C64::ZERO;
        let mut cur = f0;
        for n in 0..d - k {
            m[(n + k) * d + n] = cur;
            let nf = n as f64;
            let a = (2.0 * nf + 1.0 + kf - x) * ((nf + 1.0) / (nf + kf + 1.0)).sqrt();
            let b = if n == 0 { 0.0 } else { (nf + kf) * ((nf + 1.0) * nf / ((nf + kf + 1.0) * (nf + kf))).sqrt() };
            let next = cur.scale(a).add(prev.scale(-b)).scale(1.0 / (nf + 1.0));
            prev = cur;
            cur = next;
        }
    }
    for row in 0..d {
        for col in row + 1..d {
            let z = m[col * d + row];
            let sign = if (col - row) % 2 == 0 { 1.0 } else { -1.0 };
            m[row * d + col] = C64::new(sign * z.re, -sign * z.im);
        }
    }
    m
}

/// <m|S(z)|n> by the recurrences of the Bargmann generating function (Miatto and Quesada, 2020):
/// G(a*, b) = sech(r)^(1/2) exp(A11 a*^2 / 2 + A12 a* b + A22 b^2 / 2), A11 = -e^(i phi) tanh r, A12 = sech r,
/// A22 = e^(-i phi) tanh r, so
///   sqrt(m+1) S[m+1][n] = A11 sqrt(m) S[m-1][n] + A12 sqrt(n) S[m][n-1]
///   sqrt(n+1) S[0][n+1] = A22 sqrt(n) S[0][n-1]
/// Every coefficient has modulus at most 1, unlike the cosh/sinh column recurrence of AN 7.4.
pub fn squeezing_matrix_bargmann(r: f64, phi: f64, c: usize) -> Vec<C64> {
    let d = c + 1;
    let mut m = vec![C64::ZERO; d * d];
    let (th, sech) = (r.tanh(), 1.0 / r.cosh());
    let a11 = C64::new(-phi.cos() * th, -phi.sin() * th);
    let a22 = C64::new(phi.cos() * th, -phi.sin() * th);
    m[0] = C64::new(sech.sqrt(), 0.0);
    for n in 1..d {
        let prev2 = if n >= 2 { m[n - 2] } else { C64::ZERO };
        m[n] = a22.mul(prev2).scale(((n - 1) as f64).sqrt() / (n as f64).sqrt());
    }
    for row in 0..c {
        let inv = 1.0 / ((row + 1) as f64).sqrt();
        for col in 0..d {
            let up2 = if row >= 1 { a11.mul(m[(row - 1) * d + col]).scale((row as f64).sqrt()) } else { C64::ZERO };
            let left = if col >= 1 { m[row * d + col - 1].scale(sech * (col as f64).sqrt()) } else { C64::ZERO };
            m[(row + 1) * d + col] = up2.add(left).scale(inv);
        }
    }
    m
}

/// <m|S(z)|n> as exp(r G) in a truncated space of `size` states per parity, G = (a^2 - a^dag^2)/2 real
/// antisymmetric and tridiagonal within each parity (n = 2k + p); as for the beamsplitter, D^dagger (-i G) D is real
/// symmetric tridiagonal, off-diagonal sqrt(n (n-1))/2 between n-2 and n. Then S(r, phi)[m][n] =
/// e^(i phi (m-n)/2) S(r, 0)[m][n]. Accurate when the padding holds the squeezed tails (tanh(r)^size small).
pub fn squeezing_matrix_eigen_at(r: f64, phi: f64, c: usize, size: usize) -> Vec<C64> {
    let d = c + 1;
    let mut out = vec![C64::ZERO; d * d];
    for p in 0..2usize {
        let t = DMatrix::<f64>::from_fn(size, size, |a, b| {
            if a + 1 == b {
                let n = 2 * b + p;
                ((n * (n - 1)) as f64).sqrt() / 2.0
            } else if b + 1 == a {
                let n = 2 * a + p;
                ((n * (n - 1)) as f64).sqrt() / 2.0
            } else {
                0.0
            }
        });
        let e = SymmetricEigen::new(t);
        let ks: Vec<usize> = (0..size).filter(|k| 2 * k + p <= c).collect();
        for &a in &ks {
            for &b in &ks {
                let mut acc = C64::ZERO;
                for i in 0..size {
                    let l = e.eigenvalues[i] * r;
                    acc = acc.add(C64::new(l.cos(), l.sin()).scale(e.eigenvectors[(a, i)] * e.eigenvectors[(b, i)]));
                }
                let (m, n) = (2 * a + p, 2 * b + p);
                let rot = i_pow(a as i64 - b as i64);
                let ph = phi * (m as f64 - n as f64) / 2.0;
                out[m * d + n] = acc.mul(rot).mul(C64::new(ph.cos(), ph.sin()));
            }
        }
    }
    out
}

/// Where the Bargmann recurrence meets 1e-12 against mpmath (spike S5, C up to 80, r up to 3).
pub fn bargmann_accurate(r: f64, c: usize) -> bool {
    let r = r.abs();
    c <= 20 || r <= 0.1 || r >= 2.5 || (c <= 40 && (r <= 0.25 || r >= 1.5))
}

/// The change between two paddings below which the eigendecomposition has converged: its rounding floor grows
/// with the size, about 1e-15 per state (spike S5).
fn converged(change: f64, size: usize) -> bool {
    change < f64::max(2e-13, 1e-15 * size as f64)
}

/// Largest padded size per parity for the eigendecomposition.
pub const MAX_EIGEN_SIZE: usize = 4096;

/// <m|S(z)|n>, (c+1) x (c+1): the Bargmann recurrence inside its accurate region, else the padded
/// eigendecomposition, doubling the padding until the block stops changing (above the rounding floor, see
/// `converged`); deterministic: the sizes
/// depend only on r and c). None when MAX_EIGEN_SIZE is not enough.
pub fn squeezing_matrix(r: f64, phi: f64, c: usize) -> Option<Vec<C64>> {
    if bargmann_accurate(r, c) {
        return Some(squeezing_matrix_bargmann(r, phi, c));
    }
    let tail = (46.0 / -(r.abs().tanh().ln()) / 2.0).ceil() as usize;
    let mut size = c / 2 + 1 + tail + 10;
    let mut prev = squeezing_matrix_eigen_at(r, phi, c, size);
    while 2 * size <= MAX_EIGEN_SIZE {
        size *= 2;
        let next = squeezing_matrix_eigen_at(r, phi, c, size);
        let change = prev.iter().zip(&next).map(|(a, b)| (a.re - b.re).abs() + (a.im - b.im).abs()).fold(0.0, f64::max);
        if converged(change, size) {
            return Some(next);
        }
        prev = next;
    }
    None
}

#[cfg(test)]
#[allow(clippy::identity_op, clippy::erasing_op)] // row * width + column, written out
mod tests {
    use super::*;

    fn unitarity_error(u: &[C64], d: usize) -> f64 {
        let mut worst: f64 = 0.0;
        for a in 0..d {
            for b in 0..d {
                let mut acc = C64::ZERO;
                for k in 0..d {
                    let x = u[k * d + a];
                    acc = acc.add(C64::new(x.re, -x.im).mul(u[k * d + b]));
                }
                let want = if a == b { 1.0 } else { 0.0 };
                worst = worst.max((acc.re - want).abs() + acc.im.abs());
            }
        }
        worst
    }

    #[test]
    fn single_photon_matches_fs_3_3() {
        // N = 1: j = 1 is |1, 0>, j = 0 is |0, 1>; |1_a> -> cos t |1_a> + e^(i phi) sin t |1_b>
        let (t, p) = (0.37, 1.1);
        let u = beamsplitter_sector(1, t, p).unwrap();
        let to_b = u[0 * 2 + 1];
        let to_a = u[1 * 2 + 1];
        assert!((to_a.re - t.cos()).abs() < 1e-15 && to_a.im.abs() < 1e-15);
        assert!((to_b.re - p.cos() * t.sin()).abs() < 1e-15 && (to_b.im - p.sin() * t.sin()).abs() < 1e-15);
        // |1_b> -> -e^(-i phi) sin t |1_a> + cos t |1_b>
        let from_b = u[1 * 2];
        assert!((from_b.re + p.cos() * t.sin()).abs() < 1e-15 && (from_b.im - p.sin() * t.sin()).abs() < 1e-15);
    }

    #[test]
    fn hong_ou_mandel_in_the_sector() {
        // N = 2, 50:50: |1,1> (j = 1) never stays |1,1>
        let u = beamsplitter_sector(2, std::f64::consts::FRAC_PI_4, 0.0).unwrap();
        assert!(u[1 * 3 + 1].abs() < 1e-15);
    }

    #[test]
    fn sectors_are_unitary_up_to_100() {
        for n in [1, 2, 5, 10, 20, 50, 100] {
            let u = beamsplitter_sector(n, 0.7, -0.4).unwrap();
            let e = unitarity_error(&u, n + 1);
            assert!(e < 1e-12, "N = {n}: {e:e}");
        }
    }

    #[test]
    fn displacement_and_squeezing_column_zero() {
        let a = C64::new(0.7, 0.2);
        let d = displacement_matrix(a, 10);
        // P(0) of a coherent state = exp(-|alpha|^2) (FS C: 0.588605)
        assert!((d[0].abs().powi(2) - (-0.53f64).exp()).abs() < 1e-15);
        let s = squeezing_matrix(0.3, 0.5, 10).unwrap();
        // P(0) of squeezed vacuum = 1 / cosh r
        assert!((s[0].abs().powi(2) - 1.0 / 0.3f64.cosh()).abs() < 1e-15);
        assert_eq!(s[11].abs(), 0.0); // <1|S|0> = 0
        assert_eq!(displacement_matrix(C64::ZERO, 3)[5], C64::ONE); // identity
    }

    #[test]
    fn squeezing_methods_agree() {
        // inside the Bargmann region both methods are accurate; they share no code
        for &(r, c) in &[(0.4, 20usize), (1.0, 20), (2.0, 16)] {
            let b = squeezing_matrix_bargmann(r, 0.7, c);
            let mut size = c / 2 + 300;
            if r > 1.5 {
                size = 1200;
            }
            let e = squeezing_matrix_eigen_at(r, 0.7, c, size);
            let worst = b.iter().zip(&e).map(|(x, y)| (x.re - y.re).abs() + (x.im - y.im).abs()).fold(0.0, f64::max);
            assert!(worst < 2e-12, "r = {r}, c = {c}: {worst:e}"); // each within 1e-12 of mpmath (S5)
        }
        // outside it the dispatcher takes the eigendecomposition
        let s = squeezing_matrix(0.75, 0.7, 40).unwrap();
        assert!((s[0].abs().powi(2) - 1.0 / 0.75f64.cosh()).abs() < 1e-12);
    }

    #[test]
    fn displacement_is_unitary_where_it_must_be() {
        // columns n of the infinite matrix have unit norm; with c = 120 the first 20 columns are complete
        let d = displacement_matrix(C64::new(1.5, -0.5), 120);
        for n in 0..20 {
            let norm: f64 = (0..121).map(|m| d[m * 121 + n].abs().powi(2)).sum();
            assert!((norm - 1.0).abs() < 1e-13, "column {n}: {norm}");
        }
    }
}
