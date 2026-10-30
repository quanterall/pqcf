//! Gaussian states (FS 3.2, A.2-A.6, A.12; AN 6.4, 8): mean `mu` (2m) and covariance `sigma` (2m x 2m, row-major),
//! xxpp order: index j is x_j, index m + j is p_j. Maps are applied one instruction at a time (AN 8); nothing is
//! composed or projected.
use crate::num::{Sum, C64};
use crate::rng::Stream;
use nalgebra::{DMatrix, SymmetricEigen};

#[derive(Debug, PartialEq)]
pub enum Error {
    BadShape(&'static str),
    /// S Omega S^T - Omega, largest entry
    NotSymplectic(f64),
    /// the smallest eigenvalue of a covariance that must be positive semidefinite
    NonPhysical(f64),
    Cancelled,
}

/// the global indices of the local xxpp coordinates of `modes`: x of each mode, then p of each mode
fn coords(modes: &[usize], m: usize) -> Vec<usize> {
    modes.iter().copied().chain(modes.iter().map(|&j| j + m)).collect()
}

fn dot(a: impl Iterator<Item = f64>, len: usize) -> f64 {
    let mut s = Sum::new(len as u64);
    for x in a {
        s.add(x);
    }
    s.value()
}

/// S Omega S^T - Omega, largest absolute entry; S is 2k x 2k
pub fn symplectic_deviation(s: &[f64], k: usize) -> f64 {
    let d = 2 * k;
    let omega = |i: usize, j: usize| -> f64 {
        if i < k && j == i + k {
            1.0
        } else if i >= k && j + k == i {
            -1.0
        } else {
            0.0
        }
    };
    let mut worst: f64 = 0.0;
    for i in 0..d {
        for j in 0..d {
            // (S Omega S^T)_ij = sum_a,b S_ia Omega_ab S_jb
            let v = dot((0..d).flat_map(|a| (0..d).map(move |b| (a, b))).map(|(a, b)| s[i * d + a] * omega(a, b) * s[j * d + b]), d * d);
            worst = worst.max((v - omega(i, j)).abs());
        }
    }
    worst
}

/// mu' = S mu + D, sigma' = S sigma S^T on `modes` (A.2); `check`: reject S that is not symplectic within `tol`
pub fn apply_symplectic(
    mu: &mut [f64],
    sigma: &mut [f64],
    m: usize,
    s: &[f64],
    d: &[f64],
    modes: &[usize],
    check: Option<f64>,
) -> Result<(), Error> {
    let k = modes.len();
    let n = 2 * m;
    if mu.len() != n || sigma.len() != n * n || s.len() != 4 * k * k || d.len() != 2 * k {
        return Err(Error::BadShape("mean, covariance, S or D"));
    }
    if let Some(tol) = check {
        let dev = symplectic_deviation(s, k);
        if dev > tol {
            return Err(Error::NotSymplectic(dev));
        }
    }
    let c = coords(modes, m);
    let kk = 2 * k;
    // mean
    let old: Vec<f64> = c.iter().map(|&g| mu[g]).collect();
    for a in 0..kk {
        mu[c[a]] = dot((0..kk).map(|b| s[a * kk + b] * old[b]), kk) + d[a];
    }
    // rows: sigma[c_a, :] = sum_b S_ab sigma[c_b, :]
    let rows: Vec<Vec<f64>> = c.iter().map(|&g| sigma[g * n..(g + 1) * n].to_vec()).collect();
    for a in 0..kk {
        for col in 0..n {
            sigma[c[a] * n + col] = dot((0..kk).map(|b| s[a * kk + b] * rows[b][col]), kk);
        }
    }
    // columns: sigma[:, c_a] = sum_b sigma[:, c_b] S_ab
    let cols: Vec<Vec<f64>> = c.iter().map(|&g| (0..n).map(|r| sigma[r * n + g]).collect()).collect();
    for a in 0..kk {
        for row in 0..n {
            sigma[row * n + c[a]] = dot((0..kk).map(|b| cols[b][row] * s[a * kk + b]), kk);
        }
    }
    Ok(())
}

/// The symplectic matrix of a passive k x k unitary (A.3): [[Ur, -Ui], [Ui, Ur]]
pub fn passive_symplectic(u: &[C64], k: usize) -> Vec<f64> {
    let d = 2 * k;
    let mut s = vec![0.0; d * d];
    for i in 0..k {
        for j in 0..k {
            let z = u[i * k + j];
            s[i * d + j] = z.re;
            s[i * d + j + k] = -z.im;
            s[(i + k) * d + j] = z.im;
            s[(i + k) * d + j + k] = z.re;
        }
    }
    s
}

/// Loss on one mode (A.6): mu -> sqrt(eta) mu; sigma -> eta sigma + (1 - eta)(hbar/2) I on the mode, cross blocks
/// scaled by sqrt(eta)
pub fn apply_loss(mu: &mut [f64], sigma: &mut [f64], m: usize, mode: usize, eta: f64, hbar: f64) -> Result<(), Error> {
    let n = 2 * m;
    if mu.len() != n || sigma.len() != n * n || mode >= m {
        return Err(Error::BadShape("mean, covariance or mode"));
    }
    let r = eta.sqrt();
    let own = [mode, mode + m];
    for &g in &own {
        mu[g] *= r;
        for j in 0..n {
            sigma[g * n + j] *= r;
            sigma[j * n + g] *= r;
        }
    }
    for &g in &own {
        sigma[g * n + g] += (1.0 - eta) * hbar / 2.0;
    }
    Ok(())
}

/// The reduced state on `modes` (partial trace), xxpp in the order of `modes`
pub fn marginal(mu: &[f64], sigma: &[f64], m: usize, modes: &[usize]) -> (Vec<f64>, Vec<f64>) {
    let c = coords(modes, m);
    let n = 2 * m;
    let mu2 = c.iter().map(|&g| mu[g]).collect();
    let sigma2 = c.iter().flat_map(|&a| c.iter().map(move |&b| sigma[a * n + b])).collect();
    (mu2, sigma2)
}

/// smallest eigenvalue of the Hermitian matrix sigma + i (hbar/2) Omega (A.12)
pub fn physicality(sigma: &[f64], m: usize, hbar: f64) -> f64 {
    let n = 2 * m;
    let h = DMatrix::from_fn(n, n, |i, j| {
        let omega = if i < m && j == i + m {
            1.0
        } else if i >= m && j + m == i {
            -1.0
        } else {
            0.0
        };
        nalgebra::Complex::new(sigma[i * n + j], hbar / 2.0 * omega)
    });
    SymmetricEigen::new(h).eigenvalues.iter().cloned().fold(f64::INFINITY, f64::min)
}

/// Terminal homodyne sampling of x on every mode of the state (AN 6.4 with Phi = 0; other angles are rotated
/// into x first): y = mu_x + L z, L the Cholesky factor of sigma_xx, or with a positive semidefinite sigma_xx the
/// symmetric eigendecomposition with eigenvalues in [-1e-12 lambda_max, 0] clamped to 0. Per shot, m normals from
/// Box-Muller pairs (an odd m discards the last pair's sine). Shot i uses stream (seed, stream_start + i).
pub fn sample_homodyne(
    mu: &[f64],
    sigma: &[f64],
    m: usize,
    shots: usize,
    seed: u64,
    stream_start: u64,
    cancel: &dyn Fn() -> bool,
) -> Result<Vec<f64>, Error> {
    let n = 2 * m;
    let cxx = DMatrix::from_fn(m, m, |i, j| sigma[i * n + j]);
    let l: DMatrix<f64> = match cxx.clone().cholesky() {
        Some(ch) => ch.l(),
        None => {
            let e = SymmetricEigen::new(cxx);
            let lmax = e.eigenvalues.iter().cloned().fold(0.0, f64::max);
            let mut root = DMatrix::zeros(m, m);
            for (i, &lam) in e.eigenvalues.iter().enumerate() {
                if lam < -1e-12 * lmax.max(1e-300) {
                    return Err(Error::NonPhysical(lam));
                }
                let r = lam.max(0.0).sqrt();
                for a in 0..m {
                    for b in 0..m {
                        root[(a, b)] += e.eigenvectors[(a, i)] * r * e.eigenvectors[(b, i)];
                    }
                }
            }
            root
        }
    };
    let mut out = Vec::with_capacity(shots * m);
    let mut z = vec![0.0; m + 1];
    for i in 0..shots {
        if i % 1024 == 0 && cancel() {
            return Err(Error::Cancelled);
        }
        let mut st = Stream::new(seed, stream_start + i as u64);
        for p in 0..m.div_ceil(2) {
            let (a, b) = st.normal_pair();
            z[2 * p] = a;
            z[2 * p + 1] = b;
        }
        for a in 0..m {
            out.push(mu[a] + dot((0..m).map(|b| l[(a, b)] * z[b]), m));
        }
    }
    Ok(out)
}

/// Homodyne measurement of x_phi on `mode` with conditioning (FS A.7, FR-G08): the mode rotated by -phi, then
/// y ~ Normal(mu_x, sigma_xx) (the first, cosine-branch normal of stream (seed, stream)) or the fixed outcome, the
/// state conditioned on x = y, and the mode's two rows and columns deleted. Returns (y, mu', sigma') on m - 1 modes.
#[allow(clippy::too_many_arguments)]
pub fn measure_homodyne(
    mu: &[f64],
    sigma: &[f64],
    m: usize,
    mode: usize,
    phi: f64,
    fixed: Option<f64>,
    seed: u64,
    stream: u64,
) -> Result<(f64, Vec<f64>, Vec<f64>), Error> {
    let n = 2 * m;
    if mu.len() != n || sigma.len() != n * n || mode >= m {
        return Err(Error::BadShape("mean, covariance or mode"));
    }
    let (mut mu, mut sigma) = (mu.to_vec(), sigma.to_vec());
    if phi != 0.0 {
        let rot = passive_symplectic(&[C64::new((-phi).cos(), (-phi).sin())], 1);
        apply_symplectic(&mut mu, &mut sigma, m, &rot, &[0.0, 0.0], &[mode], None)?;
    }
    let k = mode;
    let var = sigma[k * n + k];
    if var.is_nan() || var <= 0.0 {
        return Err(Error::NonPhysical(var));
    }
    let mean = mu[k];
    let y = match fixed {
        Some(y) => y,
        None => mean + var.sqrt() * Stream::new(seed, stream).normal_pair().0,
    };
    let col: Vec<f64> = (0..n).map(|i| sigma[i * n + k]).collect();
    for i in 0..n {
        mu[i] += col[i] * (y - mean) / var;
    }
    for i in 0..n {
        for j in 0..n {
            sigma[i * n + j] -= col[i] * col[j] / var;
        }
    }
    let keep: Vec<usize> = (0..m).filter(|&j| j != k).collect();
    let (mu2, sigma2) = marginal(&mu, &sigma, m, &keep);
    Ok((y, mu2, sigma2))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fock::{displacement_matrix, squeezing_matrix};

    /// GT-G02, FS Appendix C: two-mode squeezed vacuum, r = 0.5, x measured on mode 0 with outcome 0.8:
    /// mode 1 has mean x = tanh(2r) y = 0.609275 (up to the sign of the correlation), Var x = 1/cosh 2r = 0.648054,
    /// Var p = cosh 2r = 1.543081
    #[test]
    fn gt_g02() {
        let r = 0.5f64;
        let (c, sh) = ((2.0 * r).cosh(), (2.0 * r).sinh());
        // xxpp: Cov(x0, x1) = sinh 2r, Cov(p0, p1) = -sinh 2r
        let sigma = vec![c, sh, 0.0, 0.0, sh, c, 0.0, 0.0, 0.0, 0.0, c, -sh, 0.0, 0.0, -sh, c];
        let (y, mu, sig) = measure_homodyne(&[0.0; 4], &sigma, 2, 0, 0.0, Some(0.8), 0, 0).unwrap();
        assert_eq!(y, 0.8);
        assert!((mu[0] - 0.609275).abs() < 1e-6 && mu[1] == 0.0);
        assert!((sig[0] - 0.648054).abs() < 1e-6 && (sig[3] - 1.543081).abs() < 1e-6 && sig[1].abs() < 1e-15);
        // measured along p (phi = pi/2): p0 = y gives p1 = -tanh(2r) y
        let (_, mu_p, _) = measure_homodyne(&[0.0; 4], &sigma, 2, 0, std::f64::consts::FRAC_PI_2, Some(0.8), 0, 0).unwrap();
        assert!((mu_p[1] + 0.609275).abs() < 1e-6);
        // sampled outcomes follow Normal(0, cosh 2r)
        let ys: Vec<f64> = (0..20000).map(|k| measure_homodyne(&[0.0; 4], &sigma, 2, 0, 0.0, None, 1, k).unwrap().0).collect();
        let v = ys.iter().map(|y| y * y).sum::<f64>() / ys.len() as f64;
        assert!((v / c - 1.0).abs() < 0.03);
        // a mode without x variance cannot be conditioned on
        assert_eq!(measure_homodyne(&[0.0, 0.0], &[0.0, 0.0, 0.0, 1.0], 1, 0, 0.0, Some(0.0), 0, 0).unwrap_err(), Error::NonPhysical(0.0));
    }

    fn vacuum(m: usize) -> (Vec<f64>, Vec<f64>) {
        let n = 2 * m;
        (vec![0.0; n], (0..n * n).map(|k| if k / n == k % n { 1.0 } else { 0.0 }).collect())
    }

    /// <x>, <p>, Var x, Var p, Cov(x, p) of a single-mode Fock state at hbar = 2 (x = a + a^dag, p = -i(a - a^dag))
    fn fock_moments(psi: &[C64]) -> [f64; 5] {
        let c = psi.len();
        let a = |v: &[C64]| -> Vec<C64> {
            (0..c).map(|n| if n + 1 < c { v[n + 1].scale(((n + 1) as f64).sqrt()) } else { C64::ZERO }).collect()
        };
        let ad = |v: &[C64]| -> Vec<C64> {
            (0..c).map(|n| if n > 0 { v[n - 1].scale((n as f64).sqrt()) } else { C64::ZERO }).collect()
        };
        let x = |v: &[C64]| -> Vec<C64> { a(v).iter().zip(ad(v)).map(|(p, q)| p.add(q)).collect() };
        let p = |v: &[C64]| -> Vec<C64> {
            a(v).iter().zip(ad(v)).map(|(s, t)| C64::new(s.im - t.im, -(s.re - t.re))).collect()
        };
        let ev = |l: &[C64], r: &[C64]| -> f64 { l.iter().zip(r).map(|(u, w)| u.re * w.re + u.im * w.im).sum() };
        let (xv, pv) = (x(psi), p(psi));
        let (mx, mp) = (ev(psi, &xv), ev(psi, &pv));
        let vx = ev(&xv, &xv) - mx * mx;
        let vp = ev(&pv, &pv) - mp * mp;
        // symmetrised covariance: Re <x p> - <x><p>
        let cxp = ev(&xv, &pv) - mx * mp;
        [mx, mp, vx, vp, cxp]
    }

    /// GT-C02: squeezing and displacement on the Gaussian side against the operators of FS 3.3 in Fock space
    #[test]
    fn gt_c02_against_the_fock_operators() {
        for &r in &[0.0f64, 0.5, 1.0] {
            for &phi in &[0.0, std::f64::consts::FRAC_PI_2, std::f64::consts::PI / 6.0] {
                let (c, s) = (r.cosh(), r.sinh());
                let sm = [c - phi.cos() * s, -phi.sin() * s, -phi.sin() * s, c + phi.cos() * s];
                let (mut mu, mut sig) = vacuum(1);
                apply_symplectic(&mut mu, &mut sig, 1, &sm, &[0.0, 0.0], &[0], Some(1e-12)).unwrap();
                let sq = squeezing_matrix(r, phi, 120).unwrap();
                let psi: Vec<C64> = (0..121).map(|k| sq[k * 121]).collect();
                let f = fock_moments(&psi);
                assert!((sig[0] - f[2]).abs() < 1e-10 && (sig[3] - f[3]).abs() < 1e-10, "r = {r}, phi = {phi}");
                assert!((sig[1] - f[4]).abs() < 1e-10, "cov, r = {r}, phi = {phi}: {} vs {}", sig[1], f[4]);
            }
        }
        // displacement alpha = R e^(i Phi): mean (sqrt(2 hbar) Re alpha, sqrt(2 hbar) Im alpha), covariance unchanged
        let alpha = C64::new(0.6, -0.3);
        let dm = displacement_matrix(alpha, 60);
        let psi: Vec<C64> = (0..61).map(|k| dm[k * 61]).collect();
        let f = fock_moments(&psi);
        assert!((f[0] - 2.0 * 0.6).abs() < 1e-12 && (f[1] - 2.0 * -0.3).abs() < 1e-12);
        assert!((f[2] - 1.0).abs() < 1e-12 && (f[3] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn passive_maps_and_checks() {
        // phase shift: x' = cos x - sin p
        let ph = 0.4f64;
        let s = passive_symplectic(&[C64::new(ph.cos(), ph.sin())], 1);
        assert_eq!(s, vec![ph.cos(), -ph.sin(), ph.sin(), ph.cos()]);
        assert!(symplectic_deviation(&s, 1) < 1e-15);
        let (mut mu, mut sig) = vacuum(2);
        assert_eq!(
            apply_symplectic(&mut mu, &mut sig, 2, &[2.0, 0.0, 0.0, 2.0], &[0.0, 0.0], &[1], Some(1e-10)),
            Err(Error::NotSymplectic(3.0))
        );
    }

    #[test]
    fn loss_and_physicality() {
        let (mut mu, mut sig) = vacuum(1);
        mu[0] = 2.0;
        // a squeezed state, then loss: variance eta V + (1 - eta)
        let r = 0.5f64;
        apply_symplectic(&mut mu, &mut sig, 1, &[(-r).exp(), 0.0, 0.0, r.exp()], &[0.0, 0.0], &[0], None).unwrap();
        apply_loss(&mut mu, &mut sig, 1, 0, 0.3, 2.0).unwrap();
        assert!((sig[0] - (0.3 * (-2.0 * r).exp() + 0.7)).abs() < 1e-15);
        assert!((mu[0] - 2.0 * (-r).exp() * 0.3f64.sqrt()).abs() < 1e-15);
        assert!(physicality(&sig, 1, 2.0) > -1e-12);
        // a covariance below the uncertainty bound is not physical
        assert!(physicality(&[0.5, 0.0, 0.0, 0.5], 1, 2.0) < -0.4);
    }

    #[test]
    fn homodyne_samples() {
        // squeezed vacuum, Var x = e^(-2r): sample variance within 1 %
        let r = 0.5f64;
        let (mu, sig) = (vec![0.0, 0.0], vec![(-2.0 * r).exp(), 0.0, 0.0, (2.0 * r).exp()]);
        let y = sample_homodyne(&mu, &sig, 1, 200000, 3, 0, &|| false).unwrap();
        let var: f64 = y.iter().map(|v| v * v).sum::<f64>() / y.len() as f64;
        assert!((var / (-2.0 * r).exp() - 1.0).abs() < 0.01);
        // shots do not depend on the split
        let a = sample_homodyne(&mu, &sig, 1, 10, 3, 4, &|| false).unwrap();
        assert_eq!(a[2..], sample_homodyne(&mu, &sig, 1, 8, 3, 6, &|| false).unwrap()[..]);
    }
}
