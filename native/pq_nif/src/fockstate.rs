//! Pure states in the truncated Fock space (FS 3.4, 3.5, 5.5; AN 7): amplitudes in the basis of FS 3.4 (total
//! photon number ascending, then lexicographic from mode 0), D = C(cutoff + m, m).
use crate::fock::{beamsplitter_sector, displacement_matrix, squeezing_matrix};
use crate::mesh;
use crate::num::{Sum, C64};

#[derive(Debug, PartialEq)]
pub enum Error {
    BadShape(&'static str),
    /// pq_bug: a beamsplitter sector's eigenvalue snap failed
    Internal(String),
    /// the squeezing fallback needs more than its size limit
    TooLarge,
}

pub struct Basis {
    pub m: usize,
    pub cutoff: usize,
    pub patterns: Vec<Vec<u32>>,
    binom: Vec<Vec<u64>>,
}

impl Basis {
    pub fn new(m: usize, cutoff: usize) -> Basis {
        let size = cutoff + m + 1;
        let mut binom = vec![vec![0u64; size]; size];
        for n in 0..size {
            binom[n][0] = 1;
            for k in 1..=n {
                binom[n][k] = binom[n - 1][k - 1].saturating_add(binom[n - 1][k]);
            }
        }
        let mut patterns = Vec::new();
        for n in 0..=cutoff as u32 {
            patterns.extend(crate::boson::patterns(m, n));
        }
        Basis { m, cutoff, patterns, binom }
    }

    pub fn dim(m: usize, cutoff: usize) -> u64 {
        // C(cutoff + m, m), saturating
        (0..m as u64).fold(1u64, |acc, i| acc.saturating_mul(cutoff as u64 + 1 + i) / (i + 1))
    }

    fn c(&self, n: usize, k: usize) -> u64 {
        if k > n { 0 } else { self.binom[n][k] }
    }

    /// patterns of `k` modes with total `s`
    fn count(&self, k: usize, s: usize) -> u64 {
        if k == 0 {
            (s == 0) as u64
        } else {
            self.c(s + k - 1, k - 1)
        }
    }

    pub fn index(&self, p: &[u32]) -> usize {
        let n: usize = p.iter().map(|&x| x as usize).sum();
        let mut idx = if n == 0 { 0 } else { self.c(n - 1 + self.m, self.m) };
        let mut left = n;
        for (i, &pi) in p.iter().enumerate().take(self.m - 1) {
            for a in 0..pi as usize {
                idx += self.count(self.m - i - 1, left - a);
            }
            left -= pi as usize;
        }
        idx as usize
    }
}

pub fn phaseshift(b: &Basis, amps: &mut [C64], phi: f64, mode: usize) {
    for (a, p) in amps.iter_mut().zip(&b.patterns) {
        let t = phi * p[mode] as f64;
        *a = a.mul(C64::new(t.cos(), t.sin()));
    }
}

pub fn kerr(b: &Basis, amps: &mut [C64], kappa: f64, mode: usize) {
    for (a, p) in amps.iter_mut().zip(&b.patterns) {
        let n = p[mode] as f64;
        let t = kappa * n * n;
        *a = a.mul(C64::new(t.cos(), t.sin()));
    }
}

/// AN 7.2: the sector matrix U_N on each group of states with the other modes fixed and n_a + n_b = N
pub fn beamsplitter(b: &Basis, amps: &mut [C64], theta: f64, phi: f64, ma: usize, mb: usize) -> Result<(), Error> {
    let sectors: Vec<Vec<C64>> = (0..=b.cutoff)
        .map(|n| beamsplitter_sector(n, theta, phi).map_err(|e| Error::Internal(format!("{e:?}"))))
        .collect::<Result<_, _>>()?;
    for rep in b.patterns.iter().filter(|p| p[mb] == 0) {
        let n = rep[ma] as usize;
        if n == 0 {
            continue;
        }
        let mut q = rep.clone();
        let idx: Vec<usize> = (0..=n)
            .map(|j| {
                q[ma] = j as u32;
                q[mb] = (n - j) as u32;
                b.index(&q)
            })
            .collect();
        let v: Vec<C64> = idx.iter().map(|&i| amps[i]).collect();
        let u = &sectors[n];
        for (r, &i) in idx.iter().enumerate() {
            let mut acc = C64::ZERO;
            for (c, vc) in v.iter().enumerate() {
                acc = acc.add(u[r * (n + 1) + c].mul(*vc));
            }
            amps[i] = acc;
        }
    }
    Ok(())
}

/// AN 7.3: the Clements mesh of U on `modes`, then the mesh's phases
pub fn passive(b: &Basis, amps: &mut [C64], u: &[C64], modes: &[usize]) -> Result<(), Error> {
    let k = modes.len();
    let (bs, phases) = mesh::decompose(u, k);
    for (m, t, p) in bs {
        beamsplitter(b, amps, t, p, modes[m], modes[m + 1])?;
    }
    for (j, &a) in phases.iter().enumerate() {
        phaseshift(b, amps, a, modes[j]);
    }
    Ok(())
}

/// A single-mode matrix G ((cutoff+1) x (cutoff+1)) on `mode`; components beyond the cutoff are dropped; returns
/// the norm deficit (input norm minus output norm, FS 3.5)
pub fn single_mode(b: &Basis, amps: &mut [C64], g: &[C64], mode: usize) -> f64 {
    let d = b.cutoff + 1;
    let mut norm_in = Sum::new(amps.len() as u64);
    for a in amps.iter() {
        norm_in.add(a.re * a.re + a.im * a.im);
    }
    for rep in b.patterns.iter().filter(|p| p[mode] == 0) {
        let others: usize = rep.iter().map(|&x| x as usize).sum();
        let top = b.cutoff - others;
        let mut q = rep.clone();
        let idx: Vec<usize> = (0..=top)
            .map(|n| {
                q[mode] = n as u32;
                b.index(&q)
            })
            .collect();
        let v: Vec<C64> = idx.iter().map(|&i| amps[i]).collect();
        for (r, &i) in idx.iter().enumerate() {
            let mut acc = C64::ZERO;
            for (c, vc) in v.iter().enumerate() {
                acc = acc.add(g[r * d + c].mul(*vc));
            }
            amps[i] = acc;
        }
    }
    let mut norm_out = Sum::new(amps.len() as u64);
    for a in amps.iter() {
        norm_out.add(a.re * a.re + a.im * a.im);
    }
    norm_in.value() - norm_out.value()
}

pub fn displacement(b: &Basis, amps: &mut [C64], r: f64, phi: f64, mode: usize) -> f64 {
    let g = displacement_matrix(C64::new(r * phi.cos(), r * phi.sin()), b.cutoff);
    single_mode(b, amps, &g, mode)
}

pub fn squeezing(b: &Basis, amps: &mut [C64], r: f64, phi: f64, mode: usize) -> Result<f64, Error> {
    let g = squeezing_matrix(r, phi, b.cutoff).ok_or(Error::TooLarge)?;
    Ok(single_mode(b, amps, &g, mode))
}

/// FR-F03: probabilities of the patterns of `modes`, marginalised over the others, in the order of FS 3.4 on the
/// measured modes
pub fn pnr_probabilities(b: &Basis, amps: &[C64], modes: &[usize]) -> Vec<(Vec<u32>, f64)> {
    let mut acc: std::collections::BTreeMap<(u32, Vec<u32>), Sum> = std::collections::BTreeMap::new();
    for (a, p) in amps.iter().zip(&b.patterns) {
        let key: Vec<u32> = modes.iter().map(|&j| p[j]).collect();
        let total = key.iter().sum();
        acc.entry((total, key)).or_insert_with(|| Sum::new(amps.len() as u64)).add(a.re * a.re + a.im * a.im);
    }
    acc.into_iter().map(|((_, k), s)| (k, s.value())).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::perm;

    fn basis_state(b: &Basis, p: &[u32]) -> Vec<C64> {
        let mut v = vec![C64::ZERO; b.patterns.len()];
        v[b.index(p)] = C64::ONE;
        v
    }

    fn close(a: C64, b: C64, tol: f64) -> bool {
        (a.re - b.re).abs() + (a.im - b.im).abs() <= tol
    }

    #[test]
    fn index_is_the_position() {
        for (m, c) in [(1, 5), (2, 4), (3, 3), (5, 2)] {
            let b = Basis::new(m, c);
            assert_eq!(b.patterns.len() as u64, Basis::dim(m, c));
            for (i, p) in b.patterns.iter().enumerate() {
                assert_eq!(b.index(p), i);
            }
        }
    }

    /// GT-C01: beamsplitter and phase conventions of FS 3.3 on one and two photons
    #[test]
    fn gt_c01() {
        let b = Basis::new(2, 2);
        let (t, p) = (0.37, 1.1);
        let mut v = basis_state(&b, &[1, 0]);
        beamsplitter(&b, &mut v, t, p, 0, 1).unwrap();
        assert!(close(v[b.index(&[1, 0])], C64::new(t.cos(), 0.0), 1e-15));
        assert!(close(v[b.index(&[0, 1])], C64::new(p.cos() * t.sin(), p.sin() * t.sin()), 1e-15));
        // phase shift: |n> -> e^(i phi n) |n>
        let mut w = basis_state(&b, &[2, 0]);
        phaseshift(&b, &mut w, 0.3, 0);
        assert!(close(w[b.index(&[2, 0])], C64::new(0.6f64.cos(), 0.6f64.sin()), 1e-15));
        // two photons, 50:50: |1,1> -> (e^(-i phi)|2,0> ... ) with no |1,1> (Hong-Ou-Mandel)
        let mut h = basis_state(&b, &[1, 1]);
        beamsplitter(&b, &mut h, std::f64::consts::FRAC_PI_4, 0.0, 0, 1).unwrap();
        assert!(h[b.index(&[1, 1])].abs() < 1e-15);
        assert!((h[b.index(&[2, 0])].abs().powi(2) - 0.5).abs() < 1e-15);
    }

    /// GT-C04: Kerr |n> -> e^(i kappa n^2) |n>, norm preserved
    #[test]
    fn gt_c04() {
        let b = Basis::new(1, 6);
        let mut v: Vec<C64> = (0..7).map(|n| C64::new(1.0 / 7f64.sqrt(), 0.0).scale(if n % 2 == 0 { 1.0 } else { -1.0 })).collect();
        let before = v.clone();
        kerr(&b, &mut v, 0.2, 0);
        for n in 0..7 {
            let t = 0.2 * (n * n) as f64;
            assert!(close(v[n], before[n].mul(C64::new(t.cos(), t.sin())), 1e-15));
        }
        let norm: f64 = v.iter().map(|z| z.abs().powi(2)).sum();
        assert!((norm - 1.0).abs() < 1e-15);
    }

    /// a passive unitary through the mesh gives <T|U|S> = perm(U_T,S) / sqrt(S! T!)
    #[test]
    fn passive_amplitudes_are_permanents() {
        let k = 3;
        let mut s = crate::rng::Stream::new(5, 0);
        let z = nalgebra::DMatrix::from_fn(k, k, |_, _| {
            let (a, b) = s.normal_pair();
            nalgebra::Complex::new(a, b)
        });
        let q = z.qr().q();
        let u: Vec<C64> = (0..k * k).map(|i| C64::new(q[(i / k, i % k)].re, q[(i / k, i % k)].im)).collect();
        let b = Basis::new(k, 3);
        let input = [2u32, 1, 0];
        let mut v = basis_state(&b, &input);
        passive(&b, &mut v, &u, &[0, 1, 2]).unwrap();
        let fact = |x: &[u32]| x.iter().map(|&n| (1..=n).product::<u32>() as f64).product::<f64>();
        for out in crate::boson::patterns(3, 3) {
            let p = perm::permanent_repeated(&u, 3, 3, &out, &input, 1, &|| false).unwrap().value;
            let want = p.scale(1.0 / (fact(&input) * fact(&out)).sqrt());
            assert!(close(v[b.index(&out)], want, 1e-13), "{out:?}");
        }
    }

    #[test]
    fn displacement_of_vacuum_and_the_deficit() {
        let b = Basis::new(2, 30);
        let mut v = basis_state(&b, &[0, 0]);
        let d = displacement(&b, &mut v, 0.8, 0.3, 1);
        assert!(d.abs() < 1e-14); // the Poisson tail beyond 30 photons is far below rounding
        // P(n = 0) on mode 1 = exp(-|alpha|^2)
        let probs = pnr_probabilities(&b, &v, &[1]);
        assert_eq!(probs[0].0, vec![0]);
        assert!((probs[0].1 - (-0.64f64).exp()).abs() < 1e-15);
        // a small cutoff loses norm, and says so
        let b2 = Basis::new(1, 3);
        let mut w = basis_state(&b2, &[0]);
        let deficit = displacement(&b2, &mut w, 1.5, 0.0, 0);
        let kept: f64 = (0..4).map(|n| (-2.25f64).exp() * 2.25f64.powi(n) / (1..=n as u32).product::<u32>().max(1) as f64).sum();
        assert!((deficit - (1.0 - kept)).abs() < 1e-14);
    }

    #[test]
    fn squeezed_vacuum_probabilities() {
        // P(0) = 1 / cosh r, P(2) = tanh^2 r / (2 cosh r)
        let b = Basis::new(1, 40);
        let mut v = basis_state(&b, &[0]);
        squeezing(&b, &mut v, 0.6, 0.2, 0).unwrap();
        let p = pnr_probabilities(&b, &v, &[0]);
        let ch = 0.6f64.cosh();
        assert!((p[0].1 - 1.0 / ch).abs() < 1e-14);
        assert!((p[2].1 - 0.6f64.tanh().powi(2) / (2.0 * ch)).abs() < 1e-14);
        assert!(p[1].1 < 1e-30);
    }
}
