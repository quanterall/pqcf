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

/// FR-F05 trajectory mode: one Kraus operator K_l of the loss channel (A.11) drawn with probability |K_l psi|^2
/// by the uniform u (categorical over l = 0, 1, ...), the state replaced by K_l psi / |K_l psi|. Returns l.
pub fn loss_trajectory(b: &Basis, amps: &mut Vec<C64>, eta: f64, mode: usize, u: f64) -> usize {
    let binom = |n: u32, k: u32| -> f64 { (0..k).fold(1.0, |acc, i| acc * (n - i) as f64 / (i + 1) as f64) };
    let kraus = |l: u32| -> Vec<C64> {
        let mut out = vec![C64::ZERO; amps.len()];
        for (a, p) in amps.iter().zip(&b.patterns) {
            let n = p[mode];
            if n >= l {
                let mut q = p.clone();
                q[mode] = n - l;
                let w = binom(n, l).sqrt() * eta.powf((n - l) as f64 / 2.0) * (1.0 - eta).powf(l as f64 / 2.0);
                let i = b.index(&q);
                out[i] = out[i].add(a.scale(w));
            }
        }
        out
    };
    let states: Vec<Vec<C64>> = (0..=b.cutoff as u32).map(kraus).collect();
    let probs: Vec<f64> = states.iter().map(|v| v.iter().map(|z| z.re * z.re + z.im * z.im).sum()).collect();
    let l = crate::boson::categorical(&probs, u);
    let norm = probs[l].sqrt();
    *amps = states[l].iter().map(|z| if norm > 0.0 { z.scale(1.0 / norm) } else { *z }).collect();
    l
}

/// Density matrices (FS 4.5, FR-F05, AN 7.5): rho is D x D, row-major, in the basis of FS 3.4.
pub mod dm {
    use super::{Basis, Error};
    use crate::num::{Sum, C64};

    fn conj(z: C64) -> C64 {
        C64::new(z.re, -z.im)
    }

    /// rho <- G rho G^dagger for a gate G given by its action on a pure state (`apply` maps an amplitude vector
    /// in place); returns the trace deficit (input trace minus output trace) for the truncated gates
    pub fn unitary_like(b: &Basis, rho: &mut [C64], mut apply: impl FnMut(&mut [C64]) -> Result<(), Error>) -> Result<f64, Error> {
        let d = b.patterns.len();
        let trace_in = trace(rho, d);
        // M = G rho: G on every column
        let mut col = vec![C64::ZERO; d];
        for j in 0..d {
            for i in 0..d {
                col[i] = rho[i * d + j];
            }
            apply(&mut col)?;
            for i in 0..d {
                rho[i * d + j] = col[i];
            }
        }
        // rho' = M G^dagger = (G M^dagger)^dagger: G on every column of M^dagger (the conjugated rows of M)
        for i in 0..d {
            for j in 0..d {
                col[j] = conj(rho[i * d + j]);
            }
            apply(&mut col)?;
            for j in 0..d {
                rho[i * d + j] = conj(col[j]);
            }
        }
        Ok(trace_in - trace(rho, d))
    }

    pub fn trace(rho: &[C64], d: usize) -> f64 {
        let mut s = Sum::new(d as u64);
        for i in 0..d {
            s.add(rho[i * d + i].re);
        }
        s.value()
    }

    /// FS A.11: rho -> sum_l K_l rho K_l^dagger, K_l = sum_(n >= l) sqrt(C(n, l)) eta^((n-l)/2) (1-eta)^(l/2) |n-l><n|
    pub fn loss(b: &Basis, rho: &[C64], eta: f64, mode: usize) -> Vec<C64> {
        let d = b.patterns.len();
        let binom = |n: u32, k: u32| -> f64 { (0..k).fold(1.0, |acc, i| acc * (n - i) as f64 / (i + 1) as f64) };
        let mut out = vec![C64::ZERO; d * d];
        for (i, p) in b.patterns.iter().enumerate() {
            for (j, q) in b.patterns.iter().enumerate() {
                let z = rho[i * d + j];
                if z.re == 0.0 && z.im == 0.0 {
                    continue;
                }
                let (np, nq) = (p[mode], q[mode]);
                let mut pp = p.clone();
                let mut qq = q.clone();
                for l in 0..=np.min(nq) {
                    pp[mode] = np - l;
                    qq[mode] = nq - l;
                    let w = (binom(np, l) * binom(nq, l)).sqrt()
                        * eta.powf((np + nq - 2 * l) as f64 / 2.0)
                        * (1.0 - eta).powi(l as i32);
                    let (a, c) = (b.index(&pp), b.index(&qq));
                    out[a * d + c] = out[a * d + c].add(z.scale(w));
                }
            }
        }
        out
    }

    /// the diagonal as probabilities of `modes` (FR-F03), in the order of FS 3.4 on the measured modes
    pub fn pnr_probabilities(b: &Basis, rho: &[C64], modes: &[usize]) -> Vec<(Vec<u32>, f64)> {
        let d = b.patterns.len();
        let diag: Vec<C64> = (0..d).map(|i| C64::new(rho[i * d + i].re.max(0.0).sqrt(), 0.0)).collect();
        super::pnr_probabilities(b, &diag, modes)
    }

    /// |psi><psi|
    pub fn from_pure(amps: &[C64]) -> Vec<C64> {
        amps.iter().flat_map(|&a| amps.iter().map(move |&c| a.mul(conj(c)))).collect()
    }
}

#[cfg(test)]
mod dm_tests {
    use super::dm;
    use super::*;

    #[test]
    fn loss_of_a_fock_state() {
        // |2> through eta: P(2) = eta^2, P(1) = 2 eta (1 - eta), P(0) = (1 - eta)^2; trace 1
        let b = Basis::new(1, 2);
        let mut psi = vec![C64::ZERO; 3];
        psi[2] = C64::ONE;
        let rho = dm::loss(&b, &dm::from_pure(&psi), 0.7, 0);
        let p = dm::pnr_probabilities(&b, &rho, &[0]);
        assert!((p[2].1 - 0.49).abs() < 1e-15 && (p[1].1 - 0.42).abs() < 1e-15 && (p[0].1 - 0.09).abs() < 1e-15);
        assert!((dm::trace(&rho, 3) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn trajectories_average_to_the_channel() {
        // |2> through eta = 0.7: lost-photon counts 0, 1, 2 with 0.49, 0.42, 0.09
        let b = Basis::new(1, 2);
        let mut counts = [0usize; 3];
        let mut s = crate::rng::Stream::new(4, 0);
        for _ in 0..40000 {
            let mut psi = vec![C64::ZERO, C64::ZERO, C64::ONE];
            let l = loss_trajectory(&b, &mut psi, 0.7, 0, s.uniform());
            counts[l] += 1;
            assert!((psi.iter().map(|z| z.abs().powi(2)).sum::<f64>() - 1.0).abs() < 1e-14);
        }
        for (l, want) in [0.49, 0.42, 0.09].iter().enumerate() {
            assert!((counts[l] as f64 / 40000.0 - want).abs() < 0.01, "l = {l}");
        }
    }

    #[test]
    fn unitaries_agree_with_pure_states() {
        // a beamsplitter on |psi><psi| equals |G psi><G psi|
        let b = Basis::new(2, 3);
        let mut psi = vec![C64::ZERO; b.patterns.len()];
        psi[b.index(&[1, 1])] = C64::new(0.6, 0.0);
        psi[b.index(&[2, 0])] = C64::new(0.0, 0.8);
        let mut rho = dm::from_pure(&psi);
        let def = dm::unitary_like(&b, &mut rho, |v| beamsplitter(&b, v, 0.4, 0.9, 0, 1)).unwrap();
        beamsplitter(&b, &mut psi, 0.4, 0.9, 0, 1).unwrap();
        let want = dm::from_pure(&psi);
        assert!(def.abs() < 1e-15);
        assert!(rho.iter().zip(&want).all(|(x, y)| (x.re - y.re).abs() + (x.im - y.im).abs() < 1e-15));
    }
}
