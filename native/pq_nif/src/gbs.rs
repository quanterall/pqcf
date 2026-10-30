//! Gaussian boson sampling (FS 5.7, A.9, A.10; AN 6): photon-number and threshold probabilities of a Gaussian
//! state, and the chain-rule samplers.
//!
//! Conventions of thewalrus at the pinned version (FS A.9), with an xxpp covariance:
//!   Q   = [[aidaj, conj(aiaj)], [aiaj, conj(aidaj)]] + I,   x, xp, p the blocks of Sigma times 2/hbar,
//!         aidaj = (x + p + i (xp - xp^T) - 2 I) / 4,  aiaj = (x - p + i (xp + xp^T)) / 4
//!   A   = X conj(I - Q^-1),   X = [[0, I], [I, 0]]
//!   beta = (alpha, conj(alpha)),  alpha = (mu_x + i mu_p) / sqrt(2 hbar);   gamma = conj(beta) - A beta
//!   P(S) = lhaf(A_S, gamma_S) / (prod S! Norm),  Norm = sqrt(det Q) exp(1/2 beta^T Q^-1 conj(beta))
//!   threshold, zero mean: P(B) = tor(O_B) / sqrt(det Q),  O = I - Q^-1
use crate::hafnian::hafnian_repeated;
use crate::num::{Sum, C64};
use crate::rng::Stream;
use crate::toronto::torontonian;
use nalgebra::{Complex, DMatrix, DVector};

#[derive(Debug, PartialEq)]
pub enum Error {
    BadShape(&'static str),
    /// a probability below -1e-12 (FR-GB06)
    NegativeProbability(f64),
    /// Q not invertible or det Q not positive: not a physical state
    NonPhysical(&'static str),
    /// threshold detection of a displaced state needs the loop torontonian (Phase 3)
    DisplacedThreshold,
    Cancelled,
}

type Cm = DMatrix<Complex<f64>>;

/// A Gaussian state in the complex form of A.9
pub struct Amat {
    pub m: usize,
    pub a: Cm,
    pub gamma: DVector<Complex<f64>>,
    /// sqrt(det Q) exp(1/2 beta^T Q^-1 conj(beta))
    pub norm: f64,
    pub qinv: Cm,
    pub sqrt_det_q: f64,
    pub displaced: bool,
}

pub fn amat(mu: &[f64], sigma: &[f64], m: usize, hbar: f64) -> Result<Amat, Error> {
    let n = 2 * m;
    if mu.len() != n || sigma.len() != n * n {
        return Err(Error::BadShape("mean and covariance"));
    }
    let s = |i: usize, j: usize| sigma[i * n + j] * 2.0 / hbar;
    let mut q = Cm::zeros(n, n);
    for i in 0..m {
        for j in 0..m {
            let (x, p) = (s(i, j), s(i + m, j + m));
            let (xp, px) = (s(i, j + m), s(j, i + m)); // xp[i][j] and xp[j][i] = (xp^T)[i][j]
            let id = if i == j { 1.0 } else { 0.0 };
            let aidaj = Complex::new((x + p - 2.0 * id) / 4.0, (xp - px) / 4.0);
            let aiaj = Complex::new((x - p) / 4.0, (xp + px) / 4.0);
            q[(i, j)] = aidaj + Complex::new(id, 0.0);
            q[(i, j + m)] = aiaj.conj();
            q[(i + m, j)] = aiaj;
            q[(i + m, j + m)] = aidaj.conj() + Complex::new(id, 0.0);
        }
    }
    let det = q.clone().determinant();
    if det.re.is_nan() || det.re <= 0.0 {
        return Err(Error::NonPhysical("det Q"));
    }
    let qinv = q.try_inverse().ok_or(Error::NonPhysical("Q singular"))?;
    let o = Cm::identity(n, n) - &qinv;
    let mut a = Cm::zeros(n, n);
    for i in 0..n {
        let xi = if i < m { i + m } else { i - m };
        for j in 0..n {
            a[(i, j)] = o[(xi, j)].conj();
        }
    }
    let k = (2.0 * hbar).sqrt();
    let alpha: Vec<Complex<f64>> = (0..m).map(|j| Complex::new(mu[j] / k, mu[j + m] / k)).collect();
    let beta = DVector::from_iterator(n, alpha.iter().cloned().chain(alpha.iter().map(|z| z.conj())));
    let gamma = beta.map(|z| z.conj()) - &a * &beta;
    let expo = (beta.transpose() * &qinv * beta.map(|z| z.conj()))[(0, 0)].re;
    let sqrt_det_q = det.re.sqrt();
    let displaced = mu.iter().any(|&x| x != 0.0);
    Ok(Amat { m, a, gamma, norm: sqrt_det_q * (0.5 * expo).exp(), qinv, sqrt_det_q, displaced })
}

fn c64(z: Complex<f64>) -> C64 {
    C64::new(z.re, z.im)
}

/// P(S) for a photon-number pattern S (A.9); negative values below -1e-12 are errors, smaller ones clipped
pub fn pnr_probability(am: &Amat, s: &[u32], cancel: &(dyn Fn() -> bool + Sync)) -> Result<f64, Error> {
    let n = 2 * am.m;
    let a: Vec<C64> = (0..n * n).map(|k| c64(am.a[(k / n, k % n)])).collect();
    let mult: Vec<u32> = s.iter().chain(s.iter()).copied().collect();
    let loops: Option<Vec<C64>> = am.displaced.then(|| am.gamma.iter().map(|&z| c64(z)).collect());
    let h = hafnian_repeated(&a, n, loops.as_deref(), &mult, 1, cancel).map_err(|e| match e {
        crate::hafnian::Error::Cancelled => Error::Cancelled,
        _ => Error::BadShape("hafnian"),
    })?;
    let fact: f64 = s.iter().map(|&k| (1..=k).map(|x| x as f64).product::<f64>()).product();
    let p = h.value.re / (fact * am.norm);
    check_sign(p, h.scale / (fact * am.norm))
}

/// FR-GB06 with the error scale of A-1: a probability below -max(1e-12, 1e-10 S) is an error, where S is the
/// kernel's scale normalised like the probability; anything above is rounding and is clipped to 0.
fn check_sign(p: f64, scale: f64) -> Result<f64, Error> {
    if p < -f64::max(1e-12, 1e-10 * scale) {
        return Err(Error::NegativeProbability(p));
    }
    Ok(p.max(0.0))
}

/// P(B) for a click pattern B, zero-mean states (A.10)
pub fn threshold_probability(am: &Amat, b: &[u32], cancel: &(dyn Fn() -> bool + Sync)) -> Result<f64, Error> {
    if am.displaced {
        return Err(Error::DisplacedThreshold);
    }
    let m = am.m;
    let n = 2 * m;
    let clicked: Vec<usize> = (0..m).filter(|&j| b[j] > 0).collect();
    let idx: Vec<usize> = clicked.iter().copied().chain(clicked.iter().map(|&j| j + m)).collect();
    let d = idx.len();
    let o: Vec<C64> = idx
        .iter()
        .flat_map(|&i| idx.iter().map(move |&j| (i, j)))
        .map(|(i, j)| {
            let id = if i == j { 1.0 } else { 0.0 };
            C64::new(id - am.qinv[(i, j)].re, -am.qinv[(i, j)].im)
        })
        .collect();
    let _ = n;
    let t = torontonian(&o, d, 1, cancel).map_err(|e| match e {
        crate::toronto::Error::Cancelled => Error::Cancelled,
        crate::toronto::Error::NotPositiveDefinite => Error::NonPhysical("I - O"),
        _ => Error::BadShape("torontonian"),
    })?;
    check_sign(t.value / am.sqrt_det_q, t.scale / am.sqrt_det_q)
}

/// the marginal state of modes 0..k (xxpp)
fn prefix(mu: &[f64], sigma: &[f64], m: usize, k: usize) -> (Vec<f64>, Vec<f64>) {
    let modes: Vec<usize> = (0..k).collect();
    crate::gaussian::marginal(mu, sigma, m, &modes)
}

pub struct Caps {
    pub max_photons: u32,
    pub max_per_mode: u32,
    pub max_clicks: u32,
}

/// The chain-rule samplers (AN 6.2, 6.3), mode by mode with one uniform u per mode. The candidates (photon
/// counts within the caps, or click / no click) are evaluated lazily in order until the cumulative probability
/// exceeds u p_(k-1), p_(k-1) the probability of the prefix (the untruncated conditional). If the caps end the
/// candidates first (probability t_k = 1 - Z / p_(k-1)), u is rescaled to (u p_(k-1) - Z) / (p_(k-1) - Z) and drawn
/// over the allowed candidates by Z: every allowed outcome j then has probability p_j / Z, the truncated and
/// renormalised conditional. The returned lost mass is the mean over shots of the number of such steps, an unbiased
/// estimate of the mean of sum_k t_k (FS 5.7 truncated_mass). Row-major shots x m.
#[allow(clippy::too_many_arguments)]
pub fn sample_chain_rule(
    mu: &[f64],
    sigma: &[f64],
    m: usize,
    hbar: f64,
    threshold: bool,
    caps: &Caps,
    shots: usize,
    seed: u64,
    stream_start: u64,
    cancel: &(dyn Fn() -> bool + Sync),
) -> Result<(Vec<u32>, f64), Error> {
    let states: Vec<Amat> = (1..=m)
        .map(|k| {
            let (mk, sk) = prefix(mu, sigma, m, k);
            amat(&mk, &sk, k, hbar)
        })
        .collect::<Result<_, _>>()?;
    if threshold && states[m - 1].displaced {
        return Err(Error::DisplacedThreshold);
    }
    let mut out = Vec::with_capacity(shots * m);
    let mut beyond_total = 0u64;
    for i in 0..shots {
        if cancel() {
            return Err(Error::Cancelled);
        }
        let mut st = Stream::new(seed, stream_start + i as u64);
        let mut pattern: Vec<u32> = Vec::with_capacity(m);
        let mut p_prev = 1.0;
        for state in &states {
            let u = st.uniform();
            if p_prev <= 0.0 {
                pattern.push(0);
                continue;
            }
            let used: u32 = pattern.iter().sum();
            let clicks: u32 = pattern.iter().filter(|&&x| x > 0).count() as u32;
            let candidates: Vec<u32> = if threshold {
                if clicks >= caps.max_clicks { vec![0] } else { vec![0, 1] }
            } else {
                (0..=caps.max_per_mode.min(caps.max_photons.saturating_sub(used))).collect()
            };
            let target = u * p_prev;
            let mut probs: Vec<f64> = Vec::with_capacity(candidates.len());
            let mut cum = Sum::new(candidates.len() as u64);
            let mut chosen = None;
            for (j, &c) in candidates.iter().enumerate() {
                pattern.push(c);
                let p = if threshold {
                    threshold_probability(state, &pattern, cancel)?
                } else {
                    pnr_probability(state, &pattern, cancel)?
                };
                pattern.pop();
                probs.push(p);
                cum.add(p);
                if cum.value() > target {
                    chosen = Some(j);
                    break;
                }
            }
            let j = match chosen {
                Some(j) => j,
                None => {
                    let z = cum.value();
                    if z <= 0.0 {
                        pattern.push(0);
                        p_prev = 0.0;
                        continue;
                    }
                    beyond_total += 1;
                    let u2 = ((target - z) / (p_prev - z)).clamp(0.0, 1.0 - f64::EPSILON);
                    crate::boson::categorical(&probs, u2)
                }
            };
            pattern.push(candidates[j]);
            p_prev = probs[j];
        }
        out.extend(pattern);
    }
    Ok((out, if shots == 0 { 0.0 } else { beyond_total as f64 / shots as f64 }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn never() -> bool {
        false
    }

    #[test]
    fn coherent_state_gt_h02() {
        // |alpha>, hbar = 2: P(0) = exp(-|alpha|^2), P(1) = |alpha|^2 exp(-|alpha|^2)
        let (ar, ai) = (0.7, 0.2);
        let mu = [2.0 * ar, 2.0 * ai];
        let sigma = [1.0, 0.0, 0.0, 1.0];
        let am = amat(&mu, &sigma, 1, 2.0).unwrap();
        let x = ar * ar + ai * ai;
        assert!((pnr_probability(&am, &[0], &never).unwrap() - (-x).exp()).abs() < 1e-15);
        assert!((pnr_probability(&am, &[1], &never).unwrap() - x * (-x).exp()).abs() < 1e-15);
        assert!((pnr_probability(&am, &[3], &never).unwrap() - x.powi(3) / 6.0 * (-x).exp()).abs() < 1e-15);
        assert_eq!(threshold_probability(&am, &[1], &never), Err(Error::DisplacedThreshold));
    }

    #[test]
    fn squeezed_vacuum() {
        // P(click) = 1 - 1/cosh r; P(2) = tanh^2 r / (2 cosh r)
        let r = 0.6f64;
        let sigma = [(-2.0 * r).exp(), 0.0, 0.0, (2.0 * r).exp()];
        let am = amat(&[0.0, 0.0], &sigma, 1, 2.0).unwrap();
        assert!((threshold_probability(&am, &[1], &never).unwrap() - (1.0 - 1.0 / r.cosh())).abs() < 1e-14);
        assert!((pnr_probability(&am, &[2], &never).unwrap() - r.tanh().powi(2) / (2.0 * r.cosh())).abs() < 1e-14);
        assert!(pnr_probability(&am, &[1], &never).unwrap() < 1e-15);
    }

    #[test]
    fn chain_rule_sampler() {
        // squeezed vacuum: photon counts are even, P(0) = 1/cosh r
        let r = 0.6f64;
        let sigma = [(-2.0 * r).exp(), 0.0, 0.0, (2.0 * r).exp()];
        let caps = Caps { max_photons: 16, max_per_mode: 16, max_clicks: 20 };
        let (s, lost) = sample_chain_rule(&[0.0, 0.0], &sigma, 1, 2.0, false, &caps, 20000, 1, 0, &never).unwrap();
        assert!(s.iter().all(|&n| n % 2 == 0));
        let f0 = s.iter().filter(|&&n| n == 0).count() as f64 / 20000.0;
        assert!((f0 - 1.0 / r.cosh()).abs() < 0.01);
        assert!(lost < 1e-3);
        // tight caps: the lost-mass estimate is near the true 1 - P(0) - P(2), and the samples follow
        // the renormalised distribution P(0) / (P(0) + P(2))
        let tight = Caps { max_photons: 2, max_per_mode: 2, max_clicks: 20 };
        let (s2, lost2) = sample_chain_rule(&[0.0, 0.0], &sigma, 1, 2.0, false, &tight, 20000, 1, 0, &never).unwrap();
        let (p0, p2) = (1.0 / r.cosh(), r.tanh().powi(2) / (2.0 * r.cosh()));
        assert!((lost2 - (1.0 - p0 - p2)).abs() < 0.01, "{lost2}");
        let g0 = s2.iter().filter(|&&n| n == 0).count() as f64 / 20000.0;
        assert!((g0 - p0 / (p0 + p2)).abs() < 0.01, "{g0}");
    }
}
