//! Boson-sampling samplers (AN 5, FS 5.6, FS 8.5). Shot i of a call uses its own stream (Seed, StreamStart + i),
//! so the samples do not depend on how shots are split into calls (D1).
//!
//! clifford_clifford (AN 5.1, Clifford and Clifford 2018, algorithm B), with the weights corrected: for the k-th
//! photon, with output modes y_1..y_(k-1) chosen and input columns c_1..c_k,
//!   w_j = | sum_l U[j][c_l] perm(U[y_1..y_(k-1); c_1..c_k without c_l]) |^2,
//! the Laplace expansion of the k x k permanent along the new row. AN v0.1 has sum_l |U[j][c_l]|^2 |perm(...)|^2,
//! which samples distinguishable photons (it gives p(1,1) = 1/2 for Hong-Ou-Mandel instead of 0). The draw is
//! normalised by the sum of the weights, which also covers bunched inputs (repeated columns).
//!
//! Stream order per shot: the Fisher-Yates integers for the photons' input columns (listed in mode order), then one
//! uniform per photon for the output draws. brute_force: one uniform.
use crate::num::{Sum, C64};
use crate::perm;
use crate::rng::Stream;

#[derive(Debug, PartialEq)]
pub enum Error {
    BadShape(&'static str),
    /// brute_force: more patterns than max_enumeration
    TooManyPatterns(u64),
    /// pq_bug: the brute-force probabilities do not sum to 1 within 1e-9
    ProbabilitySum(f64),
    Cancelled,
}

fn never() -> bool {
    false
}

/// The categorical draw of AN 5.1: cumulative sums in index order (Neumaier above 64 terms), the first index
/// whose cumulative value exceeds u * total; if rounding leaves none, the last index with a positive weight.
pub fn categorical(w: &[f64], u: f64) -> usize {
    let mut total = Sum::new(w.len() as u64);
    for &x in w {
        total.add(x);
    }
    let target = u * total.value();
    let mut cum = Sum::new(w.len() as u64);
    for (j, &x) in w.iter().enumerate() {
        cum.add(x);
        if cum.value() > target {
            return j;
        }
    }
    w.iter().rposition(|&x| x > 0.0).unwrap_or(0)
}

/// One Clifford-Clifford shot: the output pattern (length m).
pub fn clifford_clifford_shot(u: &[C64], m: usize, input: &[u32], s: &mut Stream) -> Vec<u32> {
    let mut cols: Vec<usize> = (0..m).flat_map(|i| std::iter::repeat_n(i, input[i] as usize)).collect();
    let n = cols.len();
    s.shuffle(&mut cols);
    let mut ys: Vec<usize> = Vec::with_capacity(n);
    let mut w = vec![0.0f64; m];
    for k in 1..=n {
        // p_l = perm of the (k-1) x (k-1) matrix: rows ys, columns c_1..c_k without c_l
        let p: Vec<C64> = (0..k)
            .map(|l| {
                if k == 1 {
                    return C64::ONE;
                }
                let sub: Vec<C64> = ys
                    .iter()
                    .flat_map(|&y| (0..k).filter(move |&c| c != l).map(move |c| (y, c)))
                    .map(|(y, c)| u[y * m + cols[c]])
                    .collect();
                perm::permanent(&sub, k - 1, 1, &never).expect("permanent").value
            })
            .collect();
        for (j, wj) in w.iter_mut().enumerate() {
            let mut acc = C64::ZERO;
            for l in 0..k {
                acc = acc.add(u[j * m + cols[l]].mul(p[l]));
            }
            *wj = acc.re * acc.re + acc.im * acc.im;
        }
        ys.push(categorical(&w, s.uniform()));
    }
    let mut out = vec![0u32; m];
    for y in ys {
        out[y] += 1;
    }
    out
}

/// Every pattern with n photons in m modes, lexicographically ascending (mode 0 most significant).
pub fn patterns(m: usize, n: u32) -> Vec<Vec<u32>> {
    fn rec(m: usize, left: u32, cur: &mut Vec<u32>, out: &mut Vec<Vec<u32>>) {
        if cur.len() == m - 1 {
            cur.push(left);
            out.push(cur.clone());
            cur.pop();
            return;
        }
        for k in 0..=left {
            cur.push(k);
            rec(m, left - k, cur, out);
            cur.pop();
        }
    }
    let mut out = Vec::new();
    rec(m, n, &mut Vec::with_capacity(m), &mut out);
    out
}

fn binomial(n: u64, k: u64) -> u64 {
    (0..k).fold(1u64, |acc, i| acc.saturating_mul(n - i) / (i + 1))
}

/// The exact distribution: every pattern with its probability (FR-B02), checked to sum to 1 within 1e-9.
pub fn distribution(u: &[C64], m: usize, input: &[u32], max_patterns: u64) -> Result<(Vec<Vec<u32>>, Vec<f64>), Error> {
    let n: u32 = input.iter().sum();
    let count = binomial(n as u64 + m as u64 - 1, n as u64);
    if count > max_patterns {
        return Err(Error::TooManyPatterns(count));
    }
    let pats = patterns(m, n);
    let probs: Vec<f64> =
        pats.iter().map(|t| perm::bs_probability(u, m, input, t, 1, &never).expect("probability").0).collect();
    let mut total = Sum::new(probs.len() as u64);
    for &p in &probs {
        total.add(p);
    }
    if (total.value() - 1.0).abs() > 1e-9 {
        return Err(Error::ProbabilitySum(total.value()));
    }
    Ok((pats, probs))
}

pub enum Algorithm {
    CliffordClifford,
    BruteForce { max_enumeration: u64 },
}

/// `shots` samples, shot i from stream (seed, stream_start + i); row-major shots x m.
#[allow(clippy::too_many_arguments)] // the arguments of the NIF
pub fn sample_boson(
    u: &[C64],
    m: usize,
    input: &[u32],
    shots: usize,
    seed: u64,
    stream_start: u64,
    algorithm: Algorithm,
    cancel: &dyn Fn() -> bool,
) -> Result<Vec<u32>, Error> {
    if u.len() != m * m || input.len() != m {
        return Err(Error::BadShape("unitary and input"));
    }
    let mut out = Vec::with_capacity(shots * m);
    match algorithm {
        Algorithm::CliffordClifford => {
            for i in 0..shots {
                if i % 64 == 0 && cancel() {
                    return Err(Error::Cancelled);
                }
                let mut s = Stream::new(seed, stream_start + i as u64);
                out.extend(clifford_clifford_shot(u, m, input, &mut s));
            }
        }
        Algorithm::BruteForce { max_enumeration } => {
            let (pats, probs) = distribution(u, m, input, max_enumeration)?;
            for i in 0..shots {
                if i % 4096 == 0 && cancel() {
                    return Err(Error::Cancelled);
                }
                let mut s = Stream::new(seed, stream_start + i as u64);
                out.extend_from_slice(&pats[categorical(&probs, s.uniform())]);
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hom() -> Vec<C64> {
        let h = std::f64::consts::FRAC_1_SQRT_2;
        vec![C64::new(h, 0.0), C64::new(-h, 0.0), C64::new(h, 0.0), C64::new(h, 0.0)]
    }

    #[test]
    fn hong_ou_mandel_never_gives_one_one() {
        let s = sample_boson(&hom(), 2, &[1, 1], 20000, 7, 0, Algorithm::CliffordClifford, &never).unwrap();
        let ones = s.chunks(2).filter(|t| t == &[1, 1]).count();
        assert_eq!(ones, 0);
        let twos = s.chunks(2).filter(|t| t == &[2, 0]).count();
        assert!((twos as f64 / 20000.0 - 0.5).abs() < 0.02);
    }

    #[test]
    fn patterns_in_order() {
        assert_eq!(patterns(3, 2), vec![
            vec![0, 0, 2], vec![0, 1, 1], vec![0, 2, 0], vec![1, 0, 1], vec![1, 1, 0], vec![2, 0, 0]
        ]);
    }

    #[test]
    fn shots_do_not_depend_on_the_split() {
        let u = hom();
        for alg in [|| Algorithm::CliffordClifford, || Algorithm::BruteForce { max_enumeration: 100 }] {
            let all = sample_boson(&u, 2, &[2, 1], 30, 99, 5, alg(), &never).unwrap();
            let mut parts = sample_boson(&u, 2, &[2, 1], 12, 99, 5, alg(), &never).unwrap();
            parts.extend(sample_boson(&u, 2, &[2, 1], 18, 99, 17, alg(), &never).unwrap());
            assert_eq!(all, parts);
        }
    }

    #[test]
    fn categorical_edges() {
        assert_eq!(categorical(&[0.0, 1.0, 0.0], 0.0), 1);
        assert_eq!(categorical(&[0.5, 0.5, 0.0], 0.999999), 1);
        assert_eq!(categorical(&[0.0, 0.25, 0.75], 0.25), 2);
    }

    #[test]
    fn stream_draws() {
        let mut s = Stream::new(1, 2);
        for _ in 0..1000 {
            let x = s.uniform();
            assert!((0.0..1.0).contains(&x));
            assert!(s.below(3) < 3);
        }
        let mut v: Vec<usize> = (0..10).collect();
        s.shuffle(&mut v);
        let mut w = v.clone();
        w.sort();
        assert_eq!(w, (0..10).collect::<Vec<_>>());
    }
}
