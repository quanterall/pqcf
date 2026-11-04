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
#![forbid(unsafe_code)] // NFR-04: kernels are safe Rust; unsafe only in the Rustler glue (lib.rs)
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

/// Loss at the inputs or the outputs (FR-B05, FS A.8): per-mode transmissivities.
#[derive(Clone, Debug)]
pub enum Loss {
    None,
    Input(Vec<f64>),
    Output(Vec<f64>),
}

/// One lossy shot (AN 5.1 stream order): input loss, one uniform per photon in mode order; the Clifford-Clifford
/// shot on the survivors; output loss, one uniform per sampled photon in output mode order.
pub fn lossy_shot(u: &[C64], m: usize, input: &[u32], loss: &Loss, s: &mut Stream) -> Vec<u32> {
    let kept: Vec<u32> = match loss {
        Loss::Input(eta) => (0..m).map(|i| (0..input[i]).filter(|_| s.uniform() < eta[i]).count() as u32).collect(),
        _ => input.to_vec(),
    };
    let mut out = clifford_clifford_shot(u, m, &kept, s);
    if let Loss::Output(eta) = loss {
        for (j, n) in out.iter_mut().enumerate() {
            *n = (0..*n).filter(|_| s.uniform() < eta[j]).count() as u32;
        }
    }
    out
}

/// One shot with partial distinguishability (Phase 4, FS 15 item 4; docs/decisions/N1-noise-models.md): photon a
/// of input mode i has the internal state sqrt(p_i)|0> + sqrt(1 - p_i)|a>, so <phi_a|phi_b> = sqrt(p_a p_b) and two
/// photons with the same p show the Hong-Ou-Mandel visibility p^2. The internal states |a> are orthogonal and the
/// unitary does not touch them, so the output distribution is exactly the mixture over which photons are in |0>:
/// those are sampled together by Clifford-Clifford, each other one alone from |U[j][i]|^2.
/// Stream order: input loss as in `lossy_shot`; one uniform per surviving photon in mode order (u < p_i: shared);
/// the Clifford-Clifford shot of the shared photons; one uniform per private photon in mode order; output loss.
pub fn noisy_shot(u: &[C64], m: usize, input: &[u32], loss: &Loss, shared: &[f64], s: &mut Stream) -> Vec<u32> {
    let kept: Vec<u32> = match loss {
        Loss::Input(eta) => (0..m).map(|i| (0..input[i]).filter(|_| s.uniform() < eta[i]).count() as u32).collect(),
        _ => input.to_vec(),
    };
    let mut together = vec![0u32; m];
    let mut alone: Vec<usize> = Vec::new();
    for i in 0..m {
        for _ in 0..kept[i] {
            if s.uniform() < shared[i] {
                together[i] += 1;
            } else {
                alone.push(i);
            }
        }
    }
    let mut out = clifford_clifford_shot(u, m, &together, s);
    for i in alone {
        let col: Vec<f64> = (0..m).map(|j| { let z = u[j * m + i]; z.re * z.re + z.im * z.im }).collect();
        out[categorical(&col, s.uniform())] += 1;
    }
    if let Loss::Output(eta) = loss {
        for (j, n) in out.iter_mut().enumerate() {
            *n = (0..*n).filter(|_| s.uniform() < eta[j]).count() as u32;
        }
    }
    out
}

/// The exact probability of output pattern `out` for photons with the Gram matrix `gram` of their internal states
/// (n x n, row-major, photons listed in input-mode order, gram[a][b] = <phi_a|phi_b>):
///   P(T) = sum over permutations pi of (prod_k gram[pi(k)][k]) perm(M o conj(M_pi)) / (prod T! * norm),
/// M = U[rows of T with multiplicity; the photons' input modes], M_pi[r][k] = M[r][pi(k)], and norm the squared
/// norm of the input state, the product over input modes of the permanent of gram on that mode's photons (prod S!
/// for identical photons, 1 for distinguishable ones). n! permanents: an oracle for small n.
pub fn partial_probability(u: &[C64], m: usize, input: &[u32], out: &[u32], gram: &[C64]) -> Result<f64, Error> {
    let cols: Vec<usize> = (0..m).flat_map(|i| std::iter::repeat_n(i, input[i] as usize)).collect();
    let rows: Vec<usize> = (0..m).flat_map(|j| std::iter::repeat_n(j, out[j] as usize)).collect();
    let n = cols.len();
    if input.len() != m || out.len() != m || gram.len() != n * n {
        return Err(Error::BadShape("input, output and an n x n Gram matrix"));
    }
    if rows.len() != n {
        return Ok(0.0);
    }
    let mm: Vec<C64> = (0..n * n).map(|x| u[rows[x / n] * m + cols[x % n]]).collect();
    let conj = |z: C64| C64::new(z.re, -z.im);
    let mut perm_idx: Vec<usize> = (0..n).collect();
    let mut total = crate::num::CSum::new(1 << n.min(20));
    let mut add = |pi: &[usize]| {
        let mut w = C64::ONE;
        for (k, &pk) in pi.iter().enumerate() {
            w = w.mul(gram[pk * n + k]);
        }
        if w.re == 0.0 && w.im == 0.0 {
            return;
        }
        let a: Vec<C64> = (0..n * n).map(|x| mm[x].mul(conj(mm[(x / n) * n + pi[x % n]]))).collect();
        let p = perm::permanent(&a, n, 1, &never).expect("permanent").value;
        total.add(w.mul(p));
    };
    // Heap's algorithm over all n! permutations
    let mut c = vec![0usize; n];
    add(&perm_idx);
    let mut i = 0;
    while i < n {
        if c[i] < i {
            if i % 2 == 0 { perm_idx.swap(0, i) } else { perm_idx.swap(c[i], i) }
            add(&perm_idx);
            c[i] += 1;
            i = 0;
        } else {
            c[i] = 0;
            i += 1;
        }
    }
    // the input state's squared norm: per input mode, the permanent of gram on its photons
    let mut norm = 1.0f64;
    let mut start = 0;
    for &k in input {
        let k = k as usize;
        if k > 1 {
            let sub: Vec<C64> = (0..k * k).map(|x| gram[(start + x / k) * n + start + x % k]).collect();
            norm *= perm::permanent(&sub, k, 1, &never).expect("permanent").value.re;
        }
        start += k;
    }
    Ok(total.value().re / (fact(out) * norm))
}

/// The permanents of the k minors of the (k-1) x k matrix U[rows; cols] (minor l without column l), by one
/// Glynn walk over the rows in Gray-code order: per position the column sums v_j, and for every l the product of
/// the others from prefix and suffix products (no division). O(k 2^k) instead of k Glynn sums (AN 5.1, Phase 3).
pub fn minor_permanents(u: &[C64], m: usize, rows: &[usize], cols: &[usize]) -> Vec<C64> {
    let k = cols.len();
    let r = rows.len();
    if r == 0 {
        return vec![C64::ONE; k];
    }
    let at = |i: usize, j: usize| u[rows[i] * m + cols[j]];
    // delta_0 = +1; delta_i = -1 when bit i-1 of the Gray code is set
    let mut delta = vec![1.0f64; r];
    let mut v: Vec<C64> = (0..k).map(|j| (0..r).fold(C64::ZERO, |acc, i| acc.add(at(i, j)))).collect();
    let positions = 1u64 << (r - 1);
    let mut sums: Vec<crate::num::CSum> = (0..k).map(|_| crate::num::CSum::new(positions)).collect();
    let mut prefix = vec![C64::ONE; k + 1];
    let mut suffix = vec![C64::ONE; k + 1];
    for pos in 0..positions {
        if pos > 0 {
            let i = pos.trailing_zeros() as usize + 1;
            delta[i] = -delta[i];
            for (j, vj) in v.iter_mut().enumerate() {
                *vj = vj.add(at(i, j).scale(2.0 * delta[i]));
            }
        }
        let sign: f64 = delta.iter().product();
        for j in 0..k {
            prefix[j + 1] = prefix[j].mul(v[j]);
        }
        for j in (0..k).rev() {
            suffix[j] = suffix[j + 1].mul(v[j]);
        }
        for (l, sl) in sums.iter_mut().enumerate() {
            sl.add(prefix[l].mul(suffix[l + 1]).scale(sign));
        }
    }
    let norm = 0.5f64.powi((r - 1) as i32);
    sums.iter().map(|sl| sl.value().scale(norm)).collect()
}

/// One Clifford-Clifford shot: the output pattern (length m).
pub fn clifford_clifford_shot(u: &[C64], m: usize, input: &[u32], s: &mut Stream) -> Vec<u32> {
    let mut cols: Vec<usize> = (0..m).flat_map(|i| std::iter::repeat_n(i, input[i] as usize)).collect();
    let n = cols.len();
    s.shuffle(&mut cols);
    let mut ys: Vec<usize> = Vec::with_capacity(n);
    let mut w = vec![0.0f64; m];
    for k in 1..=n {
        // p_l = perm of the (k-1) x (k-1) matrix: rows ys, columns c_1..c_k without c_l, all l at once
        let p = minor_permanents(u, m, &ys, &cols[..k]);
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
    Metropolis { burn_in: usize, thin: usize, chains: usize },
}

fn fact(t: &[u32]) -> f64 {
    t.iter().map(|&k| (1..=k).map(|x| x as f64).product::<f64>()).product()
}

/// AN 5.3: an independence Metropolis sampler with the distinguishable-photon proposal: each photon of input
/// column c lands in mode j with probability |U[j][c]|^2; q(T) = perm(|U|^2_(T,S)) / prod T!. Chain c uses stream
/// (seed, c): the initial proposal (one uniform per photon), then per step the proposal and one uniform for the
/// acceptance; burn_in steps discarded, every thin-th state kept, ceil(shots / chains) per chain, the chains
/// concatenated and cut to shots. Returns the samples and the acceptance rate.
#[allow(clippy::too_many_arguments)]
pub fn metropolis(u: &[C64], m: usize, input: &[u32], shots: usize, seed: u64, burn_in: usize, thin: usize, chains: usize,
                  cancel: &dyn Fn() -> bool) -> Result<(Vec<u32>, f64), Error> {
    let cols: Vec<usize> = (0..m).flat_map(|i| std::iter::repeat_n(i, input[i] as usize)).collect();
    let w: Vec<C64> = u.iter().map(|z| C64::new(z.re * z.re + z.im * z.im, 0.0)).collect();
    let p = |t: &[u32]| perm::bs_probability(u, m, input, t, 1, &never).expect("probability").0;
    let q = |t: &[u32]| perm::permanent_repeated(&w, m, m, t, input, 1, &never).expect("permanent").value.re / fact(t);
    let propose = |s: &mut Stream| -> Vec<u32> {
        let mut t = vec![0u32; m];
        for &c in &cols {
            let col: Vec<f64> = (0..m).map(|j| w[j * m + c].re).collect();
            t[categorical(&col, s.uniform())] += 1;
        }
        t
    };
    let chains = chains.max(1);
    let per = shots.div_ceil(chains);
    let (mut out, mut accepted, mut steps) = (Vec::with_capacity(per * chains * m), 0u64, 0u64);
    for c in 0..chains {
        if cancel() {
            return Err(Error::Cancelled);
        }
        let mut s = Stream::new(seed, c as u64);
        let mut cur = propose(&mut s);
        let (mut pc, mut qc) = (p(&cur), q(&cur));
        let mut kept = 0;
        let mut step = 0usize;
        while kept < per {
            let next = propose(&mut s);
            let (pn, qn) = (p(&next), q(&next));
            let a = s.uniform();
            // min(1, p(T') q(T) / (p(T) q(T'))), written without division
            if pn * qc >= a * pc * qn {
                cur = next;
                pc = pn;
                qc = qn;
                accepted += 1;
            }
            steps += 1;
            step += 1;
            if step > burn_in && (step - burn_in).is_multiple_of(thin.max(1)) {
                out.extend_from_slice(&cur);
                kept += 1;
            }
        }
    }
    out.truncate(shots * m);
    Ok((out, accepted as f64 / steps.max(1) as f64))
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
    sample_boson_lossy(u, m, input, shots, seed, stream_start, algorithm, &Loss::None, cancel)
}

/// sample_boson with loss at the inputs or outputs (FR-B05): clifford_clifford with either, brute_force with output
/// loss (one uniform per sampled photon after the shot's uniform); Metropolis without loss.
#[allow(clippy::too_many_arguments)]
pub fn sample_boson_lossy(
    u: &[C64],
    m: usize,
    input: &[u32],
    shots: usize,
    seed: u64,
    stream_start: u64,
    algorithm: Algorithm,
    loss: &Loss,
    cancel: &dyn Fn() -> bool,
) -> Result<Vec<u32>, Error> {
    sample_boson_noisy(u, m, input, shots, seed, stream_start, algorithm, loss, None, cancel)
}

/// sample_boson_lossy with partial distinguishability (`shared`: per input mode the probability p_i of the shared
/// internal state, see `noisy_shot`); clifford_clifford only.
#[allow(clippy::too_many_arguments)]
pub fn sample_boson_noisy(
    u: &[C64],
    m: usize,
    input: &[u32],
    shots: usize,
    seed: u64,
    stream_start: u64,
    algorithm: Algorithm,
    loss: &Loss,
    shared: Option<&[f64]>,
    cancel: &dyn Fn() -> bool,
) -> Result<Vec<u32>, Error> {
    if let Some(p) = shared {
        if p.len() != m || p.iter().any(|x| !(0.0..=1.0).contains(x)) {
            return Err(Error::BadShape("one shared-state probability in [0, 1] per mode"));
        }
        if !matches!(algorithm, Algorithm::CliffordClifford) {
            return Err(Error::BadShape("partial distinguishability with this algorithm"));
        }
    }
    if u.len() != m * m || input.len() != m {
        return Err(Error::BadShape("unitary and input"));
    }
    match (loss, &algorithm) {
        (Loss::Input(e) | Loss::Output(e), _) if e.len() != m || e.iter().any(|x| !(0.0..=1.0).contains(x)) => {
            return Err(Error::BadShape("one transmissivity in [0, 1] per mode"));
        }
        (Loss::Input(_), Algorithm::BruteForce { .. }) | (Loss::Input(_) | Loss::Output(_), Algorithm::Metropolis { .. }) => {
            return Err(Error::BadShape("loss with this algorithm"));
        }
        _ => {}
    }
    let mut out = Vec::with_capacity(shots * m);
    match algorithm {
        Algorithm::CliffordClifford => {
            for i in 0..shots {
                if i % 64 == 0 && cancel() {
                    return Err(Error::Cancelled);
                }
                let mut s = Stream::new(seed, stream_start + i as u64);
                match shared {
                    Some(p) => out.extend(noisy_shot(u, m, input, loss, p, &mut s)),
                    None => out.extend(lossy_shot(u, m, input, loss, &mut s)),
                }
            }
        }
        Algorithm::Metropolis { burn_in, thin, chains } => {
            return metropolis(u, m, input, shots, seed, burn_in, thin, chains, cancel).map(|(s, _)| s);
        }
        Algorithm::BruteForce { max_enumeration } => {
            let (pats, probs) = distribution(u, m, input, max_enumeration)?;
            for i in 0..shots {
                if i % 4096 == 0 && cancel() {
                    return Err(Error::Cancelled);
                }
                let mut s = Stream::new(seed, stream_start + i as u64);
                let mut t = pats[categorical(&probs, s.uniform())].clone();
                if let Loss::Output(eta) = loss {
                    for (j, n) in t.iter_mut().enumerate() {
                        *n = (0..*n).filter(|_| s.uniform() < eta[j]).count() as u32;
                    }
                }
                out.extend(t);
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

    /// the Gram formula: identical photons give |perm|^2 / (S! T!), distinguishable ones perm(|U|^2) / T!, and
    /// Hong-Ou-Mandel with overlap x gives P(1,1) = (1 - |x|^2) / 2
    #[test]
    fn partial_probability_limits() {
        let u = hom();
        let gram = |x: f64| vec![C64::ONE, C64::new(x, 0.0), C64::new(x, 0.0), C64::ONE];
        for x in [0.0, 0.3, 0.8, 1.0] {
            let p = partial_probability(&u, 2, &[1, 1], &[1, 1], &gram(x)).unwrap();
            assert!((p - (1.0 - x * x) / 2.0).abs() < 1e-15, "{x}: {p}");
        }
        // three photons, two of them in one input mode, on a random-ish 3x3 unitary from two beamsplitters
        let (c, s) = (0.6f64, 0.8f64);
        let b1 = [C64::new(c, 0.0), C64::new(-s, 0.0), C64::ZERO, C64::new(s, 0.0), C64::new(c, 0.0), C64::ZERO, C64::ZERO, C64::ZERO, C64::ONE];
        let (c2, s2) = (0.8f64, 0.6f64);
        let b2 = [C64::ONE, C64::ZERO, C64::ZERO, C64::ZERO, C64::new(c2, 0.0), C64::new(0.0, -s2), C64::ZERO, C64::new(0.0, -s2), C64::new(c2, 0.0)];
        let mut u3 = vec![C64::ZERO; 9];
        for r in 0..3 { for k in 0..3 { for l in 0..3 { u3[r * 3 + k] = u3[r * 3 + k].add(b2[r * 3 + l].mul(b1[l * 3 + k])); } } }
        let ones = vec![C64::ONE; 9];
        let id: Vec<C64> = (0..9).map(|x| if x % 4 == 0 { C64::ONE } else { C64::ZERO }).collect();
        let w: Vec<C64> = u3.iter().map(|z| C64::new(z.re * z.re + z.im * z.im, 0.0)).collect();
        let mut tot = (0.0, 0.0);
        for t in patterns(3, 3) {
            let pi = partial_probability(&u3, 3, &[2, 1, 0], &t, &ones).unwrap();
            let want = perm::bs_probability(&u3, 3, &[2, 1, 0], &t, 1, &never).unwrap().0;
            assert!((pi - want).abs() < 1e-14, "identical {t:?}");
            let pd = partial_probability(&u3, 3, &[2, 1, 0], &t, &id).unwrap();
            let want_d = perm::permanent_repeated(&w, 3, 3, &t, &[2, 1, 0], 1, &never).unwrap().value.re / fact(&t);
            assert!((pd - want_d).abs() < 1e-14, "distinguishable {t:?}");
            tot = (tot.0 + pi, tot.1 + pd);
        }
        assert!((tot.0 - 1.0).abs() < 1e-13 && (tot.1 - 1.0).abs() < 1e-13);
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
    fn lossy_and_metropolis() {
        let u = hom();
        // input loss eta = (1, 0.5): mean photon number 1.5; output loss: the same mean on the outputs
        for loss in [Loss::Input(vec![1.0, 0.5]), Loss::Output(vec![0.75, 0.75])] {
            let s = sample_boson_lossy(&u, 2, &[1, 1], 40000, 3, 0, Algorithm::CliffordClifford, &loss, &never).unwrap();
            let mean = s.iter().sum::<u32>() as f64 / 40000.0;
            assert!((mean - 1.5).abs() < 0.02, "{loss:?}: {mean}");
        }
        // input loss of one photon of HOM: the survivor is in either output with probability 1/2, never (1, 1)
        let s = sample_boson_lossy(&u, 2, &[1, 1], 20000, 3, 0, Algorithm::CliffordClifford, &Loss::Input(vec![1.0, 0.5]), &never).unwrap();
        assert_eq!(s.chunks(2).filter(|t| t == &[1, 1]).count(), 0);
        assert!(matches!(sample_boson_lossy(&u, 2, &[1, 1], 1, 0, 0, Algorithm::BruteForce { max_enumeration: 10 },
                                            &Loss::Input(vec![1.0, 1.0]), &never), Err(Error::BadShape(_))));
        // Metropolis on HOM: the exact distribution has no (1, 1)
        let (ms, rate) = metropolis(&u, 2, &[1, 1], 4000, 7, 100, 10, 8, &never).unwrap();
        assert_eq!(ms.len(), 8000);
        assert_eq!(ms.chunks(2).filter(|t| t == &[1, 1]).count(), 0);
        assert!(rate > 0.0 && rate <= 1.0);
    }

    #[test]
    fn minor_permanents_equal_separate_glynn_sums() {
        let m = 7;
        let mut st = Stream::new(3, 3);
        let u: Vec<C64> = (0..m * m).map(|_| { let (a, b) = st.normal_pair(); C64::new(a, b) }).collect();
        for k in 1..=7 {
            let rows: Vec<usize> = (0..k - 1).map(|i| (3 * i + 1) % m).collect();
            let cols: Vec<usize> = (0..k).map(|j| (2 * j) % m).collect();
            let fast = minor_permanents(&u, m, &rows, &cols);
            let ur = &u;
            for (l, f) in fast.iter().enumerate() {
                let sub: Vec<C64> = rows.iter().flat_map(|&y| cols.iter().enumerate().filter(move |&(c, _)| c != l).map(move |(_, &cc)| ur[y * m + cc])).collect();
                let slow = if k == 1 { C64::ONE } else { perm::permanent(&sub, k - 1, 1, &never).unwrap().value };
                assert!((f.re - slow.re).abs() + (f.im - slow.im).abs() < 1e-12, "k = {k}, l = {l}");
            }
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
