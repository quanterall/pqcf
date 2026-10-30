//! Hafnian and loop hafnian (AN 4.1): the power-trace formula of Bjorklund, Gupt and Quesada (2019).
//!
//! For a symmetric M x M matrix A, M = 2n, with the fixed pairing (2i, 2i+1) and X the permutation that swaps the
//! two members of each pair:
//!
//!   haf(A)  = sum over subsets Z of the n pairs of (-1)^(n - |Z|) [eta^n] exp( sum_k g_k eta^k ),
//!   g_k     = tr(E^k) / (2k)                                  E = (A X) restricted to the indices of Z
//!   lhaf(A) : g_k += (1/2) v^T X E^(k-1) v                    v = the loop weights restricted to Z
//!
//! tr(E^k) for k = 1..n from the characteristic polynomial (Hessenberg reduction by Householder reflections, then
//! La Budde's recurrence) and Newton's identities: O(|Z|^3) per subset, O(M^3 2^(M/2)) in all. The coefficient of
//! eta^n by the power series of the exponential: h_0 = 1, h_m = (1/m) sum_(k=1..m) k g_k h_(m-k).
//!
//! D1: subsets are visited in Gray-code order of their index p (subset = bits of p ^ (p >> 1)), in blocks of
//! 2^12, block sums added in block order (FS 8.5). Odd M: haf = 0; lhaf pads one vertex with loop weight 1.
use crate::num::{CSum, Sum, C64};

pub const BLOCK: u64 = 1 << 12;
pub const MAX_PAIRS: usize = 40;

#[derive(Debug, PartialEq)]
pub enum Error {
    BadShape(&'static str),
    TooLarge(usize),
    Cancelled,
    NonFinite(&'static str),
}

#[derive(Debug, Clone, Copy)]
pub struct HafResult {
    pub value: C64,
    /// the scale S of AN 2.4: sum of |terms|
    pub scale: f64,
    pub subsets: u64,
}

fn sub(a: C64, b: C64) -> C64 {
    C64::new(a.re - b.re, a.im - b.im)
}

fn conj(a: C64) -> C64 {
    C64::new(a.re, -a.im)
}

fn div_real(a: C64, x: f64) -> C64 {
    C64::new(a.re / x, a.im / x)
}

/// tr(E^k) for k = 1..=kmax, E a d x d matrix (row-major), via the characteristic polynomial.
pub fn power_traces(e: &[C64], d: usize, kmax: usize) -> Vec<C64> {
    if d == 0 {
        return vec![C64::ZERO; kmax];
    }
    let h = hessenberg(e, d);
    let c = charpoly_la_budde(&h, d); // c[0] = 1, p(x) = x^d + c[1] x^(d-1) + ... + c[d]
    let mut s = vec![C64::ZERO; kmax + 1];
    for k in 1..=kmax {
        let mut acc = C64::ZERO;
        if k <= d {
            acc = c[k].scale(-(k as f64));
        }
        // Newton: s_k = -k c_k - sum_(i<k) c_i s_(k-i) for k <= d; s_k = -sum_(i<=d) c_i s_(k-i) beyond
        for i in 1..k.min(d + 1) {
            acc = sub(acc, c[i].mul(s[k - i]));
        }
        s[k] = acc;
    }
    s[1..].to_vec()
}

/// Reduction to upper Hessenberg form by Householder reflections (similarity, so the spectrum is kept).
fn hessenberg(e: &[C64], d: usize) -> Vec<C64> {
    let mut h = e.to_vec();
    for k in 0..d.saturating_sub(2) {
        let col: Vec<C64> = (k + 1..d).map(|i| h[i * d + k]).collect();
        let norm = col.iter().map(|z| z.re * z.re + z.im * z.im).sum::<f64>().sqrt();
        if norm == 0.0 {
            continue;
        }
        let x0 = col[0];
        let a0 = x0.abs();
        let phase = if a0 == 0.0 { C64::ONE } else { div_real(x0, a0) };
        // v = x + e^(i arg x0) |x| e_1, normalised
        let mut v = col.clone();
        v[0] = x0.add(phase.scale(norm));
        let vn = v.iter().map(|z| z.re * z.re + z.im * z.im).sum::<f64>().sqrt();
        if vn == 0.0 {
            continue;
        }
        for z in v.iter_mut() {
            *z = div_real(*z, vn);
        }
        // H <- (I - 2 v v^dagger) H on rows k+1.., then H <- H (I - 2 v v^dagger) on columns k+1..
        for j in 0..d {
            let mut dot = C64::ZERO;
            for (t, vt) in v.iter().enumerate() {
                dot = dot.add(conj(*vt).mul(h[(k + 1 + t) * d + j]));
            }
            for (t, vt) in v.iter().enumerate() {
                let idx = (k + 1 + t) * d + j;
                h[idx] = sub(h[idx], vt.mul(dot).scale(2.0));
            }
        }
        for i in 0..d {
            let mut dot = C64::ZERO;
            for (t, vt) in v.iter().enumerate() {
                dot = dot.add(h[i * d + k + 1 + t].mul(*vt));
            }
            for (t, vt) in v.iter().enumerate() {
                let idx = i * d + k + 1 + t;
                h[idx] = sub(h[idx], dot.mul(conj(*vt)).scale(2.0));
            }
        }
    }
    h
}

/// La Budde's recurrence for the characteristic polynomial of an upper Hessenberg matrix; the coefficients
/// c[0..=d] of x^d, x^(d-1), ..., x^0 (c[0] = 1).
fn charpoly_la_budde(h: &[C64], d: usize) -> Vec<C64> {
    let at = |i: usize, j: usize| h[i * d + j];
    // p[i] = characteristic polynomial of the leading i x i block, as coefficients in powers of x (p[i][t] of x^t)
    let mut p: Vec<Vec<C64>> = vec![vec![C64::ONE]];
    for i in 1..=d {
        let mut next = vec![C64::ZERO; i + 1];
        // (x - h_ii) p_(i-1)
        for (t, &q) in p[i - 1].iter().enumerate() {
            next[t + 1] = next[t + 1].add(q);
            next[t] = sub(next[t], at(i - 1, i - 1).mul(q));
        }
        // - sum_m h_(i-m, i) beta_i ... beta_(i-m+1) p_(i-m-1)   (1-based: beta_j = h_(j, j-1))
        let mut prod = C64::ONE;
        for m in 1..i {
            prod = prod.mul(at(i - m, i - m - 1));
            let coef = at(i - m - 1, i - 1).mul(prod);
            for (t, &q) in p[i - m - 1].iter().enumerate() {
                next[t] = sub(next[t], coef.mul(q));
            }
        }
        p.push(next);
    }
    let pd = &p[d];
    (0..=d).map(|k| pd[d - k]).collect()
}

/// the subset term without its sign: [eta^n] exp(sum_k g_k eta^k)
fn subset_term(a: &[C64], m: usize, pairs: &[usize], loops: Option<&[C64]>, n: usize) -> C64 {
    let idx: Vec<usize> = pairs.iter().flat_map(|&i| [2 * i, 2 * i + 1]).collect();
    let d = idx.len();
    // E = (A X)_Z: column c of A X is column c ^ 1 of A
    let e: Vec<C64> = idx.iter().flat_map(|&r| idx.iter().map(move |&c| a[r * m + (c ^ 1)])).collect();
    let tr = power_traces(&e, d, n);
    let mut g: Vec<C64> = (1..=n).map(|k| div_real(tr[k - 1], 2.0 * k as f64)).collect();
    if let Some(l) = loops {
        // (1/2) v^T X E^(k-1) v; (X w)_r = w_(r ^ 1) in the local index
        let v: Vec<C64> = idx.iter().map(|&r| l[r]).collect();
        let mut w = v.clone();
        for gk in g.iter_mut() {
            let mut acc = C64::ZERO;
            for (r, vr) in v.iter().enumerate() {
                acc = acc.add(vr.mul(w[r ^ 1]));
            }
            *gk = gk.add(acc.scale(0.5));
            let mut nw = vec![C64::ZERO; d];
            for (r, nr) in nw.iter_mut().enumerate() {
                let mut s = C64::ZERO;
                for (c, wc) in w.iter().enumerate() {
                    s = s.add(e[r * d + c].mul(*wc));
                }
                *nr = s;
            }
            w = nw;
        }
    }
    let mut h = vec![C64::ZERO; n + 1];
    h[0] = C64::ONE;
    for mm in 1..=n {
        let mut acc = C64::ZERO;
        for k in 1..=mm {
            acc = acc.add(g[k - 1].scale(k as f64).mul(h[mm - k]));
        }
        h[mm] = div_real(acc, mm as f64);
    }
    h[n]
}

/// hafnian (loops = None) or loop hafnian (loops = Some(weights on the diagonal)) of a symmetric m x m matrix.
pub fn hafnian(
    a: &[C64],
    m: usize,
    loops: Option<&[C64]>,
    threads: usize,
    cancel: &(dyn Fn() -> bool + Sync),
) -> Result<HafResult, Error> {
    if a.len() != m * m || loops.is_some_and(|l| l.len() != m) {
        return Err(Error::BadShape("matrix and loop weights"));
    }
    if a.iter().chain(loops.unwrap_or(&[])).any(|z| !z.is_finite()) {
        return Err(Error::NonFinite("input"));
    }
    if m % 2 == 1 {
        return match loops {
            None => Ok(HafResult { value: C64::ZERO, scale: 0.0, subsets: 0 }),
            Some(l) => {
                // pad one vertex with loop weight 1 and no edges
                let mp = m + 1;
                let mut ap = vec![C64::ZERO; mp * mp];
                for i in 0..m {
                    for j in 0..m {
                        ap[i * mp + j] = a[i * m + j];
                    }
                }
                let mut lp = l.to_vec();
                lp.push(C64::ONE);
                hafnian(&ap, mp, Some(&lp), threads, cancel)
            }
        };
    }
    let n = m / 2;
    if n == 0 {
        return Ok(HafResult { value: C64::ONE, scale: 1.0, subsets: 0 });
    }
    if n > MAX_PAIRS {
        return Err(Error::TooLarge(m));
    }
    let total = 1u64 << n;
    let blocks = total.div_ceil(BLOCK);
    let block = |b: u64| -> Result<(C64, f64), Error> {
        let (start, end) = (b * BLOCK, ((b + 1) * BLOCK).min(total));
        let mut s = CSum::new(end - start);
        let mut sc = Sum::new(end - start);
        for p in start..end {
            if (p - start).is_multiple_of(256) && cancel() {
                return Err(Error::Cancelled);
            }
            let gray = p ^ (p >> 1);
            let pairs: Vec<usize> = (0..n).filter(|&i| gray >> i & 1 == 1).collect();
            let mut t = subset_term(a, m, &pairs, loops, n);
            if (n - pairs.len()) % 2 == 1 {
                t = t.scale(-1.0);
            }
            s.add(t);
            sc.add(t.abs());
        }
        Ok((s.value(), sc.value()))
    };
    let results: Vec<Result<(C64, f64), Error>> = if threads <= 1 || blocks == 1 {
        (0..blocks).map(block).collect()
    } else {
        let per = blocks.div_ceil(threads as u64);
        std::thread::scope(|sc| {
            let block = &block;
            let hs: Vec<_> = (0..threads as u64)
                .map(|k| sc.spawn(move || (k * per..((k + 1) * per).min(blocks)).map(block).collect::<Vec<_>>()))
                .collect();
            hs.into_iter().flat_map(|h| h.join().expect("hafnian worker panicked")).collect()
        })
    };
    let mut sum = CSum::new(blocks);
    let mut scale = Sum::new(blocks);
    for r in results {
        let (z, s) = r?;
        sum.add(z);
        scale.add(s);
    }
    let value = sum.value();
    if !value.is_finite() {
        return Err(Error::NonFinite("output"));
    }
    Ok(HafResult { value, scale: scale.value(), subsets: total })
}

/// Repeated indices (FS 8.4 hafnian_repeated, loop_hafnian with Mult): index i taken mult[i] times; the
/// expanded matrix goes through the formula above. Loop weights, when given, are repeated with their index.
pub fn hafnian_repeated(
    a: &[C64],
    m: usize,
    loops: Option<&[C64]>,
    mult: &[u32],
    threads: usize,
    cancel: &(dyn Fn() -> bool + Sync),
) -> Result<HafResult, Error> {
    if mult.len() != m || a.len() != m * m {
        return Err(Error::BadShape("multiplicities"));
    }
    let idx: Vec<usize> = (0..m).flat_map(|i| std::iter::repeat_n(i, mult[i] as usize)).collect();
    let k = idx.len();
    if k / 2 > MAX_PAIRS {
        return Err(Error::TooLarge(k));
    }
    let b: Vec<C64> = idx.iter().flat_map(|&i| idx.iter().map(move |&j| a[i * m + j])).collect();
    let l: Option<Vec<C64>> = loops.map(|l| idx.iter().map(|&i| l[i]).collect());
    hafnian(&b, k, l.as_deref(), threads, cancel)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn never() -> bool {
        false
    }

    fn sym(m: usize, seed: u64) -> Vec<C64> {
        let mut s = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        let mut next = || {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((s >> 11) as f64) / (1u64 << 53) as f64 - 0.5
        };
        let mut a = vec![C64::ZERO; m * m];
        for i in 0..m {
            for j in i..m {
                let z = C64::new(next(), next());
                a[i * m + j] = z;
                a[j * m + i] = z;
            }
        }
        a
    }

    /// the definitions: sum over perfect matchings (loops: over matchings where unmatched vertices take their loop)
    fn by_matchings(a: &[C64], m: usize, loops: Option<&[C64]>) -> C64 {
        fn rec(a: &[C64], m: usize, loops: Option<&[C64]>, free: &mut Vec<usize>) -> C64 {
            if free.is_empty() {
                return C64::ONE;
            }
            let i = free.remove(0);
            let mut acc = C64::ZERO;
            if let Some(l) = loops {
                acc = acc.add(l[i].mul(rec(a, m, loops, free)));
            }
            for t in 0..free.len() {
                let j = free.remove(t);
                acc = acc.add(a[i * m + j].mul(rec(a, m, loops, free)));
                free.insert(t, j);
            }
            free.insert(0, i);
            acc
        }
        rec(a, m, loops, &mut (0..m).collect())
    }

    #[test]
    fn power_traces_match_repeated_products() {
        let d = 7;
        let e = sym(d, 3);
        let tr = power_traces(&e, d, 12);
        let mut p = e.clone();
        for k in 1..=12 {
            let t = (0..d).fold(C64::ZERO, |acc, i| acc.add(p[i * d + i]));
            assert!((t.re - tr[k - 1].re).abs() + (t.im - tr[k - 1].im).abs() < 1e-10 * t.abs().max(1.0), "k = {k}");
            let mut q = vec![C64::ZERO; d * d];
            for i in 0..d {
                for j in 0..d {
                    for l in 0..d {
                        q[i * d + j] = q[i * d + j].add(p[i * d + l].mul(e[l * d + j]));
                    }
                }
            }
            p = q;
        }
    }

    #[test]
    fn against_the_matching_recursion() {
        for m in 0..=12 {
            for seed in 0..3 {
                let a = sym(m, 10 * m as u64 + seed);
                let l: Vec<C64> = (0..m).map(|i| a[i * m + i]).collect();
                for loops in [None, Some(&l[..])] {
                    let got = hafnian(&a, m, loops, 1, &never).unwrap();
                    let want = if m % 2 == 1 && loops.is_none() { C64::ZERO } else { by_matchings(&a, m, loops) };
                    let err = (got.value.re - want.re).abs() + (got.value.im - want.im).abs();
                    assert!(err <= 1e-12 * got.scale.max(1.0), "m = {m}, loops = {}: {err:e}", loops.is_some());
                }
            }
        }
    }

    #[test]
    fn known_values() {
        // haf of the all-ones 2n x 2n matrix is (2n-1)!!; lhaf of the all-ones n x n matrix is the number of
        // involutions (telephone numbers)
        let ones = |m: usize| vec![C64::ONE; m * m];
        assert_eq!(hafnian(&ones(8), 8, None, 1, &never).unwrap().value.re.round(), 105.0);
        let l = vec![C64::ONE; 6];
        assert_eq!(hafnian(&ones(6), 6, Some(&l), 1, &never).unwrap().value.re.round(), 76.0);
        assert_eq!(hafnian(&[], 0, None, 1, &never).unwrap().value, C64::ONE);
    }

    #[test]
    fn threads_change_no_bit() {
        let a = sym(28, 5); // 2^14 subsets = 4 blocks
        let one = hafnian(&a, 28, None, 1, &never).unwrap();
        let four = hafnian(&a, 28, None, 4, &never).unwrap();
        assert_eq!(one.value.re.to_bits(), four.value.re.to_bits());
        assert_eq!(one.value.im.to_bits(), four.value.im.to_bits());
    }

    #[test]
    fn repeated_is_the_expanded_matrix() {
        let a = sym(3, 9);
        let got = hafnian_repeated(&a, 3, None, &[2, 1, 1], 1, &never).unwrap().value;
        let idx = [0usize, 0, 1, 2];
        let ar = &a;
        let b: Vec<C64> = idx.iter().flat_map(|&i| idx.iter().map(move |&j| ar[i * 3 + j])).collect();
        let want = by_matchings(&b, 4, None);
        assert!((got.re - want.re).abs() + (got.im - want.im).abs() < 1e-13);
    }
}
