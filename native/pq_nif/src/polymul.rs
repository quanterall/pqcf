//! The product of two plain polynomials of the symbolic engine (pq_sym_poly): monomials of symbols with integer
//! exponents, exact Gaussian-integer coefficients (pq_sym_poly scales away the denominators first). The exponents
//! of a monomial are packed into one u64, each variable in the bits its range in the product needs, so a product of
//! monomials is one addition; the coefficients are i128 with checked arithmetic. The result is exact or refused --
//! an overflow or exponents that don't pack -- and pq_sym_poly then multiplies in Erlang: the same polynomial
//! either way (D1).

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

#[derive(Debug, PartialEq, Eq)]
pub enum Refused {
    /// a coefficient of the product, or a partial sum of one, outside i128
    Overflow,
    /// the exponent ranges of the product need more than 64 bits
    Packing,
}

/// A polynomial in `nvars` variables: term t has exponents `exps[t * nvars..(t + 1) * nvars]` and coefficient
/// `coef[2 t] + i coef[2 t + 1]`.
pub struct Poly {
    pub exps: Vec<i64>,
    pub coef: Vec<i128>,
}

impl Poly {
    fn len(&self) -> usize {
        self.coef.len() / 2
    }
}

/// Fibonacci hashing of a packed monomial: a multiply and a rotate, deterministic (no random state)
#[derive(Default)]
struct KeyHasher(u64);

impl Hasher for KeyHasher {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0.rotate_left(8) ^ b as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
        }
    }
    fn write_u64(&mut self, k: u64) {
        self.0 = (k ^ (k >> 29)).wrapping_mul(0x9E37_79B9_7F4A_7C15).rotate_left(26);
    }
}

type KeyMap<V> = HashMap<u64, V, BuildHasherDefault<KeyHasher>>;

/// a * b, its terms sorted by packed exponents, zero coefficients dropped
pub fn mul(nvars: usize, a: &Poly, b: &Poly) -> Result<Poly, Refused> {
    let (na, nb) = (a.len(), b.len());
    if na == 0 || nb == 0 {
        return Ok(Poly { exps: vec![], coef: vec![] });
    }
    // per variable: the smallest exponent of each factor and the bits the product's range needs
    let range = |p: &Poly, n: usize, j: usize| {
        (0..n).map(|t| p.exps[t * nvars + j]).fold((i64::MAX, i64::MIN), |(lo, hi), e| (lo.min(e), hi.max(e)))
    };
    let mut lo_a = vec![0i64; nvars];
    let mut lo_b = vec![0i64; nvars];
    let mut shift = vec![0u32; nvars];
    let mut bits = vec![0u32; nvars];
    let mut used = 0u32;
    for j in 0..nvars {
        let (la, ha) = range(a, na, j);
        let (lb, hb) = range(b, nb, j);
        let span = (ha as i128 - la as i128) + (hb as i128 - lb as i128);
        if span >= 1i128 << 63 {
            return Err(Refused::Packing);
        }
        let w = 64 - (span as u64).leading_zeros();
        lo_a[j] = la;
        lo_b[j] = lb;
        shift[j] = used;
        bits[j] = w;
        used += w;
        if used > 64 {
            return Err(Refused::Packing);
        }
    }
    let pack = |p: &Poly, t: usize, lo: &[i64]| -> u64 {
        (0..nvars).fold(0u64, |k, j| k | (((p.exps[t * nvars + j] - lo[j]) as u64) << (shift[j] % 64)))
    };
    let ka: Vec<u64> = (0..na).map(|t| pack(a, t, &lo_a)).collect();
    let kb: Vec<u64> = (0..nb).map(|t| pack(b, t, &lo_b)).collect();
    let real = a.coef.iter().skip(1).step_by(2).all(|&x| x == 0) && b.coef.iter().skip(1).step_by(2).all(|&x| x == 0);

    let mut acc: KeyMap<(i128, i128)> = KeyMap::with_capacity_and_hasher((na * nb).min(1 << 22), Default::default());
    for (&ki, ca) in ka.iter().zip(a.coef.chunks_exact(2)) {
        let (ar, ai) = (ca[0], ca[1]);
        for (&kj, cb) in kb.iter().zip(b.coef.chunks_exact(2)) {
            let (br, bi) = (cb[0], cb[1]);
            let e = acc.entry(ki + kj).or_insert((0, 0));
            if real {
                let p = ar.checked_mul(br).ok_or(Refused::Overflow)?;
                e.0 = e.0.checked_add(p).ok_or(Refused::Overflow)?;
            } else {
                let m = |x: i128, y: i128| x.checked_mul(y).ok_or(Refused::Overflow);
                let re = m(ar, br)?.checked_sub(m(ai, bi)?).ok_or(Refused::Overflow)?;
                let im = m(ar, bi)?.checked_add(m(ai, br)?).ok_or(Refused::Overflow)?;
                e.0 = e.0.checked_add(re).ok_or(Refused::Overflow)?;
                e.1 = e.1.checked_add(im).ok_or(Refused::Overflow)?;
            }
        }
    }
    let mut terms: Vec<(u64, (i128, i128))> = acc.into_iter().filter(|&(_, c)| c != (0, 0)).collect();
    terms.sort_unstable_by_key(|t| t.0);
    let mut out = Poly { exps: Vec::with_capacity(terms.len() * nvars), coef: Vec::with_capacity(terms.len() * 2) };
    for (k, (re, im)) in terms {
        for j in 0..nvars {
            let mask = if bits[j] == 64 { u64::MAX } else { (1u64 << bits[j]) - 1 };
            let off = if bits[j] == 0 { 0 } else { (k >> shift[j]) & mask };
            out.exps.push(off as i64 + lo_a[j] + lo_b[j]);
        }
        out.coef.push(re);
        out.coef.push(im);
    }
    Ok(out)
}

/// the greatest common divisor, gcd(n, 0) = n
pub fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

#[cfg(test)]
mod tests {
    use super::*;

    /// the schoolbook product over a BTreeMap of exponent vectors, as the reference
    fn naive(nvars: usize, a: &Poly, b: &Poly) -> Vec<(Vec<i64>, (i128, i128))> {
        let mut m = std::collections::BTreeMap::new();
        for i in 0..a.len() {
            for j in 0..b.len() {
                let e: Vec<i64> = (0..nvars).map(|v| a.exps[i * nvars + v] + b.exps[j * nvars + v]).collect();
                let (ar, ai, br, bi) = (a.coef[2 * i], a.coef[2 * i + 1], b.coef[2 * j], b.coef[2 * j + 1]);
                let c = m.entry(e).or_insert((0i128, 0i128));
                c.0 += ar * br - ai * bi;
                c.1 += ar * bi + ai * br;
            }
        }
        m.into_iter().filter(|(_, c)| *c != (0, 0)).collect()
    }

    fn as_sorted(nvars: usize, p: &Poly) -> Vec<(Vec<i64>, (i128, i128))> {
        let mut v: Vec<_> = (0..p.len())
            .map(|t| (p.exps[t * nvars..(t + 1) * nvars].to_vec(), (p.coef[2 * t], p.coef[2 * t + 1])))
            .collect();
        v.sort();
        v
    }

    fn random(seed: &mut u64, nvars: usize, n: usize, lo: i64, hi: i64, complex: bool) -> Poly {
        let mut next = || {
            *seed ^= *seed << 13;
            *seed ^= *seed >> 7;
            *seed ^= *seed << 17;
            *seed
        };
        let mut p = Poly { exps: vec![], coef: vec![] };
        for _ in 0..n {
            for _ in 0..nvars {
                p.exps.push(lo + (next() % (hi - lo + 1) as u64) as i64);
            }
            p.coef.push((next() % 21) as i128 - 10);
            p.coef.push(if complex { (next() % 7) as i128 - 3 } else { 0 });
        }
        p
    }

    #[test]
    fn matches_the_schoolbook_product() {
        let mut s = 0x2545_F491_4F6C_DD1D;
        for case in 0..200 {
            let nvars = case % 6;
            let (lo, hi) = if case % 3 == 0 { (-4, 4) } else { (0, 9) };
            let a = random(&mut s, nvars, 1 + case % 17, lo, hi, case % 2 == 0);
            let b = random(&mut s, nvars, 1 + case % 13, lo, hi, case % 4 == 1);
            let got = mul(nvars, &a, &b).unwrap();
            assert_eq!(as_sorted(nvars, &got), naive(nvars, &a, &b), "case {case}");
        }
    }

    #[test]
    fn refuses_overflow_and_wide_exponents() {
        let big = Poly { exps: vec![0], coef: vec![1i128 << 100, 0] };
        assert_eq!(mul(1, &big, &big).err(), Some(Refused::Overflow));
        let wide = Poly { exps: vec![0, 0, 1 << 40, 1 << 40], coef: vec![1, 0, 1, 0] };
        assert_eq!(mul(2, &wide, &wide).err(), Some(Refused::Packing));
        let cancel = Poly { exps: vec![0, 1], coef: vec![1, 0, 1, 0] };
        let minus = Poly { exps: vec![0, 1], coef: vec![1, 0, -1, 0] };
        let sq = mul(1, &cancel, &minus).unwrap();
        assert_eq!(as_sorted(1, &sq), vec![(vec![0], (1, 0)), (vec![2], (-1, 0))]);
    }
}
