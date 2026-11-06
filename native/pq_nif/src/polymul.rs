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
    pub fn len(&self) -> usize {
        self.coef.len() / 2
    }

    pub fn is_empty(&self) -> bool {
        self.coef.is_empty()
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

/// Terms in packed form: (exponents packed, (real, imaginary)).
type Packed = Vec<(u64, (i128, i128))>;

/// Per variable, where its exponent sits in a packed u64: the bits of its range, from `shift`.
struct Layout {
    shift: Vec<u32>,
    bits: Vec<u32>,
}

/// the layout for exponent ranges of the given spans (largest minus smallest), if they fit in 64 bits
fn layout(spans: &[i128]) -> Result<Layout, Refused> {
    let mut l = Layout { shift: Vec::with_capacity(spans.len()), bits: Vec::with_capacity(spans.len()) };
    let mut used = 0u32;
    for &span in spans {
        if span >= 1i128 << 63 {
            return Err(Refused::Packing);
        }
        let w = 64 - (span as u64).leading_zeros();
        l.shift.push(used);
        l.bits.push(w);
        used += w;
        if used > 64 {
            return Err(Refused::Packing);
        }
    }
    Ok(l)
}

/// the smallest and largest exponent of variable j over the terms
fn range(p: &Poly, nvars: usize, j: usize) -> (i64, i64) {
    (0..p.len()).map(|t| p.exps[t * nvars + j]).fold((i64::MAX, i64::MIN), |(lo, hi), e| (lo.min(e), hi.max(e)))
}

/// p's terms packed, every exponent offset by -lo
fn pack(l: &Layout, p: &Poly, lo: &[i64]) -> Packed {
    let nvars = lo.len();
    (0..p.len())
        .map(|t| {
            let k = (0..nvars).fold(0u64, |k, j| k | (((p.exps[t * nvars + j] - lo[j]) as u64) << (l.shift[j] % 64)));
            (k, (p.coef[2 * t], p.coef[2 * t + 1]))
        })
        .collect()
}

/// packed terms back to a polynomial, sorted by packed exponents, every exponent offset by +lo
fn unpack(l: &Layout, mut terms: Packed, lo: &[i64]) -> Poly {
    terms.sort_unstable_by_key(|t| t.0);
    let nvars = lo.len();
    let mut out = Poly { exps: Vec::with_capacity(terms.len() * nvars), coef: Vec::with_capacity(terms.len() * 2) };
    for (k, (re, im)) in terms {
        for ((&bits, &shift), &lo_j) in l.bits.iter().zip(&l.shift).zip(lo) {
            let mask = if bits == 64 { u64::MAX } else { (1u64 << bits) - 1 };
            let off = if bits == 0 { 0 } else { (k >> shift) & mask };
            out.exps.push(off as i64 + lo_j);
        }
        out.coef.push(re);
        out.coef.push(im);
    }
    out
}

/// the product of packed terms (the layout must hold the sums), zero coefficients dropped, unsorted
fn mul_packed(a: &Packed, b: &Packed, real: bool) -> Result<Packed, Refused> {
    let mut acc: KeyMap<(i128, i128)> = KeyMap::with_capacity_and_hasher((a.len() * b.len()).min(1 << 22), Default::default());
    for &(ki, (ar, ai)) in a {
        for &(kj, (br, bi)) in b {
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
    Ok(acc.into_iter().filter(|&(_, c)| c != (0, 0)).collect())
}

fn is_real(p: &Poly) -> bool {
    p.coef.iter().skip(1).step_by(2).all(|&x| x == 0)
}

/// a * b, its terms sorted by packed exponents, zero coefficients dropped
pub fn mul(nvars: usize, a: &Poly, b: &Poly) -> Result<Poly, Refused> {
    if a.is_empty() || b.is_empty() {
        return Ok(Poly { exps: vec![], coef: vec![] });
    }
    // per variable: the smallest exponent of each factor and the bits the product's range needs
    let (ra, rb): (Vec<_>, Vec<_>) = (0..nvars).map(|j| (range(a, nvars, j), range(b, nvars, j))).unzip();
    let spans: Vec<i128> = ra.iter().zip(&rb).map(|(&(la, ha), &(lb, hb))| (ha as i128 - la as i128) + (hb as i128 - lb as i128)).collect();
    let l = layout(&spans)?;
    let lo_a: Vec<i64> = ra.iter().map(|r| r.0).collect();
    let lo_b: Vec<i64> = rb.iter().map(|r| r.0).collect();
    let prod = mul_packed(&pack(&l, a, &lo_a), &pack(&l, b, &lo_b), is_real(a) && is_real(b))?;
    let lo: Vec<i64> = lo_a.iter().zip(&lo_b).map(|(x, y)| x + y).collect();
    Ok(unpack(&l, prod, &lo))
}

/// a + b, sorted, zero coefficients dropped
pub fn add(nvars: usize, a: &Poly, b: &Poly) -> Result<Poly, Refused> {
    if a.is_empty() && b.is_empty() {
        return Ok(Poly { exps: vec![], coef: vec![] });
    }
    let both = |j: usize| {
        let (la, ha) = if !a.is_empty() { range(a, nvars, j) } else { range(b, nvars, j) };
        let (lb, hb) = if !b.is_empty() { range(b, nvars, j) } else { (la, ha) };
        (la.min(lb), ha.max(hb))
    };
    let r: Vec<(i64, i64)> = (0..nvars).map(both).collect();
    let l = layout(&r.iter().map(|&(lo, hi)| hi as i128 - lo as i128).collect::<Vec<_>>())?;
    let lo: Vec<i64> = r.iter().map(|x| x.0).collect();
    let mut acc: KeyMap<(i128, i128)> = KeyMap::with_capacity_and_hasher(a.len() + b.len(), Default::default());
    for (k, (re, im)) in pack(&l, a, &lo).into_iter().chain(pack(&l, b, &lo)) {
        let e = acc.entry(k).or_insert((0, 0));
        e.0 = e.0.checked_add(re).ok_or(Refused::Overflow)?;
        e.1 = e.1.checked_add(im).ok_or(Refused::Overflow)?;
    }
    Ok(unpack(&l, acc.into_iter().filter(|&(_, c)| c != (0, 0)).collect(), &lo))
}

/// a^n by repeated multiplication with a, in the layout of a^n's exponent ranges (n times a's), sorted, zero
/// coefficients dropped; a^0 = 1
pub fn pow(nvars: usize, a: &Poly, n: u32) -> Result<Poly, Refused> {
    if n == 0 {
        return Ok(Poly { exps: vec![0; nvars], coef: vec![1, 0] });
    }
    if a.is_empty() {
        return Ok(Poly { exps: vec![], coef: vec![] });
    }
    let ra: Vec<(i64, i64)> = (0..nvars).map(|j| range(a, nvars, j)).collect();
    let spans: Vec<i128> = ra.iter().map(|&(lo, hi)| (hi as i128 - lo as i128) * n as i128).collect();
    let l = layout(&spans)?;
    let lo_a: Vec<i64> = ra.iter().map(|r| r.0).collect();
    let pa = pack(&l, a, &lo_a);
    let real = is_real(a);
    let mut acc: Packed = vec![(0, (1, 0))];
    for _ in 0..n {
        acc = mul_packed(&acc, &pa, real)?;
    }
    let lo: Vec<i64> = lo_a.iter().map(|&x| x.checked_mul(n as i64).ok_or(Refused::Packing)).collect::<Result<_, _>>()?;
    Ok(unpack(&l, acc, &lo))
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
    fn powers_match_repeated_schoolbook_products() {
        let mut s = 0x9E37_79B9_7F4A_7C15;
        for case in 0..60 {
            let nvars = 1 + case % 4;
            let (lo, hi) = if case % 3 == 0 { (-2, 2) } else { (0, 3) };
            let a = random(&mut s, nvars, 1 + case % 6, lo, hi, case % 2 == 0);
            let n = (case % 7) as u32;
            let got = pow(nvars, &a, n).unwrap();
            let mut want = vec![(vec![0i64; nvars], (1i128, 0i128))];
            for _ in 0..n {
                let w = Poly {
                    exps: want.iter().flat_map(|t| t.0.clone()).collect(),
                    coef: want.iter().flat_map(|t| [t.1 .0, t.1 .1]).collect(),
                };
                want = naive(nvars, &w, &a);
            }
            assert_eq!(as_sorted(nvars, &got), want, "case {case}");
        }
        let big = Poly { exps: vec![0, 1], coef: vec![1i128 << 40, 0, 1, 0] };
        assert_eq!(pow(1, &big, 4).err(), Some(Refused::Overflow));
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
