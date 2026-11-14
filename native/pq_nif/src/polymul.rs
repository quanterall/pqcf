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

/// bases of at most this many terms are multiplied by merging (mul_merge)
const MERGE_MAX: usize = 64;

/// acc * b, acc sorted by packed key (b in any order): b's shifted copies of acc, each sorted, merged by key with
/// a heap of |b| heads, equal keys combined as they meet; sorted, zero coefficients dropped. The layout must hold
/// the sums, as for mul_packed.
fn mul_merge(acc: &Packed, b: &Packed, real: bool) -> Result<Packed, Refused> {
    use std::cmp::Reverse;
    use std::collections::BinaryHeap;
    if acc.is_empty() || b.is_empty() {
        return Ok(vec![]);
    }
    let m = |x: i128, y: i128| x.checked_mul(y).ok_or(Refused::Overflow);
    let mut heap: BinaryHeap<Reverse<(u64, usize, usize)>> = b.iter().enumerate().map(|(j, &(kb, _))| Reverse((acc[0].0 + kb, j, 0))).collect();
    let mut out: Packed = Vec::with_capacity(acc.len() + acc.len() / 2 + b.len());
    while let Some(Reverse((k, j, i))) = heap.pop() {
        let ((_, (ar, ai)), (kb, (br, bi))) = (acc[i], b[j]);
        let (re, im) = if real { (m(ar, br)?, 0) } else { (m(ar, br)?.checked_sub(m(ai, bi)?).ok_or(Refused::Overflow)?, m(ar, bi)?.checked_add(m(ai, br)?).ok_or(Refused::Overflow)?) };
        match out.last_mut() {
            Some((lk, (lr, li))) if *lk == k => {
                *lr = lr.checked_add(re).ok_or(Refused::Overflow)?;
                *li = li.checked_add(im).ok_or(Refused::Overflow)?;
            }
            _ => out.push((k, (re, im))),
        }
        if i + 1 < acc.len() {
            heap.push(Reverse((acc[i + 1].0 + kb, j, i + 1)));
        }
    }
    out.retain(|&(_, c)| c != (0, 0));
    Ok(out)
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
        // a base of few terms (a sum of symbols, say): its shifted copies of the sorted accumulator merged, no
        // hashing; a large base: the hashed product
        acc = if pa.len() <= MERGE_MAX { mul_merge(&acc, &pa, real)? } else { mul_packed(&acc, &pa, real)? };
    }
    let lo: Vec<i64> = lo_a.iter().map(|&x| x.checked_mul(n as i64).ok_or(Refused::Packing)).collect::<Result<_, _>>()?;
    Ok(unpack(&l, acc, &lo))
}

/// The permanent of a k x k matrix of polynomials (row-major), by the subset expansion pq_symbolic:glynn/3 does in
/// Erlang: row by row, for every set of columns used so far, the sum over the assignments of the rows to them. The
/// entries are packed once in the layout of a product of k of them; each partial sum is a hashed accumulator.
/// Sorted, zero coefficients dropped; k = 0 gives 1.
pub fn permanent(nvars: usize, k: usize, e: &[Poly]) -> Result<Poly, Refused> {
    if e.len() != k * k {
        return Err(Refused::Packing);
    }
    if k == 0 {
        return Ok(Poly { exps: vec![0; nvars], coef: vec![1, 0] });
    }
    if k > 24 {
        return Err(Refused::Packing);
    }
    // per variable: the lowest exponent over every entry, and what a product of one entry from each row can add to
    // k times it -- the sum over the rows of the row's highest exponent above it (a matrix of distinct symbols: one
    // bit a symbol, so 8 x 8 fits a 64-bit key)
    let lo: Vec<i64> = (0..nvars).map(|j| e.iter().filter(|p| !p.is_empty()).map(|p| range(p, nvars, j).0).min().unwrap_or(0)).collect();
    let spans: Vec<i128> = (0..nvars)
        .map(|j| {
            (0..k)
                .map(|i| e[i * k..(i + 1) * k].iter().filter(|p| !p.is_empty()).map(|p| (range(p, nvars, j).1 - lo[j]) as i128).max().unwrap_or(0))
                .sum()
        })
        .collect();
    let l = layout(&spans)?;
    let packed: Vec<Packed> = e.iter().map(|p| if p.is_empty() { vec![] } else { pack(&l, p, &lo) }).collect();
    let real = e.iter().all(is_real);
    let mut layer: HashMap<u32, Packed> = HashMap::from([(0u32, vec![(0u64, (1i128, 0i128))])]);
    for i in 0..k {
        let mut next: HashMap<u32, KeyMap<(i128, i128)>> = HashMap::new();
        for (&used, p) in &layer {
            for j in 0..k {
                let bit = 1u32 << j;
                let entry = &packed[i * k + j];
                if used & bit != 0 || entry.is_empty() {
                    continue;
                }
                let acc = next.entry(used | bit).or_insert_with(|| KeyMap::with_capacity_and_hasher(p.len(), Default::default()));
                for (key, (re, im)) in mul_packed(p, entry, real)? {
                    let c = acc.entry(key).or_insert((0, 0));
                    c.0 = c.0.checked_add(re).ok_or(Refused::Overflow)?;
                    c.1 = c.1.checked_add(im).ok_or(Refused::Overflow)?;
                }
            }
        }
        layer = next.into_iter().map(|(m, acc)| (m, acc.into_iter().filter(|&(_, c)| c != (0, 0)).collect())).collect();
    }
    let full = layer.remove(&((1u32 << k) - 1)).unwrap_or_default();
    let lo_k: Vec<i64> = lo.iter().map(|&x| x.checked_mul(k as i64).ok_or(Refused::Packing)).collect::<Result<_, _>>()?;
    Ok(unpack(&l, full, &lo_k))
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
    fn permanent_is_the_sum_over_permutations() {
        fn perms(k: usize) -> Vec<Vec<usize>> {
            if k == 0 {
                return vec![vec![]];
            }
            perms(k - 1).into_iter().flat_map(|p| (0..k).map(move |i| { let mut q = p.clone(); q.insert(i, k - 1); q })).collect()
        }
        let mut s = 0x9E37_79B9_7F4A_7C15u64;
        for case in 0..24 {
            let k = 1 + case % 4;
            let nvars = 1 + case % 3;
            let complex = case % 2 == 0;
            let e: Vec<Poly> = (0..k * k).map(|i| if (i + case) % 7 == 3 { Poly { exps: vec![], coef: vec![] } } else { random(&mut s, nvars, 1 + i % 3, 0, 2, complex) }).collect();
            let got = permanent(nvars, k, &e).unwrap();
            // the sum over permutations of the products, each product by repeated schoolbook multiplication
            let mut want: Vec<(Vec<i64>, (i128, i128))> = vec![];
            for p in perms(k) {
                let mut prod = vec![(vec![0i64; nvars], (1i128, 0i128))];
                for (row, &col) in p.iter().enumerate() {
                    let w = Poly { exps: prod.iter().flat_map(|t| t.0.clone()).collect(), coef: prod.iter().flat_map(|t| [t.1 .0, t.1 .1]).collect() };
                    prod = naive(nvars, &w, &e[row * k + col]);
                }
                let sum = Poly { exps: want.iter().chain(&prod).flat_map(|t| t.0.clone()).collect(), coef: want.iter().chain(&prod).flat_map(|t| [t.1 .0, t.1 .1]).collect() };
                let one = Poly { exps: vec![0; nvars], coef: vec![1, 0] };
                want = naive(nvars, &sum, &one);
            }
            assert_eq!(as_sorted(nvars, &got), want, "case {case}");
        }
        assert_eq!(as_sorted(2, &permanent(2, 0, &[]).unwrap()), vec![(vec![0, 0], (1, 0))]);
    }

    #[test]
    fn merged_and_hashed_products_agree() {
        // mul_merge (small bases) against mul_packed (hashing), on sorted accumulators, real and complex
        let mut s = 0x2545_F491_4F6C_DD1D;
        for case in 0..40 {
            let nvars = 1 + case % 3;
            let complex = case % 2 == 1;
            let acc = random(&mut s, nvars, 5 + case * 7, 0, 6, complex);
            let b = random(&mut s, nvars, 1 + case % 9, 0, 4, complex);
            let spans: Vec<i128> = (0..nvars).map(|j| (range(&acc, nvars, j).1 + range(&b, nvars, j).1) as i128).collect();
            let l = layout(&spans).unwrap();
            let zero = vec![0i64; nvars];
            let mut pa = pack(&l, &acc, &zero);
            pa.sort_unstable_by_key(|t| t.0);
            let pb = pack(&l, &b, &zero);
            let mut hashed = mul_packed(&pa, &pb, !complex).unwrap();
            hashed.sort_unstable_by_key(|t| t.0);
            assert_eq!(mul_merge(&pa, &pb, !complex).unwrap(), hashed, "case {case}");
        }
        // a base of more than MERGE_MAX terms takes the hashed path; both paths give the schoolbook power
        let a = random(&mut s, 2, MERGE_MAX + 6, 0, 9, false);
        let sq = pow(2, &a, 2).unwrap();
        assert_eq!(as_sorted(2, &sq), naive(2, &a, &a));
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
