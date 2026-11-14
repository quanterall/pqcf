//! The greatest common divisor of two multivariate polynomials modulo a prime p < 2^63 (decision record SY2), by
//! Brown's dense modular algorithm (W. S. Brown, J. ACM 18 (1971) 478; Geddes, Czapor, Labahn, "Algorithms for
//! Computer Algebra", 7.4): the last variable is evaluated at points 1, 2, 3, ... (the same points on every node,
//! D1), the gcds of the images are found recursively, scaled to the gcd of the leading coefficients, and
//! interpolated back (Newton). A candidate is accepted only when it divides both inputs exactly mod p. pq_sym_gcd
//! brings the results for several primes back to the integers.
//!
//! A polynomial is a list of terms (exponents, one per variable, and a coefficient in [1, p)), sorted by the
//! exponents in descending lexicographic order, the first variable most significant.

use crate::ffslp::{invmod, Mont};
use std::collections::BTreeMap;

pub type Term = (Vec<u32>, u64);

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MPoly {
    pub terms: Vec<Term>,
}

/// Univariate polynomials mod p, dense, coefficient i of x^i, no trailing zero (the zero polynomial is empty).
pub(crate) type UPoly = Vec<u64>;

#[inline]
pub(crate) fn addm(a: u64, b: u64, p: &Fp) -> u64 {
    let s = a + b;
    if s >= p.p { s - p.p } else { s }
}

#[inline]
pub(crate) fn subm(a: u64, b: u64, p: &Fp) -> u64 {
    if a >= b { a - b } else { a + p.p - b }
}

pub(crate) fn inv(a: u64, p: &Fp) -> u64 {
    p.m.to(invmod(p.m.from(a), p.p).expect("inverse of 0 mod p"))
}

/// The field mod p, its elements in Montgomery form (x stands for x 2^64 mod p): a product without a 128-bit
/// division, a library call on aarch64. Sums and differences are the usual ones.
pub struct Fp {
    pub p: u64,
    pub(crate) m: Mont,
    pub(crate) one: u64,
}

impl Fp {
    pub fn new(p: u64) -> Fp {
        let m = Mont::new(p);
        Fp { p, m, one: m.to(1) }
    }

    #[inline]
    pub(crate) fn mul(&self, a: u64, b: u64) -> u64 {
        self.m.mul(a, b)
    }
}

pub(crate) fn u_trim(mut a: UPoly) -> UPoly {
    while a.last() == Some(&0) {
        a.pop();
    }
    a
}

pub(crate) fn u_eval(a: &[u64], x: u64, p: &Fp) -> u64 {
    a.iter().rev().fold(0, |acc, &c| addm(p.mul(acc, x), c, p))
}

pub(crate) fn u_scale(a: &[u64], k: u64, p: &Fp) -> UPoly {
    if k == 0 {
        return vec![];
    }
    a.iter().map(|&c| p.mul(c, k)).collect()
}

pub(crate) fn u_monic(a: &[u64], p: &Fp) -> UPoly {
    match a.last() {
        None => vec![],
        Some(&l) => u_scale(a, inv(l, p), p),
    }
}

pub(crate) fn u_mul(a: &[u64], b: &[u64], p: &Fp) -> UPoly {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut r = vec![0u64; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        if x == 0 {
            continue;
        }
        for (j, &y) in b.iter().enumerate() {
            r[i + j] = addm(r[i + j], p.mul(x, y), p);
        }
    }
    u_trim(r)
}

pub(crate) fn u_add(a: &[u64], b: &[u64], p: &Fp) -> UPoly {
    let mut r: UPoly = (0..a.len().max(b.len()))
        .map(|i| addm(*a.get(i).unwrap_or(&0), *b.get(i).unwrap_or(&0), p))
        .collect();
    r = u_trim(r);
    r
}

/// (quotient, remainder) of a by b (b nonzero)
pub(crate) fn u_divrem(a: &[u64], b: &[u64], p: &Fp) -> (UPoly, UPoly) {
    let mut r = a.to_vec();
    if r.len() < b.len() {
        return (vec![], r);
    }
    let lb = inv(*b.last().unwrap(), p);
    let mut q = vec![0u64; r.len() - b.len() + 1];
    for k in (0..q.len()).rev() {
        let c = p.mul(r[k + b.len() - 1], lb);
        q[k] = c;
        if c != 0 {
            for (j, &y) in b.iter().enumerate() {
                r[k + j] = subm(r[k + j], p.mul(c, y), p);
            }
        }
    }
    r.truncate(b.len() - 1);
    (u_trim(q), u_trim(r))
}

/// the monic gcd (0 when both are 0)
pub(crate) fn u_gcd(a: &[u64], b: &[u64], p: &Fp) -> UPoly {
    let (mut x, mut y) = (a.to_vec(), b.to_vec());
    while !y.is_empty() {
        let r = u_divrem(&x, &y, p).1;
        x = y;
        y = r;
    }
    u_monic(&x, p)
}

impl MPoly {
    fn is_const(&self) -> bool {
        self.terms.len() == 1 && self.terms[0].0.iter().all(|&e| e == 0)
    }

    fn monic(mut self, p: &Fp) -> MPoly {
        if let Some(&(_, l)) = self.terms.first() {
            let k = inv(l, p);
            for t in &mut self.terms {
                t.1 = p.mul(t.1, k);
            }
        }
        self
    }

    /// the terms sorted (descending) with equal exponents merged and zero coefficients removed
    pub fn normalise(mut terms: Vec<Term>, p: &Fp) -> MPoly {
        terms.sort_unstable_by(|a, b| b.0.cmp(&a.0));
        let mut out: Vec<Term> = Vec::with_capacity(terms.len());
        for (e, c) in terms {
            let c = c % p.p;
            match out.last_mut() {
                Some(last) if last.0 == e => last.1 = addm(last.1, c, p),
                _ => out.push((e, c)),
            }
        }
        out.retain(|t| t.1 != 0);
        MPoly { terms: out }
    }
}

/// a polynomial as one in the variables before v, with coefficients univariate in variable v (exponent v zeroed in
/// the keys); the keys sort as the terms do
type Rec = BTreeMap<Vec<u32>, UPoly>;

fn to_rec(a: &MPoly, v: usize) -> Rec {
    let mut r: Rec = BTreeMap::new();
    for (e, c) in &a.terms {
        let mut k = e.clone();
        let d = k[v] as usize;
        k[v] = 0;
        let u = r.entry(k).or_default();
        if u.len() <= d {
            u.resize(d + 1, 0);
        }
        u[d] = *c;
    }
    r
}

fn from_rec(r: &Rec, v: usize) -> MPoly {
    let mut terms = Vec::new();
    for (k, u) in r.iter().rev() {
        for (d, &c) in u.iter().enumerate().rev() {
            if c != 0 {
                let mut e = k.clone();
                e[v] = d as u32;
                terms.push((e, c));
            }
        }
    }
    // keys descending, then the degree in v descending: the descending order, since exponent v is 0 in the keys
    // and only the variables before v are nonzero there
    MPoly { terms }
}

/// the image at x_v = x (variable v then 0), terms in order
fn eval_rec(r: &Rec, x: u64, p: &Fp) -> MPoly {
    MPoly { terms: r.iter().rev().filter_map(|(k, u)| {
        let c = u_eval(u, x, p);
        (c != 0).then(|| (k.clone(), c))
    }).collect() }
}

/// does g divide a exactly mod p (g nonzero)?
fn divides(a: &MPoly, g: &MPoly, p: &Fp) -> bool {
    let (lm, lc) = &g.terms[0];
    let li = inv(*lc, p);
    let mut r: BTreeMap<Vec<u32>, u64> = a.terms.iter().cloned().collect();
    while let Some((e, c)) = r.iter().next_back().map(|(e, &c)| (e.clone(), c)) {
        if e.iter().zip(lm).any(|(x, y)| x < y) {
            return false;
        }
        let shift: Vec<u32> = e.iter().zip(lm).map(|(x, y)| x - y).collect();
        let q = p.mul(c, li);
        for (ge, gc) in &g.terms {
            let k: Vec<u32> = ge.iter().zip(&shift).map(|(x, y)| x + y).collect();
            let t = p.mul(q, *gc);
            let ent = r.entry(k).or_insert(0);
            *ent = subm(*ent, t, p);
            if *ent == 0 {
                let k: Vec<u32> = ge.iter().zip(&shift).map(|(x, y)| x + y).collect();
                r.remove(&k);
            }
        }
    }
    true
}

fn constant(n: usize, p: &Fp) -> MPoly {
    MPoly { terms: vec![(vec![0; n], p.one)] }
}

/// The monic gcd of a and b (both nonzero) mod p, in variables 0..k (the others absent). The variables present are
/// moved to the front first; then Zippel's sparse method, and Brown's dense one if that fails (an unlucky random
/// point). None only if both fail.
pub fn gcd(a: &MPoly, b: &MPoly, k: usize, p: &Fp) -> Option<MPoly> {
    let n = a.terms[0].0.len();
    if a.is_const() || b.is_const() {
        return Some(constant(n, p));
    }
    let used: Vec<usize> = (0..k).filter(|&i| a.terms.iter().chain(&b.terms).any(|t| t.0[i] > 0)).collect();
    if used.len() < k || used.iter().enumerate().any(|(j, &i)| i != j) {
        let m = used.len();
        let fwd = |x: &MPoly| MPoly::normalise(x.terms.iter().map(|(e, c)| {
            let mut f = vec![0; n];
            for (j, &i) in used.iter().enumerate() {
                f[j] = e[i];
            }
            (f, *c)
        }).collect(), p);
        let g = gcd(&fwd(a), &fwd(b), m, p)?;
        return Some(MPoly::normalise(g.terms.iter().map(|(f, c)| {
            let mut e = vec![0; n];
            for (j, &i) in used.iter().enumerate() {
                e[i] = f[j];
            }
            (e, *c)
        }).collect(), p));
    }
    // the monomial contents: their gcd times the gcd of what is left
    let low = |x: &MPoly| -> Vec<u32> { (0..n).map(|i| x.terms.iter().map(|t| t.0[i]).min().unwrap_or(0)).collect() };
    let (ma, mb) = (low(a), low(b));
    if ma.iter().chain(&mb).any(|&e| e > 0) {
        let shift = |x: &MPoly, m: &[u32]| MPoly { terms: x.terms.iter().map(|(e, c)| (e.iter().zip(m).map(|(x, y)| x - y).collect(), *c)).collect() };
        let g = gcd(&shift(a, &ma), &shift(b, &mb), k, p)?;
        let m: Vec<u32> = ma.iter().zip(&mb).map(|(x, y)| *x.min(y)).collect();
        return Some(MPoly { terms: g.terms.into_iter().map(|(e, c)| (e.iter().zip(&m).map(|(x, y)| x + y).collect(), c)).collect() });
    }
    if k == 1 {
        return brown(a, b, 1, p);
    }
    if coprime(a, b, k, p) {
        return Some(constant(n, p));
    }
    zippel(a, b, k, p).or_else(|| brown(a, b, k, p))
}

/// Brown's dense algorithm: the last variable evaluated, the images' gcds recursively, Newton interpolation
fn brown(a: &MPoly, b: &MPoly, k: usize, p: &Fp) -> Option<MPoly> {
    let n = a.terms[0].0.len();
    if a.is_const() || b.is_const() {
        return Some(constant(n, p));
    }
    let v = k - 1;
    if k == 1 {
        let (ua, ub) = (to_rec(a, 0), to_rec(b, 0));
        let g = u_gcd(&ua[&vec![0; n]], &ub[&vec![0; n]], p);
        let mut r = Rec::new();
        r.insert(vec![0; n], g);
        return Some(from_rec(&r, 0));
    }
    let (mut ra, mut rb) = (to_rec(a, v), to_rec(b, v));
    // the contents in x_v, divided out
    let ca = ra.values().fold(vec![], |g, u| u_gcd(&g, u, p));
    let cb = rb.values().fold(vec![], |g, u| u_gcd(&g, u, p));
    let c = u_gcd(&ca, &cb, p);
    for u in ra.values_mut() {
        *u = u_divrem(u, &ca, p).0;
    }
    for u in rb.values_mut() {
        *u = u_divrem(u, &cb, p).0;
    }
    let (la, lb) = (ra.values().next_back().unwrap().clone(), rb.values().next_back().unwrap().clone());
    let g = u_gcd(&la, &lb, p);
    let dva = ra.values().map(|u| u.len()).max().unwrap_or(1) - 1;
    let dvb = rb.values().map(|u| u.len()).max().unwrap_or(1) - 1;
    let bound = (g.len() - 1) + dva.min(dvb) + 1;
    let (pa, pb) = (from_rec(&ra, v), from_rec(&rb, v));
    let with_content = |h: &Rec| -> MPoly {
        let mut r = h.clone();
        for u in r.values_mut() {
            *u = u_mul(u, &c, p);
        }
        from_rec(&r, v).monic(p)
    };

    // the interpolation so far: H (lc g), its leading monomial, the points' product q(x_v), the number of points
    let mut state: Option<(Rec, Vec<u32>, UPoly, usize)> = None;
    let limit = 4 * bound as u64 + 64;
    for x in 1..=limit {
        if x >= p.p {
            break;
        }
        if u_eval(&la, x, p) == 0 || u_eval(&lb, x, p) == 0 {
            continue;
        }
        let ci = gcd(&eval_rec(&ra, x, p), &eval_rec(&rb, x, p), k - 1, p)?;
        if ci.is_const() {
            // the primitive parts are coprime: the gcd is the gcd of the contents
            let mut r = Rec::new();
            r.insert(vec![0; n], c.clone());
            return Some(from_rec(&r, v).monic(p));
        }
        let gx = u_eval(&g, x, p);
        let lm = ci.terms[0].0.clone();
        let fresh = |ci: &MPoly| -> (Rec, Vec<u32>, UPoly, usize) {
            let h: Rec = ci.terms.iter().map(|(e, cf)| (e.clone(), vec![p.mul(*cf, gx)])).collect();
            (h, lm.clone(), vec![subm(0, x, p), p.one], 1)
        };
        let changed = match &mut state {
            None => {
                state = Some(fresh(&ci));
                true
            }
            Some((_, cur, _, _)) if lm < *cur => {
                state = Some(fresh(&ci));
                true
            }
            Some((_, cur, _, _)) if lm > *cur => continue,
            Some((h, _, q, npts)) => {
                let qi = inv(u_eval(q, x, p), p);
                let mut changed = false;
                let img: BTreeMap<&Vec<u32>, u64> = ci.terms.iter().map(|(e, cf)| (e, p.mul(*cf, gx))).collect();
                let keys: Vec<Vec<u32>> = h.keys().cloned().chain(img.keys().map(|e| (*e).clone())).collect();
                for key in keys {
                    let hv = h.get(&key).map_or(0, |u| u_eval(u, x, p));
                    let want = img.get(&key).copied().unwrap_or(0);
                    if hv != want {
                        changed = true;
                        let corr = u_scale(q, p.mul(subm(want, hv, p), qi), p);
                        let u = h.entry(key.clone()).or_default();
                        *u = u_add(u, &corr, p);
                        if u.is_empty() {
                            h.remove(&key);
                        }
                    }
                }
                *q = u_mul(q, &[subm(0, x, p), p.one], p);
                *npts += 1;
                changed
            }
        };
        let (h, _, _, npts) = state.as_ref().unwrap();
        if !changed || *npts >= bound {
            // the primitive part of H in x_v, tested by division
            let cont = h.values().fold(vec![], |acc, u| u_gcd(&acc, u, p));
            let mut hp = h.clone();
            for u in hp.values_mut() {
                *u = u_divrem(u, &cont, p).0;
            }
            let cand = from_rec(&hp, v).monic(p);
            if divides(&pa, &cand, p) && divides(&pb, &cand, p) {
                return Some(with_content(&hp));
            }
            if *npts >= bound {
                state = None;
            }
        }
    }
    None
}

/// a * b mod p
fn mul(a: &MPoly, b: &MPoly, p: &Fp) -> MPoly {
    let mut t = Vec::with_capacity(a.terms.len() * b.terms.len());
    for (ea, ca) in &a.terms {
        for (eb, cb) in &b.terms {
            t.push((ea.iter().zip(eb).map(|(x, y)| x + y).collect(), p.mul(*ca, *cb)));
        }
    }
    MPoly::normalise(t, p)
}

/// a / g mod p when g divides a exactly
fn div_exact(a: &MPoly, g: &MPoly, p: &Fp) -> Option<MPoly> {
    let (lm, lc) = &g.terms[0];
    let li = inv(*lc, p);
    let mut r: BTreeMap<Vec<u32>, u64> = a.terms.iter().cloned().collect();
    let mut q = Vec::new();
    while let Some((e, c)) = r.iter().next_back().map(|(e, &c)| (e.clone(), c)) {
        if e.iter().zip(lm).any(|(x, y)| x < y) {
            return None;
        }
        let shift: Vec<u32> = e.iter().zip(lm).map(|(x, y)| x - y).collect();
        let qc = p.mul(c, li);
        for (ge, gc) in &g.terms {
            let k: Vec<u32> = ge.iter().zip(&shift).map(|(x, y)| x + y).collect();
            let ent = r.entry(k.clone()).or_insert(0);
            *ent = subm(*ent, p.mul(qc, *gc), p);
            if *ent == 0 {
                r.remove(&k);
            }
        }
        q.push((shift, qc));
    }
    Some(MPoly { terms: q })
}

/// the coefficients of a in x0, as polynomials free of x0, from the highest power down: (power, coefficient)
fn coeffs_x0(a: &MPoly) -> Vec<(u32, MPoly)> {
    let mut out: Vec<(u32, MPoly)> = Vec::new();
    for (e, c) in &a.terms {
        let mut f = e.clone();
        f[0] = 0;
        match out.last_mut() {
            Some((d, m)) if *d == e[0] => m.terms.push((f, *c)),
            _ => out.push((e[0], MPoly { terms: vec![(f, *c)] })),
        }
    }
    out
}

/// the gcd of a's coefficients in x0 (monic)
fn content_x0(a: &MPoly, k: usize, p: &Fp) -> Option<MPoly> {
    let n = a.terms[0].0.len();
    let mut cs = coeffs_x0(a);
    cs.sort_by_key(|(_, m)| m.terms.len());
    let mut g: Option<MPoly> = None;
    for (_, m) in cs {
        let next = match g {
            None => m.monic(p),
            Some(g) => gcd(&g, &m, k, p)?,
        };
        if next.is_const() {
            return Some(constant(n, p));
        }
        g = Some(next);
    }
    g
}

/// the gcd of every coefficient in x0 of a and of b (monic), the smallest first, stopping at a constant
fn content_pair(a: &MPoly, b: &MPoly, k: usize, p: &Fp) -> Option<MPoly> {
    let n = a.terms[0].0.len();
    let mut cs: Vec<MPoly> = coeffs_x0(a).into_iter().chain(coeffs_x0(b)).map(|(_, m)| m).collect();
    cs.sort_by_key(|m| m.terms.len());
    // a single-term coefficient makes the gcd a monomial, which divides a and b: 1, as gcd() has divided their
    // monomial contents out before calling zippel
    if cs[0].terms.len() == 1 {
        return Some(constant(n, p));
    }
    let mut g: Option<MPoly> = None;
    for m in cs {
        let next = match g {
            None => m.monic(p),
            Some(g) => gcd(&g, &m, k, p)?,
        };
        if next.is_const() {
            return Some(constant(n, p));
        }
        g = Some(next);
    }
    g
}

/// Is gcd(a, b) = 1? Decided, when it can be, from one univariate image in each variable x_i (the others at a
/// random point where the leading coefficients in x_i do not vanish): the gcd's leading coefficient in x_i divides
/// theirs, so its degree in x_i is the image gcd's degree or less. A gcd of 1 in every variable proves a gcd of 1.
/// false means "not proven".
fn coprime(a: &MPoly, b: &MPoly, k: usize, p: &Fp) -> bool {
    let mut rng = Rng(p.p ^ 0x5DEE_CE66_D1CE_4E5B ^ ((a.terms.len() as u64) << 20) ^ b.terms.len() as u64);
    let z: Vec<u64> = (0..k).map(|_| rng.next(p.p)).collect();
    // z_j^e by tables, per variable
    let maxe: Vec<u32> = (0..k).map(|j| maxdeg(a, j).max(maxdeg(b, j))).collect();
    let pw: Vec<Vec<u64>> = (0..k).map(|j| {
        let mut v = vec![p.one];
        for e in 0..maxe[j] as usize {
            let next = p.mul(v[e], z[j]);
            v.push(next);
        }
        v
    }).collect();
    let image = |x: &MPoly, i: usize| -> UPoly {
        let mut u = vec![0u64; maxe[i] as usize + 1];
        for (e, c) in &x.terms {
            let mut m = *c;
            for j in 0..k {
                if j != i && e[j] > 0 {
                    m = p.mul(m, pw[j][e[j] as usize]);
                }
            }
            u[e[i] as usize] = addm(u[e[i] as usize], m, p);
        }
        u
    };
    (0..k).all(|i| {
        let (da, db) = (maxdeg(a, i) as usize, maxdeg(b, i) as usize);
        if da == 0 || db == 0 {
            // free of x_i in one of them: so is the gcd
            return true;
        }
        let (ua, ub) = (image(a, i), image(b, i));
        ua[da] != 0 && ub[db] != 0 && u_gcd(&ua[..=da], &ub[..=db], p).len() == 1
    })
}

/// splitmix64: the random points, the same on every node for the same input (D1)
struct Rng(u64);

impl Rng {
    fn next(&mut self, p: u64) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        1 + (z ^ (z >> 31)) % (p - 1)
    }
}

fn maxdeg(a: &MPoly, i: usize) -> u32 {
    a.terms.iter().map(|t| t.0[i]).max().unwrap_or(0)
}

/// Evaluation of the parameters x1.. at a point, by tables of powers.
struct Evaluator {
    k: usize,
    maxe: Vec<u32>,
}

impl Evaluator {
    fn powers(&self, y: &[u64], p: &Fp) -> Vec<Vec<u64>> {
        (0..self.k).map(|i| {
            let mut v = Vec::with_capacity(self.maxe[i] as usize + 1);
            let mut x = p.one;
            v.push(p.one);
            for _ in 0..self.maxe[i] {
                x = p.mul(x, y[i]);
                v.push(x);
            }
            v
        }).collect()
    }

    fn mono(&self, e: &[u32], pw: &[Vec<u64>], p: &Fp) -> u64 {
        let mut m = p.one;
        for i in 1..self.k {
            if e[i] > 0 {
                m = p.mul(m, pw[i][e[i] as usize]);
            }
        }
        m
    }

    /// a at the point, as a dense univariate polynomial in x0 (untrimmed, length deg + 1)
    fn univariate(&self, a: &MPoly, deg: usize, pw: &[Vec<u64>], p: &Fp) -> UPoly {
        let mut u = vec![0u64; deg + 1];
        for (e, c) in &a.terms {
            let d = e[0] as usize;
            u[d] = addm(u[d], p.mul(*c, self.mono(e, pw, p)), p);
        }
        u
    }
}

/// A polynomial prepared for the points of one stage of zippel (parameter x_j interpolated): per term its power of x0
/// and of x_j, its coefficient times the later parameters' anchors (base), and its monomial in the earlier parameters
/// at r (step). The stage's points put the earlier parameters at r^1, r^2, ..., so each next point costs one product
/// per term.
struct Staged {
    d0: Vec<usize>,
    ej: Vec<u32>,
    base: Vec<u64>,
    step: Vec<u64>,
    deg: usize,
    maxej: u32,
}

impl Staged {
    fn new(a: &MPoly, deg: usize, j: usize, k: usize, r: &[u64], z: &[u64], p: &Fp) -> Staged {
        let mut st = Staged { d0: vec![], ej: vec![], base: vec![], step: vec![], deg, maxej: 0 };
        for (e, c) in &a.terms {
            let mut base = *c;
            for i in j + 1..k {
                if e[i] > 0 {
                    base = p.mul(base, p.m.pow(z[i], e[i] as u64));
                }
            }
            let mut step = p.one;
            for i in 1..j {
                if e[i] > 0 {
                    step = p.mul(step, p.m.pow(r[i], e[i] as u64));
                }
            }
            st.d0.push(e[0] as usize);
            st.ej.push(e[j]);
            st.base.push(base);
            st.step.push(step);
            st.maxej = st.maxej.max(e[j]);
        }
        st
    }

    /// the univariates in x0 at x_j = t and the earlier parameters at r^1 .. r^count
    fn images(&self, t: u64, count: usize, p: &Fp) -> Vec<UPoly> {
        let mut tp = vec![p.one];
        for i in 0..self.maxej as usize {
            let next = p.mul(tp[i], t);
            tp.push(next);
        }
        let mut cur: Vec<u64> = self.base.iter().zip(&self.ej).map(|(&b, &e)| p.mul(b, tp[e as usize])).collect();
        (0..count).map(|_| {
            let mut u = vec![0u64; self.deg + 1];
            for (i, c) in cur.iter_mut().enumerate() {
                *c = p.mul(*c, self.step[i]);
                u[self.d0[i]] = addm(u[self.d0[i]], *c, p);
            }
            u
        }).collect()
    }
}

/// Zippel's sparse interpolation (R. Zippel, EUROSAM 1979; Geddes, Czapor, Labahn 7.5), x0 the main variable: the
/// univariate gcds in x0 at points of the parameters x1 .. x(k-1), scaled to gamma = gcd(lc(a), lc(b)) (leading
/// coefficients in x0), are the values of H = gamma / lc(G) G. H is found one parameter at a time: densely in x1,
/// then for each further parameter from the support found so far by transposed Vandermonde systems, each
/// interpolation stopping when one more point agrees. G is the primitive part of H. None when a random point
/// proves unlucky (a degree that drops, equal nodes, a candidate that does not divide).
fn zippel(a: &MPoly, b: &MPoly, k: usize, p: &Fp) -> Option<MPoly> {
    let n = a.terms[0].0.len();
    // G = c Gp: c the gcd of every coefficient in x0 of a and of b, Gp the gcd of the primitive parts. The
    // univariate images do not see c (a nonzero constant there), so a and b are not divided by their contents.
    let c = content_pair(a, b, k, p)?;
    let (pa, pb) = (a, b);
    let (da, db) = (maxdeg(pa, 0) as usize, maxdeg(pb, 0) as usize);
    if da == 0 || db == 0 {
        return Some(c.monic(p));
    }
    let la = coeffs_x0(pa).swap_remove(0).1;
    let lb = coeffs_x0(pb).swap_remove(0).1;
    let gamma = gcd(&la, &lb, k, p)?;
    let ev = Evaluator {
        k,
        maxe: (0..k).map(|i| maxdeg(pa, i).max(maxdeg(pb, i)).max(maxdeg(&gamma, i))).collect(),
    };
    let seed = p.p ^ (k as u64).wrapping_mul(0x2545_F491_4F6C_DD1D) ^ ((a.terms.len() as u64) << 32) ^ b.terms.len() as u64;
    let mut rng = Rng(seed);
    // the image at a point y (y[0] unused): gamma(y) times the monic gcd in x0; None at a point where a leading
    // coefficient vanishes
    let image = |y: &[u64]| -> Option<UPoly> {
        let pw = ev.powers(y, p);
        let ua = ev.univariate(pa, da, &pw, p);
        let ub = ev.univariate(pb, db, &pw, p);
        if ua[da] == 0 || ub[db] == 0 {
            return None;
        }
        let gm = gamma.terms.iter().fold(0, |s, (e, cf)| addm(s, p.mul(*cf, ev.mono(e, &pw, p)), p));
        if gm == 0 {
            return None;
        }
        Some(u_scale(&u_gcd(&ua, &ub, p), gm, p))
    };
    // the anchors of the parameters not yet interpolated, and the degree of G in x0
    let z: Vec<u64> = (0..k).map(|_| rng.next(p.p)).collect();
    let d0 = image(&z)?.len().checked_sub(1)?;
    if d0 == 0 {
        return Some(c.monic(p));
    }
    // the image from the univariates of a, b and gamma's value: Some(None) at an unlucky point (a leading
    // coefficient vanishes, or a degree above d0), None for a degree below d0 (d0 itself was unlucky)
    let image_of = |ua: &[u64], ub: &[u64], gm: u64| -> Option<Option<UPoly>> {
        if ua[da] == 0 || ub[db] == 0 || gm == 0 {
            return Some(None);
        }
        let mut u = u_scale(&u_gcd(ua, ub, p), gm, p);
        match u.len() - 1 {
            d if d > d0 => Some(None),
            d if d < d0 => None,
            _ => {
                u.resize(d0 + 1, 0);
                Some(Some(u))
            }
        }
    };
    let bound = |i: usize| (maxdeg(&gamma, i) + maxdeg(pa, i).min(maxdeg(pb, i))) as usize + 1;

    // H_d for d = 0 .. d0, each a list of (exponents, coefficient) in the parameters interpolated so far
    let mut h: Vec<Vec<Term>> = vec![vec![(vec![0; n], 0)]; d0 + 1];
    for j in 1..k {
        // the coefficients' supports so far; every (d, s) slot gets a univariate interpolant in x_j
        let supports: Vec<Vec<Vec<u32>>> = h.iter().map(|hd| hd.iter().map(|t| t.0.clone()).collect()).collect();
        let r: Vec<u64> = (0..k).map(|_| rng.next(p.p)).collect();
        let nodes: Vec<Vec<u64>> = supports.iter().map(|sd| sd.iter().map(|e| {
            let mut m = p.one;
            for i in 1..j {
                for _ in 0..e[i] {
                    m = p.mul(m, r[i]);
                }
            }
            m
        }).collect()).collect();
        let t_max = supports.iter().map(|s| s.len()).max().unwrap_or(1);
        let nslots: usize = supports.iter().map(|s| s.len()).sum();
        let mut interp: Vec<UPoly> = vec![vec![]; nslots];
        let mut q: UPoly = vec![p.one];
        let mut npts = 0usize;
        let cap = bound(j) + 1;
        let (sa, sb, sg) = (Staged::new(pa, da, j, k, &r, &z, p), Staged::new(pb, db, j, k, &r, &z, p), Staged::new(&gamma, 0, j, k, &r, &z, p));
        loop {
            if npts > cap {
                return None;
            }
            let t = rng.next(p.p);
            // the values of every slot at x_j = t
            let mut vals: Vec<u64> = Vec::with_capacity(nslots);
            let count = if j == 1 { 1 } else { t_max };
            let (ia, ib, ig) = (sa.images(t, count, p), sb.images(t, count, p), sg.images(t, count, p));
            let mut imgs: Vec<UPoly> = Vec::with_capacity(count);
            for kk in 0..count {
                match image_of(&ia[kk], &ib[kk], ig[kk][0])? {
                    None => break,
                    Some(u) => imgs.push(u),
                }
            }
            if imgs.len() < count {
                continue;
            }
            if j == 1 {
                vals.extend(imgs.swap_remove(0));
            } else {
                for (d, sd) in supports.iter().enumerate() {
                    let values: Vec<u64> = (0..sd.len()).map(|kk| imgs[kk][d]).collect();
                    let nd: Vec<u64> = nodes[d].iter().map(|&x| p.m.from(x)).collect();
                    let vs: Vec<u64> = values.iter().map(|&x| p.m.from(x)).collect();
                    let sol: Vec<u64> = crate::ffslp::tvandermonde(p.p, &nd, &vs)?.into_iter().map(|x| p.m.to(x)).collect();
                    for (s, &v) in sol.iter().enumerate() {
                        vals.push(p.mul(v, inv(nodes[d][s], p)));
                    }
                }
            }
            // Newton: stop when the interpolants already give these values
            let qt = u_eval(&q, t, p);
            let agree = interp.iter().zip(&vals).all(|(u, &v)| u_eval(u, t, p) == v);
            if agree && npts > 0 {
                break;
            }
            let qi = inv(qt, p);
            for (u, &v) in interp.iter_mut().zip(&vals) {
                let diff = subm(v, u_eval(u, t, p), p);
                if diff != 0 {
                    *u = u_add(u, &u_scale(&q, p.mul(diff, qi), p), p);
                }
            }
            q = u_mul(&q, &[subm(0, t, p), p.one], p);
            npts += 1;
        }
        // the new supports, with x_j
        let mut slot = 0;
        for (d, sd) in supports.iter().enumerate() {
            let mut next = Vec::new();
            for e in sd {
                for (deg, &cf) in interp[slot].iter().enumerate() {
                    if cf != 0 {
                        let mut f = e.clone();
                        f[j] = deg as u32;
                        next.push((f, cf));
                    }
                }
                slot += 1;
            }
            h[d] = next;
        }
    }
    // H, then its primitive part in x0
    let mut terms = Vec::new();
    for (d, hd) in h.iter().enumerate() {
        for (e, cf) in hd {
            let mut f = e.clone();
            f[0] = d as u32;
            terms.push((f, *cf));
        }
    }
    let hpoly = MPoly::normalise(terms, p);
    if hpoly.terms.is_empty() {
        return None;
    }
    let g = if gamma.is_const() {
        hpoly.monic(p)
    } else {
        let cont = content_x0(&hpoly, k, p)?;
        div_exact(&hpoly, &cont, p)?.monic(p)
    };
    if maxdeg(&g, 0) as usize != d0 || !divides(pa, &g, p) || !divides(pb, &g, p) {
        return None;
    }
    Some(mul(&c, &g, p).monic(p))
}

/// the gcd of a and b mod p, monic; the zero polynomial's gcd with b is b made monic
pub fn gcd_terms(n: usize, a: Vec<Term>, b: Vec<Term>, p: u64) -> Option<MPoly> {
    let f = Fp::new(p);
    let into = |t: Vec<Term>| MPoly::normalise(t.into_iter().map(|(e, c)| (e, f.m.to(c % p))).collect(), &f);
    let (a, b) = (into(a), into(b));
    let g = match (a.terms.is_empty(), b.terms.is_empty()) {
        (true, true) => Some(MPoly { terms: vec![] }),
        (true, false) => Some(b.monic(&f)),
        (false, true) => Some(a.monic(&f)),
        _ => gcd(&a, &b, n.max(1), &f),
    }?;
    Some(MPoly { terms: g.terms.into_iter().map(|(e, c)| (e, f.m.from(c))).collect() })
}

#[cfg(test)]
mod tests {
    use super::*;

    const P: u64 = 9223372036854775783; // the largest prime below 2^63

    fn poly(terms: &[(&[u32], i64)]) -> MPoly {
        MPoly::normalise(terms.iter().map(|(e, c)| (e.to_vec(), c.rem_euclid(P as i64) as u64)).collect(), &Fp::new(P))
    }

    fn mul(a: &MPoly, b: &MPoly) -> MPoly {
        let mut t = Vec::new();
        for (ea, ca) in &a.terms {
            for (eb, cb) in &b.terms {
                t.push((ea.iter().zip(eb).map(|(x, y)| x + y).collect(), crate::ffslp::mulmod(*ca, *cb, P)));
            }
        }
        MPoly::normalise(t, &Fp::new(P))
    }

    #[test]
    fn univariate() {
        // (x - 1)(x + 2) and (x - 1)(x + 3)
        let a = poly(&[(&[2], 1), (&[1], 1), (&[0], -2)]);
        let b = poly(&[(&[2], 1), (&[1], 2), (&[0], -3)]);
        assert_eq!(gcd_terms(1, a.terms, b.terms, P).unwrap(), poly(&[(&[1], 1), (&[0], -1)]));
    }

    #[test]
    fn multivariate_with_contents() {
        // g = x y + z^2 + 3, a = g (x + y z + 1) (z + 1), b = g (x^2 - z) (z + 1)^2
        let g = poly(&[(&[1, 1, 0], 1), (&[0, 0, 2], 1), (&[0, 0, 0], 3)]);
        let z1 = poly(&[(&[0, 0, 1], 1), (&[0, 0, 0], 1)]);
        let a = mul(&mul(&g, &poly(&[(&[1, 0, 0], 1), (&[0, 1, 1], 1), (&[0, 0, 0], 1)])), &z1);
        let b = mul(&mul(&g, &poly(&[(&[2, 0, 0], 1), (&[0, 0, 1], -1)])), &mul(&z1, &z1));
        assert_eq!(gcd_terms(3, a.terms, b.terms, P).unwrap(), mul(&g, &z1));
    }

    #[test]
    fn coprime_and_monomial_factors() {
        // x^2 y and x y^3 (x + 1): x y
        let a = poly(&[(&[2, 1], 1)]);
        let b = poly(&[(&[2, 3], 1), (&[1, 3], 1)]);
        assert_eq!(gcd_terms(2, a.terms, b.terms, P).unwrap(), poly(&[(&[1, 1], 1)]));
        let c = poly(&[(&[1, 0], 1), (&[0, 1], 1)]);
        let d = poly(&[(&[1, 0], 1), (&[0, 1], -1)]);
        assert_eq!(gcd_terms(2, c.terms, d.terms, P).unwrap(), poly(&[(&[0, 0], 1)]));
    }

    #[test]
    fn zippel_agrees_with_brown() {
        // random sparse polynomials in 4 variables: g a and g b
        let f = Fp::new(P);
        let mut rng = Rng(12345);
        let mut rnd = |terms: usize, deg: u64| {
            let t: Vec<Term> = (0..terms).map(|_| ((0..4).map(|_| (rng.next(P) % (deg + 1)) as u32).collect(), rng.next(P) % 1000 + 1)).collect();
            MPoly::normalise(t, &f)
        };
        for _ in 0..20 {
            let (g, a, b) = (rnd(4, 4), rnd(5, 4), rnd(5, 4));
            // zippel's inputs have no monomial content (gcd divides it out first)
            let strip = |x: MPoly| {
                let m: Vec<u32> = (0..4).map(|i| x.terms.iter().map(|t| t.0[i]).min().unwrap()).collect();
                MPoly { terms: x.terms.into_iter().map(|(e, c)| (e.iter().zip(&m).map(|(u, v)| u - v).collect(), c)).collect() }
            };
            let (x, y) = (strip(super::mul(&a, &g, &f)), strip(super::mul(&b, &g, &f)));
            let z = zippel(&x, &y, 4, &f).expect("zippel");
            assert_eq!(z, brown(&x, &y, 4, &f).unwrap());
            let (x, y) = (super::mul(&a, &g, &f), super::mul(&b, &g, &f));
            assert_eq!(gcd(&x, &y, 4, &f).unwrap(), brown(&x, &y, 4, &f).unwrap());
            assert!(divides(&gcd(&x, &y, 4, &f).unwrap(), &g.clone().monic(&f), &f));
        }
    }
}
