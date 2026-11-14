//! Wang's multivariate Hensel lifting modulo a prime p < 2^63 (decision record SY2; P. S. Wang, Math. Comp. 32
//! (1978) 1215; Geddes, Czapor, Labahn, "Algorithms for Computer Algebra", algorithms 6.2-6.4). pq_sym_factor
//! factors the univariate image f(x1, a2, .., an) over the integers, gives each factor the leading coefficient L(a)
//! (L the leading coefficient of f in x1, the polynomial lifted being L^(r-1) f), and calls this for several primes;
//! the lifted factors are unique integer polynomials, brought back by Chinese remaindering.
//!
//! Polynomials are sparse: a map from packed exponents (an equal share of the 64 bits per variable, x1 lowest) to a
//! coefficient in Montgomery form.

use crate::mgcd::{addm, inv, subm, u_divrem, u_mul, u_trim, Fp, UPoly};
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

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

type SPoly = HashMap<u64, u64, BuildHasherDefault<KeyHasher>>;

/// terms as the NIF passes them: (exponents, coefficient)
pub type Terms = Vec<(Vec<u32>, u64)>;

/// Where each variable's exponent sits in a packed key.
struct Layout {
    bits: u32,
    n: usize,
}

impl Layout {
    fn get(&self, k: u64, j: usize) -> u32 {
        ((k >> (self.bits * j as u32)) & ((1u64 << self.bits) - 1)) as u32
    }
    fn without(&self, k: u64, j: usize) -> u64 {
        k & !(((1u64 << self.bits) - 1) << (self.bits * j as u32))
    }
    fn unit(&self, j: usize, e: u32) -> u64 {
        (e as u64) << (self.bits * j as u32)
    }
    fn pack(&self, e: &[u32]) -> u64 {
        e.iter().enumerate().fold(0, |k, (j, &x)| k | self.unit(j, x))
    }
    fn unpack(&self, k: u64) -> Vec<u32> {
        (0..self.n).map(|j| self.get(k, j)).collect()
    }
}

struct Ctx<'a> {
    p: &'a Fp,
    l: Layout,
    pts: Vec<u64>,
    /// pw[j][e] = pts[j]^e, binom[n][k]
    pw: Vec<Vec<u64>>,
    binom: Vec<Vec<u64>>,
    /// the univariate factors and s_i, the inverse of prod_(k /= i) u_k modulo u_i
    us: Vec<UPoly>,
    ss: Vec<UPoly>,
    dmax: u32,
}

fn acc(a: &mut SPoly, k: u64, v: u64, p: &Fp) {
    if v == 0 {
        return;
    }
    match a.get_mut(&k) {
        Some(w) => {
            *w = addm(*w, v, p);
            if *w == 0 {
                a.remove(&k);
            }
        }
        None => {
            a.insert(k, v);
        }
    }
}

impl Ctx<'_> {
    fn add(&self, a: &SPoly, b: &SPoly) -> SPoly {
        let mut r = a.clone();
        for (&k, &v) in b {
            acc(&mut r, k, v, self.p);
        }
        r
    }

    fn sub(&self, a: &SPoly, b: &SPoly) -> SPoly {
        let mut r = a.clone();
        for (&k, &v) in b {
            acc(&mut r, k, subm(0, v, self.p), self.p);
        }
        r
    }

    /// the products of every pair of terms, sorted by key and merged (faster than inserting each into a hash map)
    fn mul(&self, a: &SPoly, b: &SPoly) -> SPoly {
        let mut t: Vec<(u64, u64)> = Vec::with_capacity(a.len() * b.len());
        for (&ka, &va) in a {
            for (&kb, &vb) in b {
                t.push((ka + kb, self.p.mul(va, vb)));
            }
        }
        t.sort_unstable_by_key(|x| x.0);
        let mut r = SPoly::with_capacity_and_hasher(t.len().min(1 << 16), Default::default());
        let mut i = 0;
        while i < t.len() {
            let (k, mut v) = t[i];
            i += 1;
            while i < t.len() && t[i].0 == k {
                v = addm(v, t[i].1, self.p);
                i += 1;
            }
            if v != 0 {
                r.insert(k, v);
            }
        }
        r
    }

    fn prod(&self, fs: &[SPoly]) -> SPoly {
        let mut r = fs[0].clone();
        for f in &fs[1..] {
            r = self.mul(&r, f);
        }
        r
    }

    fn deg(&self, a: &SPoly, j: usize) -> u32 {
        a.keys().map(|&k| self.l.get(k, j)).max().unwrap_or(0)
    }

    /// x_j = pts[j], exponent j then 0
    fn eval(&self, a: &SPoly, j: usize) -> SPoly {
        let mut r = SPoly::default();
        for (&k, &v) in a {
            let e = self.l.get(k, j) as usize;
            acc(&mut r, self.l.without(k, j), self.p.mul(v, self.pw[j][e]), self.p);
        }
        r
    }

    /// the coefficient of (x_j - pts[j])^m
    fn taylor(&self, a: &SPoly, j: usize, m: u32) -> SPoly {
        let mut r = SPoly::default();
        for (&k, &v) in a {
            let e = self.l.get(k, j);
            if e >= m {
                let c = self.p.mul(self.p.mul(v, self.binom[e as usize][m as usize]), self.pw[j][(e - m) as usize]);
                acc(&mut r, self.l.without(k, j), c, self.p);
            }
        }
        r
    }

    /// (x_j - pts[j])^m
    fn shiftpow(&self, j: usize, m: u32) -> SPoly {
        let mut r = SPoly::default();
        let neg = subm(0, self.pts[j], self.p);
        let mut negpow = vec![self.p.one];
        for i in 0..m as usize {
            let next = self.p.mul(negpow[i], neg);
            negpow.push(next);
        }
        for i in 0..=m {
            acc(&mut r, self.l.unit(j, i), self.p.mul(self.binom[m as usize][i as usize], negpow[(m - i) as usize]), self.p);
        }
        r
    }

    fn to_dense(&self, a: &SPoly) -> UPoly {
        let d = self.deg(a, 0) as usize;
        let mut u = vec![0u64; d + 1];
        for (&k, &v) in a {
            u[self.l.get(k, 0) as usize] = v;
        }
        u_trim(u)
    }

    fn sparse_of(&self, u: &[u64]) -> SPoly {
        u.iter().enumerate().filter(|(_, &c)| c != 0).map(|(e, &c)| (self.l.unit(0, e as u32), c)).collect()
    }

    fn cofactors(&self, a: &[SPoly]) -> Vec<SPoly> {
        (0..a.len()).map(|i| {
            let others: Vec<SPoly> = a.iter().enumerate().filter(|(k, _)| *k != i).map(|(_, x)| x.clone()).collect();
            self.prod(&others)
        }).collect()
    }

    /// sigma_i with sum sigma_i prod_(k /= i) a_k = c, in x1 and x2 .. x_top (algorithm 6.2), in powers of
    /// (x_top - a_top): coefficient m of sigma_i from the residual's coefficient m, which needs only the coefficients
    /// below m of every sigma_i and of the cofactors
    fn diophant(&self, a: &[SPoly], c: &SPoly, top: usize) -> Vec<SPoly> {
        if top == 0 {
            let cd = self.to_dense(c);
            return self.us.iter().zip(&self.ss).map(|(u, s)| self.sparse_of(&u_divrem(&u_mul(&cd, s, self.p), u, self.p).1)).collect();
        }
        let anew: Vec<SPoly> = a.iter().map(|x| self.eval(x, top)).collect();
        let bs = self.cofactors(a);
        let dc = self.deg(c, top);
        let db = bs.iter().map(|b| self.deg(b, top)).max().unwrap_or(0);
        let dmax = self.dmax.max(dc);
        // the cofactors' coefficients, as needed
        let bco: Vec<Vec<SPoly>> = bs.iter().map(|b| (0..=dmax.min(db)).map(|t| self.taylor(b, top, t)).collect()).collect();
        let mut sig: Vec<Vec<SPoly>> = vec![Vec::new(); a.len()];
        for m in 0..=dmax {
            let mut r = self.taylor(c, top, m);
            for (si, bi) in sig.iter().zip(&bco) {
                for (t, st) in si.iter().enumerate() {
                    let k = m as usize - t;
                    if k < bi.len() && !st.is_empty() && !bi[k].is_empty() {
                        r = self.sub(&r, &self.mul(st, &bi[k]));
                    }
                }
            }
            let ds = if r.is_empty() { vec![SPoly::default(); a.len()] } else { self.diophant(&anew, &r, top - 1) };
            for (si, d) in sig.iter_mut().zip(ds) {
                si.push(d);
            }
        }
        sig.iter().map(|si| self.untaylor(si, top)).collect()
    }

    /// sum_m co[m] (x_j - a_j)^m
    fn untaylor(&self, co: &[SPoly], j: usize) -> SPoly {
        let mut r = SPoly::default();
        for (m, c) in co.iter().enumerate() {
            if !c.is_empty() {
                r = self.add(&r, &self.mul(c, &self.shiftpow(j, m as u32)));
            }
        }
        r
    }
}

/// {G, S} with S a + T b = G (G not normalised)
fn u_egcd(a: &[u64], b: &[u64], p: &Fp) -> (UPoly, UPoly) {
    let (mut r0, mut r1) = (a.to_vec(), b.to_vec());
    let (mut s0, mut s1): (UPoly, UPoly) = (vec![p.one], vec![]);
    while !r1.is_empty() {
        let (q, r) = u_divrem(&r0, &r1, p);
        let qs = u_mul(&q, &s1, p);
        let s2 = u_trim((0..s0.len().max(qs.len())).map(|i| subm(*s0.get(i).unwrap_or(&0), *qs.get(i).unwrap_or(&0), p)).collect());
        r0 = std::mem::replace(&mut r1, r);
        s0 = std::mem::replace(&mut s1, s2);
    }
    (r0, s0)
}

/// The lifted factors of f (terms over n variables, coefficients in [0, p)): f = prod U_i, U_i = us[i] at
/// x_j = pts[j] (j >= 1), every U_i's leading coefficient in x1 lc. None when the factors do not lift (a bad prime
/// or point) or the exponents do not pack.
pub fn lift(p: u64, n: usize, f: Terms, pts: Vec<u64>, us: Vec<Vec<u64>>, lc: Terms) -> Option<Vec<Terms>> {
    let fp = Fp::new(p);
    let maxe = f.iter().flat_map(|(e, _)| e.iter().copied()).max().unwrap_or(0);
    let bits = 64 / n.max(1) as u32;
    // products of the factors reach r times f's degree in a variable before they cancel
    let reach = maxe as u64 * us.len().max(1) as u64 + 1;
    if n == 0 || n > 16 || bits < 64 && reach >= 1u64 << bits {
        return None;
    }
    let l = Layout { bits, n };
    let to = |t: &[(Vec<u32>, u64)]| -> SPoly {
        let mut r = SPoly::default();
        for (e, c) in t {
            acc(&mut r, l.pack(e), fp.m.to(*c % p), &fp);
        }
        r
    };
    let (f, lc) = (to(&f), to(&lc));
    let us: Vec<UPoly> = us.into_iter().map(|u| u_trim(u.into_iter().map(|c| fp.m.to(c % p)).collect())).collect();
    let pts: Vec<u64> = (0..n).map(|j| if j == 0 { 0 } else { fp.m.to(pts[j - 1] % p) }).collect();
    let dmax = (1..n).map(|j| f.keys().map(|&k| l.get(k, j)).max().unwrap_or(0)).max().unwrap_or(0);
    let top = reach as usize;
    let pw: Vec<Vec<u64>> = pts.iter().map(|&x| {
        let mut v = vec![fp.one];
        for i in 0..top {
            let next = fp.mul(v[i], x);
            v.push(next);
        }
        v
    }).collect();
    let mut binom = vec![vec![fp.one]];
    for i in 1..=top {
        let prev = &binom[i - 1];
        let row: Vec<u64> = (0..=i).map(|k| {
            let a = if k > 0 { prev[k - 1] } else { 0 };
            let b = if k < i { prev[k] } else { 0 };
            addm(a, b, &fp)
        }).collect();
        binom.push(row);
    }
    // s_i: the inverse of prod_(k /= i) u_k modulo u_i
    let mut ss = Vec::with_capacity(us.len());
    for i in 0..us.len() {
        let b = us.iter().enumerate().filter(|(k, _)| *k != i).fold(vec![fp.one], |acc, (_, u)| u_mul(&acc, u, &fp));
        let b = u_divrem(&b, &us[i], &fp).1;
        let (g, s) = u_egcd(&b, &us[i], &fp);
        if g.len() != 1 {
            return None;
        }
        ss.push(u_divrem(&s.iter().map(|&c| fp.mul(c, inv(g[0], &fp))).collect::<Vec<_>>(), &us[i], &fp).1);
    }
    let cx = Ctx { p: &fp, l, pts, pw, binom, us, ss, dmax };

    // A_j: f with x_(j+1) .. at the point
    let mut as_: Vec<SPoly> = vec![f.clone()];
    for j in (1..n).rev() {
        let next = cx.eval(&as_[0], j);
        as_.insert(0, next);
    }
    let mut u: Vec<SPoly> = cx.us.iter().map(|x| cx.sparse_of(x)).collect();
    let r = u.len();
    for (j, a) in as_.iter().enumerate().skip(1) {
        let mut lcj = lc.clone();
        for k in j + 1..n {
            lcj = cx.eval(&lcj, k);
        }
        let dj = cx.deg(a, j);
        // every factor in powers of (x_j - a_j): coefficient 0 the factor so far, the others its forced leading
        // coefficient's (lc's coefficient m times x1^deg)
        let degs: Vec<u32> = u.iter().map(|x| cx.deg(x, 0)).collect();
        let lco: Vec<SPoly> = (0..=dj).map(|m| cx.taylor(&lcj, j, m)).collect();
        let mut co: Vec<Vec<SPoly>> = (0..r).map(|i| {
            (0..=dj as usize).map(|m| if m == 0 { u[i].clone() } else {
                lco[m].iter().map(|(&k, &v)| (k + cx.l.unit(0, degs[i]), v)).collect()
            }).collect()
        }).collect();
        let uold = u.clone();
        // pr[l][m]: coefficient m of the product of the first l + 1 factors
        let mut pr: Vec<Vec<SPoly>> = vec![Vec::new(); r];
        let coef = |pr: &Vec<Vec<SPoly>>, co: &Vec<Vec<SPoly>>, l: usize, m: usize| -> SPoly {
            if l == 0 {
                return co[0][m].clone();
            }
            let mut s = SPoly::default();
            for t in 0..=m {
                if !pr[l - 1][t].is_empty() && !co[l][m - t].is_empty() {
                    s = cx.add(&s, &cx.mul(&pr[l - 1][t], &co[l][m - t]));
                }
            }
            s
        };
        for m in 0..=dj as usize {
            for l in 0..r {
                let v = coef(&pr, &co, l, m);
                pr[l].push(v);
            }
            if m == 0 {
                continue;
            }
            let c = cx.sub(&cx.taylor(a, j, m as u32), &pr[r - 1][m]);
            if c.is_empty() {
                continue;
            }
            let ds = cx.diophant(&uold, &c, j - 1);
            for (i, d) in ds.into_iter().enumerate() {
                co[i][m] = cx.add(&co[i][m], &d);
            }
            for l in 0..r {
                pr[l][m] = coef(&pr, &co, l, m);
            }
        }
        u = co.iter().map(|ci| cx.untaylor(ci, j)).collect();
    }
    if !cx.sub(&f, &cx.prod(&u)).is_empty() {
        return None;
    }
    Some(u.iter().map(|x| {
        let mut t: Vec<(Vec<u32>, u64)> = x.iter().map(|(&k, &v)| (cx.l.unpack(k), fp.m.from(v))).collect();
        t.sort_unstable_by(|a, b| b.0.cmp(&a.0));
        t
    }).collect())
}

/// The monic irreducible factors of a monic square-free f modulo a prime p (odd, < 2^63): distinct-degree
/// factorisation with the Frobenius map as a matrix (h^p mod f = sum h_i (x^(i p) mod f)), then Cantor and
/// Zassenhaus's equal-degree splitting; with count_only, the number of factors without splitting. A count of 0: f
/// is not square-free modulo p.
pub fn factor_mod(p: u64, f: Vec<u64>, count_only: bool) -> (usize, Vec<Vec<u64>>) {
    let fp = Fp::new(p);
    let mut f: UPoly = u_trim(f.into_iter().map(|c| fp.m.to(c % p)).collect());
    let n = f.len() - 1;
    // square-free modulo p? (else the prime is unsuitable: count 0)
    let df: UPoly = u_trim(f.iter().enumerate().skip(1).map(|(i, &c)| fp.mul(c, fp.m.to(i as u64 % p))).collect());
    if df.is_empty() || crate::mgcd::u_gcd(&f, &df, &fp).len() > 1 {
        return (0, vec![]);
    }
    // q[i] = x^(i p) mod f
    let xp = upow(&[0, fp.one], p, &f, &fp);
    let mut q: Vec<UPoly> = vec![vec![fp.one]];
    for i in 1..n {
        let next = u_divrem(&u_mul(&q[i - 1], &xp, &fp), &f, &fp).1;
        q.push(next);
    }
    let frob = |h: &[u64], q: &[UPoly], f: &[u64]| -> UPoly {
        // q's rows are reduced modulo the original f, so they may be longer than the current one
        let mut r: UPoly = vec![0; n + 1];
        for (i, &c) in h.iter().enumerate() {
            if c != 0 {
                for (k, &v) in q[i].iter().enumerate() {
                    r[k] = addm(r[k], fp.mul(c, v), &fp);
                }
            }
        }
        u_divrem(&u_trim(r), f, &fp).1
    };
    let mut parts: Vec<(UPoly, usize)> = Vec::new();
    let mut h = xp.clone();
    let mut d = 1;
    while f.len() > 2 * d {
        let mut hx = h.clone();
        if hx.len() < 2 {
            hx.resize(2, 0);
        }
        hx[1] = subm(hx[1], fp.one, &fp);
        let g = crate::mgcd::u_gcd(&u_trim(hx), &f, &fp);
        if g.len() > 1 {
            f = u_divrem(&f, &g, &fp).0;
            h = u_divrem(&h, &f, &fp).1;
            parts.push((g, d));
        }
        d += 1;
        if f.len() > 2 * d {
            // the Frobenius matrix of the current f: q reduced mod f still maps h to h^p mod f
            h = frob(&h, &q, &f);
        }
    }
    if f.len() > 1 {
        let df = f.len() - 1;
        parts.push((f, df));
    }
    let count: usize = parts.iter().map(|(g, d)| (g.len() - 1) / d).sum();
    if count_only {
        return (count, vec![]);
    }
    let mut out = Vec::new();
    let mut seed = p ^ (n as u64) << 32;
    for (g, d) in parts {
        edf(g, d, p, &fp, &mut seed, &mut out);
    }
    let mut res: Vec<Vec<u64>> = out.into_iter().map(|g| g.into_iter().map(|c| fp.m.from(c)).collect()).collect();
    res.sort();
    (count, res)
}

fn upow(b: &[u64], mut e: u64, f: &[u64], fp: &Fp) -> UPoly {
    let mut base = u_divrem(b, f, fp).1;
    let mut acc: UPoly = vec![fp.one];
    while e > 0 {
        if e & 1 == 1 {
            acc = u_divrem(&u_mul(&acc, &base, fp), f, fp).1;
        }
        base = u_divrem(&u_mul(&base, &base, fp), f, fp).1;
        e >>= 1;
    }
    acc
}

fn edf(g: UPoly, d: usize, p: u64, fp: &Fp, seed: &mut u64, out: &mut Vec<UPoly>) {
    let n = g.len() - 1;
    if n == d {
        out.push(g);
        return;
    }
    let e = {
        // (p^d - 1) / 2, as a sequence of exponentiations: a^((p^d - 1)/2) = (a^((p^d-1)/(p-1)))^((p-1)/2), and
        // a^((p^d-1)/(p-1)) = a a^p a^(p^2) .. a^(p^(d-1))
        (p - 1) / 2
    };
    loop {
        let a: UPoly = u_trim((0..n).map(|_| {
            *seed = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = *seed;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            fp.m.to((z ^ (z >> 31)) % p)
        }).collect());
        if a.len() < 2 {
            continue;
        }
        // the norm-like product a a^p .. a^(p^(d-1)) mod g, then its ((p-1)/2)-th power
        let mut t = a.clone();
        let mut ap = a.clone();
        for _ in 1..d {
            ap = upow(&ap, p, &g, fp);
            t = u_divrem(&u_mul(&t, &ap, fp), &g, fp).1;
        }
        let mut b = upow(&t, e, &g, fp);
        if b.is_empty() {
            continue;
        }
        b[0] = subm(b[0], fp.one, fp);
        let c = crate::mgcd::u_gcd(&u_trim(b), &g, fp);
        if c.len() > 1 && c.len() < g.len() {
            let q = crate::mgcd::u_monic(&u_divrem(&g, &c, fp).0, fp);
            edf(c, d, p, fp, seed, out);
            edf(q, d, p, fp, seed, out);
            return;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lifts_two_factors() {
        // f = (x y + 1)(x + y + 2) = x^2 y + x y^2 + 2 x y + x + y + 2; at y = 3: (3x + 1)(x + 5); L = y, r = 2:
        // lifting y f with factors (y x + 1) and (y x + y^2 + 2 y): images at 3 are 3x + 1 and 3x + 15
        let p = 9223372036854775783u64;
        let f0 = [(vec![2u32, 1], 1u64), (vec![1, 2], 1), (vec![1, 1], 2), (vec![1, 0], 1), (vec![0, 1], 1), (vec![0, 0], 2)];
        let f: Vec<(Vec<u32>, u64)> = f0.iter().map(|(e, c)| (vec![e[0], e[1] + 1], *c)).collect();
        let got = lift(p, 2, f, vec![3], vec![vec![1, 3], vec![15, 3]], vec![(vec![0, 1], 1)]).expect("lifts");
        let mut want1 = vec![(vec![1u32, 1], 1u64), (vec![0, 0], 1)];
        let mut want2 = vec![(vec![1u32, 1], 1u64), (vec![0, 2], 1), (vec![0, 1], 2)];
        want1.sort_unstable_by(|a, b| b.0.cmp(&a.0));
        want2.sort_unstable_by(|a, b| b.0.cmp(&a.0));
        assert_eq!(got, vec![want1, want2]);
    }

    #[test]
    fn factors_mod_p() {
        // x^4 + 1 mod 17 splits into 4 linear factors (17 = 1 mod 8); mod 3 into two quadratics
        let (n, fs) = factor_mod(17, vec![1, 0, 0, 0, 1], false);
        assert_eq!((n, fs.len()), (4, 4));
        assert!(fs.iter().all(|g| g.len() == 2 && g[1] == 1));
        let (n, fs) = factor_mod(3, vec![1, 0, 0, 0, 1], false);
        assert_eq!((n, fs.len()), (2, 2));
        assert_eq!(factor_mod(3, vec![1, 0, 0, 0, 1], true).0, 2);
        // (x + 1)(x^2 + x + 2)(x^3 + 2x + 1) mod 7, the last two irreducible mod 7? checked by the count only
        let (n, _) = factor_mod(101, vec![1, 0, 1], false);
        assert_eq!(n, 2); // x^2 + 1 splits mod 101 (101 = 1 mod 4)
    }
}
