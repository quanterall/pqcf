//! The torontonian (AN 4.2, FS A.10):
//!
//!   tor(O) = sum over subsets Z of the c clicked modes of (-1)^(c - |Z|) / sqrt(det(I - O_Z)),
//!
//! O is 2c x 2c, mode j owning the indices j and j + c. I - O = Q^-1 is Hermitian positive definite, so each
//! I - O_Z has a Cholesky factor L (I - O_Z = L L^dagger) with a real positive diagonal and the term is
//! sign / prod(L_ii). The subsets are visited by depth-first search over the modes in increasing order (exclude,
//! then include); including a mode borders L with its two rows in O(|Z|^2). O(c^2 2^c) time, O(c^2) memory.
//!
//! The loop torontonian (FS 15 item 2; thewalrus's `ltor`, Bulmer, Quesada and Paesani) carries a displacement
//! vector gamma: each term is multiplied by exp(1/2 gamma_Z^T (I - O_Z)^-1 conj(gamma_Z)). The factor below is of
//! conj(I - O_Z) = L L^dagger, so that quadratic form is |w|^2 for w = L^-1 gamma_Z, which grows by one
//! forward-substitution entry per index pushed, O(|Z|): every term stays real and the cost O(c^2 2^c).
//!
//! D1: the first min(6, c) decisions split the tree into at most 64 subtrees, the blocks of FS 8.5, processed
//! and summed in the order of their index (bit t = mode t included).
#![forbid(unsafe_code)] // NFR-04: kernels are safe Rust; unsafe only in the Rustler glue (lib.rs)
use crate::num::{Sum, C64};

pub const SUBTREE_DEPTH: usize = 6;
pub const MAX_CLICKS: usize = 40;

#[derive(Debug, PartialEq)]
pub enum Error {
    BadShape(&'static str),
    TooLarge(usize),
    Cancelled,
    NonFinite(&'static str),
    /// I - O_Z is not positive definite: O is not of the form I - Q^-1 of a physical state
    NotPositiveDefinite,
}

#[derive(Debug, Clone, Copy)]
pub struct TorResult {
    pub value: f64,
    pub scale: f64,
    pub terms: u64,
}

fn conj(a: C64) -> C64 {
    C64::new(a.re, -a.im)
}

/// The Cholesky factor of I - O on an ordered set of indices, grown and shrunk one index at a time.
struct Factor<'a> {
    o: &'a [C64],
    dim: usize,          // 2c
    idx: Vec<usize>,     // the indices in the factor, in order
    l: Vec<Vec<C64>>,    // row k has k + 1 entries
    diag_prod: Vec<f64>, // prod of L_ii up to row k
    gamma: Option<&'a [C64]>,
    w: Vec<C64>,         // L^-1 conj(gamma) on idx (the loop torontonian only)
    w_norm: Vec<f64>,    // |w|^2 up to row k
}

impl<'a> Factor<'a> {
    fn m(&self, i: usize, j: usize) -> C64 {
        let v = self.o[i * self.dim + j];
        let id = if i == j { 1.0 } else { 0.0 };
        C64::new(id - v.re, -v.im)
    }

    fn push(&mut self, i: usize) -> Result<(), Error> {
        let k = self.idx.len();
        let mut row = vec![C64::ZERO; k + 1];
        for j in 0..k {
            let mut s = self.m(i, self.idx[j]);
            for (rt, lt) in row.iter().zip(&self.l[j]).take(j) {
                s = s.add(rt.mul(conj(*lt)).scale(-1.0));
            }
            let d = self.l[j][j].re;
            row[j] = C64::new(s.re / d, s.im / d);
        }
        let mut dd = self.m(i, i).re;
        for x in row.iter().take(k) {
            dd -= x.re * x.re + x.im * x.im;
        }
        if dd.is_nan() || dd <= 0.0 {
            return Err(Error::NotPositiveDefinite);
        }
        let lii = dd.sqrt();
        row[k] = C64::new(lii, 0.0);
        let prev = if k == 0 { 1.0 } else { self.diag_prod[k - 1] };
        self.diag_prod.push(prev * lii);
        if let Some(g) = self.gamma {
            // the factor is of conj(I - O_Z) (m() conjugates; the torontonian only needs its real determinant), so
            // w = L^-1 gamma: |w|^2 = gamma^dagger conj(I - O_Z)^-1 gamma = gamma^T (I - O_Z)^-1 conj(gamma), real.
            // Forward substitution: w_k = (gamma_i - sum_j L_kj w_j) / L_kk
            let mut s = g[i];
            for (lkj, wj) in row.iter().zip(&self.w).take(k) {
                s = s.add(lkj.mul(*wj).scale(-1.0));
            }
            let wk = C64::new(s.re / lii, s.im / lii);
            let prev = if k == 0 { 0.0 } else { self.w_norm[k - 1] };
            self.w_norm.push(prev + wk.re * wk.re + wk.im * wk.im);
            self.w.push(wk);
        }
        self.l.push(row);
        self.idx.push(i);
        Ok(())
    }

    fn pop(&mut self) {
        self.idx.pop();
        self.l.pop();
        self.diag_prod.pop();
        if self.gamma.is_some() {
            self.w.pop();
            self.w_norm.pop();
        }
    }

    /// exp(1/2 |L^-1 conj(gamma)|^2): 1 for the plain torontonian
    fn loop_factor(&self) -> f64 {
        match self.gamma {
            Some(_) => (0.5 * self.w_norm.last().copied().unwrap_or(0.0)).exp(),
            None => 1.0,
        }
    }

    fn det_sqrt(&self) -> f64 {
        self.diag_prod.last().copied().unwrap_or(1.0)
    }
}

struct Walk<'b, 'a> {
    f: Factor<'a>,
    c: usize,
    sum: &'b mut Sum,
    scale: &'b mut Sum,
    steps: u64,
    cancel: &'b (dyn Fn() -> bool + Sync),
}

impl Walk<'_, '_> {
    /// every subset of the modes from `t` on, added to the current one
    fn dfs(&mut self, t: usize, size: usize) -> Result<(), Error> {
        if t == self.c {
            let term = self.f.loop_factor() / self.f.det_sqrt();
            let term = if (self.c - size) % 2 == 1 { -term } else { term };
            self.sum.add(term);
            self.scale.add(term.abs());
            self.steps += 1;
            if self.steps.is_multiple_of(4096) && (self.cancel)() {
                return Err(Error::Cancelled);
            }
            return Ok(());
        }
        self.dfs(t + 1, size)?;
        self.f.push(t)?;
        self.f.push(t + self.c)?;
        let r = self.dfs(t + 1, size + 1);
        self.f.pop();
        self.f.pop();
        r
    }
}

/// tor(O) for a 2c x 2c matrix O (row-major).
pub fn torontonian(
    o: &[C64],
    dim: usize,
    threads: usize,
    cancel: &(dyn Fn() -> bool + Sync),
) -> Result<TorResult, Error> {
    walk(o, None, dim, threads, cancel)
}

/// ltor(O, gamma) for a 2c x 2c matrix O and a 2c-vector gamma (thewalrus's ordering: mode j at j and j + c).
pub fn loop_torontonian(
    o: &[C64],
    gamma: &[C64],
    dim: usize,
    threads: usize,
    cancel: &(dyn Fn() -> bool + Sync),
) -> Result<TorResult, Error> {
    if gamma.len() != dim {
        return Err(Error::BadShape("gamma must have 2c entries"));
    }
    if gamma.iter().any(|z| !z.is_finite()) {
        return Err(Error::NonFinite("input"));
    }
    walk(o, Some(gamma), dim, threads, cancel)
}

fn walk(
    o: &[C64],
    gamma: Option<&[C64]>,
    dim: usize,
    threads: usize,
    cancel: &(dyn Fn() -> bool + Sync),
) -> Result<TorResult, Error> {
    if dim % 2 == 1 || o.len() != dim * dim {
        return Err(Error::BadShape("O must be 2c x 2c"));
    }
    if o.iter().any(|z| !z.is_finite()) {
        return Err(Error::NonFinite("input"));
    }
    let c = dim / 2;
    if c == 0 {
        return Ok(TorResult { value: 1.0, scale: 1.0, terms: 1 });
    }
    if c > MAX_CLICKS {
        return Err(Error::TooLarge(c));
    }
    let depth = c.min(SUBTREE_DEPTH);
    let subtrees = 1u64 << depth;
    let per_tree = 1u64 << (c - depth);
    let run = |b: u64| -> Result<(f64, f64), Error> {
        let mut f = Factor { o, dim, idx: Vec::new(), l: Vec::new(), diag_prod: Vec::new(), gamma, w: Vec::new(),
                             w_norm: Vec::new() };
        let mut size = 0;
        for t in 0..depth {
            if b >> t & 1 == 1 {
                f.push(t)?;
                f.push(t + c)?;
                size += 1;
            }
        }
        let mut sum = Sum::new(per_tree);
        let mut scale = Sum::new(per_tree);
        let mut w = Walk { f, c, sum: &mut sum, scale: &mut scale, steps: 0, cancel };
        w.dfs(depth, size)?;
        Ok((sum.value(), scale.value()))
    };
    let results: Vec<Result<(f64, f64), Error>> = if threads <= 1 {
        (0..subtrees).map(run).collect()
    } else {
        let per = subtrees.div_ceil(threads as u64);
        std::thread::scope(|sc| {
            let run = &run;
            let hs: Vec<_> = (0..threads as u64)
                .map(|k| sc.spawn(move || (k * per..((k + 1) * per).min(subtrees)).map(run).collect::<Vec<_>>()))
                .collect();
            hs.into_iter().flat_map(|h| h.join().expect("torontonian worker panicked")).collect()
        })
    };
    let mut sum = Sum::new(subtrees);
    let mut scale = Sum::new(subtrees);
    for r in results {
        let (v, s) = r?;
        sum.add(v);
        scale.add(s);
    }
    if !sum.value().is_finite() {
        return Err(Error::NonFinite("output"));
    }
    Ok(TorResult { value: sum.value(), scale: scale.value(), terms: 1u64 << c })
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{Complex, DMatrix};

    fn never() -> bool {
        false
    }

    /// O = I - Q^-1 for a random Hermitian positive definite Q = I + B B^dagger / 4
    fn random_o(c: usize, seed: u64) -> Vec<C64> {
        let d = 2 * c;
        let mut s = seed.wrapping_mul(6364136223846793005).wrapping_add(7);
        let mut next = || {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((s >> 11) as f64) / (1u64 << 53) as f64 - 0.5
        };
        let b = DMatrix::from_fn(d, d, |_, _| Complex::new(next(), next()));
        let q = DMatrix::<Complex<f64>>::identity(d, d) + &b * b.adjoint() * Complex::new(0.25, 0.0);
        let qi = q.try_inverse().unwrap();
        let o = DMatrix::<Complex<f64>>::identity(d, d) - qi;
        (0..d * d).map(|k| C64::new(o[(k / d, k % d)].re, o[(k / d, k % d)].im)).collect()
    }

    /// the definition, with determinants by LU
    fn by_definition(o: &[C64], d: usize) -> f64 {
        let c = d / 2;
        let mut total = 0.0;
        for z in 0u64..(1 << c) {
            let idx: Vec<usize> = (0..c).filter(|&j| z >> j & 1 == 1).flat_map(|j| [j, j + c]).collect();
            let k = idx.len();
            let m = DMatrix::from_fn(k, k, |a, b| {
                let v = o[idx[a] * d + idx[b]];
                Complex::new(if a == b { 1.0 } else { 0.0 } - v.re, -v.im)
            });
            let det = if k == 0 { Complex::new(1.0, 0.0) } else { m.determinant() };
            let sign = if (c - z.count_ones() as usize) % 2 == 1 { -1.0 } else { 1.0 };
            total += sign / det.re.sqrt();
        }
        total
    }

    #[test]
    fn single_mode_squeezed_vacuum() {
        // Q = [[cosh^2 r, -sinh r cosh r], [-sinh r cosh r, cosh^2 r]]: P(click) = 1 - 1/cosh r (FS A.10)
        for r in [0.1f64, 0.5, 1.0, 2.0] {
            let (ch, sh) = (r.cosh(), r.sinh());
            let det = ch.powi(4) - (sh * ch).powi(2);
            let (a, b) = (ch * ch / det, sh * ch / det); // Q^-1 = [[a, b], [b, a]]
            let o = vec![C64::new(1.0 - a, 0.0), C64::new(-b, 0.0), C64::new(-b, 0.0), C64::new(1.0 - a, 0.0)];
            let t = torontonian(&o, 2, 1, &never).unwrap().value;
            let p = t / det.sqrt();
            assert!((p - (1.0 - 1.0 / ch)).abs() < 1e-14, "r = {r}: {p}");
        }
    }

    #[test]
    fn against_the_definition() {
        for c in 0..=10 {
            for seed in 0..3 {
                let o = random_o(c, 100 * c as u64 + seed);
                let r = torontonian(&o, 2 * c, 1, &never).unwrap();
                let want = by_definition(&o, 2 * c);
                assert!((r.value - want).abs() <= 1e-12 * r.scale.max(1.0), "c = {c}: {} vs {want}", r.value);
            }
        }
    }

    #[test]
    fn threads_change_no_bit() {
        let o = random_o(12, 4);
        let one = torontonian(&o, 24, 1, &never).unwrap();
        let five = torontonian(&o, 24, 5, &never).unwrap();
        assert_eq!(one.value.to_bits(), five.value.to_bits());
    }

    #[test]
    fn not_positive_definite() {
        // O = 2 I: I - O = -I
        let o = vec![C64::new(2.0, 0.0), C64::ZERO, C64::ZERO, C64::new(2.0, 0.0)];
        assert_eq!(torontonian(&o, 2, 1, &never).unwrap_err(), Error::NotPositiveDefinite);
    }
    /// the loop torontonian by its definition: every subset, an explicit inverse and determinant
    fn ltor_by_definition(o: &[C64], g: &[C64], c: usize) -> f64 {
        let n = 2 * c;
        let mut total = 0.0;
        for mask in 0u64..(1 << c) {
            let modes: Vec<usize> = (0..c).filter(|&j| mask >> j & 1 == 1).collect();
            let idx: Vec<usize> = modes.iter().copied().chain(modes.iter().map(|&j| j + c)).collect();
            let k = idx.len();
            let sign = if (c - modes.len()) % 2 == 1 { -1.0 } else { 1.0 };
            if k == 0 {
                total += sign;
                continue;
            }
            // I - O_Z itself, not its conjugate (the kernel factors the conjugate: an independent check)
            let mz = DMatrix::from_fn(k, k, |a, b| {
                let v = o[idx[a] * n + idx[b]];
                Complex::new(if a == b { 1.0 } else { 0.0 } - v.re, -v.im).conj()
            });
            let gz = nalgebra::DVector::from_fn(k, |a, _| Complex::new(g[idx[a]].re, g[idx[a]].im));
            let inv = mz.clone().try_inverse().unwrap();
            let q = (gz.transpose() * inv * gz.map(|z| z.conj()))[(0, 0)];
            let det = mz.determinant();
            total += sign * (0.5 * q.re).exp() / det.re.sqrt();
        }
        total
    }

    #[test]
    fn loop_torontonian_against_its_definition() {
        for c in 1..=7 {
            for seed in 0..3u64 {
                let o = random_o(c, 100 + seed);
                let mut s = crate::rng::Stream::new(7 + seed, c as u64);
                let g: Vec<C64> = (0..2 * c).map(|_| { let (a, b) = s.normal_pair(); C64::new(0.4 * a, 0.4 * b) }).collect();
                let got = loop_torontonian(&o, &g, 2 * c, 1, &never).unwrap();
                let want = ltor_by_definition(&o, &g, c);
                assert!((got.value - want).abs() <= 1e-12 * got.scale, "c = {c}: {} vs {}", got.value, want);
                // threads change nothing (D1)
                assert_eq!(loop_torontonian(&o, &g, 2 * c, 4, &never).unwrap().value.to_bits(), got.value.to_bits());
            }
        }
    }

    #[test]
    fn loop_torontonian_with_zero_gamma_is_the_torontonian() {
        for c in 1..=6 {
            let o = random_o(c, 42 + c as u64);
            let z = vec![C64::ZERO; 2 * c];
            assert_eq!(loop_torontonian(&o, &z, 2 * c, 1, &never).unwrap().value.to_bits(),
                       torontonian(&o, 2 * c, 1, &never).unwrap().value.to_bits());
        }
    }
}
