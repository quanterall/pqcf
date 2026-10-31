//! The permanent family (AN 3): Glynn's formula with repeated rows and columns,
//!
//!   perm(B) = 2^-n * sum over k_i in 0..r_i of (-1)^(sum k_i) * prod_i C(r_i, k_i) * prod_j w_j(k)^c_j,
//!   w_j(k)  = sum_i (r_i - 2 k_i) * U[i][j],
//!
//! enumerated in a reflected mixed-radix Gray code. With every multiplicity 1 this is Glynn's formula in
//! binary Gray-code order, term for term and bit for bit, so `permanent` is this kernel with r = c = 1.
//!
//! The terms for k and r - k are equal; the first row's digit runs only over 0..=floor(r_0/2), and a term counts
//! twice unless 2 k_0 = r_0. That digit is the most significant one.
//!
//! D1 (AN 2.1, FS 8.5): the positions are cut into blocks of BLOCK steps; each block starts from w computed
//! directly, refreshes it directly every REFRESH steps, and the block sums are added in block order. The thread
//! count only decides which thread computes which block, so results are identical for any thread count.
#![forbid(unsafe_code)] // NFR-04: kernels are safe Rust; unsafe only in the Rustler glue (lib.rs)
use crate::num::{CSum, Sum, C64};

pub const BLOCK: u64 = 1 << 16;
pub const REFRESH: u64 = 4096;
/// Largest number of Gray-code positions accepted (2^50: about a quadrillion terms).
pub const MAX_POSITIONS: u64 = 1 << 50;

#[derive(Debug, PartialEq)]
pub enum Error {
    BadShape(&'static str),
    TooLarge(u64),
    Cancelled,
    NonFinite(&'static str),
}

#[derive(Debug, Clone, Copy)]
pub struct PermResult {
    pub value: C64,
    /// The scale S of AN 2.4: the sum of the terms' absolute values, normalised like the result.
    pub scale: f64,
    pub positions: u64,
    pub blocks: u64,
}

/// The reduced problem: rows and columns with non-zero multiplicity, rows in significance order (fastest first;
/// the halved first row last).
struct Plan<'a> {
    u: &'a [C64],
    cols: usize,
    rows_sig: Vec<usize>,
    radix: Vec<u64>,
    r_sig: Vec<u32>,
    col_idx: Vec<usize>,
    col_mult: Vec<u32>,
    binom: Vec<Vec<f64>>, // per significance digit: C(r, k) for k = 0..=r
    n: u32,
    positions: u64,
}

fn binomials(r: u32) -> Vec<f64> {
    let mut row = vec![1.0f64; r as usize + 1];
    for k in 1..=r as usize {
        // C(r, k) = C(r, k-1) * (r - k + 1) / k, exact while it fits in 53 bits
        row[k] = row[k - 1] * (r as usize - k + 1) as f64 / k as f64;
    }
    row
}

impl<'a> Plan<'a> {
    fn new(u: &'a [C64], rows: usize, cols: usize, r: &[u32], c: &[u32]) -> Result<Plan<'a>, Error> {
        if u.len() != rows * cols || r.len() != rows || c.len() != cols {
            return Err(Error::BadShape("multiplicities"));
        }
        let n: u64 = r.iter().map(|&x| x as u64).sum();
        if n != c.iter().map(|&x| x as u64).sum::<u64>() {
            return Err(Error::BadShape("row and column multiplicities differ in total"));
        }
        if u.iter().any(|z| !z.is_finite()) {
            return Err(Error::NonFinite("input"));
        }
        let nz: Vec<usize> = (0..rows).filter(|&i| r[i] > 0).collect();
        let mut rows_sig: Vec<usize> = nz.iter().skip(1).copied().collect();
        if let Some(&first) = nz.first() {
            rows_sig.push(first);
        }
        let last = rows_sig.len().saturating_sub(1);
        let radix: Vec<u64> = rows_sig
            .iter()
            .enumerate()
            .map(|(t, &i)| if t == last { r[i] as u64 / 2 + 1 } else { r[i] as u64 + 1 })
            .collect();
        let mut positions: u64 = 1;
        for &d in &radix {
            positions = positions.checked_mul(d).filter(|&p| p <= MAX_POSITIONS).ok_or(Error::TooLarge(u64::MAX))?;
        }
        let r_sig: Vec<u32> = rows_sig.iter().map(|&i| r[i]).collect();
        let col_idx: Vec<usize> = (0..cols).filter(|&j| c[j] > 0).collect();
        let col_mult = col_idx.iter().map(|&j| c[j]).collect();
        let binom = r_sig.iter().map(|&x| binomials(x)).collect();
        Ok(Plan { u, cols, rows_sig, radix, r_sig, col_idx, col_mult, binom, n: n as u32, positions })
    }

    /// Gray digits and directions at position p (reflected mixed radix): digit t is reflected when the number
    /// formed by the digits above it is odd.
    fn digits_at(&self, p: u64, g: &mut [u64], dir: &mut [i64]) {
        let mut below: u64 = 1;
        for t in 0..self.radix.len() {
            let a = (p / below) % self.radix[t];
            let high = p / (below * self.radix[t]);
            if high & 1 == 0 {
                g[t] = a;
                dir[t] = 1;
            } else {
                g[t] = self.radix[t] - 1 - a;
                dir[t] = -1;
            }
            below *= self.radix[t];
        }
    }

    /// w_j = sum_i (r_i - 2 k_i) U[i][j], over the reduced rows in significance order.
    fn w_direct(&self, g: &[u64], w: &mut [C64]) {
        for (jj, &j) in self.col_idx.iter().enumerate() {
            let mut acc = C64::ZERO;
            for (t, &i) in self.rows_sig.iter().enumerate() {
                let f = self.r_sig[t] as f64 - 2.0 * g[t] as f64;
                acc = acc.add(self.u[i * self.cols + j].scale(f));
            }
            w[jj] = acc;
        }
    }

    /// The term at digits g, without the 2^-n factor; and its modulus.
    #[inline]
    fn term(&self, g: &[u64], w: &[C64]) -> C64 {
        let mut prod = C64::ONE;
        for (jj, &m) in self.col_mult.iter().enumerate() {
            for _ in 0..m {
                prod = prod.mul(w[jj]);
            }
        }
        let last = self.radix.len() - 1;
        let mut coef: f64 = if 2 * g[last] == self.r_sig[last] as u64 { 1.0 } else { 2.0 };
        let mut parity = 0u64;
        for (binom, &k) in self.binom.iter().zip(g) {
            coef *= binom[k as usize];
            parity += k;
        }
        if parity & 1 == 1 {
            coef = -coef;
        }
        prod.scale(coef)
    }

    fn block(&self, b: u64, cancel: &dyn Fn() -> bool) -> Result<(C64, f64), Error> {
        let start = b * BLOCK;
        let end = (start + BLOCK).min(self.positions);
        let len = end - start;
        let d = self.radix.len();
        let mut g = vec![0u64; d];
        let mut dir = vec![0i64; d];
        let mut w = vec![C64::ZERO; self.col_idx.len()];
        self.digits_at(start, &mut g, &mut dir);
        self.w_direct(&g, &mut w);
        let mut sum = CSum::new(len);
        let mut scale = Sum::new(len);
        let t0 = self.term(&g, &w);
        sum.add(t0);
        scale.add(t0.abs());
        for p in start + 1..end {
            // the lowest digit that can move in its direction moves; the digits below it turn around
            let mut t = 0;
            loop {
                let next = g[t] as i64 + dir[t];
                if next >= 0 && (next as u64) < self.radix[t] {
                    break;
                }
                dir[t] = -dir[t];
                t += 1;
            }
            g[t] = (g[t] as i64 + dir[t]) as u64;
            if (p - start).is_multiple_of(REFRESH) {
                if cancel() {
                    return Err(Error::Cancelled);
                }
                self.w_direct(&g, &mut w);
            } else {
                // k_i moved by dir: w_j changes by -2 dir U[i][j]
                let i = self.rows_sig[t];
                let f = -2.0 * dir[t] as f64;
                for (jj, &j) in self.col_idx.iter().enumerate() {
                    w[jj] = w[jj].add(self.u[i * self.cols + j].scale(f));
                }
            }
            let term = self.term(&g, &w);
            sum.add(term);
            scale.add(term.abs());
        }
        Ok((sum.value(), scale.value()))
    }
}

/// The permanent of the matrix U (rows x cols, row-major) with row multiplicities r and column multiplicities c.
pub fn permanent_repeated(
    u: &[C64],
    rows: usize,
    cols: usize,
    r: &[u32],
    c: &[u32],
    threads: usize,
    cancel: &(dyn Fn() -> bool + Sync),
) -> Result<PermResult, Error> {
    let plan = Plan::new(u, rows, cols, r, c)?;
    if plan.n == 0 {
        return Ok(PermResult { value: C64::ONE, scale: 1.0, positions: 0, blocks: 0 });
    }
    let blocks = plan.positions.div_ceil(BLOCK);
    let results: Vec<Result<(C64, f64), Error>> = if threads <= 1 || blocks == 1 {
        (0..blocks).map(|b| plan.block(b, cancel)).collect()
    } else {
        let per = blocks.div_ceil(threads as u64);
        std::thread::scope(|s| {
            let plan = &plan;
            let handles: Vec<_> = (0..threads as u64)
                .map(|k| {
                    s.spawn(move || {
                        (k * per..((k + 1) * per).min(blocks)).map(|b| plan.block(b, cancel)).collect::<Vec<_>>()
                    })
                })
                .collect();
            handles.into_iter().flat_map(|h| h.join().expect("permanent worker panicked")).collect()
        })
    };
    let mut sum = CSum::new(blocks);
    let mut scale = Sum::new(blocks);
    for res in results {
        let (z, s) = res?;
        sum.add(z);
        scale.add(s);
    }
    let norm = 0.5f64.powi(plan.n as i32);
    let value = sum.value().scale(norm);
    let scale = scale.value() * norm;
    if !value.is_finite() || !scale.is_finite() {
        return Err(Error::NonFinite("output"));
    }
    Ok(PermResult { value, scale, positions: plan.positions, blocks })
}

/// The permanent of a square matrix: Glynn's formula (multiplicities all 1).
pub fn permanent(
    a: &[C64],
    n: usize,
    threads: usize,
    cancel: &(dyn Fn() -> bool + Sync),
) -> Result<PermResult, Error> {
    let ones = vec![1u32; n];
    permanent_repeated(a, n, n, &ones, &ones, threads, cancel)
}

/// FR-B02: probability of output pattern `out` (rows) from input pattern `inp` (columns) through the unitary U:
/// |perm(U_{T,S})|^2 / (prod S_i! prod T_j!). A pattern with a different photon number has probability 0.
pub fn bs_probability(
    u: &[C64],
    m: usize,
    inp: &[u32],
    out: &[u32],
    threads: usize,
    cancel: &(dyn Fn() -> bool + Sync),
) -> Result<(f64, PermResult), Error> {
    if inp.len() != m || out.len() != m || u.len() != m * m {
        return Err(Error::BadShape("occupations"));
    }
    if inp.iter().map(|&x| x as u64).sum::<u64>() != out.iter().map(|&x| x as u64).sum::<u64>() {
        let zero = PermResult { value: C64::ZERO, scale: 0.0, positions: 0, blocks: 0 };
        return Ok((0.0, zero));
    }
    let p = permanent_repeated(u, m, m, out, inp, threads, cancel)?;
    let mut fact = 1.0f64;
    for &k in inp.iter().chain(out.iter()) {
        for x in 2..=k {
            fact *= x as f64;
        }
    }
    let a = p.value.abs();
    Ok((a * a / fact, p))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn never() -> bool {
        false
    }

    /// A deterministic pseudo-random complex matrix (an LCG, independent of the kernel).
    fn matrix(n: usize, m: usize, seed: u64) -> Vec<C64> {
        let mut s = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        let mut next = || {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((s >> 11) as f64) / (1u64 << 53) as f64 - 0.5
        };
        (0..n * m).map(|_| C64::new(next(), next())).collect()
    }

    /// The definition: the sum over all permutations, O(n!).
    fn perm_definition(a: &[C64], n: usize) -> C64 {
        fn rec(a: &[C64], n: usize, row: usize, used: &mut Vec<bool>, acc: C64, out: &mut C64) {
            if row == n {
                *out = out.add(acc);
                return;
            }
            for j in 0..n {
                if !used[j] {
                    used[j] = true;
                    rec(a, n, row + 1, used, acc.mul(a[row * n + j]), out);
                    used[j] = false;
                }
            }
        }
        let mut out = C64::ZERO;
        rec(a, n, 0, &mut vec![false; n], C64::ONE, &mut out);
        out
    }

    fn close(a: C64, b: C64, tol: f64) -> bool {
        (a.re - b.re).abs() <= tol && (a.im - b.im).abs() <= tol
    }

    #[test]
    fn identity_and_all_ones() {
        for n in 1..=10 {
            let mut id = vec![C64::ZERO; n * n];
            for i in 0..n {
                id[i * n + i] = C64::ONE;
            }
            assert_eq!(permanent(&id, n, 1, &never).unwrap().value, C64::ONE);
            let ones = vec![C64::ONE; n * n];
            let fact: f64 = (1..=n).map(|x| x as f64).product();
            assert_eq!(permanent(&ones, n, 1, &never).unwrap().value, C64::new(fact, 0.0), "n = {n}");
        }
    }

    #[test]
    fn empty_and_one_by_one() {
        assert_eq!(permanent(&[], 0, 1, &never).unwrap().value, C64::ONE);
        let z = C64::new(0.3, -0.7);
        assert_eq!(permanent(&[z], 1, 1, &never).unwrap().value, z);
    }

    #[test]
    fn hong_ou_mandel() {
        // 50:50 beamsplitter (FS 3.3, theta = pi/4, phi = 0): one photon in each input never gives (1, 1)
        let h = std::f64::consts::FRAC_1_SQRT_2;
        let u = [C64::new(h, 0.0), C64::new(-h, 0.0), C64::new(h, 0.0), C64::new(h, 0.0)];
        let (p11, _) = bs_probability(&u, 2, &[1, 1], &[1, 1], 1, &never).unwrap();
        let (p20, _) = bs_probability(&u, 2, &[1, 1], &[2, 0], 1, &never).unwrap();
        let (p02, _) = bs_probability(&u, 2, &[1, 1], &[0, 2], 1, &never).unwrap();
        assert!(p11.abs() < 1e-15);
        assert!((p20 - 0.5).abs() < 1e-15 && (p02 - 0.5).abs() < 1e-15);
        assert_eq!(bs_probability(&u, 2, &[1, 1], &[1, 0], 1, &never).unwrap().0, 0.0);
    }

    #[test]
    fn against_the_definition() {
        for n in 1..=8 {
            for seed in 0..5 {
                let a = matrix(n, n, seed * 100 + n as u64);
                let r = permanent(&a, n, 1, &never).unwrap();
                let d = perm_definition(&a, n);
                assert!(close(r.value, d, 1e-14 * r.scale.max(1.0)), "n = {n}: {:?} vs {:?}", r.value, d);
            }
        }
    }

    /// Expands multiplicities into the full matrix and takes the definition's permanent.
    fn expanded(u: &[C64], rows: usize, cols: usize, r: &[u32], c: &[u32]) -> C64 {
        let ri: Vec<usize> = (0..rows).flat_map(|i| std::iter::repeat_n(i, r[i] as usize)).collect();
        let cj: Vec<usize> = (0..cols).flat_map(|j| std::iter::repeat_n(j, c[j] as usize)).collect();
        let n = ri.len();
        let b: Vec<C64> = ri.iter().flat_map(|&i| cj.iter().map(move |&j| u[i * cols + j])).collect();
        perm_definition(&b, n)
    }

    #[test]
    fn repeated_against_expansion() {
        let cases: &[(&[u32], &[u32])] = &[
            (&[2, 0, 1], &[1, 1, 1]),
            (&[3, 1], &[2, 2]),
            (&[0, 4, 0], &[1, 1, 2]),
            (&[2, 2, 2], &[3, 0, 3]),
            (&[1, 1, 1, 1], &[1, 1, 1, 1]),
            (&[5], &[5]),
            (&[2, 3], &[1, 4]),
        ];
        for (k, (r, c)) in cases.iter().enumerate() {
            let u = matrix(r.len(), c.len(), 7 + k as u64);
            let got = permanent_repeated(&u, r.len(), c.len(), r, c, 1, &never).unwrap();
            let want = expanded(&u, r.len(), c.len(), r, c);
            assert!(close(got.value, want, 1e-13 * got.scale.max(1.0)), "{r:?} {c:?}: {:?} vs {want:?}", got.value);
        }
    }

    #[test]
    fn gray_code_steps_match_the_direct_digits() {
        let u = matrix(3, 3, 1);
        let plan = Plan::new(&u, 3, 3, &[3, 2, 4], &[3, 3, 3]).unwrap();
        let d = plan.radix.len();
        let (mut g, mut dir) = (vec![0; d], vec![0; d]);
        let (mut g2, mut dir2) = (vec![0; d], vec![0; d]);
        plan.digits_at(0, &mut g, &mut dir);
        for p in 1..plan.positions {
            let mut t = 0;
            loop {
                let next = g[t] as i64 + dir[t];
                if next >= 0 && (next as u64) < plan.radix[t] {
                    break;
                }
                dir[t] = -dir[t];
                t += 1;
            }
            let before = g.clone();
            g[t] = (g[t] as i64 + dir[t]) as u64;
            plan.digits_at(p, &mut g2, &mut dir2);
            assert_eq!(g, g2, "position {p}");
            assert_eq!(before.iter().zip(&g).filter(|(a, b)| a != b).count(), 1);
        }
    }

    #[test]
    fn blocks_and_threads_do_not_change_a_bit() {
        // n = 19: 2^18 positions = 4 blocks
        let n = 19;
        let a = matrix(n, n, 42);
        let one = permanent(&a, n, 1, &never).unwrap();
        assert_eq!(one.blocks, 4);
        for t in [2, 3, 4, 8] {
            let r = permanent(&a, n, t, &never).unwrap();
            assert_eq!(r.value.re.to_bits(), one.value.re.to_bits());
            assert_eq!(r.value.im.to_bits(), one.value.im.to_bits());
            assert_eq!(r.scale.to_bits(), one.scale.to_bits());
        }
    }

    #[test]
    fn cancellation() {
        let a = matrix(20, 20, 3);
        assert_eq!(permanent(&a, 20, 1, &|| true).unwrap_err(), Error::Cancelled);
    }

    #[test]
    fn errors() {
        let u = matrix(2, 2, 0);
        assert!(matches!(permanent_repeated(&u, 2, 2, &[1, 1], &[2, 1], 1, &never), Err(Error::BadShape(_))));
        let mut bad = u.clone();
        bad[1] = C64::new(f64::NAN, 0.0);
        assert_eq!(permanent(&bad, 2, 1, &never).unwrap_err(), Error::NonFinite("input"));
        assert!(permanent_repeated(&u, 2, 2, &[60, 0], &[30, 30], 1, &never).is_ok());
        assert!(matches!(permanent(&matrix(52, 52, 0), 52, 1, &never), Err(Error::TooLarge(_))));
    }
}
