//! Straight-line programs evaluated over a prime field (decision record R1): the black box of the reconstruction.
//! A program is prepared once (its operations, exact rational constants), then evaluated modulo a prime p < 2^63
//! at a batch of points. A point where the program divides by zero mod p is reported as bad (every output
//! `BAD`), for the reconstruction to skip; a constant whose denominator is 0 mod p makes the prime unusable.

/// The value of every output at a bad point (never a residue, as p < 2^63).
pub const BAD: u64 = u64::MAX;

#[derive(Clone)]
pub enum Op {
    /// the point's coordinate i
    Input(usize),
    /// sign * num / den, num and den as 64-bit limbs, least significant first (den > 0)
    Const(i8, Vec<u64>, Vec<u64>),
    Mul(usize, usize),
    Add(usize, usize),
    Sub(usize, usize),
    /// an integer power; negative ones divide
    Pow(usize, i64),
}

#[derive(Clone)]
pub struct Program {
    pub ops: Vec<Op>,
    pub outputs: Vec<usize>,
    pub ninputs: usize,
}

#[derive(Debug, PartialEq, Eq)]
pub enum FfError {
    /// a constant's denominator is 0 mod p
    BadPrime,
    /// the program refers to a slot it hasn't computed yet, or a point has the wrong length
    Shape,
}

#[inline]
pub fn mulmod(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}

/// Shoup's multiplication by a fixed w < p (p < 2^63): w' = floor(w 2^64 / p) once, then w x mod p for any x < p
/// with two 64-bit products and one high product -- no 128-bit division, which on aarch64 is a library call
#[inline]
pub fn shoup_pre(w: u64, p: u64) -> u64 {
    (((w as u128) << 64) / p as u128) as u64
}

#[inline]
pub fn shoup_mul(x: u64, w: u64, wp: u64, p: u64) -> u64 {
    let q = ((x as u128 * wp as u128) >> 64) as u64;
    let r = x.wrapping_mul(w).wrapping_sub(q.wrapping_mul(p));
    if r >= p { r - p } else { r }
}

pub fn powmod(mut a: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1 % p;
    a %= p;
    while e > 0 {
        if e & 1 == 1 {
            r = mulmod(r, a, p);
        }
        a = mulmod(a, a, p);
        e >>= 1;
    }
    r
}

/// the inverse of a mod p (p prime), None for 0
/// the inverse mod p (prime) by the extended Euclidean algorithm: 64-bit divisions only (Fermat's a^(p-2) costs
/// about 126 multiplications, each a 128-bit remainder, a library call on aarch64)
pub fn invmod(a: u64, p: u64) -> Option<u64> {
    let a = a % p;
    if a == 0 {
        return None;
    }
    let (mut r0, mut r1) = (p, a);
    let (mut t0, mut t1) = (0i128, 1i128);
    while r1 != 0 {
        let q = r0 / r1;
        (r0, r1) = (r1, r0 - q * r1);
        (t0, t1) = (t1, t0 - q as i128 * t1);
    }
    (r0 == 1).then(|| t0.rem_euclid(p as i128) as u64)
}

/// a number given as 64-bit limbs (least significant first) mod p, by Horner's rule from the top limb
fn limbs_mod(limbs: &[u64], p: u64) -> u64 {
    let base = ((1u128 << 64) % p as u128) as u64;
    limbs.iter().rev().fold(0u64, |acc, &l| ((mulmod(acc, base, p) as u128 + (l % p) as u128) % p as u128) as u64)
}

impl Program {
    pub fn check(&self) -> Result<(), FfError> {
        for (s, op) in self.ops.iter().enumerate() {
            let ok = match *op {
                Op::Input(i) => i < self.ninputs,
                Op::Const(_, _, ref d) => d.iter().any(|&x| x != 0),
                Op::Mul(a, b) | Op::Add(a, b) | Op::Sub(a, b) => a < s && b < s,
                Op::Pow(a, _) => a < s,
            };
            if !ok {
                return Err(FfError::Shape);
            }
        }
        if self.outputs.iter().all(|&o| o < self.ops.len()) { Ok(()) } else { Err(FfError::Shape) }
    }

    /// The outputs at every point (`points`: ninputs residues per point), mod p: `outputs.len()` values per point,
    /// all `BAD` at a point where the program divides by zero.
    pub fn eval(&self, p: u64, points: &[u64]) -> Result<Vec<u64>, FfError> {
        if self.ninputs == 0 && !points.is_empty() || self.ninputs > 0 && !points.len().is_multiple_of(self.ninputs) {
            return Err(FfError::Shape);
        }
        // the constants once per prime
        let mut consts = vec![0u64; self.ops.len()];
        for (s, op) in self.ops.iter().enumerate() {
            if let Op::Const(sign, ref n, ref d) = *op {
                let di = invmod(limbs_mod(d, p), p).ok_or(FfError::BadPrime)?;
                let v = mulmod(limbs_mod(n, p), di, p);
                consts[s] = if sign < 0 { (p - v) % p } else { v };
            }
        }
        let npoints = if self.ninputs == 0 { 1 } else { points.len() / self.ninputs };
        let mut out = Vec::with_capacity(npoints * self.outputs.len());
        let mut v = vec![0u64; self.ops.len()];
        'points: for k in 0..npoints {
            let x = &points[k * self.ninputs..(k + 1) * self.ninputs];
            for (s, op) in self.ops.iter().enumerate() {
                v[s] = match *op {
                    Op::Input(i) => x[i] % p,
                    Op::Const(..) => consts[s],
                    Op::Mul(a, b) => mulmod(v[a], v[b], p),
                    Op::Add(a, b) => {
                        let t = v[a] + v[b];
                        if t >= p { t - p } else { t }
                    }
                    Op::Sub(a, b) => {
                        if v[a] >= v[b] { v[a] - v[b] } else { v[a] + p - v[b] }
                    }
                    Op::Pow(a, e) => {
                        let base = if e < 0 {
                            match invmod(v[a], p) {
                                Some(i) => i,
                                None => {
                                    out.extend(std::iter::repeat_n(BAD, self.outputs.len()));
                                    continue 'points;
                                }
                            }
                        } else {
                            v[a]
                        };
                        powmod(base, e.unsigned_abs(), p)
                    }
                };
            }
            out.extend(self.outputs.iter().map(|&o| v[o]));
        }
        Ok(out)
    }
}

/// The transposed Vandermonde system of Zippel's interpolation: the c with values[i] = sum_j c_j nodes_j^i
/// (i = 0 .. n-1) mod p, for distinct nodes, in O(n^2): with P(z) = prod_k (z - v_k) and Q_j = P / (z - v_j),
/// sum_i q_{j,i} values_i = c_j Q_j(v_j). None if two nodes are equal.
pub fn tvandermonde(p: u64, nodes: &[u64], values: &[u64]) -> Option<Vec<u64>> {
    let n = nodes.len();
    if values.len() != n {
        return None;
    }
    // P(z), coefficients from the constant term up, degree n
    let mut master = vec![0u64; n + 1];
    master[0] = 1;
    for (k, &v) in nodes.iter().enumerate() {
        let neg = (p - v % p) % p;
        for t in (0..=k + 1).rev() {
            let shifted = if t > 0 { master[t - 1] } else { 0 };
            master[t] = (mulmod(master[t], neg, p) + shifted) % p;
        }
    }
    let mut out = Vec::with_capacity(n);
    let mut q = vec![0u64; n];
    for &v in nodes {
        // Q = P / (z - v) by synthetic division, from the top
        let v = v % p;
        let mut carry = master[n];
        for t in (0..n).rev() {
            q[t] = carry;
            carry = (master[t] + mulmod(carry, v, p)) % p;
        }
        // Q(v), and sum_i q_i values_i
        let (mut qv, mut num) = (0u64, 0u64);
        for t in (0..n).rev() {
            qv = (mulmod(qv, v, p) + q[t]) % p;
        }
        for (t, &w) in values.iter().enumerate() {
            num = (num + mulmod(q[t], w % p, p)) % p;
        }
        out.push(mulmod(num, invmod(qv, p)?, p));
    }
    Some(out)
}

/// What the reconstruction samples (R1): a program's first output, or the solution of a linear system whose
/// entries are a program's outputs (phase 5): A u = B, both sparse: row r of A holds the columns `pattern[r]`
/// (increasing), its entries the first outputs row by row; row r of B the right-hand sides `bpattern[r]`
/// (increasing, < m), its entries the next outputs row by row. Solved by sparse elimination mod p at each point
/// (`solve_aug`); the value is component `want` of the solution for right-hand side `col` (bad where A is singular).
pub enum BlackBox {
    Prog(Program),
    Lin { prog: Program, n: usize, m: usize, want: usize, col: usize, pattern: Patterns, cache: LinCache, scale: Scale, plan: PlanSlot },
}

/// The learned elimination of a linear-system black box: learned at its first point, replayed at the others
/// (None if the first point was singular: the search-based elimination then serves every point).
pub type PlanSlot = std::sync::Arc<std::sync::OnceLock<Option<Plan>>>;

/// A learned elimination (R1). The elimination's sequence of operations depends only on the system's structure, not
/// on the point (but for accidental zeros, probability about its size / p, which the reconstruction's final checks
/// catch): learned once, it is replayed on a flat array of slots, one for every (row, column) the elimination ever
/// touches, with no search for pivots, no sets and no merging. The wanted component then comes from the transposed
/// triangular system: z^T U = e_want (one scalar per unknown), u_want = sum_k z_k y_k over the eliminated right-hand
/// sides -- instead of every unknown as a vector over the right-hand sides.
pub struct Plan {
    m: usize,
    slots: usize,
    /// the slot of each program output (A's entries row by row, then B's)
    init: Vec<u32>,
    ops: Vec<PlanOp>,
    pairs: Vec<(u32, u32)>,
    /// for each unknown k >= want: the slot of U_kk, the entries U_jk (j - want, slot) for want <= j < k, and the
    /// right-hand sides (t, slot) of its pivot row
    diag: Vec<u32>,
    upper: Vec<Vec<(u32, u32)>>,
    rhs: Vec<Vec<(u32, u32)>>,
}

enum PlanOp {
    /// the pivot of the next column: its slot (its inverse is used by the eliminations that follow)
    Pivot(u32),
    /// row -= f * pivot row: the slot of f (the row's entry in the pivot column), and the range of (destination,
    /// source) slot pairs
    Elim(u32, u32, u32),
}

impl Plan {
    /// learn from a point's program outputs (the search-based elimination of solve_aug, recorded)
    pub fn learn(p: u64, vals: &[u64], pattern: &(Vec<Vec<u32>>, Vec<Vec<u32>>), n: usize, m: usize, want: usize) -> Option<Plan> {
        use std::collections::{BTreeSet, HashMap};
        let (pa, pb) = pattern;
        let mut slot: HashMap<(usize, usize), u32> = HashMap::new();
        let alloc = |r: usize, c: usize, slot: &mut HashMap<(usize, usize), u32>| -> u32 {
            let next = slot.len() as u32;
            *slot.entry((r, c)).or_insert(next)
        };
        // the rows with every entry the pattern gives (zeros too, so that each output has a slot)
        let mut init = Vec::with_capacity(vals.len());
        let mut rows: Vec<Vec<(usize, u64)>> = Vec::with_capacity(n);
        let nnz: usize = pa.iter().map(Vec::len).sum();
        let (mut at, mut bt) = (0, nnz);
        let mut binit = Vec::new();
        for (r, (ca, cb)) in pa.iter().zip(pb).enumerate() {
            let mut row = Vec::with_capacity(ca.len() + cb.len());
            for (&c, &v) in ca.iter().zip(&vals[at..at + ca.len()]) {
                init.push(alloc(r, c as usize, &mut slot));
                row.push((c as usize, v % p));
            }
            for (&t, &v) in cb.iter().zip(&vals[bt..bt + cb.len()]) {
                binit.push(alloc(r, n + t as usize, &mut slot));
                row.push((n + t as usize, v % p));
            }
            at += ca.len();
            bt += cb.len();
            rows.push(row);
        }
        init.extend(binit);
        let mut col_rows: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); n];
        for (r, row) in rows.iter().enumerate() {
            for &(c, _) in row {
                if c < n {
                    col_rows[c].insert(r);
                }
            }
        }
        let (mut ops, mut pairs) = (Vec::new(), Vec::new());
        let mut pivot = vec![usize::MAX; n];
        let mut done = vec![false; rows.len()];
        for k in 0..n {
            let r = if k < rows.len() && col_rows[k].contains(&k) && !done[k] { k } else { *col_rows[k].iter().find(|&&q| !done[q])? };
            done[r] = true;
            pivot[k] = r;
            for &(c, _) in &rows[r] {
                if c < n {
                    col_rows[c].remove(&r);
                }
            }
            let ark = rows[r].iter().find(|e| e.0 == k).map(|e| e.1)?;
            let inv = invmod(ark, p)?;
            ops.push(PlanOp::Pivot(slot[&(r, k)]));
            let others: Vec<usize> = col_rows[k].iter().copied().filter(|&q| !done[q]).collect();
            for q in others {
                let aqk = rows[q].iter().find(|e| e.0 == k).map(|e| e.1).unwrap_or(0);
                let start = pairs.len() as u32;
                for &(c, _) in rows[r].iter().filter(|e| e.0 != k) {
                    let d = alloc(q, c, &mut slot);
                    pairs.push((d, slot[&(r, c)]));
                }
                ops.push(PlanOp::Elim(slot[&(q, k)], start, pairs.len() as u32));
                let f = mulmod(aqk, inv, p);
                let (rq, rr) = (std::mem::take(&mut rows[q]), &rows[r]);
                let mut merged = Vec::with_capacity(rq.len() + rr.len());
                let (mut i, mut j) = (0, 0);
                while i < rq.len() || j < rr.len() {
                    let (c, v) = if j >= rr.len() || (i < rq.len() && rq[i].0 < rr[j].0) {
                        i += 1;
                        rq[i - 1]
                    } else if i >= rq.len() || rr[j].0 < rq[i].0 {
                        j += 1;
                        (rr[j - 1].0, (p - mulmod(f, rr[j - 1].1, p)) % p)
                    } else {
                        i += 1;
                        j += 1;
                        (rq[i - 1].0, (rq[i - 1].1 + p - mulmod(f, rr[j - 1].1, p)) % p)
                    };
                    if v != 0 || c != k {
                        merged.push((c, v));
                    }
                    if v == 0 && c < n {
                        col_rows[c].remove(&q);
                    }
                }
                for &(c, v) in &merged {
                    if c < n && v != 0 {
                        col_rows[c].insert(q);
                    }
                }
                rows[q] = merged.into_iter().filter(|e| e.1 != 0).collect();
            }
        }
        let w = n - want;
        let mut diag = vec![0u32; w];
        let mut upper: Vec<Vec<(u32, u32)>> = vec![Vec::new(); w];
        let mut rhs: Vec<Vec<(u32, u32)>> = vec![Vec::new(); w];
        for k in want..n {
            let r = pivot[k];
            for &(c, _) in &rows[r] {
                let sl = slot[&(r, c)];
                if c == k {
                    diag[k - want] = sl;
                } else if c < n {
                    upper[c - want].push(((k - want) as u32, sl));
                } else {
                    rhs[k - want].push(((c - n) as u32, sl));
                }
            }
        }
        Some(Plan { m, slots: slot.len(), init, ops, pairs, diag, upper, rhs })
    }

    /// the wanted component for every right-hand side at a point (its program outputs); None if a pivot vanishes
    pub fn replay(&self, p: u64, vals: &[u64]) -> Option<Vec<u64>> {
        debug_assert!(self.diag.iter().all(|&d| (d as usize) < self.slots));
        let mut v = vec![0u64; self.slots];
        for (&s, &x) in self.init.iter().zip(vals) {
            v[s as usize] = x % p;
        }
        let mut inv = 0;
        // each pivot's inverse, kept: the triangular factor's diagonal is the pivots, in column order
        let mut invs = Vec::with_capacity(self.diag.len());
        let first = self.ops.iter().filter(|o| matches!(o, PlanOp::Pivot(_))).count() - self.diag.len();
        let mut k = 0;
        for op in &self.ops {
            match *op {
                PlanOp::Pivot(s) => {
                    inv = invmod(v[s as usize], p)?;
                    if k >= first {
                        invs.push(inv);
                    }
                    k += 1;
                }
                PlanOp::Elim(f, a, b) => {
                    let fac = mulmod(v[f as usize], inv, p);
                    if fac != 0 {
                        let fp = shoup_pre(fac, p);
                        for &(d, s) in &self.pairs[a as usize..b as usize] {
                            let t = shoup_mul(v[s as usize], fac, fp, p);
                            let x = v[d as usize];
                            v[d as usize] = if x >= t { x - t } else { x + p - t };
                        }
                    }
                    v[f as usize] = 0;
                }
            }
        }
        // z^T U = e_want, then u_want = sum_k z_k y_k
        let w = self.diag.len();
        let mut z = vec![0u64; w];
        let mut out = vec![0u64; self.m];
        for c in 0..w {
            let mut acc = if c == 0 { 1 } else { 0 };
            for &(j, s) in &self.upper[c] {
                acc = (acc + p - mulmod(z[j as usize], v[s as usize], p)) % p;
            }
            z[c] = mulmod(acc, invs[c], p);
            if z[c] != 0 {
                let zp = shoup_pre(z[c], p);
                for &(t, s) in &self.rhs[c] {
                    let x = out[t as usize] + shoup_mul(v[s as usize], z[c], zp, p);
                    out[t as usize] = if x >= p { x - p } else { x };
                }
            }
        }
        Some(out)
    }
}

/// A product of linear factors (coefficients per input, constant, power) that eval1 multiplies a linear system's
/// value by: F * D for a guessed denominator D is a polynomial, reconstructed as one (R1). Empty: no scaling.
pub type Scale = std::sync::Arc<Vec<(Vec<i64>, i64, u32)>>;

/// the scale's value at a point, or None where a factor vanishes (the function has a pole there)
fn scale_at(p: u64, scale: &[(Vec<i64>, i64, u32)], pt: &[u64]) -> Option<u64> {
    let m = |c: i64| -> u64 { c.rem_euclid(p as i64) as u64 };
    let mut acc = 1u64;
    for (cs, c0, e) in scale {
        let mut v = m(*c0);
        for (c, x) in cs.iter().zip(pt) {
            v = (v + mulmod(m(*c), *x % p, p)) % p;
        }
        if v == 0 {
            return None;
        }
        acc = mulmod(acc, powmod(v, *e as u64, p), p);
    }
    Some(acc)
}

/// The solutions at the points already solved, shared by the black boxes of one program (one per right-hand side):
/// keyed by (n, m, want, p, point), every right-hand side's non-zero values (None where A is singular). One
/// elimination gives every right-hand side, and the reconstructions of different right-hand sides sample mostly the
/// same points, so each point is eliminated once. The values are those a fresh elimination gives: the cache changes
/// no result.
pub type LinCache = std::sync::Arc<std::sync::Mutex<std::collections::HashMap<LinKey, Solved>>>;

/// (n, m, want, p, point)
pub type LinKey = (usize, usize, usize, u64, Vec<u64>);

/// every right-hand side's non-zero values at a point, None where A is singular
pub type Solved = Option<Vec<(u32, u64)>>;

/// the sparsity patterns of A and B (shared by a system's black boxes)
pub type Patterns = std::sync::Arc<(Vec<Vec<u32>>, Vec<Vec<u32>>)>;

/// the rows of [A | B] from one point's outputs and the patterns, zeros dropped, B's right-hand side t as column
/// n + t; with `only`, that right-hand side alone, as column n
fn aug_rows(p: u64, vals: &[u64], pattern: &(Vec<Vec<u32>>, Vec<Vec<u32>>), n: usize, only: Option<usize>) -> Vec<Vec<(usize, u64)>> {
    let (pa, pb) = pattern;
    let nnz: usize = pa.iter().map(Vec::len).sum();
    let (mut at, mut bt) = (0, nnz);
    pa.iter()
        .zip(pb)
        .map(|(ca, cb)| {
            let mut row: Vec<(usize, u64)> =
                ca.iter().zip(&vals[at..at + ca.len()]).filter_map(|(&c, &v)| (v % p != 0).then_some((c as usize, v % p))).collect();
            row.extend(cb.iter().zip(&vals[bt..bt + cb.len()]).filter_map(|(&t, &v)| match only {
                _ if v % p == 0 => None,
                None => Some((n + t as usize, v % p)),
                Some(o) => (t as usize == o).then_some((n, v % p)),
            }));
            at += ca.len();
            bt += cb.len();
            row
        })
        .collect()
}

/// Component `want` of the solutions of A u = B (A n x n, B n x m, both row-major) mod p, for every right-hand side:
/// sparse Gaussian elimination, column k = 0 up to n-1 (the order of the columns: put the most complex first, as
/// Laporta's reductions do, to keep the fill-in low), its pivot row k if that row has the column, else the first
/// remaining row that does (a deterministic choice); rows indexed by column, so only the rows that hold a column are
/// touched; then back substitution from column n-1 down to `want` (one step for the last column). None if A is
/// singular.
pub fn solve_sparse(p: u64, a: &[u64], b: &[u64], n: usize, m: usize, want: usize) -> Option<Vec<u64>> {
    let rows: Vec<Vec<(usize, u64)>> = (0..n)
        .map(|r| (0..n).filter_map(|c| { let v = a[r * n + c] % p; (v != 0).then_some((c, v)) }).collect())
        .collect();
    solve_rows(p, rows, b, n, m, want)
}

/// solve_sparse/6 on A given as sparse rows (column, value), each sorted by column, values reduced mod p, and B dense
pub fn solve_rows(p: u64, mut rows: Vec<Vec<(usize, u64)>>, b: &[u64], n: usize, m: usize, want: usize) -> Option<Vec<u64>> {
    for (r, row) in rows.iter_mut().enumerate() {
        row.extend((0..m).filter_map(|t| { let v = b[r * m + t] % p; (v != 0).then_some((n + t, v)) }));
    }
    solve_aug(p, rows, n, m, want)
}

/// The augmented form: row r of [A | B] as sparse (column, value) pairs, A's columns 0..n-1 and right-hand side t as
/// column n + t, sorted, values reduced mod p. The elimination treats the right-hand sides as columns it never
/// pivots on, so each costs only its non-zero entries: thousands of them (the masters of a large reduction) are
/// cheap where most rows hold a few. Component `want` of every right-hand side's solution; None if A is singular.
pub fn solve_aug(p: u64, mut rows: Vec<Vec<(usize, u64)>>, n: usize, m: usize, want: usize) -> Option<Vec<u64>> {
    use std::collections::{BTreeMap, BTreeSet};
    let mut col_rows: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); n];
    for (r, row) in rows.iter().enumerate() {
        for &(c, _) in row {
            if c < n {
                col_rows[c].insert(r);
            }
        }
    }
    let mut pivot = vec![usize::MAX; n];
    let mut done = vec![false; rows.len()];
    for k in 0..n {
        let r = if k < rows.len() && col_rows[k].contains(&k) && !done[k] { k } else { *col_rows[k].iter().find(|&&q| !done[q])? };
        done[r] = true;
        pivot[k] = r;
        for &(c, _) in &rows[r] {
            if c < n {
                col_rows[c].remove(&r);
            }
        }
        let ark = rows[r].iter().find(|e| e.0 == k).map(|e| e.1)?;
        let inv = invmod(ark, p)?;
        let others: Vec<usize> = col_rows[k].iter().copied().filter(|&q| !done[q]).collect();
        for q in others {
            let aqk = rows[q].iter().find(|e| e.0 == k).map(|e| e.1).unwrap_or(0);
            if aqk == 0 {
                continue;
            }
            let f = mulmod(aqk, inv, p);
            // row_q -= f row_r, merging the sorted rows
            let (rq, rr) = (std::mem::take(&mut rows[q]), &rows[r]);
            let mut merged = Vec::with_capacity(rq.len() + rr.len());
            let (mut i, mut j) = (0, 0);
            while i < rq.len() || j < rr.len() {
                let (c, v) = if j >= rr.len() || (i < rq.len() && rq[i].0 < rr[j].0) {
                    i += 1;
                    rq[i - 1]
                } else if i >= rq.len() || rr[j].0 < rq[i].0 {
                    j += 1;
                    (rr[j - 1].0, (p - mulmod(f, rr[j - 1].1, p)) % p)
                } else {
                    i += 1;
                    j += 1;
                    (rq[i - 1].0, (rq[i - 1].1 + p - mulmod(f, rr[j - 1].1, p)) % p)
                };
                if v != 0 {
                    merged.push((c, v));
                } else if c < n {
                    col_rows[c].remove(&q);
                }
            }
            for &(c, _) in &merged {
                if c < n {
                    col_rows[c].insert(q);
                }
            }
            rows[q] = merged;
        }
    }
    // back substitution, each u_k a sparse vector over the right-hand sides: the pivot row of column k holds only
    // columns >= k (with the wanted component last, one step)
    let mut u: Vec<BTreeMap<usize, u64>> = vec![BTreeMap::new(); n];
    for k in (want..n).rev() {
        let r = pivot[k];
        let inv = invmod(rows[r].iter().find(|e| e.0 == k).map(|e| e.1)?, p)?;
        let mut acc: BTreeMap<usize, u64> = BTreeMap::new();
        for &(c, v) in &rows[r] {
            if c >= n {
                *acc.entry(c - n).or_insert(0) = v;
            }
        }
        for &(c, v) in &rows[r] {
            if c > k && c < n {
                for (&t, &x) in &u[c] {
                    let e = acc.entry(t).or_insert(0);
                    *e = (*e + p - mulmod(v, x, p)) % p;
                }
            }
        }
        u[k] = acc.into_iter().filter(|e| e.1 != 0).map(|(t, x)| (t, mulmod(x, inv, p))).collect();
    }
    let mut out = vec![0u64; m];
    for (&t, &x) in &u[want] {
        out[t] = x;
    }
    Some(out)
}

/// The pivots of a sparse m x n matrix mod p (rows of (column, value), sorted by column): Gaussian elimination
/// column by column, 0 to n-1, each column's pivot the shortest remaining row that holds it, the lowest of equals (a
/// deterministic choice that keeps the fill-in low);
/// (pivot rows, pivot columns) in column order, and the rows that `targets` need. The columns without a pivot are the
/// dependent ones: in a Laporta reduction, columns ordered most complex first, they are the master integrals.
///
/// Mark and sweep (as in FiniteFlow): a pivot row is final once chosen, and every row records the pivot rows it was
/// combined with. A target column needs its pivot row; a needed row needs the rows it was combined with (its forward
/// elimination) and the pivot row of every other non-master column left in it (back substitution). The other rows
/// never touch the targets: the needed rows, in the pivot order of their columns, solve for them alone, by the same
/// elimination. Sorted by pivot column; empty without targets.
pub fn pivots(p: u64, mut rows: Vec<Vec<(usize, u64)>>, n: usize, targets: &[usize]) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    use std::collections::BTreeSet;
    let mut col_rows: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); n];
    for (r, row) in rows.iter_mut().enumerate() {
        row.retain(|e| e.1 % p != 0);
        for e in row.iter_mut() {
            e.1 %= p;
            col_rows[e.0].insert(r);
        }
    }
    let (mut prow, mut pcol) = (Vec::new(), Vec::new());
    let mut deps: Vec<Vec<u32>> = vec![Vec::new(); rows.len()];
    let mut pivot_of = vec![usize::MAX; n];
    for k in 0..n {
        // the shortest row holding the column (the lowest index among equals): the least fill-in
        let Some(r) = col_rows[k].iter().copied().min_by_key(|&q| (rows[q].len(), q)) else { continue };
        prow.push(r);
        pcol.push(k);
        pivot_of[k] = r;
        for &(c, _) in &rows[r] {
            col_rows[c].remove(&r);
        }
        let inv = match rows[r].iter().find(|e| e.0 == k).and_then(|e| invmod(e.1, p)) {
            Some(i) => i,
            None => continue,
        };
        let others: Vec<usize> = col_rows[k].iter().copied().collect();
        for q in others {
            let aqk = rows[q].iter().find(|e| e.0 == k).map_or(0, |e| e.1);
            let f = mulmod(aqk, inv, p);
            deps[q].push(r as u32);
            let (rq, rr) = (std::mem::take(&mut rows[q]), &rows[r]);
            let mut merged = Vec::with_capacity(rq.len() + rr.len());
            let (mut i, mut j) = (0, 0);
            while i < rq.len() || j < rr.len() {
                let (c, v) = if j >= rr.len() || (i < rq.len() && rq[i].0 < rr[j].0) {
                    i += 1;
                    rq[i - 1]
                } else if i >= rq.len() || rr[j].0 < rq[i].0 {
                    j += 1;
                    (rr[j - 1].0, (p - mulmod(f, rr[j - 1].1, p)) % p)
                } else {
                    i += 1;
                    j += 1;
                    (rq[i - 1].0, (rq[i - 1].1 + p - mulmod(f, rr[j - 1].1, p)) % p)
                };
                if v != 0 {
                    merged.push((c, v));
                } else {
                    col_rows[c].remove(&q);
                }
            }
            for &(c, _) in &merged {
                col_rows[c].insert(q);
            }
            rows[q] = merged;
        }
    }
    // the sweep: the rows the targets need
    let mut needed = vec![false; rows.len()];
    let mut stack: Vec<usize> = targets.iter().filter(|&&t| t < n && pivot_of[t] != usize::MAX).map(|&t| pivot_of[t]).collect();
    while let Some(r) = stack.pop() {
        if needed[r] {
            continue;
        }
        needed[r] = true;
        stack.extend(deps[r].iter().map(|&d| d as usize).filter(|&d| !needed[d]));
        stack.extend(rows[r].iter().map(|e| pivot_of[e.0]).filter(|&q| q != usize::MAX && !needed[q]));
    }
    let need: Vec<usize> = pcol.iter().map(|&c| pivot_of[c]).filter(|&r| needed[r]).collect();
    (prow, pcol, need)
}

impl BlackBox {
    /// for a linear system: per point, component `want` of the solution for every right-hand side (None where A
    /// is singular or an entry is undefined)
    pub fn eval_all(&self, p: u64, points: &[u64]) -> Result<Vec<Option<Vec<u64>>>, FfError> {
        match self {
            BlackBox::Prog(_) => Err(FfError::Shape),
            BlackBox::Lin { m, .. } => Ok(self
                .lin_sparse(p, points)?
                .into_iter()
                .map(|r| {
                    r.map(|sparse| {
                        let mut u = vec![0u64; *m];
                        for (t, v) in sparse {
                            u[t as usize] = v;
                        }
                        u
                    })
                })
                .collect()),
        }
    }

    /// a linear system's solution at each point, every right-hand side's non-zero values, through the cache
    fn lin_sparse(&self, p: u64, points: &[u64]) -> Result<Vec<Solved>, FfError> {
        let BlackBox::Lin { prog, n, m, want, pattern, cache, plan, .. } = self else { return Err(FfError::Shape) };
        let k = prog.ninputs.max(1);
        let key = |pt: &[u64]| (*n, *m, *want, p, pt.to_vec());
        let mut out: Vec<Option<Solved>> = {
            let c = cache.lock().map_err(|_| FfError::Shape)?;
            points.chunks(k).map(|pt| c.get(&key(pt)).cloned()).collect()
        };
        let miss: Vec<usize> = out.iter().enumerate().filter(|(_, o)| o.is_none()).map(|(i, _)| i).collect();
        if !miss.is_empty() {
            let pts: Vec<u64> = miss.iter().flat_map(|&i| points[i * k..(i + 1) * k].iter().copied()).collect();
            let vals = prog.eval(p, &pts)?;
            let solved: Vec<Solved> = vals
                .chunks(prog.outputs.len())
                .map(|row| {
                    if row.contains(&BAD) {
                        return None;
                    }
                    // the learned elimination, learned at the first point (the search-based one if that failed)
                    let learned = plan.get_or_init(|| Plan::learn(p, row, pattern, *n, *m, *want));
                    match learned {
                        Some(pl) => pl.replay(p, row),
                        None => solve_aug(p, aug_rows(p, row, pattern, *n, None), *n, *m, *want),
                    }
                    .map(|u| u.into_iter().enumerate().filter(|e| e.1 != 0).map(|(t, v)| (t as u32, v)).collect())
                })
                .collect();
            let mut c = cache.lock().map_err(|_| FfError::Shape)?;
            for (&i, r) in miss.iter().zip(solved) {
                c.insert(key(&points[i * k..(i + 1) * k]), r.clone());
                out[i] = Some(r);
            }
        }
        Ok(out.into_iter().map(|o| o.flatten()).collect())
    }

    pub fn ninputs(&self) -> usize {
        match self {
            BlackBox::Prog(p) => p.ninputs,
            BlackBox::Lin { prog, .. } => prog.ninputs,
        }
    }

    /// one value per point (`points`: ninputs residues per point), BAD where the black box is undefined
    pub fn eval1(&self, p: u64, points: &[u64]) -> Result<Vec<u64>, FfError> {
        match self {
            BlackBox::Prog(prog) => {
                let vals = prog.eval(p, points)?;
                Ok(vals.iter().step_by(prog.outputs.len().max(1)).copied().collect())
            }
            // one right-hand side: every one is solved (and cached for the others), this one kept, times the scale
            BlackBox::Lin { col, scale, prog, .. } => {
                let k = prog.ninputs.max(1);
                Ok(self
                    .lin_sparse(p, points)?
                    .into_iter()
                    .zip(points.chunks(k))
                    .map(|(r, pt)| match (r, scale_at(p, scale, pt)) {
                        (None, _) | (_, None) => BAD,
                        (Some(sparse), Some(f)) => {
                            mulmod(sparse.binary_search_by_key(&(*col as u32), |e| e.0).map_or(0, |i| sparse[i].1), f, p)
                        }
                    })
                    .collect())
            }
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
pub enum ColumnError {
    /// a point where the program divides by zero
    Bad,
    /// two monomials take the same value at r (or one is 0): another r is needed
    Singular,
    BadPrime,
    Shape,
}

impl BlackBox {
    /// One column of Zippel's interpolation, mod p: the coefficients c_m of the support's n monomials m (exponent
    /// rows of r.len() each; with no variable yet, one empty monomial) in the program's first output at the points (r^i, suffix), i = 1 .. |support| --
    /// sum_m c_m m(r)^i -- by a transposed Vandermonde system on the nodes m(r), divided by the nodes.
    pub fn column(&self, p: u64, r: &[u64], n: usize, support: &[u32], suffix: &[u64]) -> Result<Vec<u64>, ColumnError> {
        let k = r.len();
        if k + suffix.len() != self.ninputs() || support.len() != n * k {
            return Err(ColumnError::Shape);
        }
        // the points r^1 .. r^n, each followed by the suffix
        let mut points = Vec::with_capacity(n * self.ninputs());
        let mut cur: Vec<u64> = r.iter().map(|&x| x % p).collect();
        for _ in 0..n {
            points.extend_from_slice(&cur);
            points.extend(suffix.iter().map(|&x| x % p));
            for (c, &x) in cur.iter_mut().zip(r) {
                *c = mulmod(*c, x, p);
            }
        }
        let first = self.eval1(p, &points).map_err(|e| match e {
            FfError::BadPrime => ColumnError::BadPrime,
            FfError::Shape => ColumnError::Shape,
        })?;
        if first.contains(&BAD) {
            return Err(ColumnError::Bad);
        }
        let nodes: Vec<u64> = (0..n)
            .map(|m| (0..k).fold(1u64, |a, l| mulmod(a, powmod(r[l], support[m * k + l] as u64, p), p)))
            .collect();
        if nodes.contains(&0) {
            return Err(ColumnError::Singular);
        }
        let cv = tvandermonde(p, &nodes, &first).ok_or(ColumnError::Singular)?;
        Ok(cv.iter().zip(&nodes).map(|(&c, &v)| mulmod(c, invmod(v, p).unwrap_or(0), p)).collect())
    }
}

/// The solution x of the n x n system a x = b mod p (a row-major), or None if it is singular.
pub fn solve_mod(p: u64, mut a: Vec<u64>, mut b: Vec<u64>, n: usize) -> Option<Vec<u64>> {
    for col in 0..n {
        let piv = (col..n).find(|&r| !a[r * n + col].is_multiple_of(p))?;
        if piv != col {
            for c in 0..n {
                a.swap(piv * n + c, col * n + c);
            }
            b.swap(piv, col);
        }
        let inv = invmod(a[col * n + col], p)?;
        for c in col..n {
            a[col * n + c] = mulmod(a[col * n + c], inv, p);
        }
        b[col] = mulmod(b[col], inv, p);
        for r in 0..n {
            if r != col && a[r * n + col] != 0 {
                let f = a[r * n + col];
                for c in col..n {
                    a[r * n + c] = (a[r * n + c] + p - mulmod(f, a[col * n + c], p)) % p;
                }
                b[r] = (b[r] + p - mulmod(f, b[col], p)) % p;
            }
        }
    }
    Some(b)
}

impl BlackBox {
    /// The coefficients of the univariate rational function g(t) = f(shift + t (1, zz)) = n(t) / d(t), with
    /// deg n = dn, deg d = dd and d(0) = 1, from its values at the dn + dd + 1 points ts: [n_0 .. n_dn, d_1 .. d_dd].
    /// Err(Bad) if f is undefined at one of the points, or the system is singular there.
    pub fn line(&self, p: u64, shift: &[u64], zz: &[u64], ts: &[u64], dn: usize, dd: usize) -> Result<Vec<u64>, ColumnError> {
        let n = shift.len();
        if zz.len() + 1 != n || n != self.ninputs() || ts.len() != dn + dd + 1 {
            return Err(ColumnError::Shape);
        }
        let mut points = Vec::with_capacity(ts.len() * n);
        for &t in ts {
            points.push((shift[0] + t) % p);
            for (i, &z) in zz.iter().enumerate() {
                points.push((shift[i + 1] + mulmod(t, z, p)) % p);
            }
        }
        let g = self.eval1(p, &points).map_err(|e| match e {
            FfError::BadPrime => ColumnError::BadPrime,
            FfError::Shape => ColumnError::Shape,
        })?;
        if g.contains(&BAD) {
            return Err(ColumnError::Bad);
        }
        // sum_r n_r t^r - g sum_{r >= 1} d_r t^r = g
        let m = dn + dd + 1;
        let mut a = vec![0u64; m * m];
        for (j, (&t, &gj)) in ts.iter().zip(&g).enumerate() {
            let mut tp = 1u64;
            for r in 0..=dn.max(dd) {
                if r <= dn {
                    a[j * m + r] = tp;
                }
                if r >= 1 && r <= dd {
                    a[j * m + dn + r] = (p - mulmod(gj, tp, p)) % p;
                }
                tp = mulmod(tp, t, p);
            }
        }
        solve_mod(p, a, g, m).ok_or(ColumnError::Bad)
    }

    /// One column of the multi-output Zippel interpolation of the coefficients of g (`line`) as polynomials in
    /// z_2 .. z_n: for the points zz = (r^i, suffix), i = 1 .. n_mono, the coefficients of the support's monomials
    /// in each of the dn + dd + 1 coefficient functions -- dn + dd + 1 rows of n_mono values.
    #[allow(clippy::too_many_arguments)]
    pub fn line_column(&self, p: u64, shift: &[u64], r: &[u64], n_mono: usize, support: &[u32], suffix: &[u64],
                       ts: &[u64], dn: usize, dd: usize) -> Result<Vec<Vec<u64>>, ColumnError> {
        let k = r.len();
        if support.len() != n_mono * k {
            return Err(ColumnError::Shape);
        }
        let outs = dn + dd + 1;
        let mut values = vec![Vec::with_capacity(n_mono); outs];
        let mut cur: Vec<u64> = r.iter().map(|&x| x % p).collect();
        for _ in 0..n_mono {
            let mut zz = cur.clone();
            zz.extend(suffix.iter().map(|&x| x % p));
            let c = self.line(p, shift, &zz, ts, dn, dd)?;
            for (o, v) in c.into_iter().enumerate() {
                values[o].push(v);
            }
            for (c, &x) in cur.iter_mut().zip(r) {
                *c = mulmod(*c, x, p);
            }
        }
        let nodes: Vec<u64> = (0..n_mono)
            .map(|m| (0..k).fold(1u64, |a, l| mulmod(a, powmod(r[l], support[m * k + l] as u64, p), p)))
            .collect();
        if nodes.contains(&0) {
            return Err(ColumnError::Singular);
        }
        values
            .iter()
            .map(|vals| {
                let cv = tvandermonde(p, &nodes, vals).ok_or(ColumnError::Singular)?;
                Ok(cv.iter().zip(&nodes).map(|(&c, &v)| mulmod(c, invmod(v, p).unwrap_or(0), p)).collect())
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const P: u64 = 9_223_372_036_854_775_783; // the largest prime below 2^63

    #[test]
    fn evaluates_a_rational_function() {
        // f(x, y) = (x - 2/3 y)^2 / (x + y) - 1 ... as ops
        let prog = Program {
            ops: vec![
                Op::Input(0),
                Op::Input(1),
                Op::Const(-1, vec![2], vec![3]),
                Op::Mul(2, 1),
                Op::Add(0, 3),
                Op::Pow(4, 2),
                Op::Add(0, 1),
                Op::Pow(6, -1),
                Op::Mul(5, 7),
                Op::Const(1, vec![1], vec![1]),
                Op::Sub(8, 9),
            ],
            outputs: vec![10],
            ninputs: 2,
        };
        prog.check().unwrap();
        let p = 101u64;
        let vals = prog.eval(p, &[5, 7, 3, 98]).unwrap(); // the second point has x + y = 101 = 0 mod p
        let inv = |a: u64| invmod(a, p).unwrap();
        let x = 5u64;
        let y = 7u64;
        let t = (x + p - mulmod(mulmod(2, inv(3), p), y, p)) % p;
        let want = (mulmod(mulmod(t, t, p), inv(x + y), p) + p - 1) % p;
        assert_eq!(vals, vec![want, BAD]);
        // a large prime works, and a constant's denominator 0 mod p is refused
        assert!(prog.eval(P, &[1, 2]).is_ok());
        // a constant of several limbs: 2^64 + 5 mod 101
        let big = Program { ops: vec![Op::Const(1, vec![5, 1], vec![1])], outputs: vec![0], ninputs: 0 };
        assert_eq!(big.eval(p, &[]).unwrap(), vec![(((1u128 << 64) + 5) % 101) as u64]);
        assert_eq!(prog.eval(3, &[1, 2]).err(), Some(FfError::BadPrime));
        assert_eq!(Program { ops: vec![Op::Mul(0, 1)], outputs: vec![0], ninputs: 0 }.check().err(), Some(FfError::Shape));
    }

    #[test]
    fn solves_transposed_vandermonde_systems() {
        let p = P;
        let mut s = 0x1234_5678_9abc_def1u64;
        let mut next = || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            s % p
        };
        for n in [1usize, 2, 7, 40] {
            let nodes: Vec<u64> = (0..n).map(|_| next()).collect();
            let c: Vec<u64> = (0..n).map(|_| next()).collect();
            let values: Vec<u64> = (0..n)
                .map(|i| (0..n).fold(0u64, |a, j| (a + mulmod(c[j], powmod(nodes[j], i as u64, p), p)) % p))
                .collect();
            assert_eq!(tvandermonde(p, &nodes, &values), Some(c));
        }
        assert_eq!(tvandermonde(p, &[5, 5], &[1, 2]), None);
    }

    #[test]
    fn finds_a_rational_function_on_a_line() {
        // f(x, y) = (x + 2 y) / (x - y + 3), shift (1, 2), direction (1, 5): g(t) = (5 + 11 t) / (2 - 4 t),
        // normalised d(0) = 1: n = 5/2 + 11/2 t, d = 1 - 2 t
        let p = 1_000_003u64;
        let prog = Program {
            ops: vec![
                Op::Input(0),
                Op::Input(1),
                Op::Const(1, vec![2], vec![1]),
                Op::Mul(2, 1),
                Op::Add(0, 3),
                Op::Sub(0, 1),
                Op::Const(1, vec![3], vec![1]),
                Op::Add(5, 6),
                Op::Pow(7, -1),
                Op::Mul(4, 8),
            ],
            outputs: vec![9],
            ninputs: 2,
        };
        let c = BlackBox::Prog(prog).line(p, &[1, 2], &[5], &[7, 11, 13], 1, 1).unwrap();
        let half = invmod(2, p).unwrap();
        assert_eq!(c, vec![mulmod(5, half, p), mulmod(11, half, p), p - 2]);
    }

    #[test]
    fn solves_a_linear_system_per_point() {
        // [[x, 1], [1, y]] u = [1, 0]: u_0 = y / (x y - 1), u_1 = -1 / (x y - 1)
        let prog = Program {
            ops: vec![Op::Input(0), Op::Const(1, vec![1], vec![1]), Op::Input(1), Op::Const(1, vec![0], vec![1])],
            outputs: vec![0, 1, 1, 2, 1, 3],
            ninputs: 2,
        };
        let p = 101u64;
        let dense = std::sync::Arc::new((vec![vec![0, 1], vec![0, 1]], vec![vec![0], vec![0]]));
        let cache = LinCache::default();
        let bb0 = BlackBox::Lin { prog: prog.clone(), n: 2, m: 1, want: 0, col: 0, pattern: dense.clone(), cache: cache.clone(), scale: Scale::default(), plan: PlanSlot::default() };
        let bb1 = BlackBox::Lin { prog, n: 2, m: 1, want: 1, col: 0, pattern: dense, cache, scale: Scale::default(), plan: PlanSlot::default() };
        let inv = |a: u64| invmod(a, p).unwrap();
        // at (3, 5): det 14; and (1, 1): singular
        assert_eq!(bb0.eval1(p, &[3, 5, 1, 1]).unwrap(), vec![mulmod(5, inv(14), p), BAD]);
        assert_eq!(bb1.eval1(p, &[3, 5]).unwrap(), vec![(p - inv(14)) % p]);
    }

    #[test]
    fn solves_a_system_given_by_its_pattern() {
        // [[x, 1], [0, y]] u = [1, 0], the 0 not an output: u_0 = 1 / x, u_1 = 0
        let prog = Program {
            ops: vec![Op::Input(0), Op::Const(1, vec![1], vec![1]), Op::Input(1), Op::Const(1, vec![0], vec![1])],
            outputs: vec![0, 1, 2, 1, 3],
            ninputs: 2,
        };
        let p = 101u64;
        let pattern = std::sync::Arc::new((vec![vec![0, 1], vec![1]], vec![vec![0], vec![0]]));
        let bb = BlackBox::Lin { prog, n: 2, m: 1, want: 0, col: 0, pattern, cache: LinCache::default(), scale: Scale::default(), plan: PlanSlot::default() };
        assert_eq!(bb.eval1(p, &[3, 5, 0, 5]).unwrap(), vec![invmod(3, p).unwrap(), BAD]);
        assert_eq!(bb.eval_all(p, &[3, 5]).unwrap(), vec![Some(vec![invmod(3, p).unwrap()])]);
    }

    #[test]
    fn scales_a_linear_system() {
        // [[x, 1], [0, y]] u = [1, 0]: u_0 = 1 / x; times (x + 0 y + 0)^1 it is 1, and bad where x = 0
        let prog = Program {
            ops: vec![Op::Input(0), Op::Const(1, vec![1], vec![1]), Op::Input(1), Op::Const(1, vec![0], vec![1])],
            outputs: vec![0, 1, 2, 1, 3],
            ninputs: 2,
        };
        let pattern = std::sync::Arc::new((vec![vec![0, 1], vec![1]], vec![vec![0], vec![0]]));
        let scale = std::sync::Arc::new(vec![(vec![1, 0], 0, 1)]);
        let bb = BlackBox::Lin { prog, n: 2, m: 1, want: 0, col: 0, pattern, cache: LinCache::default(), scale, plan: PlanSlot::default() };
        assert_eq!(bb.eval1(101, &[3, 5, 7, 2, 0, 5]).unwrap(), vec![1, 1, BAD]);
    }

    #[test]
    fn finds_pivots() {
        // rows x+y, 2x+2y, y, z over columns x, y, z, w: x pivots on row 0 (the first of two equally short), y on
        // row 2 (row 1 is then 0), z on row 3; w has none
        let rows = vec![vec![(0, 1), (1, 1)], vec![(0, 2), (1, 2)], vec![(1, 1)], vec![(2, 5)]];
        assert_eq!(pivots(P, rows.clone(), 4, &[]), (vec![0, 2, 3], vec![0, 1, 2], vec![]));
        // y needs only its own row; x needs row 0 and, for the y left in it, row 2; z needs row 3
        assert_eq!(pivots(P, rows.clone(), 4, &[1]).2, vec![2]);
        assert_eq!(pivots(P, rows.clone(), 4, &[0]).2, vec![0, 2]);
        assert_eq!(pivots(P, rows, 4, &[2]).2, vec![3]);
        // x in a long row and a short one: the short one pivots
        let rows = vec![vec![(0, 1), (1, 1), (2, 1)], vec![(0, 1), (2, 1)]];
        // y's row was combined with x's pivot row (row 1): y needs both
        assert_eq!(pivots(P, rows, 3, &[1]), (vec![1, 0], vec![0, 1], vec![1, 0]));
        assert_eq!(pivots(P, vec![vec![(1, P)]], 2, &[0]), (vec![], vec![], vec![]));
    }

    #[test]
    fn inverts() {
        for p in [3u64, 101, 1_000_000_007, P, (1u64 << 63) - 25] {
            for a in [1u64, 2, p - 1, p / 2 + 1, 0x1234_5678 % p] {
                if a % p != 0 {
                    assert_eq!(mulmod(invmod(a, p).unwrap(), a, p), 1, "p {p} a {a}");
                }
            }
            assert_eq!(invmod(0, p), None);
            assert_eq!(invmod(p, p), None);
        }
    }

    #[test]
    fn shoup_multiplication_agrees() {
        let mut s = 0x2545_f491_4f6c_dd1du64;
        for p in [P, 3, 1_000_000_007, (1u64 << 63) - 25] {
            for _ in 0..2000 {
                s ^= s << 13;
                s ^= s >> 7;
                s ^= s << 17;
                let (w, x) = (s % p, s.rotate_left(17) % p);
                assert_eq!(shoup_mul(x, w, shoup_pre(w, p), p), mulmod(w, x, p));
            }
        }
    }

    #[test]
    fn replays_the_learned_elimination() {
        let p = P;
        let mut s = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            s
        };
        for n in [1usize, 2, 5, 12, 30, 60] {
            let m = 4usize;
            // a fixed sparsity pattern: the diagonal and about a quarter of the rest; B about half full
            let pa: Vec<Vec<u32>> = (0..n).map(|r| (0..n).filter(|&c| c == r || next() % 4 == 0).map(|c| c as u32).collect()).collect();
            let pb: Vec<Vec<u32>> = (0..n).map(|_| (0..m).filter(|_| next() % 2 == 0).map(|t| t as u32).collect()).collect();
            let pattern = (pa, pb);
            let nnz: usize = pattern.0.iter().chain(&pattern.1).map(Vec::len).sum();
            for want in [0, n / 2, n - 1] {
                let learn_vals: Vec<u64> = (0..nnz).map(|_| next() % p).collect();
                let plan = Plan::learn(p, &learn_vals, &pattern, n, m, want).unwrap();
                for _ in 0..3 {
                    let vals: Vec<u64> = (0..nnz).map(|_| next() % p).collect();
                    let expect = solve_aug(p, aug_rows(p, &vals, &pattern, n, None), n, m, want);
                    assert_eq!(plan.replay(p, &vals), expect, "n {n} want {want}");
                }
            }
        }
    }

    #[test]
    fn sparse_elimination_agrees_with_dense() {
        let p = P;
        let mut s = 0x0123_4567_89ab_cdefu64;
        let mut next = || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            s
        };
        for n in [1usize, 2, 5, 12, 30] {
            let (m, want) = (3usize, (n - 1) * (n % 2));
            // sparse: the diagonal and about a quarter of the rest
            let a: Vec<u64> = (0..n * n).map(|i| if i % (n + 1) == 0 || next() % 4 == 0 { next() % p } else { 0 }).collect();
            let b: Vec<u64> = (0..n * m).map(|_| next() % p).collect();
            let got = solve_sparse(p, &a, &b, n, m, want).unwrap();
            for (t, &g) in got.iter().enumerate() {
                let col: Vec<u64> = (0..n).map(|r| b[r * m + t]).collect();
                assert_eq!(g, solve_mod(p, a.clone(), col, n).unwrap()[want], "n {n} rhs {t}");
            }
        }
        assert_eq!(solve_sparse(P, &[1, 2, 2, 4], &[1, 1], 2, 1, 1), None);
    }

    #[test]
    fn computes_a_column() {
        // f(x, y, z) = 3 x^2 y + 5 y z - 7 at z = 11: the coefficients of x^2 y and y in x, y are 3 and 5 * 11 (the
        // constant -7 is the third monomial)
        let p = P;
        let prog = Program {
            ops: vec![
                Op::Input(0),
                Op::Input(1),
                Op::Input(2),
                Op::Pow(0, 2),
                Op::Const(1, vec![3], vec![1]),
                Op::Mul(4, 3),
                Op::Mul(5, 1),
                Op::Const(1, vec![5], vec![1]),
                Op::Mul(7, 1),
                Op::Mul(8, 2),
                Op::Add(6, 9),
                Op::Const(1, vec![7], vec![1]),
                Op::Sub(10, 11),
            ],
            outputs: vec![12],
            ninputs: 3,
        };
        let support = [2u32, 1, 0, 1, 0, 0]; // x^2 y, y, 1
        let c = BlackBox::Prog(prog).column(p, &[123_456_789, 987_654_321], 3, &support, &[11]).unwrap();
        assert_eq!(c, vec![3, 55, p - 7]);
    }
}
