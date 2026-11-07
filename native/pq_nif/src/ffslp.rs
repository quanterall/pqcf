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
pub fn invmod(a: u64, p: u64) -> Option<u64> {
    if a.is_multiple_of(p) { None } else { Some(powmod(a, p - 2, p)) }
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
/// entries are a program's outputs (phase 5): A u = B with A sparse, row r holding the columns `pattern[r]`
/// (increasing), its entries the first outputs row by row, and B the next n*m (n rows of m right-hand sides),
/// solved by sparse elimination mod p at each point; the value is component `want` of the solution for right-hand
/// side `col` (bad where A is singular).
pub enum BlackBox {
    Prog(Program),
    Lin { prog: Program, n: usize, m: usize, want: usize, col: usize, pattern: std::sync::Arc<Vec<Vec<u32>>> },
}

/// the rows of A from one point's outputs and the pattern, zeros dropped
fn pattern_rows(p: u64, vals: &[u64], pattern: &[Vec<u32>]) -> Vec<Vec<(usize, u64)>> {
    let mut at = 0;
    pattern
        .iter()
        .map(|cols| {
            let row = cols.iter().zip(&vals[at..at + cols.len()]).filter_map(|(&c, &v)| (v % p != 0).then_some((c as usize, v % p))).collect();
            at += cols.len();
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

/// solve_sparse/6 on A given as sparse rows (column, value), each sorted by column, values reduced mod p
pub fn solve_rows(p: u64, mut rows: Vec<Vec<(usize, u64)>>, b: &[u64], n: usize, m: usize, want: usize) -> Option<Vec<u64>> {
    use std::collections::BTreeSet;
    let mut rhs: Vec<Vec<u64>> = (0..n).map(|r| b[r * m..(r + 1) * m].iter().map(|&v| v % p).collect()).collect();
    let mut col_rows: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); n];
    for (r, row) in rows.iter().enumerate() {
        for &(c, _) in row {
            col_rows[c].insert(r);
        }
    }
    let mut pivot = vec![usize::MAX; n];
    let mut done = vec![false; n];
    for k in 0..n {
        let r = if col_rows[k].contains(&k) && !done[k] { k } else { *col_rows[k].iter().find(|&&q| !done[q])? };
        done[r] = true;
        pivot[k] = r;
        for &(c, _) in &rows[r] {
            col_rows[c].remove(&r);
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
                } else {
                    col_rows[c].remove(&q);
                }
            }
            for &(c, _) in &merged {
                col_rows[c].insert(q);
            }
            rows[q] = merged;
            for t in 0..m {
                rhs[q][t] = (rhs[q][t] + p - mulmod(f, rhs[r][t], p)) % p;
            }
        }
    }
    // back substitution: the pivot row of column k holds only columns >= k
    let mut u = vec![vec![0u64; m]; n];
    for k in (want..n).rev() {
        let r = pivot[k];
        let inv = invmod(rows[r].iter().find(|e| e.0 == k).map(|e| e.1)?, p)?;
        for t in 0..m {
            let mut acc = rhs[r][t];
            for &(c, v) in &rows[r] {
                if c > k {
                    acc = (acc + p - mulmod(v, u[c][t], p)) % p;
                }
            }
            u[k][t] = mulmod(acc, inv, p);
        }
    }
    Some(u[want].clone())
}

impl BlackBox {
    /// for a linear system: per point, component `want` of the solution for every right-hand side (None where A
    /// is singular or an entry is undefined)
    pub fn eval_all(&self, p: u64, points: &[u64]) -> Result<Vec<Option<Vec<u64>>>, FfError> {
        match self {
            BlackBox::Prog(_) => Err(FfError::Shape),
            BlackBox::Lin { prog, n, m, want, pattern, .. } => {
                let (n, m) = (*n, *m);
                let nnz = prog.outputs.len() - n * m;
                let vals = prog.eval(p, points)?;
                Ok(vals
                    .chunks(prog.outputs.len())
                    .map(|row| {
                        if row.contains(&BAD) {
                            None
                        } else {
                            solve_rows(p, pattern_rows(p, &row[..nnz], pattern), &row[nnz..], n, m, *want)
                        }
                    })
                    .collect())
            }
        }
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
            // one right-hand side: only its column of B goes through the elimination
            BlackBox::Lin { prog, n, m, want, col, pattern } => {
                let (n, m) = (*n, *m);
                let nnz = prog.outputs.len() - n * m;
                let vals = prog.eval(p, points)?;
                Ok(vals
                    .chunks(prog.outputs.len())
                    .map(|row| {
                        if row.contains(&BAD) {
                            return BAD;
                        }
                        let b: Vec<u64> = (0..n).map(|r| row[nnz + r * m + col]).collect();
                        solve_rows(p, pattern_rows(p, &row[..nnz], pattern), &b, n, 1, *want).map_or(BAD, |u| u[0])
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
        let dense = std::sync::Arc::new(vec![vec![0, 1], vec![0, 1]]);
        let bb0 = BlackBox::Lin { prog: prog.clone(), n: 2, m: 1, want: 0, col: 0, pattern: dense.clone() };
        let bb1 = BlackBox::Lin { prog, n: 2, m: 1, want: 1, col: 0, pattern: dense };
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
        let pattern = std::sync::Arc::new(vec![vec![0, 1], vec![1]]);
        let bb = BlackBox::Lin { prog, n: 2, m: 1, want: 0, col: 0, pattern };
        assert_eq!(bb.eval1(p, &[3, 5, 0, 5]).unwrap(), vec![invmod(3, p).unwrap(), BAD]);
        assert_eq!(bb.eval_all(p, &[3, 5]).unwrap(), vec![Some(vec![invmod(3, p).unwrap()])]);
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
