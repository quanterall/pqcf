//! Straight-line programs evaluated over a prime field (decision record R1): the black box of the reconstruction.
//! A program is prepared once (its operations, exact rational constants), then evaluated modulo a prime p < 2^63
//! at a batch of points. A point where the program divides by zero mod p is reported as bad (every output
//! `BAD`), for the reconstruction to skip; a constant whose denominator is 0 mod p makes the prime unusable.

/// The value of every output at a bad point (never a residue, as p < 2^63).
pub const BAD: u64 = u64::MAX;

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

#[derive(Debug, PartialEq, Eq)]
pub enum ColumnError {
    /// a point where the program divides by zero
    Bad,
    /// two monomials take the same value at r (or one is 0): another r is needed
    Singular,
    BadPrime,
    Shape,
}

impl Program {
    /// One column of Zippel's interpolation, mod p: the coefficients c_m of the support's n monomials m (exponent
    /// rows of r.len() each; with no variable yet, one empty monomial) in the program's first output at the points (r^i, suffix), i = 1 .. |support| --
    /// sum_m c_m m(r)^i -- by a transposed Vandermonde system on the nodes m(r), divided by the nodes.
    pub fn column(&self, p: u64, r: &[u64], n: usize, support: &[u32], suffix: &[u64]) -> Result<Vec<u64>, ColumnError> {
        let k = r.len();
        if k + suffix.len() != self.ninputs || support.len() != n * k {
            return Err(ColumnError::Shape);
        }
        // the points r^1 .. r^n, each followed by the suffix
        let mut points = Vec::with_capacity(n * self.ninputs);
        let mut cur: Vec<u64> = r.iter().map(|&x| x % p).collect();
        for _ in 0..n {
            points.extend_from_slice(&cur);
            points.extend(suffix.iter().map(|&x| x % p));
            for (c, &x) in cur.iter_mut().zip(r) {
                *c = mulmod(*c, x, p);
            }
        }
        let vals = self.eval(p, &points).map_err(|e| match e {
            FfError::BadPrime => ColumnError::BadPrime,
            FfError::Shape => ColumnError::Shape,
        })?;
        let width = self.outputs.len().max(1);
        let first: Vec<u64> = vals.iter().step_by(width).copied().collect();
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

impl Program {
    /// The coefficients of the univariate rational function g(t) = f(shift + t (1, zz)) = n(t) / d(t), with
    /// deg n = dn, deg d = dd and d(0) = 1, from its values at the dn + dd + 1 points ts: [n_0 .. n_dn, d_1 .. d_dd].
    /// Err(Bad) if f is undefined at one of the points, or the system is singular there.
    pub fn line(&self, p: u64, shift: &[u64], zz: &[u64], ts: &[u64], dn: usize, dd: usize) -> Result<Vec<u64>, ColumnError> {
        let n = shift.len();
        if zz.len() + 1 != n || n != self.ninputs || ts.len() != dn + dd + 1 {
            return Err(ColumnError::Shape);
        }
        let mut points = Vec::with_capacity(ts.len() * n);
        for &t in ts {
            points.push((shift[0] + t) % p);
            for (i, &z) in zz.iter().enumerate() {
                points.push((shift[i + 1] + mulmod(t, z, p)) % p);
            }
        }
        let vals = self.eval(p, &points).map_err(|e| match e {
            FfError::BadPrime => ColumnError::BadPrime,
            FfError::Shape => ColumnError::Shape,
        })?;
        let width = self.outputs.len().max(1);
        let g: Vec<u64> = vals.iter().step_by(width).copied().collect();
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
        let c = prog.line(p, &[1, 2], &[5], &[7, 11, 13], 1, 1).unwrap();
        let half = invmod(2, p).unwrap();
        assert_eq!(c, vec![mulmod(5, half, p), mulmod(11, half, p), p - 2]);
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
        let c = prog.column(p, &[123_456_789, 987_654_321], 3, &support, &[11]).unwrap();
        assert_eq!(c, vec![3, 55, p - 7]);
    }
}
