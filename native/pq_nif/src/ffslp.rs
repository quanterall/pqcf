//! Straight-line programs evaluated over a prime field (decision record R1): the black box of the reconstruction.
//! A program is prepared once (its operations, exact rational constants), then evaluated modulo a prime p < 2^63
//! at a batch of points. A point where the program divides by zero mod p is reported as bad (every output
//! `BAD`), for the reconstruction to skip; a constant whose denominator is 0 mod p makes the prime unusable.

/// The value of every output at a bad point (never a residue, as p < 2^63).
pub const BAD: u64 = u64::MAX;

pub enum Op {
    /// the point's coordinate i
    Input(usize),
    /// num / den, den > 0
    Const(i128, i128),
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

fn residue(x: i128, p: u64) -> u64 {
    x.rem_euclid(p as i128) as u64
}

impl Program {
    pub fn check(&self) -> Result<(), FfError> {
        for (s, op) in self.ops.iter().enumerate() {
            let ok = match *op {
                Op::Input(i) => i < self.ninputs,
                Op::Const(_, d) => d > 0,
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
            if let Op::Const(n, d) = *op {
                let di = invmod(residue(d, p), p).ok_or(FfError::BadPrime)?;
                consts[s] = mulmod(residue(n, p), di, p);
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
                Op::Const(-2, 3),
                Op::Mul(2, 1),
                Op::Add(0, 3),
                Op::Pow(4, 2),
                Op::Add(0, 1),
                Op::Pow(6, -1),
                Op::Mul(5, 7),
                Op::Const(1, 1),
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
        assert_eq!(prog.eval(3, &[1, 2]).err(), Some(FfError::BadPrime));
        assert_eq!(Program { ops: vec![Op::Mul(0, 1)], outputs: vec![0], ninputs: 0 }.check().err(), Some(FfError::Shape));
    }
}
