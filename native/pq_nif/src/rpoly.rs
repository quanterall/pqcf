//! Polynomials that stay in Rust between operations (decision record K2): the symbolic engine's plain polynomials
//! (pq_sym_poly; symbols, conj and pi with integer exponents, exact Gaussian-rational coefficients) as a NIF
//! resource, so a chain of products, powers and sums never builds its intermediate results as Erlang terms. The
//! variables live on the Erlang side ({rpoly, Vars, Resource}); here a term is its exponents, one per variable, and
//! a Gaussian-integer numerator over one common denominator, kept in lowest terms (the gcd of every numerator and
//! the denominator is 1). Every operation is exact or refused (overflow, packing); pq_sym_rpoly then falls back to
//! Erlang maps, with the same result.

use crate::polymul::{self, Poly, Refused};
use std::cmp::Ordering;

/// A polynomial over `nvars` variables: Poly's terms (sorted by exponents, no zero coefficient) over `denom`.
pub struct RPoly {
    pub nvars: usize,
    pub p: Poly,
    pub denom: i128,
}

impl RPoly {
    pub fn len(&self) -> usize {
        self.p.coef.len() / 2
    }

    pub fn is_empty(&self) -> bool {
        self.p.coef.is_empty()
    }

    /// the terms over `denom` in lowest terms: the gcd of every numerator and the denominator divided out
    pub fn new(nvars: usize, p: Poly, denom: i128) -> RPoly {
        let g = p.coef.iter().fold(denom.unsigned_abs(), |g, &c| polymul::gcd(g, c.unsigned_abs()));
        let (p, denom) = if g > 1 {
            let g = g as i128;
            (Poly { exps: p.exps, coef: p.coef.into_iter().map(|c| c / g).collect() }, denom / g)
        } else {
            (p, denom)
        };
        RPoly { nvars, p, denom }
    }

    /// the terms in `n` variables, variable j going to position `map[j]` (an injective map)
    pub fn remap(&self, map: &[usize], n: usize) -> Poly {
        let mut exps = vec![0i64; self.len() * n];
        for t in 0..self.len() {
            for (j, &to) in map.iter().enumerate() {
                exps[t * n + to] = self.p.exps[t * self.nvars + j];
            }
        }
        Poly { exps, coef: self.p.coef.clone() }
    }
}

fn scaled(p: Poly, k: i128) -> Result<Poly, Refused> {
    if k == 1 {
        return Ok(p);
    }
    let coef = p.coef.iter().map(|&c| c.checked_mul(k).ok_or(Refused::Overflow)).collect::<Result<_, _>>()?;
    Ok(Poly { exps: p.exps, coef })
}

pub fn mul(a: &RPoly, ma: &[usize], b: &RPoly, mb: &[usize], n: usize) -> Result<RPoly, Refused> {
    let denom = a.denom.checked_mul(b.denom).ok_or(Refused::Overflow)?;
    Ok(RPoly::new(n, polymul::mul(n, &a.remap(ma, n), &b.remap(mb, n))?, denom))
}

pub fn add(a: &RPoly, ma: &[usize], b: &RPoly, mb: &[usize], n: usize) -> Result<RPoly, Refused> {
    let g = polymul::gcd(a.denom as u128, b.denom as u128) as i128;
    let denom = (a.denom / g).checked_mul(b.denom).ok_or(Refused::Overflow)?;
    let pa = scaled(a.remap(ma, n), denom / a.denom)?;
    let pb = scaled(b.remap(mb, n), denom / b.denom)?;
    Ok(RPoly::new(n, polymul::add(n, &pa, &pb)?, denom))
}

pub fn pow(a: &RPoly, k: u32) -> Result<RPoly, Refused> {
    let denom = a.denom.checked_pow(k).ok_or(Refused::Overflow)?;
    Ok(RPoly::new(a.nvars, polymul::pow(a.nvars, &a.p, k)?, denom))
}

/// a times (re + i im) / den
pub fn scale(a: &RPoly, re: i128, im: i128, den: i128) -> Result<RPoly, Refused> {
    if re == 0 && im == 0 {
        return Ok(RPoly { nvars: a.nvars, p: Poly { exps: vec![], coef: vec![] }, denom: 1 });
    }
    let denom = a.denom.checked_mul(den).ok_or(Refused::Overflow)?;
    let m = |x: i128, y: i128| x.checked_mul(y).ok_or(Refused::Overflow);
    let mut coef = Vec::with_capacity(a.p.coef.len());
    for c in a.p.coef.chunks_exact(2) {
        coef.push(m(c[0], re)?.checked_sub(m(c[1], im)?).ok_or(Refused::Overflow)?);
        coef.push(m(c[0], im)?.checked_add(m(c[1], re)?).ok_or(Refused::Overflow)?);
    }
    Ok(RPoly::new(a.nvars, Poly { exps: a.p.exps.clone(), coef }, denom))
}

/// The terms in pq_sym_norm:print_key's order, for variables that are all symbols ({sym, X}, sorted as Erlang sorts
/// them): the constant first, then by total degree, then by the list of names present (a list that is a prefix of
/// another first), then by the exponents.
pub fn print_order(a: &RPoly) -> Vec<usize> {
    let n = a.nvars;
    let row = |t: usize| &a.p.exps[t * n..(t + 1) * n];
    let keys: Vec<(bool, i64)> = (0..a.len()).map(|t| (row(t).iter().any(|&e| e != 0), row(t).iter().sum())).collect();
    let mut idx: Vec<usize> = (0..a.len()).collect();
    idx.sort_unstable_by(|&x, &y| {
        keys[x].cmp(&keys[y]).then_with(|| {
            let (rx, ry) = (row(x), row(y));
            // the names: the positions of the nonzero exponents, compared as lists
            let mut px = (0..n).filter(|&j| rx[j] != 0);
            let mut py = (0..n).filter(|&j| ry[j] != 0);
            loop {
                match (px.next(), py.next()) {
                    (None, None) => break,
                    (None, Some(_)) => return Ordering::Less,
                    (Some(_), None) => return Ordering::Greater,
                    (Some(i), Some(j)) if i != j => return i.cmp(&j),
                    _ => {}
                }
            }
            // the same names: the exponents in place
            (0..n).filter(|&j| rx[j] != 0).map(|j| rx[j].cmp(&ry[j])).find(|o| o.is_ne()).unwrap_or(Ordering::Equal)
        })
    });
    idx
}
