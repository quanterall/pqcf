//! The operations of pq_sym_slp's Horner scheme for one variable order, counted without building the program
//! (decision record H4): the score of the order search, as FORM scores its tree search by counting. It follows
//! pq_sym_slp's horner/3 step for step -- every node primitive (its content, the sign included, taken out and
//! shared), the sign folded into an addition or subtraction, a content other than +-1 one product, the leaf sums
//! adding positive terms and subtracting negative ones, every atom power, monomial, product and sum keyed by its
//! structure so that each is counted once -- for polynomials of symbols only, so that the count equals the
//! program's (tested against it). Coefficients are Gaussian rationals over a common denominator, exact in i128;
//! an overflow refuses the count.

use crate::polymul::{gcd, Poly, Refused};
use std::cmp::Ordering;
use std::collections::HashMap;

#[derive(Hash, PartialEq, Eq)]
enum Key {
    Const(i128, i128, i128),
    Pow(usize, i64),
    Mul(u32, u32),
    Add(u32, u32),
    Sub(u32, u32),
    Mono(Vec<i64>),
    HPoly(Vec<i64>, Vec<i128>),
}

/// A polynomial with Gaussian-integer numerators over `den` (> 0): terms (exponents, (re, im)).
#[derive(Clone)]
struct Pol {
    terms: Vec<(Vec<i64>, (i128, i128))>,
    den: i128,
}

struct Counter {
    nvars: usize,
    /// the order's rank of every variable (lower first)
    rank: Vec<usize>,
    keys: HashMap<Key, u32>,
    next: u32,
    ops: u64,
}

/// Erlang's term order of two monomials' atom lists [{{sym, X}, {E, 1}}, ...] (variables sorted as Erlang sorts
/// them): element by element, the variable first, then the exponent; a prefix first
fn mono_cmp(a: &[i64], b: &[i64]) -> Ordering {
    let mut pa = (0..a.len()).filter(|&j| a[j] != 0);
    let mut pb = (0..b.len()).filter(|&j| b[j] != 0);
    loop {
        match (pa.next(), pb.next()) {
            (None, None) => return Ordering::Equal,
            (None, Some(_)) => return Ordering::Less,
            (Some(_), None) => return Ordering::Greater,
            (Some(i), Some(j)) => {
                let o = i.cmp(&j).then(a[i].cmp(&b[j]));
                if o.is_ne() {
                    return o;
                }
            }
        }
    }
}

fn ck(x: Option<i128>) -> Result<i128, Refused> {
    x.ok_or(Refused::Overflow)
}

impl Counter {
    fn slot(&mut self, k: Key, counted: bool) -> u32 {
        if let Some(&s) = self.keys.get(&k) {
            return s;
        }
        let s = self.next;
        self.next += 1;
        if counted {
            self.ops += 1;
        }
        self.keys.insert(k, s);
        s
    }

    fn konst(&mut self, re: i128, im: i128, den: i128) -> u32 {
        let g = gcd(gcd(re.unsigned_abs(), im.unsigned_abs()), den as u128).max(1) as i128;
        self.slot(Key::Const(re / g, im / g, den / g), false)
    }

    /// pq_sym_poly:content(P, true): {G, P / G}, G = sign gcd(numerators) / (den), the sign making the first term's
    /// first non-zero part positive (the first term in Erlang's order of monomials)
    fn content(&self, p: &Pol) -> ((i128, i128), Pol) {
        let g = p.terms.iter().fold(0u128, |g, (_, (re, im))| gcd(gcd(g, re.unsigned_abs()), im.unsigned_abs())).max(1) as i128;
        let first = p.terms.iter().min_by(|x, y| mono_cmp(&x.0, &y.0)).map(|t| t.1).unwrap_or((1, 0));
        let neg = first.0 < 0 || (first.0 == 0 && first.1 < 0);
        let s = if neg { -g } else { g };
        let terms = p.terms.iter().map(|(e, (re, im))| (e.clone(), (re / s, im / s))).collect();
        let d = gcd(g as u128, p.den as u128) as i128;
        ((s / d, p.den / d), Pol { terms, den: 1 })
    }

    fn atom_pow(&mut self, x: usize, k: i64) -> u32 {
        if k == 1 { x as u32 } else { self.slot(Key::Pow(x, k), true) }
    }

    fn mono(&mut self, e: &[i64]) -> u32 {
        if let Some(&s) = self.keys.get(&Key::Mono(e.to_vec())) {
            return s;
        }
        let mut prev: Option<u32> = None;
        for (j, &ej) in e.iter().enumerate() {
            if ej != 0 {
                let f = self.atom_pow(j, ej);
                prev = Some(match prev {
                    None => f,
                    Some(p) => self.slot(Key::Mul(p, f), true),
                });
            }
        }
        let s = prev.expect("a monomial with an atom");
        self.keys.insert(Key::Mono(e.to_vec()), s);
        s
    }

    /// term(M, C): a constant is free, a monomial times 1 is the monomial, otherwise one product
    fn term(&mut self, e: &[i64], (re, im): (i128, i128), den: i128) -> u32 {
        if e.iter().all(|&x| x == 0) {
            return self.konst(re, im, den);
        }
        let m = self.mono(e);
        if re == den && im == 0 {
            m
        } else {
            let c = self.konst(re, im, den);
            self.slot(Key::Mul(c, m), true)
        }
    }

    fn scaled(&mut self, (n, d): (i128, i128), s: u32) -> u32 {
        if n == 1 && d == 1 {
            s
        } else {
            let c = self.konst(n, 0, d);
            self.slot(Key::Mul(c, s), true)
        }
    }

    fn horner(&mut self, p: &Pol) -> Result<u32, Refused> {
        let (s, g) = self.hc(p)?;
        Ok(self.scaled(g, s))
    }

    /// (slot of P's primitive part, P's content)
    fn hc(&mut self, p: &Pol) -> Result<(u32, (i128, i128)), Refused> {
        if p.terms.is_empty() {
            return Ok((self.konst(0, 0, 1), (1, 1)));
        }
        let (g, prim) = self.content(p);
        let mut sorted = prim.terms.clone();
        sorted.sort_unstable_by(|x, y| x.0.cmp(&y.0));
        let key = Key::HPoly(sorted.iter().flat_map(|t| t.0.iter().copied()).collect(), sorted.iter().flat_map(|t| [t.1 .0, t.1 .1]).collect());
        if let Some(&s) = self.keys.get(&key) {
            return Ok((s, g));
        }
        let s = self.hprim(&prim)?;
        self.keys.insert(key, s);
        Ok((s, g))
    }

    fn hprim(&mut self, p: &Pol) -> Result<u32, Refused> {
        if p.terms.len() == 1 {
            let (e, c) = &p.terms[0];
            return Ok(self.term(e, *c, p.den));
        }
        // the first variable of the order that at least two terms contain
        let mut best: Option<usize> = None;
        for j in 0..self.nvars {
            if p.terms.iter().filter(|t| t.0[j] > 0).count() >= 2 && best.is_none_or(|b| self.rank[j] < self.rank[b]) {
                best = Some(j);
            }
        }
        let Some(x) = best else { return self.sum_terms(p) };
        let (with, without): (Vec<_>, Vec<_>) = p.terms.iter().cloned().partition(|t| t.0[x] > 0);
        let k = with.iter().map(|t| t.0[x]).min().unwrap_or(1);
        let q = Pol {
            terms: with.into_iter().map(|(mut e, c)| { e[x] -= k; (e, c) }).collect(),
            den: p.den,
        };
        let sx = self.atom_pow(x, k);
        let (sq, gq) = self.hc(&q)?;
        let prod = self.slot(Key::Mul(sx, sq), true);
        if without.is_empty() {
            return Ok(self.scaled(gq, prod));
        }
        let (sr, gr) = self.hc(&Pol { terms: without, den: p.den })?;
        Ok(self.combine(gr, sr, gq, prod))
    }

    fn magnitude(&mut self, (n, d): (i128, i128), s: u32) -> (i32, u32) {
        (if n < 0 { -1 } else { 1 }, self.scaled((n.abs(), d), s))
    }

    fn combine(&mut self, gr: (i128, i128), sr: u32, gq: (i128, i128), sp: u32) -> u32 {
        let (sign_r, ar) = self.magnitude(gr, sr);
        let (sign_q, aq) = self.magnitude(gq, sp);
        match (sign_r, sign_q) {
            (1, 1) => self.slot(Key::Add(ar, aq), true),
            (1, _) => self.slot(Key::Sub(ar, aq), true),
            (_, 1) => self.slot(Key::Sub(aq, ar), true),
            _ => {
                let t = self.slot(Key::Add(ar, aq), true);
                self.scaled((-1, 1), t)
            }
        }
    }

    /// positive terms added in Erlang's order of monomials, then the negative real ones subtracted by magnitude
    fn sum_terms(&mut self, p: &Pol) -> Result<u32, Refused> {
        let mut terms = p.terms.clone();
        terms.sort_unstable_by(|x, y| mono_cmp(&x.0, &y.0));
        let mut pos = Vec::new();
        let mut neg = Vec::new();
        for (e, (re, im)) in &terms {
            if im == &0 && *re < 0 {
                neg.push(self.term(e, (ck(re.checked_neg())?, 0), p.den));
            } else {
                pos.push(self.term(e, (*re, *im), p.den));
            }
        }
        let add = |c: &mut Counter, v: &[u32]| v[1..].iter().fold(v[0], |a, &b| c.slot(Key::Add(a, b), true));
        Ok(if pos.is_empty() {
            let t = add(self, &neg);
            self.scaled((-1, 1), t)
        } else {
            let a = add(self, &pos);
            neg.iter().fold(a, |a, &b| self.slot(Key::Sub(a, b), true))
        })
    }
}

/// The operations of the Horner programs of `polys` (one after the other, sharing every slot) for the order
/// `order` (variable indices, first to last), each polynomial as Poly over `den`.
pub fn count(nvars: usize, polys: &[(Poly, i128)], order: &[usize]) -> Result<u64, Refused> {
    let mut rank = vec![usize::MAX; nvars];
    for (r, &j) in order.iter().enumerate() {
        rank[j] = r;
    }
    let mut c = Counter { nvars, rank, keys: HashMap::new(), next: nvars as u32, ops: 0 };
    for (p, den) in polys {
        let pol = Pol {
            terms: (0..p.len()).map(|t| (p.exps[t * nvars..(t + 1) * nvars].to_vec(), (p.coef[2 * t], p.coef[2 * t + 1]))).collect(),
            den: *den,
        };
        match pol.terms.len() {
            0 => {}
            1 => {
                let (e, co) = &pol.terms[0];
                c.term(e, *co, pol.den);
            }
            _ => {
                c.horner(&pol)?;
            }
        }
    }
    Ok(c.ops)
}
