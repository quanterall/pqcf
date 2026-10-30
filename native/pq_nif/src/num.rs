//! Complex numbers and summation with the rules of AN 2.2: plain f64 products, no fused multiply-add, Neumaier
//! compensated summation for sums of more than 64 terms (real and imaginary parts separately).

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct C64 {
    pub re: f64,
    pub im: f64,
}

// add and mul are named like the operator traits on purpose: explicit calls keep the operation order visible
#[allow(clippy::should_implement_trait)]
impl C64 {
    pub const ZERO: C64 = C64 { re: 0.0, im: 0.0 };
    pub const ONE: C64 = C64 { re: 1.0, im: 0.0 };

    pub fn new(re: f64, im: f64) -> C64 {
        C64 { re, im }
    }
    #[inline]
    pub fn mul(self, o: C64) -> C64 {
        C64 { re: self.re * o.re - self.im * o.im, im: self.re * o.im + self.im * o.re }
    }
    #[inline]
    pub fn add(self, o: C64) -> C64 {
        C64 { re: self.re + o.re, im: self.im + o.im }
    }
    #[inline]
    pub fn scale(self, s: f64) -> C64 {
        C64 { re: self.re * s, im: self.im * s }
    }
    /// The modulus; sqrt is correctly rounded in IEEE-754, so this is the same on every target.
    #[inline]
    pub fn abs(self) -> f64 {
        (self.re * self.re + self.im * self.im).sqrt()
    }
    pub fn is_finite(self) -> bool {
        self.re.is_finite() && self.im.is_finite()
    }
}

/// Sums more than 64 terms: Neumaier summation.
pub const NEUMAIER_ABOVE: u64 = 64;

/// A real sum whose term count is known in advance, so that the rule of AN 2.2 picks the algorithm.
#[derive(Clone, Copy, Debug)]
pub struct Sum {
    s: f64,
    c: f64,
    compensated: bool,
}

impl Sum {
    pub fn new(terms: u64) -> Sum {
        Sum { s: 0.0, c: 0.0, compensated: terms > NEUMAIER_ABOVE }
    }
    #[inline]
    pub fn add(&mut self, x: f64) {
        if self.compensated {
            let t = self.s + x;
            if self.s.abs() >= x.abs() {
                self.c += (self.s - t) + x;
            } else {
                self.c += (x - t) + self.s;
            }
            self.s = t;
        } else {
            self.s += x;
        }
    }
    pub fn value(&self) -> f64 {
        self.s + self.c
    }
}

#[derive(Clone, Copy, Debug)]
pub struct CSum {
    re: Sum,
    im: Sum,
}

impl CSum {
    pub fn new(terms: u64) -> CSum {
        CSum { re: Sum::new(terms), im: Sum::new(terms) }
    }
    #[inline]
    pub fn add(&mut self, z: C64) {
        self.re.add(z.re);
        self.im.add(z.im);
    }
    pub fn value(&self) -> C64 {
        C64::new(self.re.value(), self.im.value())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn neumaier_recovers_what_plain_summation_loses() {
        // 1, then 100 times 1e-16, then -1: the exact sum is 1e-14
        let mut xs = vec![1.0];
        xs.extend(std::iter::repeat_n(1e-16, 100));
        xs.push(-1.0);
        let mut n = Sum::new(xs.len() as u64);
        let mut p = Sum::new(2);
        for &x in &xs {
            n.add(x);
            p.add(x);
        }
        assert!((n.value() - 1e-14).abs() < 1e-28);
        assert_eq!(p.value(), 0.0); // the plain sum loses every small term
    }
}
