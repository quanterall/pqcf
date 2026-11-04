//! Adjoint (reverse-mode) gradients of a pure Fock state's expectation value (Phase 4, FS 15 item 3; decision
//! record docs/decisions/AD1-adjoint-gradients.md).
//!
//! f = <psi_K| O |psi_K>, psi_k = G_k psi_(k-1), O real diagonal in the basis of FS 3.4 (a weighted sum of
//! marginal PNR projectors and photon numbers). With lambda_K = O psi_K and lambda_(k-1) = G_k^dagger lambda_k,
//! df/dp = 2 Re <lambda_k| (dG_k/dp) psi_(k-1)> for every parameter p of gate k: one forward and one backward pass.
//!
//! The derivatives are those of the truncated gates the backend applies, so the gradient is exact for the
//! truncated model:
//! - phase shift e^(i phi n), Kerr e^(i kappa n^2): i n psi_k, i n^2 psi_k;
//! - beamsplitter (fock.rs): d/dtheta = G_phi B, G_phi the sector generator with the phase phi; d/dphi =
//!   i (B n_a - n_a B), since U_N(theta, phi)[j][k] = e^(i phi (k - j)) U_N(theta, 0)[j][k];
//! - displacement D(r e^(i phi)): d/dr = (e^(i phi) a^dagger - e^(-i phi) a) D, rows of D one level beyond the
//!   cutoff; d/dphi [m][n] = i (m - n) D[m][n];
//! - squeezing S(r e^(i phi)) = exp((z* a^2 - z a^dagger^2)/2): d/dr = (e^(-i phi) a^2 - e^(i phi) a^dagger^2)/2 S,
//!   rows two levels beyond; d/dphi [m][n] = i (m - n)/2 S[m][n].
//!
//! Earlier states come back by the inverse of each exact gate (B(-theta, phi), the phase and Kerr with the
//! opposite sign, the mesh of U^dagger); only the state before each truncated gate is kept.
#![forbid(unsafe_code)] // NFR-04: kernels are safe Rust; unsafe only in the Rustler glue (lib.rs)
use crate::fock::{displacement_matrix, squeezing_matrix};
use crate::fockstate::{self, Basis, Error};
use crate::num::{Sum, C64};

pub enum Gate {
    Beamsplitter { theta: f64, phi: f64, ma: usize, mb: usize },
    Phaseshift { phi: f64, mode: usize },
    Kerr { kappa: f64, mode: usize },
    Passive { u: Vec<C64>, modes: Vec<usize> },
    Displacement { r: f64, phi: f64, mode: usize },
    Squeezing { r: f64, phi: f64, mode: usize },
}

impl Gate {
    /// the number of differentiable parameters, in the order of the instruction
    pub fn arity(&self) -> usize {
        match self {
            Gate::Beamsplitter { .. } | Gate::Displacement { .. } | Gate::Squeezing { .. } => 2,
            Gate::Phaseshift { .. } | Gate::Kerr { .. } => 1,
            Gate::Passive { .. } => 0,
        }
    }

    fn truncated(&self) -> bool {
        matches!(self, Gate::Displacement { .. } | Gate::Squeezing { .. })
    }
}

/// An observable term: weight times the marginal probability of `pattern` on `modes`, or times the photon number
/// of `modes`.
pub enum Term {
    Probability { weight: f64, modes: Vec<usize>, pattern: Vec<u32> },
    Photons { weight: f64, modes: Vec<usize> },
}

impl Term {
    /// the term's diagonal entry for a basis pattern, without its weight
    fn diag(&self, p: &[u32]) -> f64 {
        match self {
            Term::Probability { modes, pattern, .. } => modes.iter().zip(pattern).all(|(&m, &n)| p[m] == n) as u8 as f64,
            Term::Photons { modes, .. } => modes.iter().map(|&m| p[m] as f64).sum(),
        }
    }
    fn weight(&self) -> f64 {
        match self {
            Term::Probability { weight, .. } | Term::Photons { weight, .. } => *weight,
        }
    }
}

/// a truncated gate's single-mode matrix and its derivatives with respect to its two parameters
type Mats = (Vec<C64>, Vec<C64>, Vec<C64>);

fn conj(z: C64) -> C64 {
    C64::new(z.re, -z.im)
}

/// Re <a|b>
fn re_inner(a: &[C64], b: &[C64]) -> f64 {
    let mut s = Sum::new(a.len() as u64);
    for (x, y) in a.iter().zip(b) {
        s.add(x.re * y.re + x.im * y.im);
    }
    s.value()
}

fn conj_transpose(g: &[C64], d: usize) -> Vec<C64> {
    let mut t = vec![C64::ZERO; d * d];
    for r in 0..d {
        for c in 0..d {
            t[c * d + r] = conj(g[r * d + c]);
        }
    }
    t
}

/// rows 0..=c of G' (c + 1 + extra square, row-major) times the ladder operator L, as a (c+1) x (c+1) matrix:
/// (L G')[m][n] = sum over (offset, coefficient(m)) of coefficient(m) G'[m + offset][n]
fn ladder_times(gp: &[C64], c: usize, extra: usize, terms: &[(isize, &dyn Fn(usize) -> C64)]) -> Vec<C64> {
    let (d, dp) = (c + 1, c + 1 + extra);
    let mut out = vec![C64::ZERO; d * d];
    for m in 0..d {
        for &(off, coef) in terms {
            let row = m as isize + off;
            if row < 0 {
                continue;
            }
            let k = coef(m);
            for n in 0..d {
                out[m * d + n] = out[m * d + n].add(k.mul(gp[row as usize * dp + n]));
            }
        }
    }
    out
}

/// the matrix times i (m - n) s
fn phase_derivative(g: &[C64], d: usize, s: f64) -> Vec<C64> {
    let mut out = vec![C64::ZERO; d * d];
    for m in 0..d {
        for n in 0..d {
            let z = g[m * d + n].scale(s * (m as f64 - n as f64));
            out[m * d + n] = C64::new(-z.im, z.re);
        }
    }
    out
}

/// The single-mode matrices of a truncated gate: the gate and its two derivatives.
fn single_mode_matrices(g: &Gate, c: usize) -> Result<Mats, Error> {
    let d = c + 1;
    match *g {
        Gate::Displacement { r, phi, .. } => {
            let alpha = C64::new(r * phi.cos(), r * phi.sin());
            let m = displacement_matrix(alpha, c);
            let mp = displacement_matrix(alpha, c + 1);
            let (e, ec) = (C64::new(phi.cos(), phi.sin()), C64::new(phi.cos(), -phi.sin()));
            let up = move |m: usize| e.scale((m as f64).sqrt());
            let down = move |m: usize| ec.scale(-((m + 1) as f64).sqrt());
            let dr = ladder_times(&mp, c, 1, &[(-1, &up), (1, &down)]);
            let dphi = phase_derivative(&m, d, 1.0);
            Ok((m, dr, dphi))
        }
        Gate::Squeezing { r, phi, .. } => {
            let m = squeezing_matrix(r, phi, c).ok_or(Error::TooLarge)?;
            let mp = squeezing_matrix(r, phi, c + 2).ok_or(Error::TooLarge)?;
            let (e, ec) = (C64::new(phi.cos(), phi.sin()), C64::new(phi.cos(), -phi.sin()));
            let lower = move |m: usize| ec.scale((((m + 1) * (m + 2)) as f64).sqrt() / 2.0);
            let raise = move |m: usize| e.scale(-(((m * m.saturating_sub(1)) as f64).sqrt()) / 2.0);
            let dr = ladder_times(&mp, c, 2, &[(2, &lower), (-2, &raise)]);
            let dphi = phase_derivative(&m, d, 0.5);
            Ok((m, dr, dphi))
        }
        _ => Err(Error::Internal("not a single-mode truncated gate".into())),
    }
}

fn apply(b: &Basis, v: &mut [C64], g: &Gate, mats: Option<&Mats>) -> Result<(), Error> {
    match g {
        Gate::Beamsplitter { theta, phi, ma, mb } => fockstate::beamsplitter(b, v, *theta, *phi, *ma, *mb),
        Gate::Phaseshift { phi, mode } => {
            fockstate::phaseshift(b, v, *phi, *mode);
            Ok(())
        }
        Gate::Kerr { kappa, mode } => {
            fockstate::kerr(b, v, *kappa, *mode);
            Ok(())
        }
        Gate::Passive { u, modes } => fockstate::passive(b, v, u, modes),
        Gate::Displacement { mode, .. } | Gate::Squeezing { mode, .. } => {
            let m = mats.ok_or(Error::Internal("matrices".into()))?;
            fockstate::single_mode(b, v, &m.0, *mode);
            Ok(())
        }
    }
}

/// v <- G^dagger v
fn apply_adjoint(b: &Basis, v: &mut [C64], g: &Gate, mats: Option<&Mats>) -> Result<(), Error> {
    match g {
        Gate::Beamsplitter { theta, phi, ma, mb } => fockstate::beamsplitter(b, v, -theta, *phi, *ma, *mb),
        Gate::Phaseshift { phi, mode } => {
            fockstate::phaseshift(b, v, -phi, *mode);
            Ok(())
        }
        Gate::Kerr { kappa, mode } => {
            fockstate::kerr(b, v, -kappa, *mode);
            Ok(())
        }
        Gate::Passive { u, modes } => {
            let k = modes.len();
            fockstate::passive(b, v, &conj_transpose(u, k), modes)
        }
        Gate::Displacement { mode, .. } | Gate::Squeezing { mode, .. } => {
            let m = mats.ok_or(Error::Internal("matrices".into()))?;
            fockstate::single_mode(b, v, &conj_transpose(&m.0, b.cutoff + 1), *mode);
            Ok(())
        }
    }
}

/// v <- n_mode^power v (times i when `imag`)
fn number_times(b: &Basis, v: &[C64], mode: usize, power: i32, imag: bool) -> Vec<C64> {
    v.iter()
        .zip(&b.patterns)
        .map(|(a, p)| {
            let z = a.scale((p[mode] as f64).powi(power));
            if imag { C64::new(-z.im, z.re) } else { z }
        })
        .collect()
}

/// G_phi v on every sector of (ma, mb): (G_phi v)[j] = e^(i phi) sqrt((j+1)(N-j)) v[j+1] - e^(-i phi) sqrt(j(N-j+1)) v[j-1]
fn beamsplitter_generator(b: &Basis, v: &[C64], phi: f64, ma: usize, mb: usize) -> Vec<C64> {
    let mut out = vec![C64::ZERO; v.len()];
    let (e, ec) = (C64::new(phi.cos(), phi.sin()), C64::new(phi.cos(), -phi.sin()));
    for rep in b.patterns.iter().filter(|p| p[mb] == 0) {
        let n = rep[ma] as usize;
        if n == 0 {
            continue;
        }
        let mut q = rep.clone();
        let idx: Vec<usize> = (0..=n)
            .map(|j| {
                q[ma] = j as u32;
                q[mb] = (n - j) as u32;
                b.index(&q)
            })
            .collect();
        for j in 0..=n {
            let mut acc = C64::ZERO;
            if j < n {
                acc = acc.add(e.mul(v[idx[j + 1]]).scale((((j + 1) * (n - j)) as f64).sqrt()));
            }
            if j > 0 {
                acc = acc.add(ec.mul(v[idx[j - 1]]).scale(-(((j * (n - j + 1)) as f64).sqrt())));
            }
            out[idx[j]] = acc;
        }
    }
    out
}

/// df/dp for the parameters of gate g, given the state before it (x), after it (y) and the adjoint after it
fn gate_gradient(b: &Basis, g: &Gate, x: &[C64], y: &[C64], lambda: &[C64], mats: Option<&Mats>) -> Result<Vec<f64>, Error> {
    let two_re = |z: &[C64]| 2.0 * re_inner(lambda, z);
    Ok(match g {
        Gate::Phaseshift { mode, .. } => vec![two_re(&number_times(b, y, *mode, 1, true))],
        Gate::Kerr { mode, .. } => vec![two_re(&number_times(b, y, *mode, 2, true))],
        Gate::Beamsplitter { theta, phi, ma, mb } => {
            let dtheta = beamsplitter_generator(b, y, *phi, *ma, *mb);
            // i (B n_a x - n_a y)
            let mut bnx = number_times(b, x, *ma, 1, false);
            fockstate::beamsplitter(b, &mut bnx, *theta, *phi, *ma, *mb)?;
            let nay = number_times(b, y, *ma, 1, false);
            let dphi: Vec<C64> = bnx.iter().zip(&nay).map(|(p, q)| C64::new(-(p.im - q.im), p.re - q.re)).collect();
            vec![two_re(&dtheta), two_re(&dphi)]
        }
        Gate::Passive { .. } => vec![],
        Gate::Displacement { mode, .. } | Gate::Squeezing { mode, .. } => {
            let m = mats.ok_or(Error::Internal("matrices".into()))?;
            let mut zr = x.to_vec();
            fockstate::single_mode(b, &mut zr, &m.1, *mode);
            let mut zp = x.to_vec();
            fockstate::single_mode(b, &mut zp, &m.2, *mode);
            vec![two_re(&zr), two_re(&zp)]
        }
    })
}

/// The value of every term at the output state and the gradient of the weighted sum with respect to every gate
/// parameter (one vector per gate, `Gate::arity` entries).
pub fn adjoint(b: &Basis, x0: &[C64], gates: &[Gate], terms: &[Term]) -> Result<(Vec<f64>, Vec<Vec<f64>>), Error> {
    let c = b.cutoff;
    let mats: Vec<Option<Mats>> =
        gates.iter().map(|g| if g.truncated() { single_mode_matrices(g, c).map(Some) } else { Ok(None) }).collect::<Result<_, _>>()?;
    // forward, keeping the state before each truncated gate
    let mut psi = x0.to_vec();
    let mut kept: Vec<Option<Vec<C64>>> = Vec::with_capacity(gates.len());
    for (g, m) in gates.iter().zip(&mats) {
        kept.push(if g.truncated() { Some(psi.clone()) } else { None });
        apply(b, &mut psi, g, m.as_ref())?;
    }
    let d = psi.len() as u64;
    let values: Vec<f64> = terms
        .iter()
        .map(|t| {
            let mut s = Sum::new(d);
            for (a, p) in psi.iter().zip(&b.patterns) {
                s.add(t.diag(p) * (a.re * a.re + a.im * a.im));
            }
            s.value()
        })
        .collect();
    let mut lambda: Vec<C64> = psi
        .iter()
        .zip(&b.patterns)
        .map(|(a, p)| a.scale(terms.iter().map(|t| t.weight() * t.diag(p)).sum::<f64>()))
        .collect();
    // backward
    let mut grads = vec![Vec::new(); gates.len()];
    for k in (0..gates.len()).rev() {
        let (g, m) = (&gates[k], mats[k].as_ref());
        let x = match kept[k].take() {
            Some(x) => x,
            None => {
                let mut x = psi.clone();
                apply_adjoint(b, &mut x, g, m)?;
                x
            }
        };
        grads[k] = gate_gradient(b, g, &x, &psi, &lambda, m)?;
        apply_adjoint(b, &mut lambda, g, m)?;
        psi = x;
    }
    Ok((values, grads))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn value(b: &Basis, x0: &[C64], gates: &[Gate], terms: &[Term]) -> f64 {
        let mut psi = x0.to_vec();
        for g in gates {
            let mats = if g.truncated() { Some(single_mode_matrices(g, b.cutoff).unwrap()) } else { None };
            apply(b, &mut psi, g, mats.as_ref()).unwrap();
        }
        terms
            .iter()
            .map(|t| t.weight() * psi.iter().zip(&b.patterns).map(|(a, p)| t.diag(p) * (a.re * a.re + a.im * a.im)).sum::<f64>())
            .sum()
    }

    fn with_param(g: &Gate, j: usize, h: f64) -> Gate {
        match *g {
            Gate::Beamsplitter { theta, phi, ma, mb } => {
                if j == 0 { Gate::Beamsplitter { theta: theta + h, phi, ma, mb } } else { Gate::Beamsplitter { theta, phi: phi + h, ma, mb } }
            }
            Gate::Phaseshift { phi, mode } => Gate::Phaseshift { phi: phi + h, mode },
            Gate::Kerr { kappa, mode } => Gate::Kerr { kappa: kappa + h, mode },
            Gate::Displacement { r, phi, mode } => {
                if j == 0 { Gate::Displacement { r: r + h, phi, mode } } else { Gate::Displacement { r, phi: phi + h, mode } }
            }
            Gate::Squeezing { r, phi, mode } => {
                if j == 0 { Gate::Squeezing { r: r + h, phi, mode } } else { Gate::Squeezing { r, phi: phi + h, mode } }
            }
            Gate::Passive { .. } => unreachable!(),
        }
    }

    fn clone_gate(g: &Gate) -> Gate {
        match g {
            Gate::Passive { u, modes } => Gate::Passive { u: u.clone(), modes: modes.clone() },
            _ => with_param(g, 0, 0.0),
        }
    }

    /// every adjoint derivative against a fourth-order central difference
    fn check(m: usize, c: usize, occ: &[u32], gates: Vec<Gate>, terms: Vec<Term>, tol: f64) {
        let b = Basis::new(m, c);
        let mut x0 = vec![C64::ZERO; b.patterns.len()];
        x0[b.index(occ)] = C64::ONE;
        let (_, grads) = adjoint(&b, &x0, &gates, &terms).unwrap();
        let h = 1e-3;
        for (k, g) in gates.iter().enumerate() {
            for j in 0..g.arity() {
                let f = |s: f64| {
                    let gs: Vec<Gate> = gates.iter().enumerate().map(|(i, x)| if i == k { with_param(x, j, s) } else { clone_gate(x) }).collect();
                    value(&b, &x0, &gs, &terms)
                };
                let fd = (-f(2.0 * h) + 8.0 * f(h) - 8.0 * f(-h) + f(-2.0 * h)) / (12.0 * h);
                assert!((grads[k][j] - fd).abs() < tol, "gate {k} param {j}: adjoint {} vs fd {fd}", grads[k][j]);
            }
        }
    }

    #[test]
    fn passive_and_kerr_against_differences() {
        let u = vec![
            C64::new(0.6, 0.0), C64::new(0.0, 0.8),
            C64::new(0.0, 0.8), C64::new(0.6, 0.0),
        ];
        let gates = vec![
            Gate::Beamsplitter { theta: 0.4, phi: 0.3, ma: 0, mb: 1 },
            Gate::Phaseshift { phi: 0.7, mode: 1 },
            Gate::Kerr { kappa: 0.2, mode: 0 },
            Gate::Passive { u, modes: vec![1, 2] },
            Gate::Beamsplitter { theta: 1.1, phi: -0.5, ma: 1, mb: 2 },
        ];
        let terms = vec![
            Term::Probability { weight: 1.0, modes: vec![0, 2], pattern: vec![1, 1] },
            Term::Photons { weight: 0.5, modes: vec![1] },
        ];
        check(3, 3, &[2, 1, 0], gates, terms, 1e-9);
    }

    #[test]
    fn truncated_gates_against_differences() {
        let gates = vec![
            Gate::Displacement { r: 0.4, phi: 0.3, mode: 0 },
            Gate::Squeezing { r: 0.3, phi: -0.6, mode: 1 },
            Gate::Beamsplitter { theta: 0.8, phi: 0.2, ma: 0, mb: 1 },
            Gate::Kerr { kappa: 0.3, mode: 1 },
            Gate::Displacement { r: 0.2, phi: 1.2, mode: 1 },
        ];
        let terms = vec![
            Term::Probability { weight: 1.0, modes: vec![0], pattern: vec![1] },
            Term::Photons { weight: -0.3, modes: vec![0, 1] },
        ];
        check(2, 8, &[0, 1], gates, terms, 1e-9);
    }
}
