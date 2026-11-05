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

/// d/dp of gate g applied to v, parameter j in the gate's order (the pure-state derivatives above, per vector)
fn derivative_vector(b: &Basis, g: &Gate, j: usize, v: &[C64], mats: Option<&Mats>) -> Result<Vec<C64>, Error> {
    let mut y = v.to_vec();
    apply(b, &mut y, g, mats)?;
    Ok(match g {
        Gate::Phaseshift { mode, .. } => number_times(b, &y, *mode, 1, true),
        Gate::Kerr { mode, .. } => number_times(b, &y, *mode, 2, true),
        Gate::Beamsplitter { theta, phi, ma, mb } => {
            if j == 0 {
                beamsplitter_generator(b, &y, *phi, *ma, *mb)
            } else {
                let mut bnx = number_times(b, v, *ma, 1, false);
                fockstate::beamsplitter(b, &mut bnx, *theta, *phi, *ma, *mb)?;
                let nay = number_times(b, &y, *ma, 1, false);
                bnx.iter().zip(&nay).map(|(p, q)| C64::new(-(p.im - q.im), p.re - q.re)).collect()
            }
        }
        Gate::Displacement { mode, .. } | Gate::Squeezing { mode, .. } => {
            let m = mats.ok_or(Error::Internal("matrices".into()))?;
            let mut z = v.to_vec();
            fockstate::single_mode(b, &mut z, if j == 0 { &m.1 } else { &m.2 }, *mode);
            z
        }
        Gate::Passive { .. } => return Err(Error::Internal("no parameters".into())),
    })
}

/// A channel of the density-matrix adjoint: a gate, pure loss (its transmissivity differentiable), or thermal loss
/// (passed through; its parameters are not differentiated).
pub enum Channel {
    Gate(Gate),
    Loss { eta: f64, mode: usize },
    Thermal { eta: f64, nth: f64, mode: usize },
}

impl Channel {
    pub fn arity(&self) -> usize {
        match self {
            Channel::Gate(g) => g.arity(),
            Channel::Loss { .. } => 1,
            Channel::Thermal { .. } => 0,
        }
    }
    fn kept(&self) -> bool {
        !matches!(self, Channel::Gate(g) if !g.truncated())
    }
}

/// G on every column of a D x D matrix (row-major)
fn on_columns(m: &mut [C64], d: usize, mut f: impl FnMut(&mut [C64]) -> Result<(), Error>) -> Result<(), Error> {
    let mut col = vec![C64::ZERO; d];
    for j in 0..d {
        for i in 0..d {
            col[i] = m[i * d + j];
        }
        f(&mut col)?;
        for i in 0..d {
            m[i * d + j] = col[i];
        }
    }
    Ok(())
}

/// Re tr(A B)
fn re_trace_product(a: &[C64], bm: &[C64], d: usize) -> f64 {
    let mut s = Sum::new((d * d) as u64);
    for i in 0..d {
        for j in 0..d {
            let (x, y) = (a[i * d + j], bm[j * d + i]);
            s.add(x.re * y.re - x.im * y.im);
        }
    }
    s.value()
}

fn binom(n: u32, k: u32) -> f64 {
    (0..k).fold(1.0, |acc, i| acc * (n - i) as f64 / (i + 1) as f64)
}

/// the adjoint of FS A.11 loss: Lambda'[p][q] = sum_l w(np, nq, l) Lambda[p - l][q - l]
fn loss_adjoint(b: &Basis, lam: &[C64], eta: f64, mode: usize) -> Vec<C64> {
    let d = b.patterns.len();
    let mut out = vec![C64::ZERO; d * d];
    for (i, p) in b.patterns.iter().enumerate() {
        for (j, q) in b.patterns.iter().enumerate() {
            let (np, nq) = (p[mode], q[mode]);
            let (mut pp, mut qq) = (p.clone(), q.clone());
            let mut acc = C64::ZERO;
            for l in 0..=np.min(nq) {
                pp[mode] = np - l;
                qq[mode] = nq - l;
                let w = (binom(np, l) * binom(nq, l)).sqrt() * eta.powf((np + nq - 2 * l) as f64 / 2.0) * (1.0 - eta).powi(l as i32);
                acc = acc.add(lam[b.index(&pp) * d + b.index(&qq)].scale(w));
            }
            out[i * d + j] = acc;
        }
    }
    out
}

/// d/d eta of eta^a (1 - eta)^l, without 0 * inf at the ends
fn dpow(eta: f64, a: f64, l: i32) -> f64 {
    let first = if a == 0.0 { 0.0 } else { a * eta.powf(a - 1.0) * (1.0 - eta).powi(l) };
    let second = if l == 0 { 0.0 } else { l as f64 * eta.powf(a) * (1.0 - eta).powi(l - 1) };
    first - second
}

/// tr(Lambda dL/deta (rho)) for loss on `mode`
fn loss_gradient(b: &Basis, rho: &[C64], lam: &[C64], eta: f64, mode: usize) -> f64 {
    let d = b.patterns.len();
    let mut s = Sum::new((d * d) as u64);
    for (i, p) in b.patterns.iter().enumerate() {
        for (j, q) in b.patterns.iter().enumerate() {
            let z = rho[i * d + j];
            if z.re == 0.0 && z.im == 0.0 {
                continue;
            }
            let (np, nq) = (p[mode], q[mode]);
            let (mut pp, mut qq) = (p.clone(), q.clone());
            for l in 0..=np.min(nq) {
                pp[mode] = np - l;
                qq[mode] = nq - l;
                let w = (binom(np, l) * binom(nq, l)).sqrt() * dpow(eta, (np + nq - 2 * l) as f64 / 2.0, l as i32);
                // tr(Lambda X) picks Lambda[target col][target row] X[target row][target col]
                let lv = lam[b.index(&qq) * d + b.index(&pp)];
                s.add(w * (z.re * lv.re - z.im * lv.im));
            }
        }
    }
    s.value()
}

/// the adjoint of the truncated amplifier (fockstate::dm::amplifier): Lambda'[p][q] = sum_k w_k Lambda[p + k][q + k]
/// over the same k (the room under the cutoff)
fn amplifier_adjoint(b: &Basis, lam: &[C64], g: f64, mode: usize) -> Vec<C64> {
    let d = b.patterns.len();
    let binom_up = |n: u32, k: u32| -> f64 { (0..k).fold(1.0, |acc, i| acc * (n + k - i) as f64 / (i + 1) as f64) };
    let ratio = (g - 1.0) / g;
    let total: Vec<u32> = b.patterns.iter().map(|p| p.iter().sum()).collect();
    let mut out = vec![C64::ZERO; d * d];
    for (i, p) in b.patterns.iter().enumerate() {
        for (j, q) in b.patterns.iter().enumerate() {
            let (np, nq) = (p[mode], q[mode]);
            let room = (b.cutoff as u32).saturating_sub(total[i].max(total[j]));
            let base = g.powf(-((np + nq + 2) as f64) / 2.0);
            let (mut pp, mut qq) = (p.clone(), q.clone());
            let mut acc = C64::ZERO;
            for k in 0..=room {
                if k > 0 && ratio == 0.0 {
                    break;
                }
                pp[mode] = np + k;
                qq[mode] = nq + k;
                let w = (binom_up(np, k) * binom_up(nq, k)).sqrt() * base * ratio.powi(k as i32);
                acc = acc.add(lam[b.index(&pp) * d + b.index(&qq)].scale(w));
            }
            out[i * d + j] = acc;
        }
    }
    out
}

/// The density-matrix adjoint (Phase 4 follow-up; docs/decisions/AD1-adjoint-gradients.md): f = tr(O rho_K),
/// rho_k = Phi_k(rho_(k-1)); Lambda_K = O, Lambda_(k-1) = Phi_k^dagger(Lambda_k); df/dp = tr(Lambda_k dPhi_k/dp
/// (rho_(k-1))): for a gate 2 Re tr(G^dagger Lambda (dG/dp) rho), for loss the derivative of its Kraus weights.
/// The state before each loss, thermal loss and truncated gate is kept; exact gates are inverted.
pub fn adjoint_dm(b: &Basis, x0: &[C64], chans: &[Channel], terms: &[Term]) -> Result<(Vec<f64>, Vec<Vec<f64>>), Error> {
    use fockstate::dm;
    let d = b.patterns.len();
    let c = b.cutoff;
    let mats: Vec<Option<Mats>> = chans
        .iter()
        .map(|ch| match ch {
            Channel::Gate(g) if g.truncated() => single_mode_matrices(g, c).map(Some),
            _ => Ok(None),
        })
        .collect::<Result<_, _>>()?;
    let mut rho = dm::from_pure(x0);
    let mut kept: Vec<Option<Vec<C64>>> = Vec::with_capacity(chans.len());
    for (ch, m) in chans.iter().zip(&mats) {
        kept.push(if ch.kept() { Some(rho.clone()) } else { None });
        match ch {
            Channel::Gate(g) => {
                dm::unitary_like(b, &mut rho, |v| apply(b, v, g, m.as_ref()))?;
            }
            Channel::Loss { eta, mode } => rho = dm::loss(b, &rho, *eta, *mode),
            Channel::Thermal { eta, nth, mode } => rho = dm::thermal_loss(b, &rho, *eta, *nth, *mode).0,
        }
    }
    let values: Vec<f64> = terms
        .iter()
        .map(|t| {
            let mut s = Sum::new(d as u64);
            for (i, p) in b.patterns.iter().enumerate() {
                s.add(t.diag(p) * rho[i * d + i].re);
            }
            s.value()
        })
        .collect();
    let mut lam = vec![C64::ZERO; d * d];
    for (i, p) in b.patterns.iter().enumerate() {
        lam[i * d + i] = C64::new(terms.iter().map(|t| t.weight() * t.diag(p)).sum::<f64>(), 0.0);
    }
    let mut grads = vec![Vec::new(); chans.len()];
    for k in (0..chans.len()).rev() {
        let m = mats[k].as_ref();
        let x = match kept[k].take() {
            Some(x) => x,
            None => {
                let Channel::Gate(g) = &chans[k] else { unreachable!() };
                let mut x = rho.clone();
                dm::unitary_like(b, &mut x, |v| apply_adjoint(b, v, g, m))?;
                x
            }
        };
        match &chans[k] {
            Channel::Gate(g) => {
                if g.arity() > 0 {
                    // W = G^dagger Lambda; grad_j = 2 Re tr(W (dG/dp_j) x)
                    let mut w = lam.clone();
                    on_columns(&mut w, d, |v| apply_adjoint(b, v, g, m))?;
                    for j in 0..g.arity() {
                        let mut dx = x.clone();
                        on_columns(&mut dx, d, |v| {
                            let r = derivative_vector(b, g, j, v, m)?;
                            v.copy_from_slice(&r);
                            Ok(())
                        })?;
                        grads[k].push(2.0 * re_trace_product(&w, &dx, d));
                    }
                }
                // Lambda <- G^dagger Lambda G
                dm::unitary_like(b, &mut lam, |v| apply_adjoint(b, v, g, m))?;
            }
            Channel::Loss { eta, mode } => {
                grads[k].push(loss_gradient(b, &x, &lam, *eta, *mode));
                lam = loss_adjoint(b, &lam, *eta, *mode);
            }
            Channel::Thermal { eta, nth, mode } => {
                let g = 1.0 + (1.0 - eta) * nth;
                lam = loss_adjoint(b, &amplifier_adjoint(b, &lam, g, *mode), eta / g, *mode);
            }
        }
        rho = x;
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

    /// the density-matrix adjoint against a fourth-order central difference, with loss (differentiated), thermal
    /// loss and truncated gates on the way
    #[test]
    fn density_matrix_against_differences() {
        let b = Basis::new(2, 6);
        let mut x0 = vec![C64::ZERO; b.patterns.len()];
        x0[b.index(&[1, 1])] = C64::ONE;
        // every differentiable parameter, in channel order: 12
        let build = |h: &[f64]| -> Vec<Channel> {
            vec![
                Channel::Gate(Gate::Beamsplitter { theta: 0.4 + h[0], phi: 0.3 + h[1], ma: 0, mb: 1 }),
                Channel::Loss { eta: 0.7 + h[2], mode: 0 },
                Channel::Gate(Gate::Kerr { kappa: 0.2 + h[3], mode: 1 }),
                Channel::Gate(Gate::Displacement { r: 0.3 + h[4], phi: 0.5 + h[5], mode: 1 }),
                Channel::Thermal { eta: 0.9, nth: 0.2, mode: 1 },
                Channel::Gate(Gate::Beamsplitter { theta: 0.8 + h[6], phi: h[7], ma: 0, mb: 1 }),
                Channel::Gate(Gate::Phaseshift { phi: 0.6 + h[8], mode: 0 }),
                Channel::Loss { eta: 0.85 + h[9], mode: 1 },
                Channel::Gate(Gate::Beamsplitter { theta: 0.5 + h[10], phi: 0.2 + h[11], ma: 0, mb: 1 }),
            ]
        };
        let terms = vec![
            Term::Probability { weight: 1.0, modes: vec![0], pattern: vec![1] },
            Term::Photons { weight: -0.4, modes: vec![1] },
        ];
        let value = |h: &[f64]| -> f64 {
            let (v, _) = adjoint_dm(&b, &x0, &build(h), &terms).unwrap();
            v[0] - 0.4 * v[1]
        };
        let (_, grads) = adjoint_dm(&b, &x0, &build(&[0.0; 12]), &terms).unwrap();
        let flat: Vec<f64> = grads.iter().flatten().copied().collect();
        assert_eq!(flat.len(), 12);
        let hh = 1e-3;
        for (p, &g) in flat.iter().enumerate() {
            let at = |s: f64| { let mut h = [0.0; 12]; h[p] = s; value(&h) };
            let fd = (-at(2.0 * hh) + 8.0 * at(hh) - 8.0 * at(-hh) + at(-2.0 * hh)) / (12.0 * hh);
            assert!((g - fd).abs() < 1e-9, "parameter {p}: adjoint {g} vs fd {fd}");
        }
    }
}
