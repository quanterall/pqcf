//! Gate plans (1.2; docs/decisions/B2-fock-gate-plans.md): a Fock-space gate prepared once and applied to many
//! vectors -- the 2D columns and rows of a density matrix (rho -> G rho G^dagger). The per-vector functions of
//! fockstate.rs recompute everything on every call: basis indices (O(D m) each), the beamsplitter's sector matrices,
//! the Clements mesh of a passive unitary, the truncated displacement and squeezing matrices. A plan holds them, and
//! `apply` performs exactly the per-vector function's arithmetic in the same order, so results are bit-identical.
//!
//! Number-conserving plans (beamsplitter, phase shift, Kerr, passive) also apply to selected photon-number sectors
//! only (`apply_sectors`): the stored blocks of a column in fockdmb.rs.
#![forbid(unsafe_code)] // NFR-04: kernels are safe Rust; unsafe only in the Rustler glue (lib.rs)
use crate::fock::{beamsplitter_sector, displacement_matrix, squeezing_matrix};
use crate::fockstate::{Basis, Error};
use crate::mesh;
use crate::num::C64;

/// a group of basis indices transformed together by one small matrix
struct Group {
    total: usize,
    idx: Vec<usize>,
    mat: usize,
}

enum Step {
    /// amplitude i times factor i (phase shift, Kerr)
    Diag(Vec<C64>),
    /// for each group: v[idx[r]] <- sum_c mats[mat][r * stride + c] v[idx[c]]
    Groups { groups: Vec<Group>, mats: Vec<Vec<C64>>, strides: Vec<usize> },
}

pub struct Plan {
    steps: Vec<Step>,
    /// maps every photon-number sector into itself
    pub conserving: bool,
}

fn totals(b: &Basis) -> Vec<usize> {
    b.patterns.iter().map(|p| p.iter().sum::<u32>() as usize).collect()
}

fn phase_step(b: &Basis, phi: f64, mode: usize) -> Step {
    // as fockstate::phaseshift
    Step::Diag(b.patterns.iter().map(|p| {
        let t = phi * p[mode] as f64;
        C64::new(t.cos(), t.sin())
    }).collect())
}

fn beamsplitter_step(b: &Basis, theta: f64, phi: f64, ma: usize, mb: usize) -> Result<Step, Error> {
    // as fockstate::beamsplitter
    let mats: Vec<Vec<C64>> = (0..=b.cutoff)
        .map(|n| beamsplitter_sector(n, theta, phi).map_err(|e| Error::Internal(format!("{e:?}"))))
        .collect::<Result<_, _>>()?;
    let strides = (0..=b.cutoff).map(|n| n + 1).collect();
    let mut groups = Vec::new();
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
        groups.push(Group { total: rep.iter().sum::<u32>() as usize, idx, mat: n });
    }
    Ok(Step::Groups { groups, mats, strides })
}

fn single_mode_step(b: &Basis, g: Vec<C64>, mode: usize) -> Step {
    // as fockstate::single_mode (without its norm bookkeeping)
    let d = b.cutoff + 1;
    let mut groups = Vec::new();
    for rep in b.patterns.iter().filter(|p| p[mode] == 0) {
        let others: usize = rep.iter().map(|&x| x as usize).sum();
        let top = b.cutoff - others;
        let mut q = rep.clone();
        let idx: Vec<usize> = (0..=top)
            .map(|n| {
                q[mode] = n as u32;
                b.index(&q)
            })
            .collect();
        groups.push(Group { total: others, idx, mat: 0 });
    }
    Step::Groups { groups, mats: vec![g], strides: vec![d] }
}

impl Plan {
    pub fn phaseshift(b: &Basis, phi: f64, mode: usize) -> Plan {
        Plan { steps: vec![phase_step(b, phi, mode)], conserving: true }
    }

    pub fn kerr(b: &Basis, kappa: f64, mode: usize) -> Plan {
        // as fockstate::kerr
        let f = b.patterns.iter().map(|p| {
            let n = p[mode] as f64;
            let t = kappa * n * n;
            C64::new(t.cos(), t.sin())
        }).collect();
        Plan { steps: vec![Step::Diag(f)], conserving: true }
    }

    pub fn beamsplitter(b: &Basis, theta: f64, phi: f64, ma: usize, mb: usize) -> Result<Plan, Error> {
        Ok(Plan { steps: vec![beamsplitter_step(b, theta, phi, ma, mb)?], conserving: true })
    }

    /// as fockstate::passive: the Clements mesh's beamsplitters, then its phases
    pub fn passive(b: &Basis, u: &[C64], modes: &[usize]) -> Result<Plan, Error> {
        let (bs, phases) = mesh::decompose(u, modes.len());
        let mut steps = Vec::with_capacity(bs.len() + phases.len());
        for (m, t, p) in bs {
            steps.push(beamsplitter_step(b, t, p, modes[m], modes[m + 1])?);
        }
        for (j, &a) in phases.iter().enumerate() {
            steps.push(phase_step(b, a, modes[j]));
        }
        Ok(Plan { steps, conserving: true })
    }

    pub fn displacement(b: &Basis, r: f64, phi: f64, mode: usize) -> Plan {
        let g = displacement_matrix(C64::new(r * phi.cos(), r * phi.sin()), b.cutoff);
        Plan { steps: vec![single_mode_step(b, g, mode)], conserving: false }
    }

    pub fn squeezing(b: &Basis, r: f64, phi: f64, mode: usize) -> Result<Plan, Error> {
        let g = squeezing_matrix(r, phi, b.cutoff).ok_or(Error::TooLarge)?;
        Ok(Plan { steps: vec![single_mode_step(b, g, mode)], conserving: false })
    }

    /// a single-mode matrix G ((cutoff+1) x (cutoff+1)) on `mode`, as fockstate::single_mode
    pub fn single_mode(b: &Basis, g: Vec<C64>, mode: usize) -> Plan {
        Plan { steps: vec![single_mode_step(b, g, mode)], conserving: false }
    }

    /// v <- G v
    pub fn apply(&self, v: &mut [C64]) {
        for s in &self.steps {
            step(s, v, None);
        }
    }

    /// v <- G v on the sectors marked in `sectors` (indexed by total photon number); the others are left alone,
    /// which is G v when they are zero and G conserves the photon number
    pub fn apply_sectors(&self, v: &mut [C64], sectors: &[bool]) {
        debug_assert!(self.conserving);
        for s in &self.steps {
            step(s, v, Some(sectors));
        }
    }
}

fn step(s: &Step, v: &mut [C64], sectors: Option<&[bool]>) {
    match s {
        Step::Diag(f) => {
            for (a, m) in v.iter_mut().zip(f) {
                *a = a.mul(*m);
            }
        }
        Step::Groups { groups, mats, strides } => {
            let mut vals: Vec<C64> = Vec::new();
            for g in groups {
                if sectors.is_some_and(|m| !m[g.total]) {
                    continue;
                }
                let (u, w) = (&mats[g.mat], strides[g.mat]);
                vals.clear();
                vals.extend(g.idx.iter().map(|&i| v[i]));
                for (r, &i) in g.idx.iter().enumerate() {
                    let mut acc = C64::ZERO;
                    for (c, vc) in vals.iter().enumerate() {
                        acc = acc.add(u[r * w + c].mul(*vc));
                    }
                    v[i] = acc;
                }
            }
        }
    }
}

/// the totals of the basis, for callers that build sector masks
pub fn sector_of_index(b: &Basis) -> Vec<usize> {
    totals(b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fockstate;

    fn vector(b: &Basis, seed: u64) -> Vec<C64> {
        let mut s = crate::rng::Stream::new(seed, 0);
        (0..b.patterns.len()).map(|_| C64::new(s.uniform() - 0.5, s.uniform() - 0.5)).collect()
    }

    /// every plan bit-identical to its per-vector function
    #[test]
    #[allow(clippy::type_complexity)]
    fn plans_equal_the_per_vector_functions() {
        let b = Basis::new(3, 6);
        let u: Vec<C64> = {
            let (c, s) = (0.6f64, 0.8f64);
            vec![C64::new(c, 0.0), C64::new(0.0, s), C64::ZERO, C64::new(0.0, s), C64::new(c, 0.0), C64::ZERO, C64::ZERO, C64::ZERO, C64::new(0.0, 1.0)]
        };
        let cases: Vec<(Plan, Box<dyn Fn(&mut [C64])>)> = vec![
            (Plan::beamsplitter(&b, 0.7, 0.3, 0, 2).unwrap(), Box::new(|v| fockstate::beamsplitter(&b, v, 0.7, 0.3, 0, 2).unwrap())),
            (Plan::phaseshift(&b, 0.4, 1), Box::new(|v| fockstate::phaseshift(&b, v, 0.4, 1))),
            (Plan::kerr(&b, 0.3, 2), Box::new(|v| fockstate::kerr(&b, v, 0.3, 2))),
            (Plan::passive(&b, &u, &[0, 1, 2]).unwrap(), Box::new(|v| fockstate::passive(&b, v, &u, &[0, 1, 2]).unwrap())),
            (Plan::displacement(&b, 0.4, 0.2, 1), Box::new(|v| { fockstate::displacement(&b, v, 0.4, 0.2, 1); })),
            (Plan::squeezing(&b, 0.3, -0.5, 0).unwrap(), Box::new(|v| { fockstate::squeezing(&b, v, 0.3, -0.5, 0).unwrap(); })),
        ];
        for (k, (plan, f)) in cases.iter().enumerate() {
            let mut a = vector(&b, k as u64);
            let mut c = a.clone();
            plan.apply(&mut a);
            f(&mut c);
            assert_eq!(a, c, "case {k}");
        }
    }

    /// on a vector supported in some sectors, applying to those sectors only is the full application
    #[test]
    fn sectors() {
        let b = Basis::new(3, 6);
        let t = sector_of_index(&b);
        let plan = Plan::beamsplitter(&b, 0.9, 0.1, 1, 2).unwrap();
        let mask: Vec<bool> = (0..=6).map(|n| n == 2 || n == 5).collect();
        let mut v = vector(&b, 9);
        for (x, &n) in v.iter_mut().zip(&t) {
            if !mask[n] {
                *x = C64::ZERO;
            }
        }
        let (mut full, mut part) = (v.clone(), v.clone());
        plan.apply(&mut full);
        plan.apply_sectors(&mut part, &mask);
        assert!(full.iter().zip(&part).all(|(x, y)| x.re == y.re && x.im == y.im));
    }

    /// the speed-up on density matrices: cargo test --release --lib fockplan::tests::speed -- --ignored --nocapture
    #[test]
    #[ignore]
    #[allow(clippy::type_complexity)]
    fn speed() {
        use crate::fockdmb::Blocks;
        use crate::fockstate::dm;
        use std::time::Instant;
        let time = |f: &mut dyn FnMut()| { let t = Instant::now(); f(); t.elapsed().as_secs_f64() };
        // dense: 3 modes, cutoff 10 (D = 286), a random state
        let b = Basis::new(3, 10);
        let psi = vector(&b, 1);
        let rho0 = dm::from_pure(&psi);
        let u: Vec<C64> = vec![C64::new(0.6, 0.0), C64::new(0.0, 0.8), C64::ZERO, C64::new(0.0, 0.8), C64::new(0.6, 0.0), C64::ZERO,
                               C64::ZERO, C64::ZERO, C64::ONE];
        let gates: Vec<(&str, Box<dyn Fn(&mut [C64])>, Plan)> = vec![
            ("beamsplitter", Box::new(|v| fockstate::beamsplitter(&b, v, 0.7, 0.3, 0, 1).unwrap()), Plan::beamsplitter(&b, 0.7, 0.3, 0, 1).unwrap()),
            ("passive 3x3", Box::new(|v| fockstate::passive(&b, v, &u, &[0, 1, 2]).unwrap()), Plan::passive(&b, &u, &[0, 1, 2]).unwrap()),
            ("displacement", Box::new(|v| { fockstate::displacement(&b, v, 0.4, 0.2, 2); }), Plan::displacement(&b, 0.4, 0.2, 2)),
            ("kerr", Box::new(|v| fockstate::kerr(&b, v, 0.3, 1)), Plan::kerr(&b, 0.3, 1)),
        ];
        for (name, f, plan) in &gates {
            let (mut r1, mut r2) = (rho0.clone(), rho0.clone());
            let t_old = time(&mut || { dm::unitary_like(&b, &mut r1, |v| { f(v); Ok(()) }).unwrap(); });
            let t_new = time(&mut || { dm::unitary_like(&b, &mut r2, |v| { plan.apply(v); Ok(()) }).unwrap(); });
            assert_eq!(r1, r2);
            eprintln!("dense D = 286, {name}: per-vector {t_old:.3} s, plan {t_new:.3} s, {:.1}x", t_old / t_new);
        }
        // blocks: 4 modes, cutoff 16 (D = 4845), the Fock input |4,4,4,4> after loss: diagonal blocks only
        let b = Basis::new(4, 16);
        let mut x = vec![C64::ZERO; b.patterns.len()];
        x[b.index(&[4, 4, 4, 4])] = C64::ONE;
        let blk = Blocks::from_pure(&b, &x).loss(&b, 0.8, 0);
        let plan = Plan::beamsplitter(&b, 0.7, 0.1, 0, 1).unwrap();
        let (mut a, mut c) = (Blocks::from_dense(&b, &blk.to_dense(b.patterns.len())), Blocks::from_dense(&b, &blk.to_dense(b.patterns.len())));
        let t_old = time(&mut || { a.unitary_like(&b, |v| fockstate::beamsplitter(&b, v, 0.7, 0.1, 0, 1)).unwrap(); });
        let t_new = time(&mut || { c.unitary_plan(&b, &plan).unwrap(); });
        let d = b.patterns.len();
        assert!(a.to_dense(d).iter().zip(&c.to_dense(d)).all(|(p, q)| p.re == q.re && p.im == q.im));
        eprintln!("blocks D = 4845 (diagonal blocks), beamsplitter: per-vector {t_old:.2} s, plan {t_new:.3} s, {:.0}x", t_old / t_new);
    }
}
