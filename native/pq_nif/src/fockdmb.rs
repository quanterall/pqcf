//! Density matrices stored by photon-number blocks (Phase 4, FS 15 item 6; Decision 3 allows it internally;
//! docs/decisions/B1-dm-blocks.md).
//!
//! The basis of FS 3.4 orders states by total photon number, so sector n is the contiguous index range
//! [off_n, off_n + d_n), d_n = C(n + m - 1, m - 1), and rho is a grid of blocks rho[n][n'] (d_n x d_n'). A block that is
//! zero is not stored. Number-conserving gates and loss keep n - n' fixed, so a Fock or number-diagonal input keeps
//! only the diagonal blocks: sum_n d_n^2 entries instead of D^2.
//!
//! Every operation is the dense kernel's (fockstate::dm) operation in the same order, entry by entry, so the
//! stored entries are bit-identical to the dense matrix's (the tests compare them); only zero blocks are skipped.
#![forbid(unsafe_code)] // NFR-04: kernels are safe Rust; unsafe only in the Rustler glue (lib.rs)
use crate::fockstate::{Basis, Error};
use crate::num::{Sum, C64};

pub struct Blocks {
    pub sizes: Vec<usize>,
    pub offsets: Vec<usize>,
    /// row-major grid, index n * (cutoff + 1) + n'
    pub blocks: Vec<Option<Vec<C64>>>,
}

fn conj(z: C64) -> C64 {
    C64::new(z.re, -z.im)
}

impl Blocks {
    pub fn empty(b: &Basis) -> Blocks {
        let k = b.cutoff + 1;
        let mut sizes = vec![0usize; k];
        for p in &b.patterns {
            sizes[p.iter().sum::<u32>() as usize] += 1;
        }
        let mut offsets = vec![0usize; k];
        for n in 1..k {
            offsets[n] = offsets[n - 1] + sizes[n - 1];
        }
        Blocks { sizes, offsets, blocks: vec![None; k * k] }
    }

    fn k(&self) -> usize {
        self.sizes.len()
    }

    pub fn entries(&self) -> usize {
        self.blocks.iter().flatten().map(|v| v.len()).sum()
    }

    fn sector_of(&self, i: usize) -> usize {
        // the last sector whose offset is at most i
        match self.offsets.binary_search(&i) {
            Ok(mut n) => {
                while n + 1 < self.k() && self.offsets[n + 1] == i {
                    n += 1;
                }
                n
            }
            Err(n) => n - 1,
        }
    }

    pub fn get(&self, i: usize, j: usize) -> C64 {
        let (a, c) = (self.sector_of(i), self.sector_of(j));
        match &self.blocks[a * self.k() + c] {
            Some(v) => v[(i - self.offsets[a]) * self.sizes[c] + j - self.offsets[c]],
            None => C64::ZERO,
        }
    }

    /// entry (i, j) of block (a, c) set, the block created when the value is not zero
    fn set(&mut self, a: usize, c: usize, i: usize, j: usize, z: C64) {
        let k = self.k();
        let (da, dc) = (self.sizes[a], self.sizes[c]);
        let blk = &mut self.blocks[a * k + c];
        if blk.is_none() {
            if z.re == 0.0 && z.im == 0.0 {
                return;
            }
            *blk = Some(vec![C64::ZERO; da * dc]);
        }
        blk.as_mut().unwrap()[(i - self.offsets[a]) * dc + j - self.offsets[c]] = z;
    }

    fn add(&mut self, a: usize, c: usize, i: usize, j: usize, z: C64) {
        let cur = self.get(i, j);
        self.set(a, c, i, j, cur.add(z));
    }

    /// blocks that became entirely zero are dropped
    fn prune(&mut self) {
        for blk in self.blocks.iter_mut() {
            if blk.as_ref().is_some_and(|v| v.iter().all(|z| z.re == 0.0 && z.im == 0.0)) {
                *blk = None;
            }
        }
    }

    pub fn from_pure(b: &Basis, amps: &[C64]) -> Blocks {
        let mut out = Blocks::empty(b);
        let k = out.k();
        for a in 0..k {
            for c in 0..k {
                for i in out.offsets[a]..out.offsets[a] + out.sizes[a] {
                    for j in out.offsets[c]..out.offsets[c] + out.sizes[c] {
                        out.set(a, c, i, j, amps[i].mul(conj(amps[j])));
                    }
                }
            }
        }
        out
    }

    pub fn from_dense(b: &Basis, rho: &[C64]) -> Blocks {
        let d = b.patterns.len();
        let mut out = Blocks::empty(b);
        let k = out.k();
        for a in 0..k {
            for c in 0..k {
                for i in out.offsets[a]..out.offsets[a] + out.sizes[a] {
                    for j in out.offsets[c]..out.offsets[c] + out.sizes[c] {
                        out.set(a, c, i, j, rho[i * d + j]);
                    }
                }
            }
        }
        out
    }

    pub fn to_dense(&self, d: usize) -> Vec<C64> {
        let mut rho = vec![C64::ZERO; d * d];
        let k = self.k();
        for a in 0..k {
            for c in 0..k {
                if let Some(v) = &self.blocks[a * k + c] {
                    for (x, z) in v.iter().enumerate() {
                        let (i, j) = (self.offsets[a] + x / self.sizes[c], self.offsets[c] + x % self.sizes[c]);
                        rho[i * d + j] = *z;
                    }
                }
            }
        }
        rho
    }

    pub fn trace(&self) -> f64 {
        let k = self.k();
        let d: usize = self.sizes.iter().sum();
        let mut s = Sum::new(d as u64);
        for a in 0..k {
            for x in 0..self.sizes[a] {
                s.add(match &self.blocks[a * k + a] {
                    Some(v) => v[x * self.sizes[a] + x].re,
                    None => 0.0,
                });
            }
        }
        s.value()
    }

    /// rho <- G rho G^dagger, as fockstate::dm::unitary_like: G on every column, then on every conjugated row; only
    /// the columns (rows) of sectors with a stored block are visited, the others being zero
    pub fn unitary_like(&mut self, b: &Basis, mut apply: impl FnMut(&mut [C64]) -> Result<(), Error>) -> Result<f64, Error> {
        let d = b.patterns.len();
        let k = self.k();
        let before = self.trace();
        let mut col = vec![C64::ZERO; d];
        for c in 0..k {
            if (0..k).all(|a| self.blocks[a * k + c].is_none()) {
                continue;
            }
            for j in self.offsets[c]..self.offsets[c] + self.sizes[c] {
                for a in 0..k {
                    let (o, da, dc, jj) = (self.offsets[a], self.sizes[a], self.sizes[c], j - self.offsets[c]);
                    match &self.blocks[a * k + c] {
                        Some(v) => (0..da).for_each(|x| col[o + x] = v[x * dc + jj]),
                        None => col[o..o + da].fill(C64::ZERO),
                    }
                }
                apply(&mut col)?;
                for a in 0..k {
                    let o = self.offsets[a];
                    for (x, &z) in col[o..o + self.sizes[a]].iter().enumerate() {
                        self.set(a, c, o + x, j, z);
                    }
                }
            }
        }
        for a in 0..k {
            if (0..k).all(|c| self.blocks[a * k + c].is_none()) {
                continue;
            }
            for i in self.offsets[a]..self.offsets[a] + self.sizes[a] {
                for c in 0..k {
                    let (o, dc, ii) = (self.offsets[c], self.sizes[c], i - self.offsets[a]);
                    match &self.blocks[a * k + c] {
                        Some(v) => (0..dc).for_each(|x| col[o + x] = conj(v[ii * dc + x])),
                        None => col[o..o + dc].fill(C64::ZERO),
                    }
                }
                apply(&mut col)?;
                for c in 0..k {
                    let o = self.offsets[c];
                    for (x, &z) in col[o..o + self.sizes[c]].iter().enumerate() {
                        self.set(a, c, i, o + x, conj(z));
                    }
                }
            }
        }
        self.prune();
        Ok(before - self.trace())
    }

    /// a Kraus-type map entry by entry, in the dense kernels' order (rows i, columns j, then the branch l): each
    /// stored entry (i, j) sends weight(np, nq, l) times itself to the patterns `target(p, l)`, `target(q, l)`
    fn map_entries(&self, b: &Basis, mode: usize, branches: impl Fn(u32, u32) -> Vec<(i64, f64)>) -> Blocks {
        let k = self.k();
        let mut out = Blocks { sizes: self.sizes.clone(), offsets: self.offsets.clone(), blocks: vec![None; k * k] };
        // the dense order: rows i ascending, then columns j ascending (sector by sector), stored blocks only
        for a in 0..k {
            for i in self.offsets[a]..self.offsets[a] + self.sizes[a] {
            for c in 0..k {
                let Some(v) = &self.blocks[a * k + c] else { continue };
            for j in self.offsets[c]..self.offsets[c] + self.sizes[c] {
                let z = v[(i - self.offsets[a]) * self.sizes[c] + j - self.offsets[c]];
                if z.re == 0.0 && z.im == 0.0 {
                    continue;
                }
                let (p, q) = (&b.patterns[i], &b.patterns[j]);
                let mut pp = p.clone();
                let mut qq = q.clone();
                for (shift, w) in branches(p[mode], q[mode]) {
                    pp[mode] = (p[mode] as i64 + shift) as u32;
                    qq[mode] = (q[mode] as i64 + shift) as u32;
                    let (ti, tj) = (b.index(&pp), b.index(&qq));
                    let (ta, tc) = ((a as i64 + shift) as usize, (c as i64 + shift) as usize);
                    out.add(ta, tc, ti, tj, z.scale(w));
                }
            }
            }
            }
        }
        out
    }

    /// FS A.11 loss, as fockstate::dm::loss
    pub fn loss(&self, b: &Basis, eta: f64, mode: usize) -> Blocks {
        let binom = |n: u32, k: u32| -> f64 { (0..k).fold(1.0, |acc, i| acc * (n - i) as f64 / (i + 1) as f64) };
        self.map_entries(b, mode, |np, nq| {
            (0..=np.min(nq))
                .map(|l| {
                    let w = (binom(np, l) * binom(nq, l)).sqrt() * eta.powf((np + nq - 2 * l) as f64 / 2.0) * (1.0 - eta).powi(l as i32);
                    (-(l as i64), w)
                })
                .collect()
        })
    }

    /// the amplifier of fockstate::dm::amplifier, room limited by the cutoff on both patterns
    fn amplifier(&self, b: &Basis, g: f64, mode: usize) -> Blocks {
        let binom = |n: u32, k: u32| -> f64 { (0..k).fold(1.0, |acc, i| acc * (n + k - i) as f64 / (i + 1) as f64) };
        let ratio = (g - 1.0) / g;
        let k = self.k();
        let mut out = Blocks { sizes: self.sizes.clone(), offsets: self.offsets.clone(), blocks: vec![None; k * k] };
        let total: Vec<u32> = b.patterns.iter().map(|p| p.iter().sum()).collect();
        for a in 0..k {
            for i in self.offsets[a]..self.offsets[a] + self.sizes[a] {
            for c in 0..k {
                let Some(v) = &self.blocks[a * k + c] else { continue };
            for j in self.offsets[c]..self.offsets[c] + self.sizes[c] {
                let z = v[(i - self.offsets[a]) * self.sizes[c] + j - self.offsets[c]];
                if z.re == 0.0 && z.im == 0.0 {
                    continue;
                }
                let (p, q) = (&b.patterns[i], &b.patterns[j]);
                let (np, nq) = (p[mode], q[mode]);
                let room = (b.cutoff as u32).saturating_sub(total[i].max(total[j]));
                let base = g.powf(-((np + nq + 2) as f64) / 2.0);
                let mut pp = p.clone();
                let mut qq = q.clone();
                for kk in 0..=room {
                    if kk > 0 && ratio == 0.0 {
                        break;
                    }
                    pp[mode] = np + kk;
                    qq[mode] = nq + kk;
                    let w = (binom(np, kk) * binom(nq, kk)).sqrt() * base * ratio.powi(kk as i32);
                    out.add(a + kk as usize, c + kk as usize, b.index(&pp), b.index(&qq), z.scale(w));
                }
            }
            }
            }
        }
        out
    }

    /// thermal loss as fockstate::dm::thermal_loss; returns the trace deficit
    pub fn thermal_loss(&self, b: &Basis, eta: f64, nth: f64, mode: usize) -> (Blocks, f64) {
        let g = 1.0 + (1.0 - eta) * nth;
        let lost = self.loss(b, eta / g, mode);
        let before = lost.trace();
        let out = lost.amplifier(b, g, mode);
        let after = out.trace();
        (out, before - after)
    }

    /// the diagonal, for the PNR probabilities (fockstate::dm::pnr_probabilities)
    pub fn diagonal(&self) -> Vec<f64> {
        let k = self.k();
        let mut out = Vec::new();
        for a in 0..k {
            for x in 0..self.sizes[a] {
                out.push(match &self.blocks[a * k + a] {
                    Some(v) => v[x * self.sizes[a] + x].re,
                    None => 0.0,
                });
            }
        }
        out
    }

    /// FR-M07 on blocks, as the dense branch of fock_measure_pnr: entries kept when both patterns match the outcome,
    /// divided by the outcome's probability
    pub fn measure(&mut self, b: &Basis, modes: &[usize], outcome: &[u32]) -> f64 {
        let keep: Vec<bool> = b.patterns.iter().map(|p| modes.iter().zip(outcome).all(|(&m, &o)| p[m] == o)).collect();
        let diag = self.diagonal();
        let prob = (0..diag.len()).filter(|&i| keep[i]).map(|i| diag[i]).sum::<f64>();
        let k = self.k();
        for a in 0..k {
            for c in 0..k {
                let (oa, oc, dc) = (self.offsets[a], self.offsets[c], self.sizes[c]);
                if let Some(v) = &mut self.blocks[a * k + c] {
                    for (x, z) in v.iter_mut().enumerate() {
                        let (i, j) = (oa + x / dc, oc + x % dc);
                        *z = if keep[i] && keep[j] && prob > 0.0 { z.scale(1.0 / prob) } else { C64::ZERO };
                    }
                }
            }
        }
        self.prune();
        prob
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fockstate::{self, dm};

    fn state(b: &Basis, occ: &[u32]) -> Vec<C64> {
        let mut v = vec![C64::ZERO; b.patterns.len()];
        v[b.index(occ)] = C64::ONE;
        v
    }

    /// every operation bit-identical to the dense kernel, and the Fock input stays block-diagonal
    #[test]
    #[allow(clippy::type_complexity)]
    fn equal_to_the_dense_kernels() {
        let b = Basis::new(3, 5);
        let d = b.patterns.len();
        let psi = state(&b, &[2, 1, 0]);
        let mut dense = dm::from_pure(&psi);
        let mut blk = Blocks::from_pure(&b, &psi);
        assert_eq!(blk.to_dense(d), dense);
        let steps: Vec<Box<dyn Fn(&mut [C64]) -> Result<(), Error>>> = vec![
            Box::new(|v| fockstate::beamsplitter(&b, v, 0.7, 0.3, 0, 1)),
            Box::new(|v| { fockstate::kerr(&b, v, 0.4, 1); Ok(()) }),
            Box::new(|v| fockstate::beamsplitter(&b, v, 1.1, -0.2, 1, 2)),
        ];
        for f in &steps {
            dm::unitary_like(&b, &mut dense, f).unwrap();
            blk.unitary_like(&b, f).unwrap();
            assert_eq!(blk.to_dense(d), dense);
        }
        dense = dm::loss(&b, &dense, 0.6, 0);
        blk = blk.loss(&b, 0.6, 0);
        assert_eq!(blk.to_dense(d), dense);
        // only the diagonal blocks are stored: n - n' is conserved
        let k = b.cutoff + 1;
        assert!((0..k).all(|a| (0..k).all(|c| a == c || blk.blocks[a * k + c].is_none())));
        assert!(blk.entries() < d * d / 3);
        // thermal loss and a displacement fill other blocks; still identical
        let (dt, _) = dm::thermal_loss(&b, &dense, 0.8, 0.3, 2);
        let (bt, _) = blk.thermal_loss(&b, 0.8, 0.3, 2);
        assert_eq!(bt.to_dense(d), dt);
        let (mut dd, mut bd) = (dt, bt);
        dm::unitary_like(&b, &mut dd, |v| { fockstate::displacement(&b, v, 0.3, 0.5, 1); Ok(()) }).unwrap();
        bd.unitary_like(&b, |v| { fockstate::displacement(&b, v, 0.3, 0.5, 1); Ok(()) }).unwrap();
        assert_eq!(bd.to_dense(d), dd);
        assert_eq!(bd.trace(), dm::trace(&dd, d));
        // a measurement
        let prob = bd.measure(&b, &[0], &[1]);
        assert!(prob > 0.0 && prob < 1.0);
        assert!((bd.trace() - 1.0).abs() < 1e-12);
    }
}
