//! Spike S6 (AN 10): Clifford-Clifford against the exact distribution, 10^6 shots per case (chi-square in
//! test/oracle/s6_chi2.py), the AN v0.1 weight formula as a mutant that must fail, and the time per sample.
//!   cargo run --release --example s6 > s6.json
use nalgebra::{Complex, DMatrix};
use pq_nif::boson::{distribution, patterns, sample_boson, Algorithm};
use pq_nif::num::C64;
use pq_nif::perm;
use pq_nif::rng::Stream;
use std::time::Instant;

/// a Haar-random m x m unitary: QR of a complex Gaussian matrix, R's diagonal phases moved into Q
fn haar(m: usize, seed: u64) -> Vec<C64> {
    let mut s = Stream::new(seed, 0);
    let z = DMatrix::from_fn(m, m, |_, _| {
        let (a, b) = s.normal_pair();
        Complex::new(a / 2f64.sqrt(), b / 2f64.sqrt())
    });
    let qr = z.qr();
    let (q, r) = (qr.q(), qr.r());
    let mut u = vec![C64::ZERO; m * m];
    for j in 0..m {
        let d = r[(j, j)];
        let ph = d / d.norm();
        for i in 0..m {
            let x = q[(i, j)] * ph;
            u[i * m + j] = C64::new(x.re, x.im);
        }
    }
    u
}

fn never() -> bool {
    false
}

/// the AN v0.1 weights, w_j = sum_l |U[j][c_l]|^2 |p_l|^2 (distinguishable photons): must fail the test
fn mutant_shot(u: &[C64], m: usize, input: &[u32], s: &mut Stream) -> Vec<u32> {
    let mut cols: Vec<usize> = (0..m).flat_map(|i| std::iter::repeat_n(i, input[i] as usize)).collect();
    let n = cols.len();
    s.shuffle(&mut cols);
    let mut ys = Vec::new();
    for k in 1..=n {
        let r: Vec<f64> = (0..k)
            .map(|l| {
                if k == 1 {
                    return 1.0;
                }
                let sub: Vec<C64> = ys
                    .iter()
                    .flat_map(|&y| (0..k).filter(move |&c| c != l).map(move |c| (y, c)))
                    .map(|(y, c): (usize, usize)| u[y * m + cols[c]])
                    .collect();
                let p = perm::permanent(&sub, k - 1, 1, &never).unwrap().value;
                p.re * p.re + p.im * p.im
            })
            .collect();
        let w: Vec<f64> = (0..m).map(|j| (0..k).map(|l| u[j * m + cols[l]].abs().powi(2) * r[l]).sum()).collect();
        let total: f64 = w.iter().sum();
        let t = s.uniform() * total;
        let mut cum = 0.0;
        let mut pick = m - 1;
        for (j, x) in w.iter().enumerate() {
            cum += x;
            if cum > t {
                pick = j;
                break;
            }
        }
        ys.push(pick);
    }
    let mut out = vec![0u32; m];
    for y in ys {
        out[y] += 1;
    }
    out
}

fn main() {
    let cases: Vec<(usize, Vec<u32>)> = vec![
        (2, vec![1, 1]),
        (4, vec![1, 1, 1, 0]),
        (5, vec![1, 1, 1, 1, 0]),
        (6, vec![1, 1, 1, 1, 1, 1]),
        (8, vec![1, 1, 1, 1, 1, 1, 0, 0]),
        (4, vec![2, 1, 0, 0]),
        (5, vec![3, 0, 1, 0, 0]),
        (8, vec![2, 2, 1, 1, 0, 0, 0, 0]),
    ];
    let shots = 1_000_000;
    let mut items = Vec::new();
    for (k, (m, input)) in cases.iter().enumerate() {
        let u = haar(*m, 1000 + k as u64);
        let (pats, probs) = distribution(&u, *m, input, 1_000_000).unwrap();
        let index = |t: &[u32]| pats.iter().position(|p| p == t).unwrap();
        let t0 = Instant::now();
        let samples = sample_boson(&u, *m, input, shots, 2026, 0, Algorithm::CliffordClifford, &never).unwrap();
        let secs = t0.elapsed().as_secs_f64();
        let mut counts = vec![0u64; pats.len()];
        for t in samples.chunks(*m) {
            counts[index(t)] += 1;
        }
        let mut mcounts = vec![0u64; pats.len()];
        for i in 0..200_000u64 {
            let mut s = Stream::new(2026, i);
            mcounts[index(&mutant_shot(&u, *m, input, &mut s))] += 1;
        }
        eprintln!("case {k}: m = {m}, {} patterns, {secs:.1} s", pats.len());
        items.push(format!(
            "{{\"m\":{m},\"input\":{:?},\"probs\":{:?},\"counts\":{:?},\"mutant_counts\":{:?},\"seconds\":{secs}}}",
            input, probs, counts, mcounts
        ));
    }
    // time per sample at n = 12, 16, 20 photons in m = 2n modes
    let mut timing = Vec::new();
    for &(n, reps) in &[(12usize, 200usize), (16, 20), (20, 3)] {
        let m = 2 * n;
        let u = haar(m, 7 + n as u64);
        let mut input = vec![0u32; m];
        input[..n].fill(1);
        let t0 = Instant::now();
        sample_boson(&u, m, &input, reps, 1, 0, Algorithm::CliffordClifford, &never).unwrap();
        let per = t0.elapsed().as_secs_f64() / reps as f64;
        eprintln!("n = {n}: {per:.4} s per sample");
        timing.push(format!("{{\"n\":{n},\"m\":{m},\"seconds_per_sample\":{per}}}"));
    }
    let _ = patterns;
    println!("{{\"cases\":[{}],\"timing\":[{}]}}", items.join(",\n"), timing.join(","));
}
