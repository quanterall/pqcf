//! Spike S2 (AN 10): nalgebra or faer for the dense complex linear algebra of the Gaussian and torontonian
//! kernels: Hermitian eigendecomposition, complex Cholesky, complex LU (log|det|, inverse), sizes 2 to 200.
//!
//!   s2_linalg gen  > matrices.json      the test matrices (deterministic)
//!   s2_linalg run  > results.json       both crates' results, timings and the bitwise-repeat check
//! then compare.py checks the results against numpy.
use faer::linalg::solvers::DenseSolveCore;
use faer::{c64, Mat, Par, Side};
use nalgebra::{Complex, DMatrix};
use std::time::Instant;

const SIZES: [usize; 9] = [2, 4, 8, 16, 32, 64, 100, 128, 200];

type Run = fn(usize, &[(f64, f64)], &[(f64, f64)], &[(f64, f64)]) -> Out;

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / (1u64 << 53) as f64 - 0.5
    }
}

/// (H Hermitian, A Hermitian positive definite, G general), row-major (re, im)
#[allow(clippy::type_complexity)]
fn matrices(n: usize) -> (Vec<(f64, f64)>, Vec<(f64, f64)>, Vec<(f64, f64)>) {
    let mut r = Lcg(1000 + n as u64);
    let g: Vec<(f64, f64)> = (0..n * n).map(|_| (r.next(), r.next())).collect();
    let b: Vec<(f64, f64)> = (0..n * n).map(|_| (r.next(), r.next())).collect();
    let at = |m: &Vec<(f64, f64)>, i: usize, j: usize| m[i * n + j];
    let mut h = vec![(0.0, 0.0); n * n];
    let mut a = vec![(0.0, 0.0); n * n];
    for i in 0..n {
        for j in 0..n {
            // H = (B + B^dagger) / 2
            let (x, y) = at(&b, i, j);
            let (u, v) = at(&b, j, i);
            h[i * n + j] = ((x + u) / 2.0, (y - v) / 2.0);
            // A = B B^dagger + n I
            let mut s = (0.0, 0.0);
            for k in 0..n {
                let (p, q) = at(&b, i, k);
                let (s2, t2) = at(&b, j, k);
                s.0 += p * s2 + q * t2;
                s.1 += q * s2 - p * t2;
            }
            if i == j {
                s.0 += n as f64;
            }
            a[i * n + j] = s;
        }
    }
    (h, a, g)
}

fn json_mat(m: &[(f64, f64)]) -> String {
    let parts: Vec<String> = m.iter().map(|(a, b)| format!("[{a:e},{b:e}]")).collect();
    format!("[{}]", parts.join(","))
}

fn json_vec(v: &[f64]) -> String {
    let parts: Vec<String> = v.iter().map(|a| format!("{a:e}")).collect();
    format!("[{}]", parts.join(","))
}

fn na(m: &[(f64, f64)], n: usize) -> DMatrix<Complex<f64>> {
    DMatrix::from_fn(n, n, |i, j| Complex::new(m[i * n + j].0, m[i * n + j].1))
}

fn fa(m: &[(f64, f64)], n: usize) -> Mat<c64> {
    Mat::from_fn(n, n, |i, j| c64::new(m[i * n + j].0, m[i * n + j].1))
}

/// eigenvalues ascending, the eigen-residual max|H V - V L|, Cholesky L (row-major), log|det G|, inverse of G
struct Out {
    eig: Vec<f64>,
    eig_residual: f64,
    chol: Vec<(f64, f64)>,
    logdet: f64,
    inv: Vec<(f64, f64)>,
}

fn run_nalgebra(n: usize, h: &[(f64, f64)], a: &[(f64, f64)], g: &[(f64, f64)]) -> Out {
    let hm = na(h, n);
    let e = nalgebra::linalg::SymmetricEigen::new(hm.clone());
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&i, &j| e.eigenvalues[i].partial_cmp(&e.eigenvalues[j]).unwrap());
    let eig: Vec<f64> = idx.iter().map(|&i| e.eigenvalues[i]).collect();
    let lam = DMatrix::from_diagonal(&e.eigenvalues.map(|x| Complex::new(x, 0.0)));
    let res = &hm * &e.eigenvectors - &e.eigenvectors * lam;
    let eig_residual = res.iter().map(|z| z.norm()).fold(0.0, f64::max);
    let l = nalgebra::linalg::Cholesky::new(na(a, n)).expect("not positive definite").l();
    let chol = (0..n * n).map(|k| (l[(k / n, k % n)].re, l[(k / n, k % n)].im)).collect();
    let lu = na(g, n).lu();
    let logdet = (0..n).map(|i| lu.u()[(i, i)].norm().ln()).sum();
    let inv_m = lu.try_inverse().expect("singular");
    let inv = (0..n * n).map(|k| (inv_m[(k / n, k % n)].re, inv_m[(k / n, k % n)].im)).collect();
    Out { eig, eig_residual, chol, logdet, inv }
}

fn run_faer(n: usize, h: &[(f64, f64)], a: &[(f64, f64)], g: &[(f64, f64)]) -> Out {
    let hm = fa(h, n);
    let e = hm.self_adjoint_eigen(Side::Lower).expect("evd");
    let s = e.S();
    let u = e.U();
    let eig: Vec<f64> = (0..n).map(|i| s[i].re).collect();
    let mut eig_residual: f64 = 0.0;
    for i in 0..n {
        for j in 0..n {
            let mut acc = c64::new(0.0, 0.0);
            for k in 0..n {
                acc += hm[(i, k)] * u[(k, j)];
            }
            eig_residual = eig_residual.max((acc - u[(i, j)] * s[j]).norm());
        }
    }
    let llt = fa(a, n).llt(Side::Lower).expect("not positive definite");
    let l = llt.L();
    let chol = (0..n * n).map(|k| (l[(k / n, k % n)].re, l[(k / n, k % n)].im)).collect();
    let lu = fa(g, n).partial_piv_lu();
    let um = lu.U();
    let logdet = (0..n).map(|i| um[(i, i)].norm().ln()).sum();
    let inv_m = lu.inverse();
    let inv = (0..n * n).map(|k| (inv_m[(k / n, k % n)].re, inv_m[(k / n, k % n)].im)).collect();
    Out { eig, eig_residual, chol, logdet, inv }
}

fn bits(o: &Out) -> Vec<u64> {
    let mut v: Vec<u64> = o.eig.iter().map(|x| x.to_bits()).collect();
    v.push(o.eig_residual.to_bits());
    v.extend(o.chol.iter().flat_map(|(a, b)| [a.to_bits(), b.to_bits()]));
    v.push(o.logdet.to_bits());
    v.extend(o.inv.iter().flat_map(|(a, b)| [a.to_bits(), b.to_bits()]));
    v
}

/// median seconds of `reps` runs of f
fn time(reps: usize, mut f: impl FnMut()) -> f64 {
    let mut t: Vec<f64> = (0..reps)
        .map(|_| {
            let s = Instant::now();
            f();
            s.elapsed().as_secs_f64()
        })
        .collect();
    t.sort_by(|a, b| a.partial_cmp(b).unwrap());
    t[reps / 2]
}

fn main() {
    faer::set_global_parallelism(Par::Seq);
    let mode = std::env::args().nth(1).unwrap_or_default();
    let mut items = Vec::new();
    for &n in &SIZES {
        let (h, a, g) = matrices(n);
        if mode == "gen" {
            items.push(format!("{{\"n\":{n},\"H\":{},\"A\":{},\"G\":{}}}", json_mat(&h), json_mat(&a), json_mat(&g)));
            continue;
        }
        let reps = if n <= 32 { 201 } else if n <= 128 { 21 } else { 7 };
        let mut parts = Vec::new();
        let runs: [(&str, Run); 2] = [("nalgebra", run_nalgebra), ("faer", run_faer)];
        for (name, f) in runs {
            let o1 = f(n, &h, &a, &g);
            let o2 = f(n, &h, &a, &g);
            let repeat = bits(&o1) == bits(&o2);
            let (hn, an, gn) = (na(&h, n), na(&a, n), na(&g, n));
            let (hf, af, gf) = (fa(&h, n), fa(&a, n), fa(&g, n));
            let (t_eig, t_chol, t_lu, t_inv) = if name == "nalgebra" {
                (
                    time(reps, || {
                        std::hint::black_box(nalgebra::linalg::SymmetricEigen::new(hn.clone()));
                    }),
                    time(reps, || {
                        std::hint::black_box(nalgebra::linalg::Cholesky::new(an.clone()));
                    }),
                    time(reps, || {
                        std::hint::black_box(gn.clone().lu());
                    }),
                    time(reps, || {
                        std::hint::black_box(gn.clone().lu().try_inverse());
                    }),
                )
            } else {
                (
                    time(reps, || {
                        std::hint::black_box(hf.self_adjoint_eigen(Side::Lower).unwrap());
                    }),
                    time(reps, || {
                        std::hint::black_box(af.llt(Side::Lower).unwrap());
                    }),
                    time(reps, || {
                        std::hint::black_box(gf.partial_piv_lu());
                    }),
                    time(reps, || {
                        std::hint::black_box(gf.partial_piv_lu().inverse());
                    }),
                )
            };
            parts.push(format!(
                "\"{name}\":{{\"eig\":{},\"eig_residual\":{:e},\"chol\":{},\"logdet\":{:e},\"inv\":{},\"repeat_bitwise\":{repeat},\"seconds\":{{\"eig\":{t_eig:e},\"chol\":{t_chol:e},\"lu\":{t_lu:e},\"inverse\":{t_inv:e}}}}}",
                json_vec(&o1.eig),
                o1.eig_residual,
                json_mat(&o1.chol),
                o1.logdet,
                json_mat(&o1.inv)
            ));
        }
        items.push(format!("{{\"n\":{n},{}}}", parts.join(",")));
        eprintln!("n = {n} done");
    }
    println!("[{}]", items.join(",\n"));
}
