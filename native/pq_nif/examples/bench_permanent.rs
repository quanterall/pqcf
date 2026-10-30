//! NFR-01 benchmark (spike S1): single-thread permanent for n = 16, 18, 20, median of repeated runs, on the
//! matrix A[i][j] = sin(i n + j) + i cos(i j + 0.5) (1-based), the same matrix bench/thewalrus_permanent.py uses.
use pq_nif::num::C64;
use pq_nif::perm::permanent;
use std::time::Instant;

fn main() {
    let never = || false;
    let mut out = Vec::new();
    for &(n, reps) in &[(16usize, 31usize), (18, 15), (20, 7)] {
        let a: Vec<C64> = (1..=n)
            .flat_map(|i| (1..=n).map(move |j| C64::new(((i * n + j) as f64).sin(), ((i * j) as f64 + 0.5).cos())))
            .collect();
        let _ = permanent(&a, n, 1, &never).unwrap();
        let mut t: Vec<f64> = (0..reps)
            .map(|_| {
                let s = Instant::now();
                std::hint::black_box(permanent(&a, n, 1, &never).unwrap());
                s.elapsed().as_secs_f64()
            })
            .collect();
        t.sort_by(|x, y| x.partial_cmp(y).unwrap());
        let v = permanent(&a, n, 1, &never).unwrap().value;
        out.push(format!("{{\"n\":{n},\"median_seconds\":{:e},\"runs\":{reps},\"value\":[{:e},{:e}]}}", t[reps / 2], v.re, v.im));
    }
    println!("[{}]", out.join(","));
}
