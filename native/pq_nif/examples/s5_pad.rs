//! S5: the padded eigendecomposition at C = 80 against the padding.
use pq_nif::fock::squeezing_matrix_eigen_at as squeezing_matrix_eigen;
fn main() {
    let mut out = Vec::new();
    for &r in &[0.5f64, 0.75, 1.0, 2.0, 3.0] {
        for &size in &[100usize, 200, 400, 800, 1600] {
            let m = squeezing_matrix_eigen(r, 0.7, 80, size);
            let z = m[80 * 81 + 80];
            let w = m[60 * 81 + 40];
            out.push(format!("{{\"r\":{r},\"size\":{size},\"s80\":[{:e},{:e}],\"s60_40\":[{:e},{:e}]}}", z.re, z.im, w.re, w.im));
        }
    }
    println!("[{}]", out.join(",\n"));
}
