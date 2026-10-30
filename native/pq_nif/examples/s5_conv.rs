//! S5: convergence of the padded eigendecomposition by doubling, at the points outside the Bargmann region.
use pq_nif::fock::{bargmann_accurate, squeezing_matrix_eigen_at};
// the stopping rule of fock::squeezing_matrix, repeated here to print the trace
fn main() {
    for &c in &[40usize, 80] {
        for &r in &[0.25f64, 0.5, 0.75, 1.0, 1.5, 2.0] {
            if bargmann_accurate(r, c) {
                continue;
            }
            let tail = (46.0 / -(r.tanh().ln()) / 2.0).ceil() as usize;
            let mut size = c / 2 + 1 + tail + 10;
            let mut prev = squeezing_matrix_eigen_at(r, 0.7, c, size);
            let mut line = format!("c={c} r={r}: {size}");
            while 2 * size <= 4096 {
                size *= 2;
                let next = squeezing_matrix_eigen_at(r, 0.7, c, size);
                let ch = prev.iter().zip(&next).map(|(a, b)| (a.re - b.re).abs() + (a.im - b.im).abs()).fold(0.0, f64::max);
                line += &format!(" -> {size}: {ch:.1e}");
                if ch < f64::max(2e-13, 1e-15 * size as f64) {
                    break;
                }
                prev = next;
            }
            println!("{line}");
        }
    }
}
