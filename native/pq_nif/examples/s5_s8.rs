//! Spikes S5 and S8 (AN 10): prints the kernels' matrices as JSON for test/oracle/s5_s8_compare.py.
//!   cargo run --release --example s5_s8 > s5_s8.json
//! The AN 7.4 column recurrences and the float sum are kept here, for the record: they are not used.
use pq_nif::fock::{beamsplitter_sector, squeezing_matrix as squeezing_final, displacement_matrix as displacement_matrix_laguerre, squeezing_matrix_bargmann,
                   squeezing_matrix_eigen_at as squeezing_matrix_eigen};
use pq_nif::num::C64;

fn ln_fact(k: usize) -> f64 {
    (2..=k).map(|x| (x as f64).ln()).sum()
}

/// <m|D(alpha)|n> for m, n = 0..=c, row-major (c+1) x (c+1). D(alpha) = exp(alpha a^dagger - conj(alpha) a).
///   D[m][0]   = exp(-|alpha|^2 / 2) alpha^m / sqrt(m!)
///   D[m][n+1] = (sqrt(m) D[m-1][n] - conj(alpha) D[m][n]) / sqrt(n+1)
fn displacement_matrix(alpha: C64, c: usize) -> Vec<C64> {
    let d = c + 1;
    let mut m = vec![C64::ZERO; d * d];
    let r = alpha.abs();
    let arg = alpha.im.atan2(alpha.re);
    for row in 0..d {
        m[row * d] = if r == 0.0 {
            if row == 0 { C64::ONE } else { C64::ZERO }
        } else {
            let mag = (-r * r / 2.0 + row as f64 * r.ln() - 0.5 * ln_fact(row)).exp();
            let ph = row as f64 * arg;
            C64::new(mag * ph.cos(), mag * ph.sin())
        };
    }
    let ac = C64::new(-alpha.re, alpha.im); // -conj(alpha)
    for col in 0..c {
        let inv = 1.0 / ((col + 1) as f64).sqrt();
        for row in 0..d {
            let up = if row > 0 { m[(row - 1) * d + col].scale((row as f64).sqrt()) } else { C64::ZERO };
            m[row * d + col + 1] = up.add(ac.mul(m[row * d + col])).scale(inv);
        }
    }
    m
}

/// <m|S(z)|n> for m, n = 0..=c, row-major (c+1) x (c+1). S(z) = exp((conj(z) a^2 - z a^dagger^2) / 2),
/// z = r e^(i phi). Built on rows 0..=2c, since column n needs rows up to 2c - n.
///   S[2k][0]  = (cosh r)^-1/2 (-e^(i phi) tanh r)^k sqrt((2k)!) / (2^k k!)
///   S[m][n+1] = (sqrt(m) cosh r S[m-1][n] + sqrt(m+1) e^(-i phi) sinh r S[m+1][n]) / sqrt(n+1)
fn squeezing_matrix(r: f64, phi: f64, c: usize) -> Vec<C64> {
    let rows = 2 * c + 1;
    let d = c + 1;
    let mut m = vec![C64::ZERO; rows * d];
    let (ch, sh) = (r.cosh(), r.sinh());
    let th = r.tanh();
    // ln cosh r, stable for large r
    let ln_cosh = r.abs() + (-2.0 * r.abs()).exp().ln_1p() - std::f64::consts::LN_2;
    for k in 0..=c {
        let row = 2 * k;
        if row >= rows {
            break;
        }
        m[row * d] = if th == 0.0 {
            if k == 0 { C64::ONE } else { C64::ZERO }
        } else {
            let mag = (-0.5 * ln_cosh + k as f64 * th.abs().ln() + 0.5 * ln_fact(2 * k)
                - k as f64 * std::f64::consts::LN_2
                - ln_fact(k))
            .exp();
            // (-e^(i phi) tanh r)^k: phase k (phi + pi), and the sign of tanh r
            let mut ph = k as f64 * (phi + std::f64::consts::PI);
            if th < 0.0 {
                ph += k as f64 * std::f64::consts::PI;
            }
            C64::new(mag * ph.cos(), mag * ph.sin())
        };
    }
    let e = C64::new(phi.cos() * sh, -phi.sin() * sh); // e^(-i phi) sinh r
    for col in 0..c {
        let inv = 1.0 / ((col + 1) as f64).sqrt();
        for row in 0..rows - col - 1 {
            let up = if row > 0 { m[(row - 1) * d + col].scale((row as f64).sqrt() * ch) } else { C64::ZERO };
            let down = e.mul(m[(row + 1) * d + col]).scale(((row + 1) as f64).sqrt());
            m[row * d + col + 1] = up.add(down).scale(inv);
        }
    }
    m.truncate(d * d);
    m
}

/// <m|S(z)|n> by the normal-ordered form S = exp(-t a^dag^2/2) sech(r)^(n+1/2) exp(conj(t) a^2/2),
/// t = e^(i phi) tanh r: <m|S|n> = sqrt(m! n! sech r) e^(i phi (m-n)/2) sum_k (-1)^i sech^k tanh^(i+j)
/// / (2^(i+j) k! i! j!), i = (m-k)/2, j = (n-k)/2; every term in logs, the sum by Neumaier.
fn squeezing_matrix_sum(r: f64, phi: f64, c: usize) -> Vec<C64> {
    let d = c + 1;
    let mut m = vec![C64::ZERO; d * d];
    if r == 0.0 {
        for i in 0..d {
            m[i * d + i] = C64::ONE;
        }
        return m;
    }
    let ln_cosh = r.abs() + (-2.0 * r.abs()).exp().ln_1p() - std::f64::consts::LN_2;
    let ln_sech = -ln_cosh;
    let th = r.tanh();
    let ln_half_th = (th.abs() / 2.0).ln();
    let lf: Vec<f64> = (0..d).map(ln_fact).collect();
    for row in 0..d {
        for col in 0..d {
            if (row + col) % 2 == 1 {
                continue;
            }
            let terms = row.min(col) / 2 + 1;
            let mut s = pq_nif::num::Sum::new(terms as u64);
            let mut k = row % 2;
            while k <= row.min(col) {
                let (i, j) = ((row - k) / 2, (col - k) / 2);
                let lm = 0.5 * (lf[row] + lf[col] + ln_sech) + k as f64 * ln_sech + (i + j) as f64 * ln_half_th
                    - lf[k]
                    - lf[i]
                    - lf[j];
                // (-1)^i, and the sign of tanh^(i+j) when r < 0
                let neg = (i % 2 == 1) ^ (th < 0.0 && (i + j) % 2 == 1);
                s.add(if neg { -lm.exp() } else { lm.exp() });
                k += 2;
            }
            let ph = phi * (row as f64 - col as f64) / 2.0;
            let v = s.value();
            m[row * d + col] = C64::new(v * ph.cos(), v * ph.sin());
        }
    }
    m
}



fn mat(m: &[C64]) -> String {
    let p: Vec<String> = m.iter().map(|z| format!("[{:e},{:e}]", z.re, z.im)).collect();
    format!("[{}]", p.join(","))
}

fn main() {
    let mut out = Vec::new();
    // S8: sectors N = 0..20 at three angle pairs; snap deviations up to N = 100 are checked by the unit tests
    for &(t, p) in &[(0.7, -0.4), (std::f64::consts::FRAC_PI_4, 0.0), (2.9, 1.3)] {
        for n in 0..=20 {
            let u = beamsplitter_sector(n, t, p).expect("snap");
            out.push(format!("{{\"kind\":\"bs\",\"n\":{n},\"theta\":{t:e},\"phi\":{p:e},\"m\":{}}}", mat(&u)));
        }
    }
    // S5: the (r, |alpha|) grid for cutoffs 10, 20, 40, 80
    let phi = 0.7;
    for &c in &[10usize, 20, 40, 80] {
        for &r in &[0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0] {
            let s = squeezing_matrix(r, phi, c);
            out.push(format!("{{\"kind\":\"sq\",\"c\":{c},\"r\":{r:e},\"phi\":{phi:e},\"m\":{}}}", mat(&s)));
            let sf = squeezing_final(r, phi, c).expect("did not converge");
            out.push(format!("{{\"kind\":\"sq_final\",\"c\":{c},\"r\":{r:e},\"phi\":{phi:e},\"m\":{}}}", mat(&sf)));
            let s3 = squeezing_matrix_bargmann(r, phi, c);
            out.push(format!("{{\"kind\":\"sq_barg\",\"c\":{c},\"r\":{r:e},\"phi\":{phi:e},\"m\":{}}}", mat(&s3)));
            // padding per parity so that tanh(r)^(2 * pad) < 1e-20, at most 2000 states
            let pad = ((46.0 / -(f64::tanh(r).ln())) / 2.0).ceil() as usize + 10;
            let size = (c / 2 + 1 + pad).min(2000);
            let s4 = squeezing_matrix_eigen(r, phi, c, size);
            out.push(format!("{{\"kind\":\"sq_eig\",\"c\":{c},\"r\":{r:e},\"phi\":{phi:e},\"size\":{size},\"m\":{}}}", mat(&s4)));
            let s2 = squeezing_matrix_sum(r, phi, c);
            out.push(format!("{{\"kind\":\"sq_sum\",\"c\":{c},\"r\":{r:e},\"phi\":{phi:e},\"m\":{}}}", mat(&s2)));
        }
        for &a in &[0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0] {
            let alpha = C64::new(a * phi.cos(), a * phi.sin());
            let d = displacement_matrix(alpha, c);
            out.push(format!("{{\"kind\":\"disp\",\"c\":{c},\"a\":{a:e},\"phi\":{phi:e},\"m\":{}}}", mat(&d)));
            let d2 = displacement_matrix_laguerre(alpha, c);
            out.push(format!("{{\"kind\":\"disp_lag\",\"c\":{c},\"a\":{a:e},\"phi\":{phi:e},\"m\":{}}}", mat(&d2)));
        }
    }
    println!("[{}]", out.join(",\n"));
}
