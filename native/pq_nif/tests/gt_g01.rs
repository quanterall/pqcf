//! GT-G01: GBS probabilities (photon counting, displaced and mixed states; threshold, zero mean) against
//! thewalrus 0.22.0 (test/oracle/gt_g01.py) within 1e-10, and the chain-rule sampler's frequencies against them.
use pq_nif::gbs::{amat, pnr_probability, sample_chain_rule, threshold_probability, Caps};
use pq_nif::rng::Source;
use serde_json::Value;

fn never() -> bool {
    false
}

fn fixture() -> Value {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../test/fixtures/gt_g01.json");
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

fn state(c: &Value) -> (Vec<f64>, Vec<f64>, usize) {
    let mu: Vec<f64> = c["mean"].as_array().unwrap().iter().map(|x| x.as_f64().unwrap()).collect();
    let sigma: Vec<f64> =
        c["cov"].as_array().unwrap().iter().flat_map(|r| r.as_array().unwrap().iter().map(|x| x.as_f64().unwrap())).collect();
    let m = mu.len() / 2;
    (mu, sigma, m)
}

fn pattern(v: &Value) -> Vec<u32> {
    v["pattern"].as_array().unwrap().iter().map(|x| x.as_u64().unwrap() as u32).collect()
}

#[test]
fn probabilities_against_thewalrus() {
    let mut worst: f64 = 0.0;
    for c in fixture()["cases"].as_array().unwrap() {
        let (mu, sigma, m) = state(c);
        let am = amat(&mu, &sigma, m, 2.0).unwrap();
        for e in c["pnr"].as_array().unwrap() {
            let got = pnr_probability(&am, &pattern(e), &never).unwrap();
            worst = worst.max((got - e["p"].as_f64().unwrap()).abs());
        }
        for e in c["threshold"].as_array().unwrap() {
            let got = threshold_probability(&am, &pattern(e), &never).unwrap();
            worst = worst.max((got - e["p"].as_f64().unwrap()).abs());
        }
    }
    println!("worst |p - thewalrus|: {worst:e}");
    assert!(worst < 1e-10, "{worst:e}");
}

/// the samplers' frequencies match the exact probabilities (4 sigma per pattern with p >= 1e-3); states of up to
/// two modes, so that the test stays short (every shot recomputes its prefix probabilities)
#[test]
fn samplers_follow_the_probabilities() {
    let shots = 5000;
    let mut checked = 0;
    for c in fixture()["cases"].as_array().unwrap().iter().filter(|c| c["n_modes"].as_u64().unwrap() <= 2) {
        let (mu, sigma, m) = state(c);
        for (key, threshold) in [("pnr", false), ("threshold", true)] {
            let t0 = std::time::Instant::now();
            let exact = c[key].as_array().unwrap();
            if exact.is_empty() {
                continue;
            }
            // caps of 8: a pattern of n photons costs a hafnian of 2^n subsets while repeated indices are expanded
            let caps = Caps { max_photons: 8, max_per_mode: 8, max_clicks: 20 };
            let (s, lost) = sample_chain_rule(&mu, &sigma, m, 2.0, threshold, &caps, shots, &mut Source::Seeded { seed: 9, start: 0 }, &never)
                .unwrap_or_else(|e| panic!("{key}, {m} modes: {e:?}"));
            // with mass beyond the caps the sampler draws the truncated, renormalised distribution (AN 6.2):
            // compare with the untruncated probabilities only where that mass is negligible
            if lost > 1e-3 {
                println!("{key}, {m} modes: {lost:.3} of the mass beyond the caps, not compared");
                continue;
            }
            checked += 1;
            for e in exact {
                let p = e["p"].as_f64().unwrap();
                if p < 1e-3 || pattern(e).iter().sum::<u32>() > 8 {
                    continue;
                }
                let pat = pattern(e);
                let f = s.chunks(m).filter(|t| *t == pat.as_slice()).count() as f64 / shots as f64;
                let sd = (p * (1.0 - p) / shots as f64).sqrt();
                assert!((f - p).abs() < 4.0 * sd + 1e-4, "{key} {pat:?}: {f} vs {p}");
            }
            println!("{key}, {m} modes: {shots} shots in {:.2} s", t0.elapsed().as_secs_f64());
        }
    }
    assert!(checked >= 2, "only {checked} cases compared");
}

/// weakly squeezed, displaced and mixed by beamsplitters: almost no mass beyond 8 photons, so the PNR sampler is
/// compared with the exact probabilities (pnr_probability, checked against thewalrus above)
#[test]
fn pnr_sampler_on_weak_states() {
    use pq_nif::gaussian::{apply_symplectic, passive_symplectic};
    use pq_nif::num::C64;
    for m in [2usize, 3] {
        let n = 2 * m;
        let mut mu = vec![0.0; n];
        let mut sigma: Vec<f64> = (0..n * n).map(|k| if k / n == k % n { 1.0 } else { 0.0 }).collect();
        for j in 0..m {
            let r = 0.25 + 0.1 * j as f64;
            apply_symplectic(&mut mu, &mut sigma, m, &[(-r).exp(), 0.0, 0.0, r.exp()], &[0.3 * j as f64, 0.2], &[j], None).unwrap();
        }
        for j in 0..m - 1 {
            let (t, p): (f64, f64) = (0.6 + 0.2 * j as f64, 0.4);
            let u = [C64::new(t.cos(), 0.0), C64::new(-p.cos() * t.sin(), p.sin() * t.sin()),
                     C64::new(p.cos() * t.sin(), p.sin() * t.sin()), C64::new(t.cos(), 0.0)];
            apply_symplectic(&mut mu, &mut sigma, m, &passive_symplectic(&u, 2), &[0.0; 4], &[j, j + 1], None).unwrap();
        }
        let am = amat(&mu, &sigma, m, 2.0).unwrap();
        let caps = Caps { max_photons: 8, max_per_mode: 8, max_clicks: 20 };
        let shots = 20000;
        let (s, lost) = sample_chain_rule(&mu, &sigma, m, 2.0, false, &caps, shots, &mut Source::Seeded { seed: 3, start: 0 }, &never).unwrap();
        assert!(lost < 1e-3, "{m} modes: lost {lost}");
        let mut seen = std::collections::BTreeMap::new();
        for t in s.chunks(m) {
            *seen.entry(t.to_vec()).or_insert(0u32) += 1;
        }
        let mut compared = 0;
        for (pat, count) in &seen {
            let p = pnr_probability(&am, pat, &never).unwrap();
            if p < 1e-3 {
                continue;
            }
            compared += 1;
            let f = *count as f64 / shots as f64;
            let sd = (p * (1.0 - p) / shots as f64).sqrt();
            assert!((f - p).abs() < 4.0 * sd + 1e-4, "{m} modes {pat:?}: {f} vs {p}");
        }
        println!("{m} modes: {compared} patterns compared, lost mass {lost:.1e}");
        assert!(compared >= 5);
    }
}
