//! GT-H01 (torontonian) and spike S4: against thewalrus 0.22.0 (test/oracle/gt_tor.py) within 1e-10 of the
//! scale S; P(B) = tor / sqrt(det Q) against thewalrus's threshold probability. c > 16 only with PQ_LONG=1.
use pq_nif::num::C64;
use pq_nif::toronto::torontonian;
use serde_json::Value;
use std::time::Instant;

fn never() -> bool {
    false
}

#[test]
fn against_thewalrus() {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../test/fixtures/gt_tor.json");
    let f: Value = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    let long = std::env::var("PQ_LONG").is_ok();
    for c in f["cases"].as_array().unwrap() {
        let rows = c["o"].as_array().unwrap();
        let d = rows.len();
        if d / 2 > 16 && !long {
            continue;
        }
        let o: Vec<C64> = rows
            .iter()
            .flat_map(|r| r.as_array().unwrap().iter().map(|z| C64::new(z[0].as_f64().unwrap(), z[1].as_f64().unwrap())))
            .collect();
        let t0 = Instant::now();
        let r = torontonian(&o, d, 1, &never).unwrap();
        let secs = t0.elapsed().as_secs_f64();
        let want = c["tor"][0].as_f64().unwrap();
        let err = (r.value - want).abs() / r.scale.max(1.0);
        let mut perr = f64::NAN;
        if let Some(p) = c["threshold_prob"].as_f64().filter(|p| p.is_finite()) {
            perr = (r.value / c["sqrt_det_q"].as_f64().unwrap() - p).abs();
            assert!(perr < 1e-10, "P(B), c = {}: {perr:e}", d / 2);
        }
        if long {
            println!("c = {}: error / S = {err:.1e}, S = {:.2e}, tor = {want:.3e}, P(B) error {perr:.1e}, ours {secs:.3} s, thewalrus {:.3} s",
                     d / 2, r.scale, c["thewalrus_seconds"].as_f64().unwrap());
        }
        assert!(err <= 1e-10, "c = {}: {err:e}", d / 2);
    }
}
