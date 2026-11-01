//! GT-H03 (FS 15 item 2): threshold detection of displaced Gaussian states -- the loop torontonian -- against
//! thewalrus 0.22.0's threshold_detection_prob (test/oracle/gt_h03.py): within 1e-10 relative (or 1e-15 absolute),
//! random mixed states with random displacements, 1 to 12 clicks.
use pq_nif::gbs::{amat, threshold_probability};
use serde_json::Value;

fn never() -> bool {
    false
}

#[test]
fn displaced_threshold_against_thewalrus() {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../test/fixtures/gt_h03.json");
    let fx: Value = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    let cases = fx["cases"].as_array().unwrap();
    assert!(cases.len() >= 15);
    for c in cases {
        let m = c["n_modes"].as_u64().unwrap() as usize;
        let mu: Vec<f64> = c["mean"].as_array().unwrap().iter().map(|x| x.as_f64().unwrap()).collect();
        let sigma: Vec<f64> = c["cov"].as_array().unwrap().iter()
            .flat_map(|r| r.as_array().unwrap().iter().map(|x| x.as_f64().unwrap()).collect::<Vec<_>>()).collect();
        let pattern: Vec<u32> = c["pattern"].as_array().unwrap().iter().map(|x| x.as_u64().unwrap() as u32).collect();
        let want = c["prob"].as_f64().unwrap();
        let am = amat(&mu, &sigma, m, 2.0).unwrap();
        assert!(am.displaced);
        let got = threshold_probability(&am, &pattern, &never).unwrap();
        let err = (got - want).abs();
        assert!(err <= 1e-15 || err <= 1e-10 * want.abs(), "m = {m}, {pattern:?}: {got:e} vs {want:e}");
    }
}
