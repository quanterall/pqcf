//! GT-H01 (hafnians) and spike S3: the hafnian family against thewalrus 0.22.0 (test/oracle/gt_h01.py), within
//! 1e-10 of the scale S (NFR-06, amendment A-1). Matrices above 28 run only with PQ_LONG=1 (the spike).
use pq_nif::hafnian::{hafnian, hafnian_repeated};
use pq_nif::num::C64;
use serde_json::Value;
use std::time::Instant;

fn never() -> bool {
    false
}

fn cv(v: &Value) -> C64 {
    C64::new(v[0].as_f64().unwrap(), v[1].as_f64().unwrap())
}

#[test]
fn against_thewalrus() {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../test/fixtures/gt_h01.json");
    let f: Value = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    let long = std::env::var("PQ_LONG").is_ok();
    let threads = std::env::var("PQ_THREADS").ok().and_then(|t| t.parse().ok()).unwrap_or(1);
    let mut worst: f64 = 0.0;
    for c in f["cases"].as_array().unwrap() {
        let rows = c["matrix"].as_array().unwrap();
        let m = rows.len();
        if m > 28 && !long {
            continue;
        }
        let a: Vec<C64> = rows.iter().flat_map(|r| r.as_array().unwrap().iter().map(cv)).collect();
        let loops: Option<Vec<C64>> = c.get("loops").map(|l| l.as_array().unwrap().iter().map(cv).collect());
        let kind = c["kind"].as_str().unwrap();
        let t0 = Instant::now();
        let r = match c.get("mult") {
            Some(mu) => {
                let mult: Vec<u32> = mu.as_array().unwrap().iter().map(|x| x.as_u64().unwrap() as u32).collect();
                hafnian_repeated(&a, m, loops.as_deref(), &mult, threads, &never).unwrap()
            }
            None => hafnian(&a, m, loops.as_deref(), threads, &never).unwrap(),
        };
        let secs = t0.elapsed().as_secs_f64();
        let w = cv(&c["value"]);
        let err = ((r.value.re - w.re).abs() + (r.value.im - w.im).abs()) / r.scale.max(1.0);
        worst = worst.max(err);
        if long {
            println!("{kind} m = {m}: error / S = {err:.1e}, S = {:.2e}, |value| = {:.2e}, ours {secs:.3} s, thewalrus {:.3} s",
                     r.scale, w.abs(), c.get("thewalrus_seconds").and_then(|t| t.as_f64()).unwrap_or(f64::NAN));
        }
        assert!(err <= 1e-10, "{kind} m = {m}: {err:e}");
    }
    println!("worst error / S: {worst:e}");
}
