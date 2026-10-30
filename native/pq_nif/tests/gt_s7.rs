//! S7: the stream draws and the boson samplers' stream order against the independent reference
//! (test/oracle/gt_s7.py): uniforms, bounded integers, shuffles and samples bit for bit, normals to 1e-15.
use pq_nif::boson::{sample_boson, Algorithm};
use pq_nif::num::C64;
use pq_nif::rng::Stream;
use serde_json::Value;

fn fixture() -> Value {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../test/fixtures/gt_s7.json");
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

fn never() -> bool {
    false
}

#[test]
fn draws() {
    for d in fixture()["draws"].as_array().unwrap() {
        let mut s = Stream::new(d["seed"].as_u64().unwrap(), d["stream"].as_u64().unwrap());
        for v in d["uniform"].as_array().unwrap() {
            assert_eq!(s.uniform().to_bits(), v.as_f64().unwrap().to_bits());
        }
        for b in d["below"].as_array().unwrap() {
            let bound = b[0].as_u64().unwrap();
            for v in b[1].as_array().unwrap() {
                assert_eq!(s.below(bound), v.as_u64().unwrap(), "below({bound})");
            }
        }
        for pair in d["normal"].as_array().unwrap() {
            let (a, b) = s.normal_pair();
            for (x, w) in [(a, &pair[0]), (b, &pair[1])] {
                let w = w.as_f64().unwrap();
                assert!((x - w).abs() <= 1e-15 * w.abs().max(1.0), "{x} vs {w}");
            }
        }
        let mut v: Vec<u64> = (0..10).collect();
        s.shuffle(&mut v);
        let want: Vec<u64> = d["shuffle"].as_array().unwrap().iter().map(|x| x.as_u64().unwrap()).collect();
        assert_eq!(v, want);
    }
}

#[test]
fn sampler_stream_order() {
    let f = fixture();
    let u: Vec<C64> = f["unitary"]
        .as_array()
        .unwrap()
        .iter()
        .flat_map(|row| row.as_array().unwrap().iter().map(|z| C64::new(z[0].as_f64().unwrap(), z[1].as_f64().unwrap())))
        .collect();
    for c in f["samplers"].as_array().unwrap() {
        let input: Vec<u32> = c["input"].as_array().unwrap().iter().map(|x| x.as_u64().unwrap() as u32).collect();
        let (seed, start, shots) =
            (c["seed"].as_u64().unwrap(), c["stream_start"].as_u64().unwrap(), c["shots"].as_u64().unwrap() as usize);
        for (key, alg) in [("brute_force", Algorithm::BruteForce { max_enumeration: 1000 }),
                           ("clifford_clifford", Algorithm::CliffordClifford)] {
            let got = sample_boson(&u, 3, &input, shots, seed, start, alg, &never).unwrap();
            let want: Vec<u32> = c[key]
                .as_array()
                .unwrap()
                .iter()
                .flat_map(|t| t.as_array().unwrap().iter().map(|x| x.as_u64().unwrap() as u32))
                .collect();
            assert_eq!(got, want, "{key} {input:?}");
        }
    }
}
