//! FS 8.3: ChaCha8, the 64-bit seed expanded to a 256-bit key with SplitMix64, the 64-bit stream selecting an
//! independent sequence. Part of the NIF contract: any change is a major NIF version bump.
#![forbid(unsafe_code)] // NFR-04: kernels are safe Rust; unsafe only in the Rustler glue (lib.rs)
use rand_chacha::ChaCha8Rng;
use rand_core::{RngCore, SeedableRng};

pub fn splitmix64_key(seed: u64) -> [u8; 32] {
    let mut state = seed;
    let mut key = [0u8; 32];
    for chunk in key.chunks_exact_mut(8) {
        state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        chunk.copy_from_slice(&(z ^ (z >> 31)).to_le_bytes());
    }
    key
}

pub fn stream(seed: u64, stream: u64) -> ChaCha8Rng {
    let mut rng = ChaCha8Rng::from_seed(splitmix64_key(seed));
    rng.set_stream(stream);
    rng
}

pub fn u64s(seed: u64, s: u64, count: usize) -> Vec<u64> {
    let mut rng = stream(seed, s);
    (0..count).map(|_| rng.next_u64()).collect()
}

/// One shot's stream with the derived draws of FS 8.5: a uniform double from one u64, bounded integers by Lemire's
/// rejection, normal pairs by Box-Muller (cosine branch first, a first uniform of 0 rejected).
pub struct Stream(ChaCha8Rng);

impl Stream {
    pub fn new(seed: u64, s: u64) -> Stream {
        Stream(stream(seed, s))
    }

    pub fn u64(&mut self) -> u64 {
        self.0.next_u64()
    }

    /// (x >> 11) * 2^-53, in [0, 1)
    pub fn uniform(&mut self) -> f64 {
        (self.u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// uniform in 0..s (s > 0), Lemire's multiply-and-reject on 64 bits
    pub fn below(&mut self, s: u64) -> u64 {
        assert!(s > 0);
        let mut m = (self.u64() as u128) * (s as u128);
        let mut low = m as u64;
        if low < s {
            let t = s.wrapping_neg() % s;
            while low < t {
                m = (self.u64() as u128) * (s as u128);
                low = m as u64;
            }
        }
        (m >> 64) as u64
    }

    /// two independent standard normals
    pub fn normal_pair(&mut self) -> (f64, f64) {
        let mut u1 = self.uniform();
        while u1 == 0.0 {
            u1 = self.uniform();
        }
        let u2 = self.uniform();
        let r = (-2.0 * u1.ln()).sqrt();
        let a = 2.0 * std::f64::consts::PI * u2;
        (r * a.cos(), r * a.sin())
    }

    /// Fisher-Yates shuffle: for i = len-1 down to 1, swap i with below(i + 1)
    pub fn shuffle<T>(&mut self, v: &mut [T]) {
        for i in (1..v.len()).rev() {
            let j = self.below(i as u64 + 1) as usize;
            v.swap(i, j);
        }
    }
}

/// Where a sampler's shots get their streams: shot i from (seed, start + i), or every draw from one shared stream
/// (a shot driven step by step from Erlang, FS 9.3: the shot's stream consumed in instruction order).
pub enum Source<'a> {
    Seeded { seed: u64, start: u64 },
    Shared(&'a mut Stream),
}

impl Source<'_> {
    pub fn with_shot<R>(&mut self, i: u64, f: impl FnOnce(&mut Stream) -> R) -> R {
        match self {
            Source::Seeded { seed, start } => f(&mut Stream::new(*seed, *start + i)),
            Source::Shared(s) => f(s),
        }
    }
}
