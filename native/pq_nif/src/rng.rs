//! FS 8.3: ChaCha8, the 64-bit seed expanded to a 256-bit key with SplitMix64, the 64-bit stream selecting an
//! independent sequence. Part of the NIF contract: any change is a major NIF version bump.
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
