"""Independent reference for GT-R01 (FS 8.3): ChaCha8 written from the ChaCha definition (Bernstein, 2008),
64-bit block counter in words 12-13 and the 64-bit stream in words 14-15, the key from SplitMix64 of the seed.
It shares no code with the Rust crate; writes test/fixtures/gt_r01.json.

    python3 test/oracle/chacha8_ref.py
"""
import json
import pathlib

M32, M64 = 0xFFFFFFFF, 0xFFFFFFFFFFFFFFFF


def splitmix64(seed, n):
    out, state = [], seed
    for _ in range(n):
        state = (state + 0x9E3779B97F4A7C15) & M64
        z = state
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & M64
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & M64
        out.append(z ^ (z >> 31))
    return out


def rotl(x, n):
    return ((x << n) | (x >> (32 - n))) & M32


def quarter(s, a, b, c, d):
    s[a] = (s[a] + s[b]) & M32; s[d] = rotl(s[d] ^ s[a], 16)
    s[c] = (s[c] + s[d]) & M32; s[b] = rotl(s[b] ^ s[c], 12)
    s[a] = (s[a] + s[b]) & M32; s[d] = rotl(s[d] ^ s[a], 8)
    s[c] = (s[c] + s[d]) & M32; s[b] = rotl(s[b] ^ s[c], 7)


def block(key_words, counter, stream, rounds=8):
    init = [0x61707865, 0x3320646E, 0x79622D32, 0x6B206574, *key_words,
            counter & M32, counter >> 32, stream & M32, stream >> 32]
    s = list(init)
    for _ in range(rounds // 2):
        quarter(s, 0, 4, 8, 12); quarter(s, 1, 5, 9, 13); quarter(s, 2, 6, 10, 14); quarter(s, 3, 7, 11, 15)
        quarter(s, 0, 5, 10, 15); quarter(s, 1, 6, 11, 12); quarter(s, 2, 7, 8, 13); quarter(s, 3, 4, 9, 14)
    return [(x + y) & M32 for x, y in zip(s, init)]


def rng_u64(seed, stream, count):
    key_words = []
    for k in splitmix64(seed, 4):            # each u64 laid out little-endian: low word first
        key_words += [k & M32, k >> 32]
    words, counter = [], 0
    while len(words) < 2 * count:
        words += block(key_words, counter, stream)
        counter += 1
    return [words[2 * i] | (words[2 * i + 1] << 32) for i in range(count)]


CASES = [(0, 0), (0, 1), (1, 0), (42, 7), (M64, 0), (12345, M64)]
COUNT = 40          # 5 ChaCha blocks: crosses block and 4-block buffer boundaries

if __name__ == "__main__":
    # the ChaCha20 block function against RFC 7539 2.3.2 (same core, 20 rounds; that test uses a 96-bit nonce,
    # which in our 64/64 layout is counter word 13 = 0x09000000 and stream = 0x4a000000 | 0 << 32)
    key = [int.from_bytes(bytes(range(i, i + 4)), "little") for i in range(0, 32, 4)]
    out = block(key, 1 | (0x09000000 << 32), 0x4A000000, rounds=20)
    assert out[0] == 0xE4E7F110 and out[15] == 0x4E3C50A2, [hex(x) for x in out]
    fixtures = [{"seed": s, "stream": t, "values": rng_u64(s, t, COUNT)} for s, t in CASES]
    path = pathlib.Path(__file__).resolve().parents[1] / "fixtures" / "gt_r01.json"
    path.write_text(json.dumps({"test": "GT-R01", "generator": "test/oracle/chacha8_ref.py",
                                "algorithm": "ChaCha8, SplitMix64 key, 64-bit counter, 64-bit stream",
                                "cases": fixtures}, indent=1) + "\n")
    print(f"wrote {path} ({len(fixtures)} cases x {COUNT})")
