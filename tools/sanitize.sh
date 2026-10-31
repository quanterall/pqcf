#!/bin/sh
# FS 12.6 sanitizer runs: the Rust crate's unit and integration tests under AddressSanitizer. Sanitizers need a
# nightly compiler, so this runs in a separate nightly container and leaves the pinned toolchain alone.
#   sh tools/sanitize.sh            (from the repository root; the tests read test/fixtures)
set -e
docker run --rm -v "$PWD":/w -w /w/native/pq_nif -e CARGO_TARGET_DIR=/tmp/t rustlang/rust:nightly sh -c '
  rustup toolchain install nightly --profile minimal -c rust-src >/dev/null 2>&1
  T=$(rustc +nightly -vV | sed -n "s/host: //p")
  RUSTFLAGS="-Zsanitizer=address" cargo +nightly test -Zbuild-std --target "$T" --release -q --lib --tests'
