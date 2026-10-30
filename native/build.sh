#!/bin/sh
# Builds the Rust NIF and puts it where pq_nif.erl loads it from (priv/pq_nif.so).
set -e
cd "$(dirname "$0")/pq_nif"
cargo build --release --locked 2>/dev/null || cargo build --release
mkdir -p ../../priv
cp target/release/libpq_nif.so ../../priv/pq_nif.so
