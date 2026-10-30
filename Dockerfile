# The build and test environment: OTP 27 and the pinned Rust toolchain, no system BLAS or LAPACK (FS 2.3).
#   docker build -t pqcf-dev .              then   docker run --rm pqcf-dev        (runs the tests)
FROM erlang:27
ENV RUSTUP_HOME=/usr/local/rustup CARGO_HOME=/usr/local/cargo PATH=/usr/local/cargo/bin:$PATH
COPY native/pq_nif/rust-toolchain.toml /tmp/rust-toolchain.toml
RUN curl -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal --default-toolchain none \
 && cd /tmp && rustup show active-toolchain || rustup toolchain install "$(sed -n 's/^channel = "\(.*\)"/\1/p' /tmp/rust-toolchain.toml)" --profile minimal -c clippy,rustfmt
WORKDIR /src
COPY . .
RUN rebar3 compile
CMD ["sh", "-c", "rebar3 eunit && rebar3 dialyzer && cd native/pq_nif && cargo test --release -q && cargo clippy --release --all-targets -- -D warnings"]
