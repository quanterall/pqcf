// The target triple is part of nif_version/0 (FS 8.1).
fn main() {
    println!("cargo:rustc-env=PQ_TARGET={}", std::env::var("TARGET").unwrap());
}
