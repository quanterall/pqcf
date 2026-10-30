//! pq_nif: the numerical kernels of the Photonic QC Framework (FS section 8).
use rustler::{Atom, Binary, Encoder, Env, NifResult, ResourceArc, Term};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::atomic::{AtomicBool, Ordering};

pub mod boson;
pub mod fock;
pub mod hafnian;
pub mod num;
pub mod perm;
pub mod rng;
pub mod toronto;

use num::C64;

mod atoms {
    rustler::atoms! {
        ok, error, dimension_limit, count, max, cmat, bad_shape, non_finite, cancelled, nif_panic,
        cancel, threads, scale, positions, blocks, input, output, permanent_scale,
        algorithm, seed, stream_start, max_enumeration, clifford_clifford, brute_force, imat, shots, patterns,
        probability_sum, photon_number_overflow, cvec, subsets, terms, non_physical_state, not_positive_definite,
        size
    }
}

/// NIF version {Major, Minor, Patch}; checked by pq_nif.erl at load time and in the node handshake.
const VERSION: (u32, u32, u32) = (0, 1, 0);
/// The largest rng_u64 count served on a normal scheduler (FS 8.1: only a small rng_u64 may run there).
const RNG_MAX_COUNT: usize = 65_536;

pub struct CancelToken(AtomicBool);

#[rustler::resource_impl]
impl rustler::Resource for CancelToken {}

#[rustler::nif]
fn nif_version() -> (u32, u32, u32, String) {
    (VERSION.0, VERSION.1, VERSION.2, env!("PQ_TARGET").to_string())
}

#[rustler::nif]
fn cancel_token_new() -> ResourceArc<CancelToken> {
    ResourceArc::new(CancelToken(AtomicBool::new(false)))
}

#[rustler::nif]
fn cancel_token_set(token: ResourceArc<CancelToken>) -> Atom {
    token.0.store(true, Ordering::Relaxed);
    atoms::ok()
}

#[rustler::nif]
fn rng_u64<'a>(env: Env<'a>, seed: u64, stream: u64, count: usize) -> NifResult<Term<'a>> {
    if count > RNG_MAX_COUNT {
        let detail = rustler::Term::map_from_pairs(
            env,
            &[(atoms::count().encode(env), count.encode(env)), (atoms::max().encode(env), RNG_MAX_COUNT.encode(env))],
        )?;
        return Ok((atoms::error(), (atoms::dimension_limit(), detail)).encode(env));
    }
    Ok((atoms::ok(), rng::u64s(seed, stream, count)).encode(env))
}

/// A packed complex matrix {cmat, Rows, Cols, Binary} (FS 8.2); a wrong size is bad_shape.
fn cmat(term: Term) -> Result<(Vec<C64>, usize, usize), String> {
    let (tag, rows, cols, bin): (Atom, usize, usize, Binary) = term.decode().map_err(|_| "not a cmat".to_string())?;
    if tag != atoms::cmat() || rows == 0 || cols == 0 {
        return Err("not a cmat".into());
    }
    if bin.len() != 16 * rows * cols {
        return Err(format!("cmat {rows}x{cols} needs {} bytes, got {}", 16 * rows * cols, bin.len()));
    }
    let f = |k: usize| f64::from_le_bytes(bin[8 * k..8 * k + 8].try_into().unwrap());
    Ok(((0..rows * cols).map(|e| C64::new(f(2 * e), f(2 * e + 1))).collect(), rows, cols))
}

struct KernelOpts {
    cancel: Option<ResourceArc<CancelToken>>,
    threads: usize,
}

fn kernel_opts(opts: Term) -> NifResult<KernelOpts> {
    let cancel = match opts.map_get(atoms::cancel()) {
        Ok(t) => Some(t.decode::<ResourceArc<CancelToken>>()?),
        Err(_) => None,
    };
    let threads = match opts.map_get(atoms::threads()) {
        Ok(t) => t.decode::<usize>()?.max(1),
        Err(_) => 1,
    };
    Ok(KernelOpts { cancel, threads })
}

fn error<'a, T: Encoder>(env: Env<'a>, detail: T) -> Term<'a> {
    (atoms::error(), detail).encode(env)
}

fn bad_shape<'a>(env: Env<'a>, msg: &str) -> Term<'a> {
    error(env, (atoms::bad_shape(), msg.to_string()))
}

fn perm_error<'a>(env: Env<'a>, e: perm::Error) -> Term<'a> {
    match e {
        perm::Error::BadShape(m) => bad_shape(env, m),
        perm::Error::TooLarge(_) => {
            let detail = Term::map_from_pairs(env, &[(atoms::positions().encode(env), perm::MAX_POSITIONS.encode(env))])
                .unwrap();
            error(env, (atoms::dimension_limit(), detail))
        }
        perm::Error::Cancelled => error(env, atoms::cancelled()),
        perm::Error::NonFinite(w) => error(env, (atoms::non_finite(), if w == "input" { atoms::input() } else { atoms::output() })),
    }
}

fn info<'a>(env: Env<'a>, r: &perm::PermResult, scale_key: Atom) -> Term<'a> {
    Term::map_from_pairs(
        env,
        &[
            (scale_key.encode(env), r.scale.encode(env)),
            (atoms::positions().encode(env), r.positions.encode(env)),
            (atoms::blocks().encode(env), r.blocks.encode(env)),
        ],
    )
    .unwrap()
}

/// Runs a kernel with panics caught at the boundary (FS 8.1: {error, {nif_panic, Message}}).
fn guarded<'a>(env: Env<'a>, f: impl FnOnce() -> Term<'a>) -> Term<'a> {
    match catch_unwind(AssertUnwindSafe(f)) {
        Ok(t) => t,
        Err(p) => {
            let msg = p
                .downcast_ref::<&str>()
                .map(|s| s.to_string())
                .or_else(|| p.downcast_ref::<String>().cloned())
                .unwrap_or_else(|| "panic".into());
            error(env, (atoms::nif_panic(), msg))
        }
    }
}

fn cancelled_fn(o: &KernelOpts) -> impl Fn() -> bool + Sync + '_ {
    move || o.cancel.as_ref().is_some_and(|t| t.0.load(Ordering::Relaxed))
}

#[rustler::nif(schedule = "DirtyCpu")]
fn permanent<'a>(env: Env<'a>, m: Term<'a>, opts: Term<'a>) -> NifResult<Term<'a>> {
    let o = kernel_opts(opts)?;
    let (a, rows, cols) = match cmat(m) {
        Ok(x) => x,
        Err(e) => return Ok(bad_shape(env, &e)),
    };
    if rows != cols {
        return Ok(bad_shape(env, "permanent needs a square matrix"));
    }
    let cancel = cancelled_fn(&o);
    Ok(guarded(env, || match perm::permanent(&a, rows, o.threads, &cancel) {
        Ok(r) => (atoms::ok(), ((r.value.re, r.value.im), info(env, &r, atoms::scale()))).encode(env),
        Err(e) => perm_error(env, e),
    }))
}

#[rustler::nif(schedule = "DirtyCpu")]
fn permanent_repeated<'a>(
    env: Env<'a>,
    m: Term<'a>,
    row_mult: Vec<u32>,
    col_mult: Vec<u32>,
    opts: Term<'a>,
) -> NifResult<Term<'a>> {
    let o = kernel_opts(opts)?;
    let (a, rows, cols) = match cmat(m) {
        Ok(x) => x,
        Err(e) => return Ok(bad_shape(env, &e)),
    };
    let cancel = cancelled_fn(&o);
    Ok(guarded(env, || match perm::permanent_repeated(&a, rows, cols, &row_mult, &col_mult, o.threads, &cancel) {
        Ok(r) => (atoms::ok(), ((r.value.re, r.value.im), info(env, &r, atoms::scale()))).encode(env),
        Err(e) => perm_error(env, e),
    }))
}

#[rustler::nif(schedule = "DirtyCpu")]
fn bs_probability<'a>(
    env: Env<'a>,
    u: Term<'a>,
    input: Vec<u32>,
    output: Vec<u32>,
    opts: Term<'a>,
) -> NifResult<Term<'a>> {
    let o = kernel_opts(opts)?;
    let (a, rows, cols) = match cmat(u) {
        Ok(x) => x,
        Err(e) => return Ok(bad_shape(env, &e)),
    };
    if rows != cols {
        return Ok(bad_shape(env, "bs_probability needs a square unitary"));
    }
    let cancel = cancelled_fn(&o);
    Ok(guarded(env, || match perm::bs_probability(&a, rows, &input, &output, o.threads, &cancel) {
        Ok((p, r)) => (atoms::ok(), (p, info(env, &r, atoms::permanent_scale()))).encode(env),
        Err(e) => perm_error(env, e),
    }))
}

fn opt<'a, T: rustler::Decoder<'a>>(opts: Term<'a>, key: Atom, default: T) -> NifResult<T> {
    match opts.map_get(key) {
        Ok(t) => t.decode(),
        Err(_) => Ok(default),
    }
}

/// sample_boson(U, InputOcc, NShots, Opts) (FS 8.4): shot i uses stream (seed, stream_start + i).
/// Opts: algorithm (clifford_clifford | brute_force), seed, stream_start, max_enumeration, cancel.
#[rustler::nif(schedule = "DirtyCpu")]
fn sample_boson_raw<'a>(env: Env<'a>, u: Term<'a>, input: Vec<u32>, shots: usize, opts: Term<'a>) -> NifResult<Term<'a>> {
    let o = kernel_opts(opts)?;
    let (a, rows, cols) = match cmat(u) {
        Ok(x) => x,
        Err(e) => return Ok(bad_shape(env, &e)),
    };
    if rows != cols || input.len() != rows {
        return Ok(bad_shape(env, "sample_boson needs a square unitary and one occupation per mode"));
    }
    let n: u64 = input.iter().map(|&x| x as u64).sum();
    if n > 65535 {
        return Ok(error(env, (atoms::photon_number_overflow(), n)));
    }
    let alg_atom: Atom = opt(opts, atoms::algorithm(), atoms::clifford_clifford())?;
    let seed: u64 = opt(opts, atoms::seed(), 0u64)?;
    let start: u64 = opt(opts, atoms::stream_start(), 0u64)?;
    let max_enum: u64 = opt(opts, atoms::max_enumeration(), 1_000_000u64)?;
    let algorithm = if alg_atom == atoms::brute_force() {
        boson::Algorithm::BruteForce { max_enumeration: max_enum }
    } else if alg_atom == atoms::clifford_clifford() {
        boson::Algorithm::CliffordClifford
    } else {
        return Err(rustler::Error::BadArg);
    };
    let cancel = cancelled_fn(&o);
    Ok(guarded(env, || match boson::sample_boson(&a, rows, &input, shots, seed, start, algorithm, &cancel) {
        Ok(samples) => {
            let mut bin = rustler::OwnedBinary::new(samples.len() * 2).expect("binary");
            for (k, x) in samples.iter().enumerate() {
                bin.as_mut_slice()[2 * k..2 * k + 2].copy_from_slice(&(*x as u16).to_le_bytes());
            }
            let packed = (atoms::imat(), shots, rows, bin.release(env));
            let info = Term::map_from_pairs(
                env,
                &[(atoms::algorithm().encode(env), alg_atom.encode(env)), (atoms::shots().encode(env), shots.encode(env))],
            )
            .unwrap();
            (atoms::ok(), (packed, info)).encode(env)
        }
        Err(boson::Error::BadShape(m)) => bad_shape(env, m),
        Err(boson::Error::TooManyPatterns(p)) => {
            let detail = Term::map_from_pairs(
                env,
                &[(atoms::patterns().encode(env), p.encode(env)), (atoms::max().encode(env), max_enum.encode(env))],
            )
            .unwrap();
            error(env, (atoms::dimension_limit(), detail))
        }
        Err(boson::Error::ProbabilitySum(t)) => error(env, (atoms::probability_sum(), t)),
        Err(boson::Error::Cancelled) => error(env, atoms::cancelled()),
    }))
}

/// A packed complex vector {cvec, Len, Binary}.
fn cvec(term: Term) -> Result<Vec<C64>, String> {
    let (tag, len, bin): (Atom, usize, Binary) = term.decode().map_err(|_| "not a cvec".to_string())?;
    if tag != atoms::cvec() || bin.len() != 16 * len {
        return Err("not a cvec of its length".into());
    }
    let f = |k: usize| f64::from_le_bytes(bin[8 * k..8 * k + 8].try_into().unwrap());
    Ok((0..len).map(|e| C64::new(f(2 * e), f(2 * e + 1))).collect())
}

fn haf_result<'a>(env: Env<'a>, r: Result<hafnian::HafResult, hafnian::Error>) -> Term<'a> {
    match r {
        Ok(h) => {
            let info = Term::map_from_pairs(
                env,
                &[(atoms::scale().encode(env), h.scale.encode(env)), (atoms::subsets().encode(env), h.subsets.encode(env))],
            )
            .unwrap();
            (atoms::ok(), ((h.value.re, h.value.im), info)).encode(env)
        }
        Err(hafnian::Error::BadShape(m)) => bad_shape(env, m),
        Err(hafnian::Error::TooLarge(m)) => {
            let d = Term::map_from_pairs(env, &[(atoms::size().encode(env), m.encode(env))]).unwrap();
            error(env, (atoms::dimension_limit(), d))
        }
        Err(hafnian::Error::Cancelled) => error(env, atoms::cancelled()),
        Err(hafnian::Error::NonFinite(w)) => {
            error(env, (atoms::non_finite(), if w == "input" { atoms::input() } else { atoms::output() }))
        }
    }
}

fn square(term: Term) -> Result<(Vec<C64>, usize), String> {
    let (a, r, c) = cmat(term)?;
    if r != c {
        return Err("a square matrix is needed".into());
    }
    Ok((a, r))
}

#[rustler::nif(schedule = "DirtyCpu")]
fn hafnian<'a>(env: Env<'a>, m: Term<'a>, opts: Term<'a>) -> NifResult<Term<'a>> {
    let o = kernel_opts(opts)?;
    let (a, n) = match square(m) {
        Ok(x) => x,
        Err(e) => return Ok(bad_shape(env, &e)),
    };
    let cancel = cancelled_fn(&o);
    Ok(guarded(env, || haf_result(env, hafnian::hafnian(&a, n, None, o.threads, &cancel))))
}

#[rustler::nif(schedule = "DirtyCpu")]
fn hafnian_repeated<'a>(env: Env<'a>, m: Term<'a>, mult: Vec<u32>, opts: Term<'a>) -> NifResult<Term<'a>> {
    let o = kernel_opts(opts)?;
    let (a, n) = match square(m) {
        Ok(x) => x,
        Err(e) => return Ok(bad_shape(env, &e)),
    };
    let cancel = cancelled_fn(&o);
    Ok(guarded(env, || haf_result(env, hafnian::hafnian_repeated(&a, n, None, &mult, o.threads, &cancel))))
}

/// loop_hafnian(M, Diag, Mult, Opts): the hafnian of M with its diagonal replaced by Diag (FS A.9)
#[rustler::nif(schedule = "DirtyCpu")]
fn loop_hafnian<'a>(env: Env<'a>, m: Term<'a>, diag: Term<'a>, mult: Vec<u32>, opts: Term<'a>) -> NifResult<Term<'a>> {
    let o = kernel_opts(opts)?;
    let (a, n) = match square(m) {
        Ok(x) => x,
        Err(e) => return Ok(bad_shape(env, &e)),
    };
    let d = match cvec(diag) {
        Ok(x) => x,
        Err(e) => return Ok(bad_shape(env, &e)),
    };
    let cancel = cancelled_fn(&o);
    Ok(guarded(env, || haf_result(env, hafnian::hafnian_repeated(&a, n, Some(&d), &mult, o.threads, &cancel))))
}

#[rustler::nif(schedule = "DirtyCpu")]
fn torontonian<'a>(env: Env<'a>, m: Term<'a>, opts: Term<'a>) -> NifResult<Term<'a>> {
    let o = kernel_opts(opts)?;
    let (a, n) = match square(m) {
        Ok(x) => x,
        Err(e) => return Ok(bad_shape(env, &e)),
    };
    let cancel = cancelled_fn(&o);
    Ok(guarded(env, || match toronto::torontonian(&a, n, o.threads, &cancel) {
        Ok(t) => {
            let info = Term::map_from_pairs(
                env,
                &[(atoms::scale().encode(env), t.scale.encode(env)), (atoms::terms().encode(env), t.terms.encode(env))],
            )
            .unwrap();
            (atoms::ok(), (t.value, info)).encode(env)
        }
        Err(toronto::Error::BadShape(m)) => bad_shape(env, m),
        Err(toronto::Error::TooLarge(c)) => {
            let d = Term::map_from_pairs(env, &[(atoms::size().encode(env), c.encode(env))]).unwrap();
            error(env, (atoms::dimension_limit(), d))
        }
        Err(toronto::Error::Cancelled) => error(env, atoms::cancelled()),
        Err(toronto::Error::NonFinite(w)) => {
            error(env, (atoms::non_finite(), if w == "input" { atoms::input() } else { atoms::output() }))
        }
        Err(toronto::Error::NotPositiveDefinite) => error(env, (atoms::non_physical_state(), atoms::not_positive_definite())),
    }))
}

rustler::init!("pq_nif");
