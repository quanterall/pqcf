//! pq_nif: the numerical kernels of the Photonic QC Framework (FS section 8).
use rustler::{Atom, Binary, Encoder, Env, NifResult, ResourceArc, Term};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::atomic::{AtomicBool, Ordering};

pub mod fock;
pub mod num;
pub mod perm;
mod rng;

use num::C64;

mod atoms {
    rustler::atoms! {
        ok, error, dimension_limit, count, max, cmat, bad_shape, non_finite, cancelled, nif_panic,
        cancel, threads, scale, positions, blocks, input, output, permanent_scale
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

rustler::init!("pq_nif");
