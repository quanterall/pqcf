//! pq_nif: the numerical kernels of the Photonic QC Framework (FS section 8).
use rustler::{Atom, Env, NifResult, ResourceArc, Term};
use std::sync::atomic::{AtomicBool, Ordering};

mod rng;

mod atoms {
    rustler::atoms! { ok, error, dimension_limit, count, max }
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
    use rustler::Encoder;
    if count > RNG_MAX_COUNT {
        let detail = rustler::Term::map_from_pairs(
            env,
            &[(atoms::count().encode(env), count.encode(env)), (atoms::max().encode(env), RNG_MAX_COUNT.encode(env))],
        )?;
        return Ok((atoms::error(), (atoms::dimension_limit(), detail)).encode(env));
    }
    Ok((atoms::ok(), rng::u64s(seed, stream, count)).encode(env))
}

rustler::init!("pq_nif");
