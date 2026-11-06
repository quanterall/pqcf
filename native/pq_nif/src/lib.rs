//! pq_nif: the numerical kernels of the Photonic QC Framework (FS section 8).
use rustler::{Atom, Binary, Encoder, Env, NifResult, ResourceArc, Term};
use std::collections::HashMap;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::atomic::{AtomicBool, Ordering};

pub mod boson;
pub mod fock;
pub mod fockdmb;
pub mod ffslp;
pub mod fockgrad;
pub mod fockplan;
pub mod fockstate;
pub mod gaussian;
pub mod gbs;
pub mod gpu;
pub mod hcount;
pub mod hafnian;
pub mod mesh;
pub mod num;
pub mod perm;
pub mod polymul;
pub mod rpoly;
pub mod rng;
pub mod toronto;

use num::C64;

mod atoms {
    rustler::atoms! {
        ok, error, dimension_limit, count, max, cmat, bad_shape, non_finite, cancelled, nif_panic,
        cancel, threads, scale, positions, blocks, input, output, permanent_scale,
        algorithm, seed, stream_start, max_enumeration, clifford_clifford, brute_force, imat, shots, patterns,
        probability_sum, photon_number_overflow, cvec, subsets, terms, non_physical_state, not_positive_definite,
        size, rvec, rmat, not_symplectic, homodyne, pnr, threshold, unsupported_input, pq_bug, beamsplitter,
        max_photons, max_per_mode, max_clicks, truncated_mass, negative_probability, nonzero_mean_threshold,
        phaseshift, kerr, passive, displacement, squeezing, loss, fix_outcome, stream, singular_matrix, stream_ref, metropolis, burn_in, thin, chains, acceptance_rate,
        probability, photons, shared, thermal_loss, extra, device, cpu, gpu, auto, reason,
        overflow, packing, q, add, mul, pow, rat, complex, expi, sym, exp, bad_prime, sub, singular, bad, const_ = "const"
    }
}

/// NIF version {Major, Minor, Patch}; checked by pq_nif.erl at load time and in the node handshake.
const VERSION: (u32, u32, u32) = (0, 2, 0);
/// The largest rng_u64 count served on a normal scheduler (FS 8.1: only a small rng_u64 may run there).
const RNG_MAX_COUNT: usize = 65_536;

pub struct CancelToken(AtomicBool);

/// One shot's stream (FS 9.3), held between calls when a shot is driven step by step from Erlang.
pub struct StreamRes(std::sync::Mutex<rng::Stream>);

#[rustler::resource_impl]
impl rustler::Resource for StreamRes {}

/// stream_new(Seed, Stream): the stream (Seed, Stream) as a resource; kernels given it as stream_ref draw from it in
/// order
#[rustler::nif]
fn stream_new(seed: u64, stream: u64) -> ResourceArc<StreamRes> {
    ResourceArc::new(StreamRes(std::sync::Mutex::new(rng::Stream::new(seed, stream))))
}

fn stream_ref(opts: Term) -> NifResult<Option<ResourceArc<StreamRes>>> {
    match opts.map_get(atoms::stream_ref()) {
        Ok(t) => Ok(Some(t.decode()?)),
        Err(_) => Ok(None),
    }
}

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
    device: Device,
}

/// Where a permanent runs (Phase 4, FS 15 item 5): device => cpu (default) | gpu | auto (the GPU from n = GPU_AUTO_N
/// when one is available)
#[derive(Clone, Copy, PartialEq)]
enum Device {
    Cpu,
    Gpu,
    Auto,
}

/// Spike G1 on the GB10: the GPU passes 20 CPU threads from n = 24
const GPU_AUTO_N: usize = 24;

fn kernel_opts(opts: Term) -> NifResult<KernelOpts> {
    let cancel = match opts.map_get(atoms::cancel()) {
        Ok(t) => Some(t.decode::<ResourceArc<CancelToken>>()?),
        Err(_) => None,
    };
    let threads = match opts.map_get(atoms::threads()) {
        Ok(t) => t.decode::<usize>()?.max(1),
        Err(_) => 1,
    };
    let device = match opts.map_get(atoms::device()) {
        Ok(t) => {
            let a: Atom = t.decode()?;
            if a == atoms::gpu() { Device::Gpu } else if a == atoms::auto() { Device::Auto } else if a == atoms::cpu() { Device::Cpu } else { return Err(rustler::Error::BadArg) }
        }
        Err(_) => Device::Cpu,
    };
    Ok(KernelOpts { cancel, threads, device })
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
    Ok(guarded(env, || match permanent_on(&a, rows, &o, &cancel) {
        Ok((r, dev)) => (atoms::ok(), ((r.value.re, r.value.im), info_on(env, &r, atoms::scale(), dev))).encode(env),
        Err(t) => t(env),
    }))
}

type ErrFn = Box<dyn for<'b> FnOnce(Env<'b>) -> Term<'b>>;

/// perm of a square matrix on the device the options ask for: the GPU kernel for device => gpu (an error when there
/// is none), or for auto from GPU_AUTO_N when there is one; else the CPU kernel
fn permanent_on(a: &[C64], n: usize, o: &KernelOpts, cancel: &(dyn Fn() -> bool + Sync)) -> Result<(perm::PermResult, Atom), ErrFn> {
    let use_gpu = match o.device {
        Device::Cpu => false,
        Device::Gpu => true,
        Device::Auto => n >= GPU_AUTO_N && gpu::device().is_ok(),
    };
    if use_gpu {
        return match gpu::permanent_cancellable(a, n, cancel) {
            Ok(r) => Ok((r, atoms::gpu())),
            Err(gpu::Error::Unavailable(why)) => Err(Box::new(move |env| {
                let d = Term::map_from_pairs(env, &[(atoms::device().encode(env), atoms::gpu().encode(env)),
                                                    (atoms::reason().encode(env), why.encode(env))]).unwrap();
                error(env, (atoms::unsupported_input(), d))
            })),
            Err(gpu::Error::Cancelled) => Err(Box::new(|env| error(env, atoms::cancelled()))),
            Err(gpu::Error::NonFinite) => Err(Box::new(|env| perm_error(env, perm::Error::NonFinite("output")))),
            Err(gpu::Error::BadShape(m)) => Err(Box::new(move |env| bad_shape(env, m))),
        };
    }
    match perm::permanent(a, n, o.threads, cancel) {
        Ok(r) => Ok((r, atoms::cpu())),
        Err(e) => Err(Box::new(move |env| perm_error(env, e))),
    }
}

/// the submatrix with rows and columns repeated by their multiplicities, when every multiplicity is 0 or 1 (the
/// GPU kernel's case); None otherwise
fn plain_submatrix(a: &[C64], rows: usize, cols: usize, r: &[u32], c: &[u32]) -> Option<(Vec<C64>, usize)> {
    if r.len() != rows || c.len() != cols || r.iter().chain(c).any(|&x| x > 1) {
        return None;
    }
    let ri: Vec<usize> = (0..rows).filter(|&i| r[i] == 1).collect();
    let ci: Vec<usize> = (0..cols).filter(|&j| c[j] == 1).collect();
    if ri.len() != ci.len() {
        return None;
    }
    Some((ri.iter().flat_map(|&i| ci.iter().map(move |&j| a[i * cols + j])).collect(), ri.len()))
}

fn info_on<'a>(env: Env<'a>, r: &perm::PermResult, scale_key: Atom, dev: Atom) -> Term<'a> {
    info(env, r, scale_key).map_put(atoms::device().encode(env), dev.encode(env)).unwrap()
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
    // the GPU kernel covers multiplicities 0 and 1; repeated rows or columns stay on the CPU
    if o.device != Device::Cpu {
        if let Some((sub, n)) = plain_submatrix(&a, rows, cols, &row_mult, &col_mult) {
            return Ok(guarded(env, || match permanent_on(&sub, n, &o, &cancel) {
                Ok((r, dev)) => (atoms::ok(), ((r.value.re, r.value.im), info_on(env, &r, atoms::scale(), dev))).encode(env),
                Err(t) => t(env),
            }));
        }
    }
    Ok(guarded(env, || match perm::permanent_repeated(&a, rows, cols, &row_mult, &col_mult, o.threads, &cancel) {
        Ok(r) => (atoms::ok(), ((r.value.re, r.value.im), info_on(env, &r, atoms::scale(), atoms::cpu()))).encode(env),
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
    // collision-free patterns on the GPU: |perm(U_(T,S))|^2, every factorial 1
    if o.device != Device::Cpu && input.iter().sum::<u32>() == output.iter().sum::<u32>() {
        if let Some((sub, n)) = plain_submatrix(&a, rows, cols, &output, &input) {
            return Ok(guarded(env, || match permanent_on(&sub, n, &o, &cancel) {
                Ok((r, dev)) => {
                    let p = r.value.re * r.value.re + r.value.im * r.value.im;
                    (atoms::ok(), (p, info_on(env, &r, atoms::permanent_scale(), dev))).encode(env)
                }
                Err(t) => t(env),
            }));
        }
    }
    Ok(guarded(env, || match perm::bs_probability(&a, rows, &input, &output, o.threads, &cancel) {
        Ok((p, r)) => (atoms::ok(), (p, info_on(env, &r, atoms::permanent_scale(), atoms::cpu()))).encode(env),
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
    } else if alg_atom == atoms::metropolis() {
        boson::Algorithm::Metropolis {
            burn_in: opt(opts, atoms::burn_in(), 100usize)?,
            thin: opt(opts, atoms::thin(), 10usize)?,
            chains: opt(opts, atoms::chains(), shots.clamp(1, 8))?,
        }
    } else {
        return Err(rustler::Error::BadArg);
    };
    // loss => {input | output, [Eta]} (FR-B05)
    let loss = match opts.map_get(atoms::loss()) {
        Ok(t) => {
            let (place, etas): (Atom, Vec<f64>) = t.decode()?;
            if place == atoms::input() { boson::Loss::Input(etas) } else if place == atoms::output() { boson::Loss::Output(etas) } else { return Err(rustler::Error::BadArg) }
        }
        Err(_) => boson::Loss::None,
    };
    // shared => [P] (Phase 4): per input mode the probability of the shared internal state
    let shared: Option<Vec<f64>> = match opts.map_get(atoms::shared()) {
        Ok(t) => Some(t.decode()?),
        Err(_) => None,
    };
    // extra => [Q] (Phase 4): per input mode the probability of one more, distinguishable photon
    let extra: Option<Vec<f64>> = match opts.map_get(atoms::extra()) {
        Ok(t) => Some(t.decode()?),
        Err(_) => None,
    };
    let cancel = cancelled_fn(&o);
    let result = match (algorithm, &loss, &shared, &extra) {
        (boson::Algorithm::Metropolis { burn_in, thin, chains }, boson::Loss::None, None, None) => {
            boson::metropolis(&a, rows, &input, shots, seed, burn_in, thin, chains, &cancel).map(|(s, r)| (s, Some(r)))
        }
        (alg, _, _, _) => boson::sample_boson_noisy(&a, rows, &input, shots, seed, start, alg, &loss, shared.as_deref(),
                                                    extra.as_deref(), &cancel)
            .map(|s| (s, None)),
    };
    Ok(guarded(env, || match result {
        Ok((samples, rate)) => {
            let mut bin = rustler::OwnedBinary::new(samples.len() * 2).expect("binary");
            for (k, x) in samples.iter().enumerate() {
                bin.as_mut_slice()[2 * k..2 * k + 2].copy_from_slice(&(*x as u16).to_le_bytes());
            }
            let packed = (atoms::imat(), shots, rows, bin.release(env));
            let mut pairs = vec![(atoms::algorithm().encode(env), alg_atom.encode(env)), (atoms::shots().encode(env), shots.encode(env))];
            if let Some(r) = rate {
                pairs.push((atoms::acceptance_rate().encode(env), r.encode(env)));
            }
            let info = Term::map_from_pairs(env, &pairs).unwrap();
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
    let ones = vec![1u32; n];
    Ok(guarded(env, || haf_on(env, &a, n, None, &ones, &o, &cancel)))
}

/// the (loop) hafnian with repeated indices on the device the options ask for: device => gpu runs the expanded
/// matrix on the GPU unless Kan's formula is cheaper (then the CPU); auto keeps the CPU (spike G1: the GB10's FP64
/// rate does not beat 20 CPU threads up to m = 40)
fn haf_on<'a>(env: Env<'a>, a: &[C64], n: usize, loops: Option<&[C64]>, mult: &[u32], o: &KernelOpts,
              cancel: &(dyn Fn() -> bool + Sync)) -> Term<'a> {
    let on_cpu = |r: Result<hafnian::HafResult, hafnian::Error>| put_haf(env, r, atoms::cpu());
    if o.device != Device::Gpu || mult.len() != n || a.len() != n * n || hafnian::kan_is_cheaper(n, mult) {
        return on_cpu(hafnian::hafnian_repeated(a, n, loops, mult, o.threads, cancel));
    }
    let (b, k, l) = match hafnian::expand(a, n, loops, mult) {
        Ok(x) => x,
        Err(e) => return on_cpu(Err(e)),
    };
    match gpu::hafnian(&b, k, l.as_deref(), cancel) {
        Ok(r) => put_haf(env, Ok(r), atoms::gpu()),
        Err(gpu::Error::Unavailable(why)) => {
            let d = Term::map_from_pairs(env, &[(atoms::device().encode(env), atoms::gpu().encode(env)),
                                                (atoms::reason().encode(env), why.encode(env))]).unwrap();
            error(env, (atoms::unsupported_input(), d))
        }
        Err(gpu::Error::Cancelled) => error(env, atoms::cancelled()),
        Err(gpu::Error::NonFinite) => on_cpu(Err(hafnian::Error::NonFinite("output"))),
        Err(gpu::Error::BadShape(m)) => bad_shape(env, m),
    }
}

/// haf_result with the device in the info map
fn put_haf<'a>(env: Env<'a>, r: Result<hafnian::HafResult, hafnian::Error>, dev: Atom) -> Term<'a> {
    let t = haf_result(env, r);
    match t.decode::<(Atom, (Term, Term))>() {
        Ok((ok, (v, info))) if ok == atoms::ok() => (ok, (v, info.map_put(atoms::device().encode(env), dev.encode(env)).unwrap())).encode(env),
        _ => t,
    }
}

#[rustler::nif(schedule = "DirtyCpu")]
fn hafnian_repeated<'a>(env: Env<'a>, m: Term<'a>, mult: Vec<u32>, opts: Term<'a>) -> NifResult<Term<'a>> {
    let o = kernel_opts(opts)?;
    let (a, n) = match square(m) {
        Ok(x) => x,
        Err(e) => return Ok(bad_shape(env, &e)),
    };
    let cancel = cancelled_fn(&o);
    Ok(guarded(env, || haf_on(env, &a, n, None, &mult, &o, &cancel)))
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
    Ok(guarded(env, || haf_on(env, &a, n, Some(&d), &mult, &o, &cancel)))
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

fn rvec(term: Term) -> Result<Vec<f64>, String> {
    let (tag, len, bin): (Atom, usize, Binary) = term.decode().map_err(|_| "not an rvec".to_string())?;
    if tag != atoms::rvec() || bin.len() != 8 * len {
        return Err("not an rvec of its length".into());
    }
    Ok((0..len).map(|k| f64::from_le_bytes(bin[8 * k..8 * k + 8].try_into().unwrap())).collect())
}

fn rmat(term: Term) -> Result<(Vec<f64>, usize, usize), String> {
    let (tag, r, c, bin): (Atom, usize, usize, Binary) = term.decode().map_err(|_| "not an rmat".to_string())?;
    if tag != atoms::rmat() || bin.len() != 8 * r * c {
        return Err("not an rmat of its size".into());
    }
    Ok(((0..r * c).map(|k| f64::from_le_bytes(bin[8 * k..8 * k + 8].try_into().unwrap())).collect(), r, c))
}

fn put_rvec<'a>(env: Env<'a>, v: &[f64]) -> Term<'a> {
    let mut bin = rustler::OwnedBinary::new(v.len() * 8).expect("binary");
    for (k, x) in v.iter().enumerate() {
        bin.as_mut_slice()[8 * k..8 * k + 8].copy_from_slice(&x.to_le_bytes());
    }
    (atoms::rvec(), v.len(), bin.release(env)).encode(env)
}

fn put_rmat<'a>(env: Env<'a>, v: &[f64], r: usize, c: usize) -> Term<'a> {
    let mut bin = rustler::OwnedBinary::new(v.len() * 8).expect("binary");
    for (k, x) in v.iter().enumerate() {
        bin.as_mut_slice()[8 * k..8 * k + 8].copy_from_slice(&x.to_le_bytes());
    }
    (atoms::rmat(), r, c, bin.release(env)).encode(env)
}

/// mean (2m) and covariance (2m x 2m) of a Gaussian state, checked against each other
fn gauss(mean: Term, cov: Term) -> Result<(Vec<f64>, Vec<f64>, usize), String> {
    let mu = rvec(mean)?;
    let (sig, r, c) = rmat(cov)?;
    if r != c || r != mu.len() || r % 2 == 1 || r == 0 {
        return Err("mean of length 2m and covariance 2m x 2m".into());
    }
    if mu.iter().chain(&sig).any(|x| !x.is_finite()) {
        return Err("non-finite".into());
    }
    Ok((mu, sig, r / 2))
}

fn modes_ok(modes: &[usize], m: usize) -> bool {
    let mut s = modes.to_vec();
    s.sort();
    s.dedup();
    s.len() == modes.len() && modes.iter().all(|&j| j < m)
}

fn gauss_result<'a>(env: Env<'a>, r: Result<(Vec<f64>, Vec<f64>), gaussian::Error>) -> Term<'a> {
    match r {
        Ok((mu, sig)) => {
            let n = mu.len();
            (atoms::ok(), (put_rvec(env, &mu), put_rmat(env, &sig, n, n))).encode(env)
        }
        Err(gaussian::Error::BadShape(m)) => bad_shape(env, m),
        Err(gaussian::Error::NotSymplectic(d)) => error(env, (atoms::not_symplectic(), d)),
        Err(gaussian::Error::NonPhysical(l)) => error(env, (atoms::non_physical_state(), l)),
        Err(gaussian::Error::Cancelled) => error(env, atoms::cancelled()),
    }
}

/// gaussian_apply_symplectic(Mean, Cov, S, D, Modes, Tol) (FS 8.4, A.2); Tol from pq_nif.erl (physicality_tol)
#[rustler::nif]
fn gaussian_apply_symplectic_raw<'a>(
    env: Env<'a>,
    mean: Term<'a>,
    cov: Term<'a>,
    s: Term<'a>,
    d: Term<'a>,
    modes: Vec<usize>,
    tol: f64,
) -> NifResult<Term<'a>> {
    let ((mut mu, mut sig, m), (sm, r, c), dv) = match (gauss(mean, cov), rmat(s), rvec(d)) {
        (Ok(g), Ok(sm), Ok(dv)) => (g, sm, dv),
        _ => return Ok(bad_shape(env, "mean, covariance, S or D")),
    };
    if r != c || r != 2 * modes.len() || dv.len() != r || !modes_ok(&modes, m) {
        return Ok(bad_shape(env, "S is 2k x 2k and D has length 2k for k distinct modes"));
    }
    let res = gaussian::apply_symplectic(&mut mu, &mut sig, m, &sm, &dv, &modes, Some(tol)).map(|_| (mu, sig));
    Ok(gauss_result(env, res))
}

/// gaussian_apply_passive(Mean, Cov, U, Modes) (FS 8.4, A.3)
#[rustler::nif]
fn gaussian_apply_passive<'a>(env: Env<'a>, mean: Term<'a>, cov: Term<'a>, u: Term<'a>, modes: Vec<usize>) -> NifResult<Term<'a>> {
    let ((mut mu, mut sig, m), (uv, r, c)) = match (gauss(mean, cov), cmat(u)) {
        (Ok(g), Ok(x)) => (g, x),
        _ => return Ok(bad_shape(env, "mean, covariance or U")),
    };
    if r != c || r != modes.len() || !modes_ok(&modes, m) {
        return Ok(bad_shape(env, "U is k x k for k distinct modes"));
    }
    let sm = gaussian::passive_symplectic(&uv, r);
    let zero = vec![0.0; 2 * r];
    let res = gaussian::apply_symplectic(&mut mu, &mut sig, m, &sm, &zero, &modes, None).map(|_| (mu, sig));
    Ok(gauss_result(env, res))
}

/// gaussian_apply_loss(Mean, Cov, Mode, Eta, Hbar) (FS 8.4, A.6)
#[rustler::nif]
fn gaussian_apply_loss<'a>(env: Env<'a>, mean: Term<'a>, cov: Term<'a>, mode: usize, eta: f64, hbar: f64) -> NifResult<Term<'a>> {
    let (mut mu, mut sig, m) = match gauss(mean, cov) {
        Ok(g) => g,
        Err(e) => return Ok(bad_shape(env, &e)),
    };
    if !(0.0..=1.0).contains(&eta) {
        return Ok(bad_shape(env, "eta outside [0, 1]"));
    }
    let res = gaussian::apply_loss(&mut mu, &mut sig, m, mode, eta, hbar).map(|_| (mu, sig));
    Ok(gauss_result(env, res))
}

/// gaussian_apply_thermal_loss(Mean, Cov, Mode, Eta, Nth, Hbar) (Phase 4): loss to a thermal environment
#[rustler::nif]
fn gaussian_apply_thermal_loss<'a>(env: Env<'a>, mean: Term<'a>, cov: Term<'a>, mode: usize, eta: f64, nth: f64, hbar: f64) -> NifResult<Term<'a>> {
    let (mut mu, mut sig, m) = match gauss(mean, cov) {
        Ok(g) => g,
        Err(e) => return Ok(bad_shape(env, &e)),
    };
    if !(0.0..=1.0).contains(&eta) || nth.is_nan() || nth < 0.0 {
        return Ok(bad_shape(env, "eta outside [0, 1] or nth negative"));
    }
    let res = gaussian::apply_thermal_loss(&mut mu, &mut sig, m, mode, eta, nth, hbar).map(|_| (mu, sig));
    Ok(gauss_result(env, res))
}

/// gaussian_marginal(Mean, Cov, Modes) (FS 8.4)
#[rustler::nif]
fn gaussian_marginal<'a>(env: Env<'a>, mean: Term<'a>, cov: Term<'a>, modes: Vec<usize>) -> NifResult<Term<'a>> {
    let (mu, sig, m) = match gauss(mean, cov) {
        Ok(g) => g,
        Err(e) => return Ok(bad_shape(env, &e)),
    };
    if modes.is_empty() || !modes_ok(&modes, m) {
        return Ok(bad_shape(env, "distinct modes in range"));
    }
    Ok(gauss_result(env, Ok(gaussian::marginal(&mu, &sig, m, &modes))))
}

/// the smallest eigenvalue of Cov + i (hbar/2) Omega (A.12), for pq_state:check_physical/1
#[rustler::nif(schedule = "DirtyCpu")]
fn gaussian_physicality<'a>(env: Env<'a>, cov: Term<'a>, hbar: f64) -> NifResult<Term<'a>> {
    match rmat(cov) {
        Ok((sig, r, c)) if r == c && r % 2 == 0 && r > 0 => {
            Ok((atoms::ok(), gaussian::physicality(&sig, r / 2, hbar)).encode(env))
        }
        _ => Ok(bad_shape(env, "covariance 2m x 2m")),
    }
}

/// sample_gaussian(Mean, Cov, NShots, Detector, Hbar, Opts) (FS 8.4): {homodyne, Phi} measures x_Phi on every mode
/// (rotated into x by a phase -Phi); photon counting comes with Gaussian boson sampling (Phase 2).
#[rustler::nif(schedule = "DirtyCpu")]
fn sample_gaussian<'a>(
    env: Env<'a>,
    mean: Term<'a>,
    cov: Term<'a>,
    shots: usize,
    detector: Term<'a>,
    _hbar: f64,
    opts: Term<'a>,
) -> NifResult<Term<'a>> {
    let o = kernel_opts(opts)?;
    let (mut mu, mut sig, m) = match gauss(mean, cov) {
        Ok(g) => g,
        Err(e) => return Ok(bad_shape(env, &e)),
    };
    let seed: u64 = opt(opts, atoms::seed(), 0u64)?;
    let start: u64 = opt(opts, atoms::stream_start(), 0u64)?;
    if let Ok(tag) = detector.decode::<Atom>() {
        if tag == atoms::pnr() || tag == atoms::threshold() {
            let caps = gbs::Caps {
                max_photons: opt(opts, atoms::max_photons(), 16u32)?,
                max_per_mode: opt(opts, atoms::max_per_mode(), opt(opts, atoms::max_photons(), 16u32)?)?,
                max_clicks: opt(opts, atoms::max_clicks(), 20u32)?,
            };
            let cancel = cancelled_fn(&o);
            let thr = tag == atoms::threshold();
            let sref = stream_ref(opts)?;
            let mut guard = sref.as_ref().map(|r| r.0.lock().expect("stream"));
            let mut src = match guard.as_deref_mut() {
                Some(st) => rng::Source::Shared(st),
                None => rng::Source::Seeded { seed, start },
            };
            return Ok(guarded(env, || match gbs::sample_chain_rule(&mu, &sig, m, _hbar, thr, &caps, shots, &mut src, &cancel) {
                Ok((samples, lost)) => {
                    let mut bin = rustler::OwnedBinary::new(samples.len() * 2).expect("binary");
                    for (k, x) in samples.iter().enumerate() {
                        bin.as_mut_slice()[2 * k..2 * k + 2].copy_from_slice(&(*x as u16).to_le_bytes());
                    }
                    let info = Term::map_from_pairs(env, &[(atoms::truncated_mass().encode(env), lost.encode(env))]).unwrap();
                    (atoms::ok(), ((atoms::imat(), shots, m, bin.release(env)), info)).encode(env)
                }
                Err(e) => gbs_error(env, e),
            }));
        }
    }
    let phi = match detector.decode::<(Atom, f64)>() {
        Ok((tag, phi)) if tag == atoms::homodyne() => phi,
        _ => return Ok(error(env, (atoms::unsupported_input(), detector))),
    };
    if phi != 0.0 {
        let rot = C64::new((-phi).cos(), (-phi).sin());
        for j in 0..m {
            let sm = gaussian::passive_symplectic(&[rot], 1);
            gaussian::apply_symplectic(&mut mu, &mut sig, m, &sm, &[0.0, 0.0], &[j], None).expect("rotation");
        }
    }
    let cancel = cancelled_fn(&o);
    let sref = stream_ref(opts)?;
    let mut guard = sref.as_ref().map(|r| r.0.lock().expect("stream"));
    let mut src = match guard.as_deref_mut() {
        Some(st) => rng::Source::Shared(st),
        None => rng::Source::Seeded { seed, start },
    };
    Ok(guarded(env, || match gaussian::sample_homodyne(&mu, &sig, m, shots, &mut src, &cancel) {
        Ok(y) => (atoms::ok(), (put_rmat(env, &y, shots, m), Term::map_new(env))).encode(env),
        Err(gaussian::Error::NonPhysical(l)) => error(env, (atoms::non_physical_state(), l)),
        Err(gaussian::Error::Cancelled) => error(env, atoms::cancelled()),
        Err(_) => bad_shape(env, "state"),
    }))
}

fn put_cvec<'a>(env: Env<'a>, v: &[C64]) -> Term<'a> {
    let mut bin = rustler::OwnedBinary::new(v.len() * 16).expect("binary");
    for (k, z) in v.iter().enumerate() {
        bin.as_mut_slice()[16 * k..16 * k + 8].copy_from_slice(&z.re.to_le_bytes());
        bin.as_mut_slice()[16 * k + 8..16 * k + 16].copy_from_slice(&z.im.to_le_bytes());
    }
    (atoms::cvec(), v.len(), bin.release(env)).encode(env)
}

/// the amplitudes of a pure state with its basis, checked against NModes and Cutoff (FS 3.4, FR-F04)
fn fock_state(amps: Term, n_modes: usize, cutoff: usize) -> Result<(Vec<C64>, fockstate::Basis), String> {
    let v = cvec(amps)?;
    if n_modes == 0 || fockstate::Basis::dim(n_modes, cutoff) != v.len() as u64 {
        return Err(format!("{} amplitudes for {n_modes} modes and cutoff {cutoff}", v.len()));
    }
    Ok((v, fockstate::Basis::new(n_modes, cutoff)))
}

fn fock_error<'a>(env: Env<'a>, e: fockstate::Error) -> Term<'a> {
    match e {
        fockstate::Error::BadShape(m) => bad_shape(env, m),
        fockstate::Error::Internal(m) => error(env, (atoms::pq_bug(), m)),
        fockstate::Error::TooLarge => error(env, (atoms::dimension_limit(), atoms::size())),
    }
}

macro_rules! fock_args {
    ($env:expr, $amps:expr, $n:expr, $c:expr, $modes:expr) => {
        match fock_state($amps, $n, $c) {
            Ok((v, b)) if $modes.iter().all(|&m: &usize| m < $n) => (v, b),
            Ok(_) => return Ok(bad_shape($env, "mode out of range")),
            Err(e) => return Ok(bad_shape($env, &e)),
        }
    };
}

#[rustler::nif(schedule = "DirtyCpu")]
#[allow(clippy::too_many_arguments)] // FS 8.4's signature
fn fock_apply_beamsplitter<'a>(env: Env<'a>, amps: Term<'a>, theta: f64, phi: f64, m1: usize, m2: usize, n: usize, c: usize) -> NifResult<Term<'a>> {
    let (mut v, b) = fock_args!(env, amps, n, c, [m1, m2]);
    if m1 == m2 {
        return Ok(bad_shape(env, "two distinct modes"));
    }
    Ok(match fockstate::beamsplitter(&b, &mut v, theta, phi, m1, m2) {
        Ok(()) => (atoms::ok(), put_cvec(env, &v)).encode(env),
        Err(e) => fock_error(env, e),
    })
}

#[rustler::nif(schedule = "DirtyCpu")]
fn fock_apply_phaseshift<'a>(env: Env<'a>, amps: Term<'a>, phi: f64, mode: usize, n: usize, c: usize) -> NifResult<Term<'a>> {
    let (mut v, b) = fock_args!(env, amps, n, c, [mode]);
    fockstate::phaseshift(&b, &mut v, phi, mode);
    Ok((atoms::ok(), put_cvec(env, &v)).encode(env))
}

#[rustler::nif(schedule = "DirtyCpu")]
fn fock_apply_kerr<'a>(env: Env<'a>, amps: Term<'a>, kappa: f64, mode: usize, n: usize, c: usize) -> NifResult<Term<'a>> {
    let (mut v, b) = fock_args!(env, amps, n, c, [mode]);
    fockstate::kerr(&b, &mut v, kappa, mode);
    Ok((atoms::ok(), put_cvec(env, &v)).encode(env))
}

#[rustler::nif(schedule = "DirtyCpu")]
fn fock_apply_passive<'a>(env: Env<'a>, amps: Term<'a>, u: Term<'a>, modes: Vec<usize>, n: usize, c: usize) -> NifResult<Term<'a>> {
    let (mut v, b) = fock_args!(env, amps, n, c, modes);
    let (uv, r, cc) = match cmat(u) {
        Ok(x) => x,
        Err(e) => return Ok(bad_shape(env, &e)),
    };
    if r != cc || r != modes.len() || !modes_ok(&modes, n) {
        return Ok(bad_shape(env, "U is k x k for k distinct modes"));
    }
    Ok(match fockstate::passive(&b, &mut v, &uv, &modes) {
        Ok(()) => (atoms::ok(), put_cvec(env, &v)).encode(env),
        Err(e) => fock_error(env, e),
    })
}

#[rustler::nif(schedule = "DirtyCpu")]
fn fock_apply_displacement<'a>(env: Env<'a>, amps: Term<'a>, r: f64, phi: f64, mode: usize, n: usize, c: usize) -> NifResult<Term<'a>> {
    let (mut v, b) = fock_args!(env, amps, n, c, [mode]);
    let deficit = fockstate::displacement(&b, &mut v, r, phi, mode);
    Ok((atoms::ok(), (put_cvec(env, &v), deficit)).encode(env))
}

#[rustler::nif(schedule = "DirtyCpu")]
fn fock_apply_squeezing<'a>(env: Env<'a>, amps: Term<'a>, r: f64, phi: f64, mode: usize, n: usize, c: usize) -> NifResult<Term<'a>> {
    let (mut v, b) = fock_args!(env, amps, n, c, [mode]);
    Ok(match fockstate::squeezing(&b, &mut v, r, phi, mode) {
        Ok(deficit) => (atoms::ok(), (put_cvec(env, &v), deficit)).encode(env),
        Err(e) => fock_error(env, e),
    })
}

/// fock_pnr_probabilities(State, Modes, NModes, Cutoff) -> [{Pattern, Prob}] (FR-F03), pure states
#[rustler::nif(schedule = "DirtyCpu")]
fn fock_pnr_probabilities<'a>(env: Env<'a>, amps: Term<'a>, modes: Vec<usize>, n: usize, c: usize) -> NifResult<Term<'a>> {
    if let Ok((rho, r, cc)) = cmat(amps) {
        if r != cc || n == 0 || fockstate::Basis::dim(n, c) != r as u64 || modes.is_empty() || !modes_ok(&modes, n) {
            return Ok(bad_shape(env, "density matrix D x D and distinct modes"));
        }
        let b = fockstate::Basis::new(n, c);
        return Ok((atoms::ok(), fockstate::dm::pnr_probabilities(&b, &rho, &modes)).encode(env));
    }
    let (v, b) = fock_args!(env, amps, n, c, modes);
    if modes.is_empty() || !modes_ok(&modes, n) {
        return Ok(bad_shape(env, "distinct modes in range"));
    }
    Ok((atoms::ok(), fockstate::pnr_probabilities(&b, &v, &modes)).encode(env))
}

/// unitary_decompose(U) -> {Mesh, Phases} (FR-C10, AN 9): Mesh = [{M, M+1, Theta, Phi}] in application order
#[rustler::nif(schedule = "DirtyCpu")]
fn unitary_decompose<'a>(env: Env<'a>, u: Term<'a>) -> NifResult<Term<'a>> {
    let (uv, r, c) = match cmat(u) {
        Ok(x) => x,
        Err(e) => return Ok(bad_shape(env, &e)),
    };
    if r != c {
        return Ok(bad_shape(env, "a square unitary"));
    }
    let (bs, phases) = mesh::decompose(&uv, r);
    let list: Vec<(usize, usize, f64, f64)> = bs.iter().map(|&(m, t, p)| (m, m + 1, t, p)).collect();
    Ok((atoms::ok(), (list, put_rvec(env, &phases))).encode(env))
}

/// sample_categorical(Probs, NShots, Opts): one uniform per shot from stream (seed, stream_start + i), the
/// categorical draw of AN 5.1 (Opts: seed, stream_start); the drawn indices
#[rustler::nif(schedule = "DirtyCpu")]
fn sample_categorical<'a>(env: Env<'a>, probs: Term<'a>, shots: usize, opts: Term<'a>) -> NifResult<Term<'a>> {
    let p = match rvec(probs) {
        Ok(p) if !p.is_empty() && p.iter().all(|x| x.is_finite() && *x >= 0.0) => p,
        _ => return Ok(bad_shape(env, "non-negative probabilities")),
    };
    let seed: u64 = opt(opts, atoms::seed(), 0u64)?;
    let start: u64 = opt(opts, atoms::stream_start(), 0u64)?;
    let sref = stream_ref(opts)?;
    let mut guard = sref.as_ref().map(|r| r.0.lock().expect("stream"));
    let mut src = match guard.as_deref_mut() {
        Some(st) => rng::Source::Shared(st),
        None => rng::Source::Seeded { seed, start },
    };
    let out: Vec<usize> = (0..shots).map(|i| src.with_shot(i as u64, |st| boson::categorical(&p, st.uniform()))).collect();
    Ok((atoms::ok(), out).encode(env))
}

fn gbs_error<'a>(env: Env<'a>, e: gbs::Error) -> Term<'a> {
    match e {
        gbs::Error::BadShape(m) => bad_shape(env, m),
        gbs::Error::NegativeProbability(p) => error(env, (atoms::negative_probability(), p)),
        gbs::Error::NonPhysical(w) => error(env, (atoms::non_physical_state(), w.to_string())),
        gbs::Error::Cancelled => error(env, atoms::cancelled()),
    }
}

/// gaussian_to_amat(Mean, Cov, Hbar) -> {A, Gamma, Norm} (FS 8.4, A.9)
#[rustler::nif(schedule = "DirtyCpu")]
fn gaussian_to_amat<'a>(env: Env<'a>, mean: Term<'a>, cov: Term<'a>, hbar: f64) -> NifResult<Term<'a>> {
    let (mu, sig, m) = match gauss(mean, cov) {
        Ok(g) => g,
        Err(e) => return Ok(bad_shape(env, &e)),
    };
    Ok(match gbs::amat(&mu, &sig, m, hbar) {
        Ok(am) => {
            let n = 2 * m;
            let a: Vec<C64> = (0..n * n).map(|k| C64::new(am.a[(k / n, k % n)].re, am.a[(k / n, k % n)].im)).collect();
            let g: Vec<C64> = am.gamma.iter().map(|z| C64::new(z.re, z.im)).collect();
            let mut bin = rustler::OwnedBinary::new(a.len() * 16).expect("binary");
            for (k, z) in a.iter().enumerate() {
                bin.as_mut_slice()[16 * k..16 * k + 8].copy_from_slice(&z.re.to_le_bytes());
                bin.as_mut_slice()[16 * k + 8..16 * k + 16].copy_from_slice(&z.im.to_le_bytes());
            }
            let amat = (atoms::cmat(), n, n, bin.release(env));
            (atoms::ok(), (amat, put_cvec(env, &g), am.norm)).encode(env)
        }
        Err(e) => gbs_error(env, e),
    })
}

/// gaussian_probability(Mean, Cov, Hbar, Pattern, Detector): P(pattern) for pnr (A.9) or threshold (A.10)
#[rustler::nif(schedule = "DirtyCpu")]
fn gaussian_probability<'a>(env: Env<'a>, mean: Term<'a>, cov: Term<'a>, hbar: f64, pattern: Vec<u32>, detector: Atom) -> NifResult<Term<'a>> {
    let (mu, sig, m) = match gauss(mean, cov) {
        Ok(g) => g,
        Err(e) => return Ok(bad_shape(env, &e)),
    };
    if pattern.len() != m {
        return Ok(bad_shape(env, "one count per mode"));
    }
    let never = || false;
    let r = gbs::amat(&mu, &sig, m, hbar).and_then(|am| {
        if detector == atoms::threshold() {
            gbs::threshold_probability(&am, &pattern, &never)
        } else {
            gbs::pnr_probability(&am, &pattern, &never)
        }
    });
    Ok(match r {
        Ok(p) => (atoms::ok(), p).encode(env),
        Err(e) => gbs_error(env, e),
    })
}

fn put_cmat<'a>(env: Env<'a>, v: &[C64], r: usize, c: usize) -> Term<'a> {
    let mut bin = rustler::OwnedBinary::new(v.len() * 16).expect("binary");
    for (k, z) in v.iter().enumerate() {
        bin.as_mut_slice()[16 * k..16 * k + 8].copy_from_slice(&z.re.to_le_bytes());
        bin.as_mut_slice()[16 * k + 8..16 * k + 16].copy_from_slice(&z.im.to_le_bytes());
    }
    (atoms::cmat(), r, c, bin.release(env)).encode(env)
}

/// fock_dm_apply_*(Rho, ...) (FS 8.4, AN 7.5): one entry point, the operation as a tuple; returns {Rho, Deficit}
/// (the trace deficit of displacement and squeezing, 0 otherwise)
#[rustler::nif(schedule = "DirtyCpu")]
fn fock_dm_apply_raw<'a>(env: Env<'a>, rho: Term<'a>, op: Term<'a>, n: usize, c: usize) -> NifResult<Term<'a>> {
    let (mut r, d, dc) = match cmat(rho) {
        Ok(x) => x,
        Err(e) => return Ok(bad_shape(env, &e)),
    };
    if d != dc || n == 0 || fockstate::Basis::dim(n, c) != d as u64 {
        return Ok(bad_shape(env, "density matrix D x D for NModes and Cutoff"));
    }
    let b = fockstate::Basis::new(n, c);
    let tag: Atom = rustler::types::tuple::get_tuple(op)?.first().ok_or(rustler::Error::BadArg)?.decode()?;
    let in_range = |ms: &[usize]| ms.iter().all(|&m| m < n);
    let res: Result<f64, fockstate::Error> = if tag == atoms::beamsplitter() {
        let (_, t, p, m1, m2): (Atom, f64, f64, usize, usize) = op.decode()?;
        if !in_range(&[m1, m2]) || m1 == m2 { return Ok(bad_shape(env, "modes")); }
        match fockplan::Plan::beamsplitter(&b, t, p, m1, m2) { Ok(pl) => fockstate::dm::unitary_like(&b, &mut r, |v| { pl.apply(v); Ok(()) }), Err(e) => Err(e) }
    } else if tag == atoms::phaseshift() {
        let (_, p, m): (Atom, f64, usize) = op.decode()?;
        if !in_range(&[m]) { return Ok(bad_shape(env, "mode")); }
        { let pl = fockplan::Plan::phaseshift(&b, p, m); fockstate::dm::unitary_like(&b, &mut r, |v| { pl.apply(v); Ok(()) }) }
    } else if tag == atoms::kerr() {
        let (_, k, m): (Atom, f64, usize) = op.decode()?;
        if !in_range(&[m]) { return Ok(bad_shape(env, "mode")); }
        { let pl = fockplan::Plan::kerr(&b, k, m); fockstate::dm::unitary_like(&b, &mut r, |v| { pl.apply(v); Ok(()) }) }
    } else if tag == atoms::passive() {
        let (_, u, ms): (Atom, Term, Vec<usize>) = op.decode()?;
        let (uv, ur, uc) = match cmat(u) { Ok(x) => x, Err(e) => return Ok(bad_shape(env, &e)) };
        if ur != uc || ur != ms.len() || !modes_ok(&ms, n) { return Ok(bad_shape(env, "U is k x k for k distinct modes")); }
        match fockplan::Plan::passive(&b, &uv, &ms) { Ok(pl) => fockstate::dm::unitary_like(&b, &mut r, |v| { pl.apply(v); Ok(()) }), Err(e) => Err(e) }
    } else if tag == atoms::displacement() {
        let (_, rr, p, m): (Atom, f64, f64, usize) = op.decode()?;
        if !in_range(&[m]) { return Ok(bad_shape(env, "mode")); }
        { let pl = fockplan::Plan::displacement(&b, rr, p, m); fockstate::dm::unitary_like(&b, &mut r, |v| { pl.apply(v); Ok(()) }) }
    } else if tag == atoms::squeezing() {
        let (_, rr, p, m): (Atom, f64, f64, usize) = op.decode()?;
        if !in_range(&[m]) { return Ok(bad_shape(env, "mode")); }
        match fockplan::Plan::squeezing(&b, rr, p, m) { Ok(pl) => fockstate::dm::unitary_like(&b, &mut r, |v| { pl.apply(v); Ok(()) }), Err(e) => Err(e) }
    } else if tag == atoms::loss() {
        let (_, eta, m): (Atom, f64, usize) = op.decode()?;
        if !in_range(&[m]) || !(0.0..=1.0).contains(&eta) { return Ok(bad_shape(env, "mode or eta")); }
        r = fockstate::dm::loss(&b, &r, eta, m);
        Ok(0.0)
    } else if tag == atoms::thermal_loss() {
        let (_, eta, nth, m): (Atom, f64, f64, usize) = op.decode()?;
        if !in_range(&[m]) || !(0.0..=1.0).contains(&eta) || nth.is_nan() || nth < 0.0 { return Ok(bad_shape(env, "mode, eta or nth")); }
        let (out, def) = fockstate::dm::thermal_loss(&b, &r, eta, nth, m);
        r = out;
        Ok(def)
    } else {
        return Err(rustler::Error::BadArg);
    };
    Ok(match res {
        Ok(def) => (atoms::ok(), (put_cmat(env, &r, d, d), def)).encode(env),
        Err(e) => fock_error(env, e),
    })
}

/// gaussian_measure_homodyne(Mean, Cov, Mode, Phi, Hbar, Opts) -> {Outcome, Mean, Cov} (FS 8.4, A.7, FR-G08).
/// Opts: fix_outcome, or seed and stream.
#[rustler::nif]
fn gaussian_measure_homodyne<'a>(env: Env<'a>, mean: Term<'a>, cov: Term<'a>, mode: usize, phi: f64, _hbar: f64, opts: Term<'a>) -> NifResult<Term<'a>> {
    let (mu, sig, m) = match gauss(mean, cov) {
        Ok(g) => g,
        Err(e) => return Ok(bad_shape(env, &e)),
    };
    if m < 2 {
        return Ok(bad_shape(env, "a mode must remain after the measurement"));
    }
    let fixed: Option<f64> = match opts.map_get(atoms::fix_outcome()) {
        Ok(t) => Some(t.decode()?),
        Err(_) => None,
    };
    let seed: u64 = opt(opts, atoms::seed(), 0u64)?;
    let stream: u64 = opt(opts, atoms::stream(), 0u64)?;
    let sref = stream_ref(opts)?;
    let mut guard = sref.as_ref().map(|r| r.0.lock().expect("stream"));
    let mut src = match guard.as_deref_mut() {
        Some(st) => rng::Source::Shared(st),
        None => rng::Source::Seeded { seed, start: stream },
    };
    Ok(match gaussian::measure_homodyne(&mu, &sig, m, mode, phi, fixed, &mut src) {
        Ok((y, mu2, sig2)) => {
            let n = mu2.len();
            (atoms::ok(), (y, put_rvec(env, &mu2), put_rmat(env, &sig2, n, n))).encode(env)
        }
        Err(gaussian::Error::NonPhysical(v)) => error(env, (atoms::singular_matrix(), v)),
        Err(_) => bad_shape(env, "mean, covariance or mode"),
    })
}

/// fock_measure_pnr(State, Modes, Outcome, NModes, Cutoff) -> {Prob, ConditionalState} (FS 8.4, FR-M07): the
/// state projected onto Outcome on Modes and renormalised; the measured modes stay in the state, holding their
/// outcome (no renumbering). A pure state (cvec) or a density matrix (cmat).
#[rustler::nif(schedule = "DirtyCpu")]
fn fock_measure_pnr<'a>(env: Env<'a>, state: Term<'a>, modes: Vec<usize>, outcome: Vec<u32>, n: usize, c: usize) -> NifResult<Term<'a>> {
    if modes.is_empty() || modes.len() != outcome.len() || !modes_ok(&modes, n) {
        return Ok(bad_shape(env, "distinct modes in range, one outcome each"));
    }
    if let Ok((mut rho, d, dc)) = cmat(state) {
        if d != dc || fockstate::Basis::dim(n, c) != d as u64 {
            return Ok(bad_shape(env, "density matrix D x D"));
        }
        let b = fockstate::Basis::new(n, c);
        let keep: Vec<bool> = b.patterns.iter().map(|p| modes.iter().zip(&outcome).all(|(&m, &o)| p[m] == o)).collect();
        let prob = (0..d).filter(|&i| keep[i]).map(|i| rho[i * d + i].re).sum::<f64>();
        for i in 0..d {
            for j in 0..d {
                rho[i * d + j] = if keep[i] && keep[j] && prob > 0.0 { rho[i * d + j].scale(1.0 / prob) } else { C64::ZERO };
            }
        }
        return Ok((atoms::ok(), (prob, put_cmat(env, &rho, d, d))).encode(env));
    }
    let (mut v, b) = fock_args!(env, state, n, c, modes);
    let keep: Vec<bool> = b.patterns.iter().map(|p| modes.iter().zip(&outcome).all(|(&m, &o)| p[m] == o)).collect();
    let prob: f64 = v.iter().zip(&keep).filter(|(_, &k)| k).map(|(a, _)| a.re * a.re + a.im * a.im).sum();
    let norm = prob.sqrt();
    for (a, &k) in v.iter_mut().zip(&keep) {
        *a = if k && norm > 0.0 { a.scale(1.0 / norm) } else { C64::ZERO };
    }
    Ok((atoms::ok(), (prob, put_cvec(env, &v))).encode(env))
}

/// fock_apply_loss_trajectory(Amps, Eta, Mode, NModes, Cutoff, Opts) -> {Amps, Lost} (FR-F05 trajectory mode):
/// one uniform from Opts stream_ref (or stream (seed, stream)) draws the lost-photon count
#[rustler::nif(schedule = "DirtyCpu")]
fn fock_apply_loss_trajectory<'a>(env: Env<'a>, amps: Term<'a>, eta: f64, mode: usize, n: usize, c: usize, opts: Term<'a>) -> NifResult<Term<'a>> {
    let (mut v, b) = fock_args!(env, amps, n, c, [mode]);
    if !(0.0..=1.0).contains(&eta) {
        return Ok(bad_shape(env, "eta outside [0, 1]"));
    }
    let seed: u64 = opt(opts, atoms::seed(), 0u64)?;
    let stream: u64 = opt(opts, atoms::stream(), 0u64)?;
    let sref = stream_ref(opts)?;
    let mut guard = sref.as_ref().map(|r| r.0.lock().expect("stream"));
    let mut src = match guard.as_deref_mut() {
        Some(st) => rng::Source::Shared(st),
        None => rng::Source::Seeded { seed, start: stream },
    };
    let u = src.with_shot(0, |st| st.uniform());
    let lost = fockstate::loss_trajectory(&b, &mut v, eta, mode, u);
    Ok((atoms::ok(), (put_cvec(env, &v), lost)).encode(env))
}

/// boson_partial_probability(U, Input, Output, Gram) -> P (Phase 4, FS 15 item 4): the probability of Output for
/// photons whose internal states have the Gram matrix Gram (n x n, photons in input-mode order); n! permanents
#[rustler::nif(schedule = "DirtyCpu")]
fn boson_partial_probability<'a>(env: Env<'a>, u: Term<'a>, input: Vec<u32>, output: Vec<u32>, gram: Term<'a>) -> NifResult<Term<'a>> {
    let (a, rows, cols) = match cmat(u) { Ok(x) => x, Err(e) => return Ok(bad_shape(env, &e)) };
    let (g, gr, gc) = match cmat(gram) { Ok(x) => x, Err(e) => return Ok(bad_shape(env, &e)) };
    let n: usize = input.iter().map(|&x| x as usize).sum();
    if rows != cols || input.len() != rows || output.len() != rows || gr != n || gc != n || n > 12 {
        return Ok(bad_shape(env, "a square unitary, one occupation per mode, an n x n Gram matrix, n <= 12"));
    }
    Ok(match boson::partial_probability(&a, rows, &input, &output, &g) {
        Ok(p) => (atoms::ok(), p).encode(env),
        Err(_) => bad_shape(env, "input, output and Gram matrix"),
    })
}

/// Photon-number blocks (Phase 4, FS 15 item 6; internal, fockdmb.rs): a density matrix as [{N1, N2, Binary}], each
/// block d_N1 x d_N2 complex row-major little-endian; zero blocks absent
fn get_blocks(term: Term, b: &fockstate::Basis) -> Result<fockdmb::Blocks, String> {
    let list: Vec<(usize, usize, Binary)> = term.decode().map_err(|_| "not a block list".to_string())?;
    let mut out = fockdmb::Blocks::empty(b);
    let k = b.cutoff + 1;
    for (a, c, bin) in list {
        if a >= k || c >= k || bin.len() != 16 * out.sizes[a] * out.sizes[c] {
            return Err(format!("block ({a}, {c}) of {} bytes", bin.len()));
        }
        let f = |x: usize| f64::from_le_bytes(bin[8 * x..8 * x + 8].try_into().unwrap());
        out.blocks[a * k + c] = Some((0..out.sizes[a] * out.sizes[c]).map(|e| C64::new(f(2 * e), f(2 * e + 1))).collect());
    }
    Ok(out)
}

fn put_blocks<'a>(env: Env<'a>, blk: &fockdmb::Blocks) -> Term<'a> {
    let k = blk.sizes.len();
    let mut items: Vec<Term<'a>> = Vec::new();
    for a in 0..k {
        for c in 0..k {
            if let Some(v) = &blk.blocks[a * k + c] {
                let mut bin = rustler::OwnedBinary::new(v.len() * 16).expect("binary");
                for (x, z) in v.iter().enumerate() {
                    bin.as_mut_slice()[16 * x..16 * x + 8].copy_from_slice(&z.re.to_le_bytes());
                    bin.as_mut_slice()[16 * x + 8..16 * x + 16].copy_from_slice(&z.im.to_le_bytes());
                }
                items.push((a, c, bin.release(env)).encode(env));
            }
        }
    }
    items.encode(env)
}

fn dmb_basis(n: usize, c: usize) -> Option<fockstate::Basis> {
    if n == 0 || fockstate::Basis::dim(n, c) > 50_000_000 { None } else { Some(fockstate::Basis::new(n, c)) }
}

/// fock_dmb_from(State, N, C): blocks from amplitudes (cvec) or a full matrix (cmat)
#[rustler::nif(schedule = "DirtyCpu")]
fn fock_dmb_from<'a>(env: Env<'a>, state: Term<'a>, n: usize, c: usize) -> NifResult<Term<'a>> {
    let Some(b) = dmb_basis(n, c) else { return Ok(bad_shape(env, "modes and cutoff")) };
    let d = b.patterns.len();
    if let Ok(v) = cvec(state) {
        if v.len() != d { return Ok(bad_shape(env, "amplitudes for NModes and Cutoff")); }
        return Ok((atoms::ok(), put_blocks(env, &fockdmb::Blocks::from_pure(&b, &v))).encode(env));
    }
    match cmat(state) {
        Ok((r, rr, cc)) if rr == d && cc == d => Ok((atoms::ok(), put_blocks(env, &fockdmb::Blocks::from_dense(&b, &r))).encode(env)),
        _ => Ok(bad_shape(env, "a cvec or a D x D cmat")),
    }
}

/// fock_dmb_to_dense(Blocks, N, C) -> cmat (D x D)
#[rustler::nif(schedule = "DirtyCpu")]
fn fock_dmb_to_dense<'a>(env: Env<'a>, blocks: Term<'a>, n: usize, c: usize) -> NifResult<Term<'a>> {
    let Some(b) = dmb_basis(n, c) else { return Ok(bad_shape(env, "modes and cutoff")) };
    let blk = match get_blocks(blocks, &b) { Ok(x) => x, Err(e) => return Ok(bad_shape(env, &e)) };
    let d = b.patterns.len();
    Ok((atoms::ok(), put_cmat(env, &blk.to_dense(d), d, d)).encode(env))
}

/// fock_dmb_apply(Blocks, Op, N, C) -> {Blocks, Deficit}: the operations of fock_dm_apply_raw on blocks
#[rustler::nif(schedule = "DirtyCpu")]
fn fock_dmb_apply<'a>(env: Env<'a>, blocks: Term<'a>, op: Term<'a>, n: usize, c: usize) -> NifResult<Term<'a>> {
    let Some(b) = dmb_basis(n, c) else { return Ok(bad_shape(env, "modes and cutoff")) };
    let mut r = match get_blocks(blocks, &b) { Ok(x) => x, Err(e) => return Ok(bad_shape(env, &e)) };
    let tag: Atom = rustler::types::tuple::get_tuple(op)?.first().ok_or(rustler::Error::BadArg)?.decode()?;
    let in_range = |ms: &[usize]| ms.iter().all(|&m| m < n);
    let res: Result<f64, fockstate::Error> = if tag == atoms::beamsplitter() {
        let (_, t, p, m1, m2): (Atom, f64, f64, usize, usize) = op.decode()?;
        if !in_range(&[m1, m2]) || m1 == m2 { return Ok(bad_shape(env, "modes")); }
        fockplan::Plan::beamsplitter(&b, t, p, m1, m2).and_then(|pl| r.unitary_plan(&b, &pl))
    } else if tag == atoms::phaseshift() {
        let (_, p, m): (Atom, f64, usize) = op.decode()?;
        if !in_range(&[m]) { return Ok(bad_shape(env, "mode")); }
        r.unitary_plan(&b, &fockplan::Plan::phaseshift(&b, p, m))
    } else if tag == atoms::kerr() {
        let (_, k, m): (Atom, f64, usize) = op.decode()?;
        if !in_range(&[m]) { return Ok(bad_shape(env, "mode")); }
        r.unitary_plan(&b, &fockplan::Plan::kerr(&b, k, m))
    } else if tag == atoms::passive() {
        let (_, u, ms): (Atom, Term, Vec<usize>) = op.decode()?;
        let (uv, ur, uc) = match cmat(u) { Ok(x) => x, Err(e) => return Ok(bad_shape(env, &e)) };
        if ur != uc || ur != ms.len() || !modes_ok(&ms, n) { return Ok(bad_shape(env, "U is k x k for k distinct modes")); }
        fockplan::Plan::passive(&b, &uv, &ms).and_then(|pl| r.unitary_plan(&b, &pl))
    } else if tag == atoms::displacement() {
        let (_, rr, p, m): (Atom, f64, f64, usize) = op.decode()?;
        if !in_range(&[m]) { return Ok(bad_shape(env, "mode")); }
        r.unitary_plan(&b, &fockplan::Plan::displacement(&b, rr, p, m))
    } else if tag == atoms::squeezing() {
        let (_, rr, p, m): (Atom, f64, f64, usize) = op.decode()?;
        if !in_range(&[m]) { return Ok(bad_shape(env, "mode")); }
        fockplan::Plan::squeezing(&b, rr, p, m).and_then(|pl| r.unitary_plan(&b, &pl))
    } else if tag == atoms::loss() {
        let (_, eta, m): (Atom, f64, usize) = op.decode()?;
        if !in_range(&[m]) || !(0.0..=1.0).contains(&eta) { return Ok(bad_shape(env, "mode or eta")); }
        r = r.loss(&b, eta, m);
        Ok(0.0)
    } else if tag == atoms::thermal_loss() {
        let (_, eta, nth, m): (Atom, f64, f64, usize) = op.decode()?;
        if !in_range(&[m]) || !(0.0..=1.0).contains(&eta) || nth.is_nan() || nth < 0.0 { return Ok(bad_shape(env, "mode, eta or nth")); }
        let (out, def) = r.thermal_loss(&b, eta, nth, m);
        r = out;
        Ok(def)
    } else {
        return Err(rustler::Error::BadArg);
    };
    Ok(match res {
        Ok(def) => (atoms::ok(), (put_blocks(env, &r), def, r.entries())).encode(env),
        Err(e) => fock_error(env, e),
    })
}

/// fock_dmb_pnr_probabilities(Blocks, Modes, N, C) -> [{Pattern, Prob}] (FR-F03 on blocks)
#[rustler::nif(schedule = "DirtyCpu")]
fn fock_dmb_pnr_probabilities<'a>(env: Env<'a>, blocks: Term<'a>, modes: Vec<usize>, n: usize, c: usize) -> NifResult<Term<'a>> {
    let Some(b) = dmb_basis(n, c) else { return Ok(bad_shape(env, "modes and cutoff")) };
    if modes.is_empty() || !modes_ok(&modes, n) { return Ok(bad_shape(env, "distinct modes in range")); }
    let blk = match get_blocks(blocks, &b) { Ok(x) => x, Err(e) => return Ok(bad_shape(env, &e)) };
    let diag: Vec<C64> = blk.diagonal().into_iter().map(|x| C64::new(x.max(0.0).sqrt(), 0.0)).collect();
    Ok((atoms::ok(), fockstate::pnr_probabilities(&b, &diag, &modes)).encode(env))
}

/// fock_dmb_measure_pnr(Blocks, Modes, Outcome, N, C) -> {Prob, Blocks} (FR-M07 on blocks)
#[rustler::nif(schedule = "DirtyCpu")]
fn fock_dmb_measure_pnr<'a>(env: Env<'a>, blocks: Term<'a>, modes: Vec<usize>, outcome: Vec<u32>, n: usize, c: usize) -> NifResult<Term<'a>> {
    let Some(b) = dmb_basis(n, c) else { return Ok(bad_shape(env, "modes and cutoff")) };
    if modes.is_empty() || modes.len() != outcome.len() || !modes_ok(&modes, n) {
        return Ok(bad_shape(env, "distinct modes in range, one outcome each"));
    }
    let mut blk = match get_blocks(blocks, &b) { Ok(x) => x, Err(e) => return Ok(bad_shape(env, &e)) };
    let prob = blk.measure(&b, &modes, &outcome);
    Ok((atoms::ok(), (prob, put_blocks(env, &blk))).encode(env))
}

/// fock_adjoint_gradient(Amps, Gates, Terms, NModes, Cutoff) -> {TermValues, GateGradients} (Phase 4, FS 15 item 3):
/// the gates as {beamsplitter, T, P, M1, M2}, {phaseshift, P, M}, {kerr, K, M}, {passive, U, Ms},
/// {displacement, R, P, M}, {squeezing, R, P, M}; the terms {probability, W, Modes, Pattern} or {photons, W, Modes};
/// the gradient of sum W f, one list per gate in its parameters' order (fockgrad.rs)
#[rustler::nif(schedule = "DirtyCpu")]
fn fock_adjoint_gradient<'a>(env: Env<'a>, amps: Term<'a>, gates: Vec<Term<'a>>, terms: Vec<Term<'a>>, n: usize, c: usize) -> NifResult<Term<'a>> {
    let (v, b) = fock_args!(env, amps, n, c, [0usize]);
    let mut gs = Vec::with_capacity(gates.len());
    for op in gates {
        match decode_gate(op, n)? {
            Ok(g) => gs.push(g),
            Err(m) => return Ok(bad_shape(env, m)),
        }
    }
    let ts = match decode_terms(terms, n)? {
        Ok(t) => t,
        Err(m) => return Ok(bad_shape(env, m)),
    };
    Ok(match fockgrad::adjoint(&b, &v, &gs, &ts) {
        Ok((values, grads)) => (atoms::ok(), (values, grads)).encode(env),
        Err(e) => fock_error(env, e),
    })
}

/// a gate tuple of the adjoint NIFs: Ok(Err(message)) for a bad shape
fn decode_gate(op: Term, n: usize) -> NifResult<Result<fockgrad::Gate, &'static str>> {
    let tag: Atom = rustler::types::tuple::get_tuple(op)?.first().ok_or(rustler::Error::BadArg)?.decode()?;
    let in_range = |ms: &[usize]| ms.iter().all(|&m| m < n);
    Ok(if tag == atoms::beamsplitter() {
        let (_, theta, phi, ma, mb): (Atom, f64, f64, usize, usize) = op.decode()?;
        if !in_range(&[ma, mb]) || ma == mb { Err("modes") } else { Ok(fockgrad::Gate::Beamsplitter { theta, phi, ma, mb }) }
    } else if tag == atoms::phaseshift() {
        let (_, phi, mode): (Atom, f64, usize) = op.decode()?;
        if !in_range(&[mode]) { Err("mode") } else { Ok(fockgrad::Gate::Phaseshift { phi, mode }) }
    } else if tag == atoms::kerr() {
        let (_, kappa, mode): (Atom, f64, usize) = op.decode()?;
        if !in_range(&[mode]) { Err("mode") } else { Ok(fockgrad::Gate::Kerr { kappa, mode }) }
    } else if tag == atoms::passive() {
        let (_, u, modes): (Atom, Term, Vec<usize>) = op.decode()?;
        match cmat(u) {
            Ok((uv, ur, uc)) if ur == uc && ur == modes.len() && modes_ok(&modes, n) => Ok(fockgrad::Gate::Passive { u: uv, modes }),
            _ => Err("U is k x k for k distinct modes"),
        }
    } else if tag == atoms::displacement() {
        let (_, r, phi, mode): (Atom, f64, f64, usize) = op.decode()?;
        if !in_range(&[mode]) { Err("mode") } else { Ok(fockgrad::Gate::Displacement { r, phi, mode }) }
    } else if tag == atoms::squeezing() {
        let (_, r, phi, mode): (Atom, f64, f64, usize) = op.decode()?;
        if !in_range(&[mode]) { Err("mode") } else { Ok(fockgrad::Gate::Squeezing { r, phi, mode }) }
    } else {
        return Err(rustler::Error::BadArg);
    })
}

fn decode_terms(terms: Vec<Term>, n: usize) -> NifResult<Result<Vec<fockgrad::Term>, &'static str>> {
    let mut ts = Vec::with_capacity(terms.len());
    for t in terms {
        let tag: Atom = rustler::types::tuple::get_tuple(t)?.first().ok_or(rustler::Error::BadArg)?.decode()?;
        if tag == atoms::probability() {
            let (_, weight, modes, pattern): (Atom, f64, Vec<usize>, Vec<u32>) = t.decode()?;
            if !modes_ok(&modes, n) || modes.len() != pattern.len() { return Ok(Err("modes and pattern")); }
            ts.push(fockgrad::Term::Probability { weight, modes, pattern });
        } else if tag == atoms::photons() {
            let (_, weight, modes): (Atom, f64, Vec<usize>) = t.decode()?;
            if !modes_ok(&modes, n) { return Ok(Err("modes")); }
            ts.push(fockgrad::Term::Photons { weight, modes });
        } else {
            return Err(rustler::Error::BadArg);
        }
    }
    Ok(Ok(ts))
}

/// fock_dm_adjoint_gradient(Amps, Channels, Terms, NModes, Cutoff) -> {TermValues, ChannelGradients}: the adjoint
/// on density matrices (fockgrad::adjoint_dm), from the pure input Amps; channels the gates of
/// fock_adjoint_gradient, {loss, Eta, M} (one gradient) and {thermal_loss, Eta, Nth, M} (none)
#[rustler::nif(schedule = "DirtyCpu")]
fn fock_dm_adjoint_gradient<'a>(env: Env<'a>, amps: Term<'a>, chans: Vec<Term<'a>>, terms: Vec<Term<'a>>, n: usize, c: usize) -> NifResult<Term<'a>> {
    let (v, b) = fock_args!(env, amps, n, c, [0usize]);
    let mut cs = Vec::with_capacity(chans.len());
    for op in chans {
        let tag: Atom = rustler::types::tuple::get_tuple(op)?.first().ok_or(rustler::Error::BadArg)?.decode()?;
        if tag == atoms::loss() {
            let (_, eta, m): (Atom, f64, usize) = op.decode()?;
            if m >= n || !(0.0..=1.0).contains(&eta) { return Ok(bad_shape(env, "mode or eta")); }
            cs.push(fockgrad::Channel::Loss { eta, mode: m });
        } else if tag == atoms::thermal_loss() {
            let (_, eta, nth, m): (Atom, f64, f64, usize) = op.decode()?;
            if m >= n || !(0.0..=1.0).contains(&eta) || nth.is_nan() || nth < 0.0 { return Ok(bad_shape(env, "mode, eta or nth")); }
            cs.push(fockgrad::Channel::Thermal { eta, nth, mode: m });
        } else {
            match decode_gate(op, n)? {
                Ok(g) => cs.push(fockgrad::Channel::Gate(g)),
                Err(m) => return Ok(bad_shape(env, m)),
            }
        }
    }
    let ts = match decode_terms(terms, n)? {
        Ok(t) => t,
        Err(m) => return Ok(bad_shape(env, m)),
    };
    Ok(match fockgrad::adjoint_dm(&b, &v, &cs, &ts) {
        Ok((values, grads)) => (atoms::ok(), (values, grads)).encode(env),
        Err(e) => fock_error(env, e),
    })
}

/// sym_poly_mul(Vars, ExpsA, CoefA, ExpsB, CoefB, L) -> {ok, Poly} | overflow | packing (the symbolic engine's plain
/// products, polymul.rs): a polynomial as its terms' exponents, one per variable of Vars (sorted atom keys), and their
/// coefficients times L, two integers (real, imaginary) per term; the product as pq_sym_poly's map, its
/// coefficients over L * L in lowest terms
#[rustler::nif(schedule = "DirtyCpu")]
fn sym_poly_mul<'a>(env: Env<'a>, vars: Vec<Term<'a>>, exps_a: Vec<i64>, coef_a: Vec<i128>, exps_b: Vec<i64>, coef_b: Vec<i128>, denom: i128) -> NifResult<Term<'a>> {
    let nvars = vars.len();
    let shape = |e: &Vec<i64>, c: &Vec<i128>| c.len().is_multiple_of(2) && e.len() == nvars * (c.len() / 2);
    if !shape(&exps_a, &coef_a) || !shape(&exps_b, &coef_b) || denom < 1 {
        return Err(rustler::Error::BadArg);
    }
    let a = polymul::Poly { exps: exps_a, coef: coef_a };
    let b = polymul::Poly { exps: exps_b, coef: coef_b };
    sym_poly_term(env, &vars, polymul::mul(nvars, &a, &b), denom)
}

/// sym_poly_pow(Vars, Exps, Coef, N, D) -> {ok, Poly} | overflow | packing: a polynomial as for sym_poly_mul, its
/// coefficients times L, to the power N; D = L^N, the coefficients of the result over D in lowest terms
#[rustler::nif(schedule = "DirtyCpu")]
fn sym_poly_pow<'a>(env: Env<'a>, vars: Vec<Term<'a>>, exps: Vec<i64>, coef: Vec<i128>, n: u32, denom: i128) -> NifResult<Term<'a>> {
    let nvars = vars.len();
    if !coef.len().is_multiple_of(2) || exps.len() != nvars * (coef.len() / 2) || denom < 1 {
        return Err(rustler::Error::BadArg);
    }
    let a = polymul::Poly { exps, coef };
    sym_poly_term(env, &vars, polymul::pow(nvars, &a, n), denom)
}

/// a kernel result as pq_sym_poly's map, or the atom of its refusal
fn sym_poly_term<'a>(env: Env<'a>, vars: &[Term<'a>], r: Result<polymul::Poly, polymul::Refused>, denom: i128) -> NifResult<Term<'a>> {
    match r {
        Ok(p) => Ok((atoms::ok(), poly_map(env, vars, &p, denom)?).encode(env)),
        Err(e) => Ok(refused(env, e)),
    }
}

fn refused(env: Env, e: polymul::Refused) -> Term {
    match e {
        polymul::Refused::Overflow => atoms::overflow().encode(env),
        polymul::Refused::Packing => atoms::packing().encode(env),
    }
}

/// n / denom in lowest terms with a positive denominator, as pq_sym_exact:rat/2 makes it
fn rat_parts(n: i128, denom: i128) -> (i128, i128) {
    let g = polymul::gcd(n.unsigned_abs(), denom as u128) as i128;
    (n / g, denom / g)
}

/// A variable of the kernel: an atom key of a monomial's atom list, or {expi, X}, the exponential factor
/// exp(k i X) of the monomial's exponential list (Base the monomial X, Sym the term {sym, X}).
enum VarKind<'a> {
    Atom(Term<'a>),
    Exp { base: Term<'a>, sym: Term<'a> },
}

fn var_kinds<'a>(env: Env<'a>, vars: &[Term<'a>]) -> NifResult<Vec<VarKind<'a>>> {
    let nil: Term<'a> = Vec::<Term<'a>>::new().encode(env);
    vars.iter()
        .map(|&v| {
            if let Ok(t) = rustler::types::tuple::get_tuple(v) {
                if t.len() == 2 && t[0].decode::<Atom>().is_ok_and(|a| a == atoms::expi()) {
                    let sym = (atoms::sym(), t[1]).encode(env);
                    let base = rustler::types::tuple::make_tuple(env, &[vec![(sym, (1i64, 1i64)).encode(env)].encode(env), nil]);
                    return Ok(VarKind::Exp { base, sym });
                }
            }
            Ok(VarKind::Atom(v))
        })
        .collect()
}

/// terms as pq_sym_poly's map: monomials {Atoms, Exponentials} from Vars ({V, {E, 1}} for an atom key, {Base,
/// {q, {0, 1}, {E, 1}}} for {expi, X}), coefficients over Denom in lowest terms
fn poly_map<'a>(env: Env<'a>, vars: &[Term<'a>], p: &polymul::Poly, denom: i128) -> NifResult<Term<'a>> {
    let nvars = vars.len();
    let kinds = var_kinds(env, vars)?;
    let q = atoms::q().encode(env);
    let zero = (0i64, 1i64).encode(env);
    let mut factor_terms: HashMap<(usize, i64), Term<'a>> = HashMap::new();
    let n = p.coef.len() / 2;
    let mut keys = Vec::with_capacity(n);
    let mut values = Vec::with_capacity(n);
    for t in 0..n {
        let mut alg = Vec::with_capacity(nvars);
        let mut exp = Vec::new();
        for (j, kind) in kinds.iter().enumerate() {
            let e = p.exps[t * nvars + j];
            if e != 0 {
                let f = *factor_terms.entry((j, e)).or_insert_with(|| match kind {
                    VarKind::Atom(v) => (*v, (e, 1i64)).encode(env),
                    VarKind::Exp { base, .. } => (*base, rustler::types::tuple::make_tuple(env, &[q, zero, (e, 1i64).encode(env)])).encode(env),
                });
                match kind {
                    VarKind::Atom(_) => alg.push(f),
                    VarKind::Exp { .. } => exp.push(f),
                }
            }
        }
        keys.push(rustler::types::tuple::make_tuple(env, &[alg.encode(env), exp.encode(env)]));
        values.push(rustler::types::tuple::make_tuple(
            env,
            &[q, rat_parts(p.coef[2 * t], denom).encode(env), rat_parts(p.coef[2 * t + 1], denom).encode(env)],
        ));
    }
    Term::map_from_term_arrays(env, &keys, &values)
}

/// A polynomial kept in Rust (rpoly.rs; pq_sym_rpoly holds its variables).
pub struct RPolyRes(rpoly::RPoly);

#[rustler::resource_impl]
impl rustler::Resource for RPolyRes {}

fn rpoly_result<'a>(env: Env<'a>, r: Result<rpoly::RPoly, polymul::Refused>) -> Term<'a> {
    match r {
        Ok(p) => (atoms::ok(), ResourceArc::new(RPolyRes(p))).encode(env),
        Err(e) => refused(env, e),
    }
}

/// rpoly_new(NVars, Exps, Coef, Denom) -> Resource: terms as for sym_poly_mul, their coefficients over Denom
#[rustler::nif(schedule = "DirtyCpu")]
fn rpoly_new(nvars: usize, exps: Vec<i64>, coef: Vec<i128>, denom: i128) -> NifResult<ResourceArc<RPolyRes>> {
    if !coef.len().is_multiple_of(2) || exps.len() != nvars * (coef.len() / 2) || denom < 1 {
        return Err(rustler::Error::BadArg);
    }
    Ok(ResourceArc::new(RPolyRes(rpoly::RPoly::new(nvars, polymul::Poly { exps, coef }, denom))))
}

fn maps_ok(a: &rpoly::RPoly, ma: &[usize], b: &rpoly::RPoly, mb: &[usize], n: usize) -> bool {
    ma.len() == a.nvars && mb.len() == b.nvars && ma.iter().chain(mb).all(|&j| j < n)
}

/// rpoly_mul(A, MapA, B, MapB, N) -> {ok, Resource} | overflow | packing: A's variable j is the result's MapA[j]
#[rustler::nif(schedule = "DirtyCpu")]
fn rpoly_mul<'a>(env: Env<'a>, a: ResourceArc<RPolyRes>, ma: Vec<usize>, b: ResourceArc<RPolyRes>, mb: Vec<usize>, n: usize) -> NifResult<Term<'a>> {
    if !maps_ok(&a.0, &ma, &b.0, &mb, n) {
        return Err(rustler::Error::BadArg);
    }
    Ok(rpoly_result(env, rpoly::mul(&a.0, &ma, &b.0, &mb, n)))
}

/// rpoly_add(A, MapA, B, MapB, N) -> {ok, Resource} | overflow | packing
#[rustler::nif(schedule = "DirtyCpu")]
fn rpoly_add<'a>(env: Env<'a>, a: ResourceArc<RPolyRes>, ma: Vec<usize>, b: ResourceArc<RPolyRes>, mb: Vec<usize>, n: usize) -> NifResult<Term<'a>> {
    if !maps_ok(&a.0, &ma, &b.0, &mb, n) {
        return Err(rustler::Error::BadArg);
    }
    Ok(rpoly_result(env, rpoly::add(&a.0, &ma, &b.0, &mb, n)))
}

/// rpoly_pow(A, K) -> {ok, Resource} | overflow | packing
#[rustler::nif(schedule = "DirtyCpu")]
fn rpoly_pow(env: Env, a: ResourceArc<RPolyRes>, k: u32) -> Term {
    rpoly_result(env, rpoly::pow(&a.0, k))
}

/// rpoly_scale(A, Re, Im, Den) -> {ok, Resource} | overflow: A times (Re + i Im) / Den
#[rustler::nif(schedule = "DirtyCpu")]
fn rpoly_scale(env: Env, a: ResourceArc<RPolyRes>, re: i128, im: i128, den: i128) -> NifResult<Term> {
    if den < 1 {
        return Err(rustler::Error::BadArg);
    }
    Ok(rpoly_result(env, rpoly::scale(&a.0, re, im, den)))
}

/// rpoly_info(A) -> {Terms, NVars}
#[rustler::nif]
fn rpoly_info(a: ResourceArc<RPolyRes>) -> (usize, usize) {
    (a.0.len(), a.0.nvars)
}

/// rpoly_to_map(A, Vars) -> pq_sym_poly's map
#[rustler::nif(schedule = "DirtyCpu")]
fn rpoly_to_map<'a>(env: Env<'a>, a: ResourceArc<RPolyRes>, vars: Vec<Term<'a>>) -> NifResult<Term<'a>> {
    if vars.len() != a.0.nvars {
        return Err(rustler::Error::BadArg);
    }
    poly_map(env, &vars, &a.0.p, a.0.denom)
}

/// rpoly_to_ast(A, Vars) -> pq_sym_norm:to_ast/1's result, for variables that are all {sym, X} or {expi, X}: the
/// terms in print order, a term {mul, [Coeff | Factors]} (a constant just its coefficient), the factors {sym, X} or
/// {pow, {sym, X}, E} and then {exp, {mul, [{complex, 0, K}, {sym, X}]}}, a coefficient N, {rat, N, D} or
/// {complex, Re, Im}; no term 0, one term itself, more {add, Terms}
#[rustler::nif(schedule = "DirtyCpu")]
fn rpoly_to_ast<'a>(env: Env<'a>, a: ResourceArc<RPolyRes>, vars: Vec<Term<'a>>) -> NifResult<Term<'a>> {
    let r = &a.0;
    if vars.len() != r.nvars {
        return Err(rustler::Error::BadArg);
    }
    let kinds = var_kinds(env, &vars)?;
    let is_exp: Vec<bool> = kinds.iter().map(|k| matches!(k, VarKind::Exp { .. })).collect();
    let rat_ast = |x: i128| -> Term<'a> {
        let (n, d) = rat_parts(x, r.denom);
        if d == 1 { n.encode(env) } else { (atoms::rat(), n, d).encode(env) }
    };
    let (mul, pow) = (atoms::mul().encode(env), atoms::pow().encode(env));
    let mut factor_terms: HashMap<(usize, i64), Term<'a>> = HashMap::new();
    let mut out = Vec::with_capacity(r.len());
    for t in rpoly::print_order(r, &is_exp) {
        let (re, im) = (r.p.coef[2 * t], r.p.coef[2 * t + 1]);
        let c = if im == 0 { rat_ast(re) } else { (atoms::complex(), rat_ast(re), rat_ast(im)).encode(env) };
        let mut fs = vec![c];
        // the atoms first, then the exponentials, each in variable order
        for pass_exp in [false, true] {
            for (j, kind) in kinds.iter().enumerate() {
                let e = r.p.exps[t * r.nvars + j];
                if e == 0 || is_exp[j] != pass_exp {
                    continue;
                }
                fs.push(*factor_terms.entry((j, e)).or_insert_with(|| match kind {
                    VarKind::Atom(v) if e == 1 => *v,
                    VarKind::Atom(v) => rustler::types::tuple::make_tuple(env, &[pow, *v, e.encode(env)]),
                    VarKind::Exp { sym, .. } => {
                        let arg = rustler::types::tuple::make_tuple(env, &[mul, vec![(atoms::complex(), 0i64, e).encode(env), *sym].encode(env)]);
                        (atoms::exp(), arg).encode(env)
                    }
                }));
            }
        }
        out.push(if fs.len() == 1 { c } else { rustler::types::tuple::make_tuple(env, &[mul, fs.encode(env)]) });
    }
    Ok(match out.len() {
        0 => 0.encode(env),
        1 => out[0],
        _ => (atoms::add(), out).encode(env),
    })
}

/// Polynomials prepared for counting their Horner operations (hcount.rs; H4).
pub struct HornerRes(hcount::Prepared);

#[rustler::resource_impl]
impl rustler::Resource for HornerRes {}

/// horner_prepare(NVars, Polys) -> {ok, Resource} | packing | overflow: plain polynomials of symbols ({Exps, Coef,
/// Den} each, as for sym_poly_mul, over the same variables), packed once for horner_count/2
#[rustler::nif(schedule = "DirtyCpu")]
fn horner_prepare<'a>(env: Env<'a>, nvars: usize, polys: Vec<(Vec<i64>, Vec<i128>, i128)>) -> NifResult<Term<'a>> {
    let mut ps = Vec::with_capacity(polys.len());
    for (exps, coef, den) in polys {
        if !coef.len().is_multiple_of(2) || exps.len() != nvars * (coef.len() / 2) || den < 1 {
            return Err(rustler::Error::BadArg);
        }
        ps.push((polymul::Poly { exps, coef }, den));
    }
    Ok(match hcount::prepare(nvars, &ps) {
        Ok(p) => (atoms::ok(), ResourceArc::new(HornerRes(p))).encode(env),
        Err(e) => refused(env, e),
    })
}

/// horner_count(Resource, Order) -> {ok, Ops} | overflow: the operations of pq_sym_slp's Horner programs of the
/// prepared polynomials for the variable order Order (variable indices, first to last), without building them
#[rustler::nif(schedule = "DirtyCpu")]
fn horner_count<'a>(env: Env<'a>, prep: ResourceArc<HornerRes>, order: Vec<usize>) -> NifResult<Term<'a>> {
    if order.iter().any(|&j| j >= prep.0.nvars_of()) {
        return Err(rustler::Error::BadArg);
    }
    Ok(match hcount::count(&prep.0, &order) {
        Ok(n) => (atoms::ok(), n).encode(env),
        Err(e) => refused(env, e),
    })
}

/// A straight-line program prepared for evaluation over prime fields (ffslp.rs; R1).
pub struct FfRes(ffslp::Program);

#[rustler::resource_impl]
impl rustler::Resource for FfRes {}

/// ff_prepare(NInputs, Ops, Outputs) -> {ok, Resource} | {error, bad_shape}: Ops in slot order, each {input, I} |
/// {const, Sign, NumLimbs, DenLimbs} (64-bit limbs, least significant first) | {mul, A, B} | {add, A, B} | {sub, A, B} | {pow, A, E} (A, B earlier slots; E an integer)
#[rustler::nif(schedule = "DirtyCpu")]
fn ff_prepare<'a>(env: Env<'a>, ninputs: usize, ops: Vec<Term<'a>>, outputs: Vec<usize>) -> NifResult<Term<'a>> {
    let mut prog = Vec::with_capacity(ops.len());
    for op in ops {
        let t = rustler::types::tuple::get_tuple(op)?;
        let tag: Atom = t.first().ok_or(rustler::Error::BadArg)?.decode()?;
        let arg = |i: usize| -> NifResult<Term<'a>> { t.get(i).copied().ok_or(rustler::Error::BadArg) };
        prog.push(if tag == atoms::input() {
            ffslp::Op::Input(arg(1)?.decode()?)
        } else if tag == atoms::const_() {
            ffslp::Op::Const(arg(1)?.decode()?, arg(2)?.decode()?, arg(3)?.decode()?)
        } else if tag == atoms::mul() {
            ffslp::Op::Mul(arg(1)?.decode()?, arg(2)?.decode()?)
        } else if tag == atoms::add() {
            ffslp::Op::Add(arg(1)?.decode()?, arg(2)?.decode()?)
        } else if tag == atoms::sub() {
            ffslp::Op::Sub(arg(1)?.decode()?, arg(2)?.decode()?)
        } else if tag == atoms::pow() {
            ffslp::Op::Pow(arg(1)?.decode()?, arg(2)?.decode()?)
        } else {
            return Ok(bad_shape(env, "ff op"));
        });
    }
    let p = ffslp::Program { ops: prog, outputs, ninputs };
    Ok(match p.check() {
        Ok(()) => (atoms::ok(), ResourceArc::new(FfRes(p))).encode(env),
        Err(_) => bad_shape(env, "ff program"),
    })
}

/// ff_eval(Resource, P, Points) -> {ok, Values} | bad_prime: the outputs mod the prime P (< 2^63) at every point
/// (Points: NInputs residues per point, flat), output by output and point by point, flat; a point where the
/// program divides by zero has every output 18446744073709551615
#[rustler::nif(schedule = "DirtyCpu")]
fn ff_eval<'a>(env: Env<'a>, prog: ResourceArc<FfRes>, p: u64, points: Vec<u64>) -> NifResult<Term<'a>> {
    if !(3..(1u64 << 63)).contains(&p) {
        return Err(rustler::Error::BadArg);
    }
    Ok(match prog.0.eval(p, &points) {
        Ok(v) => (atoms::ok(), v).encode(env),
        Err(ffslp::FfError::BadPrime) => atoms::bad_prime().encode(env),
        Err(ffslp::FfError::Shape) => bad_shape(env, "ff points"),
    })
}

/// ff_tvandermonde(P, Nodes, Values) -> {ok, Coefficients} | singular: the c with Values[i] = sum_j c_j Nodes_j^i
/// (i from 0) mod the prime P, for distinct nodes (Zippel's interpolation, R1)
#[rustler::nif(schedule = "DirtyCpu")]
fn ff_tvandermonde<'a>(env: Env<'a>, p: u64, nodes: Vec<u64>, values: Vec<u64>) -> NifResult<Term<'a>> {
    if !(3..(1u64 << 63)).contains(&p) || nodes.len() != values.len() {
        return Err(rustler::Error::BadArg);
    }
    Ok(match ffslp::tvandermonde(p, &nodes, &values) {
        Some(c) => (atoms::ok(), c).encode(env),
        None => atoms::singular().encode(env),
    })
}

/// ff_column(Resource, P, R, N, Support, Suffix) -> {ok, Coefficients} | bad | singular | bad_prime: one column of
/// Zippel's interpolation (ffslp.rs; R1): the coefficients of the N monomials of Support (exponent rows of
/// length(R), flat) in the program's first output at the points (R^i, Suffix), i = 1 .. N
#[rustler::nif(schedule = "DirtyCpu")]
fn ff_column<'a>(env: Env<'a>, prog: ResourceArc<FfRes>, p: u64, r: Vec<u64>, n: usize, support: Vec<u32>, suffix: Vec<u64>) -> NifResult<Term<'a>> {
    if !(3..(1u64 << 63)).contains(&p) {
        return Err(rustler::Error::BadArg);
    }
    Ok(match prog.0.column(p, &r, n, &support, &suffix) {
        Ok(c) => (atoms::ok(), c).encode(env),
        Err(ffslp::ColumnError::Bad) => atoms::bad().encode(env),
        Err(ffslp::ColumnError::Singular) => atoms::singular().encode(env),
        Err(ffslp::ColumnError::BadPrime) => atoms::bad_prime().encode(env),
        Err(ffslp::ColumnError::Shape) => bad_shape(env, "ff column"),
    })
}

/// ff_line_column(Resource, P, Shift, R, N, Support, Suffix, Ts, Dn, Dd) -> {ok, [Coefficients]} | bad | singular |
/// bad_prime: one column of the multi-output Zippel interpolation of a rational function's coefficients along lines
/// (ffslp.rs; R1 phase 4): Dn + Dd + 1 rows (n_0 .. n_Dn, d_1 .. d_Dd) of N values each
#[rustler::nif(schedule = "DirtyCpu")]
#[allow(clippy::too_many_arguments)]
fn ff_line_column<'a>(env: Env<'a>, prog: ResourceArc<FfRes>, p: u64, shift: Vec<u64>, r: Vec<u64>, n: usize, support: Vec<u32>,
                      suffix: Vec<u64>, ts: Vec<u64>, dn: usize, dd: usize) -> NifResult<Term<'a>> {
    if !(3..(1u64 << 63)).contains(&p) {
        return Err(rustler::Error::BadArg);
    }
    Ok(match prog.0.line_column(p, &shift, &r, n, &support, &suffix, &ts, dn, dd) {
        Ok(c) => (atoms::ok(), c).encode(env),
        Err(ffslp::ColumnError::Bad) => atoms::bad().encode(env),
        Err(ffslp::ColumnError::Singular) => atoms::singular().encode(env),
        Err(ffslp::ColumnError::BadPrime) => atoms::bad_prime().encode(env),
        Err(ffslp::ColumnError::Shape) => bad_shape(env, "ff line column"),
    })
}

rustler::init!("pq_nif");
