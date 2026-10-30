//! pq_nif: the numerical kernels of the Photonic QC Framework (FS section 8).
use rustler::{Atom, Binary, Encoder, Env, NifResult, ResourceArc, Term};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::atomic::{AtomicBool, Ordering};

pub mod boson;
pub mod fock;
pub mod gaussian;
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
        size, rvec, rmat, not_symplectic, homodyne, pnr, threshold, unsupported_input
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
    let phi = match detector.decode::<(Atom, f64)>() {
        Ok((tag, phi)) if tag == atoms::homodyne() => phi,
        _ => return Ok(error(env, (atoms::unsupported_input(), detector))),
    };
    let seed: u64 = opt(opts, atoms::seed(), 0u64)?;
    let start: u64 = opt(opts, atoms::stream_start(), 0u64)?;
    if phi != 0.0 {
        let rot = C64::new((-phi).cos(), (-phi).sin());
        for j in 0..m {
            let sm = gaussian::passive_symplectic(&[rot], 1);
            gaussian::apply_symplectic(&mut mu, &mut sig, m, &sm, &[0.0, 0.0], &[j], None).expect("rotation");
        }
    }
    let cancel = cancelled_fn(&o);
    Ok(guarded(env, || match gaussian::sample_homodyne(&mu, &sig, m, shots, seed, start, &cancel) {
        Ok(y) => (atoms::ok(), (put_rmat(env, &y, shots, m), Term::map_new(env))).encode(env),
        Err(gaussian::Error::NonPhysical(l)) => error(env, (atoms::non_physical_state(), l)),
        Err(gaussian::Error::Cancelled) => error(env, atoms::cancelled()),
        Err(_) => bad_shape(env, "state"),
    }))
}

rustler::init!("pq_nif");
