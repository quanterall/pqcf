//! GPU offload (Phase 4, FS 15 item 5; docs/decisions/G1-gpu-offload.md): Glynn's permanent and the (loop)
//! hafnian on an NVIDIA GPU through the CUDA driver API, loaded at run time.
//!
//! The kernel is gpu/perm.cu, compiled once to PTX with -fmad=false (gpu/build_ptx.sh) and embedded; the driver
//! JIT-compiles it for the device. Nothing CUDA is needed to build the crate: without libcuda.so.1 or a device every
//! call returns `Unavailable`, and the CPU kernels are used.
//!
//! Determinism: thread t computes the CPU kernel's positions [4096 t, 4096 (t+1)), starting where the CPU kernel
//! refreshes, so every term equals the CPU term bit for bit; the per-thread sums (Neumaier) are added here in thread
//! order. Results are reproducible bit for bit on any device and launch, and equal to the CPU kernel within the
//! rounding of the regrouped sum (a few ulps of the scale S).
//!
//! NFR-04: the kernel logic is the PTX and the safe code below; `unsafe` is confined to the `driver` submodule, the
//! FFI boundary to libcuda, like the Rustler glue in lib.rs.
use crate::num::{CSum, Sum, C64};
use crate::hafnian::HafResult;
use crate::perm::PermResult;
use std::sync::OnceLock;

const PTX: &str = include_str!("../gpu/perm.ptx");
const HAF_PTX: &str = include_str!("../gpu/haf.ptx");
/// Positions per GPU thread: the CPU kernel's refresh interval (perm::REFRESH).
pub const THREAD_POSITIONS: u64 = 4096;
const MAX_N: usize = 40;

#[derive(Debug, PartialEq)]
pub enum Error {
    /// no driver, no device, or the PTX did not load; the reason
    Unavailable(String),
    BadShape(&'static str),
    NonFinite,
    Cancelled,
}

/// The device's name, or why there is none.
pub fn device() -> Result<String, Error> {
    driver::get().map(|d| d.name.clone())
}

/// Threads per launch (2^28 positions): between launches the cancel token is checked.
const LAUNCH_THREADS: u64 = 1 << 16;

/// perm(A) for a square n x n matrix on the GPU, with the scale S of AN 2.4.
pub fn permanent(a: &[C64], n: usize) -> Result<PermResult, Error> {
    permanent_cancellable(a, n, &|| false)
}

/// `permanent`, checking `cancel` between launches.
pub fn permanent_cancellable(a: &[C64], n: usize, cancel: &dyn Fn() -> bool) -> Result<PermResult, Error> {
    if a.len() != n * n {
        return Err(Error::BadShape("square matrix"));
    }
    if n == 0 {
        return Ok(PermResult { value: C64::ONE, scale: 1.0, positions: 0, blocks: 0 });
    }
    if n > MAX_N {
        return Err(Error::BadShape("n at most 40"));
    }
    if a.iter().any(|z| !z.is_finite()) {
        return Err(Error::NonFinite);
    }
    let positions: u64 = 1u64 << (n - 1);
    let threads = positions.div_ceil(THREAD_POSITIONS);
    let flat: Vec<f64> = a.iter().flat_map(|z| [z.re, z.im]).collect();
    let dev = driver::get()?;
    let mut partial = Vec::with_capacity(3 * threads as usize);
    let mut first = 0;
    while first < threads {
        if cancel() {
            return Err(Error::Cancelled);
        }
        let count = LAUNCH_THREADS.min(threads - first);
        partial.extend(dev.glynn(&flat, n, positions, first, count)?);
        first += count;
    }
    let mut sum = CSum::new(threads);
    let mut scale = Sum::new(threads);
    for t in 0..threads as usize {
        sum.add(C64::new(partial[3 * t], partial[3 * t + 1]));
        scale.add(partial[3 * t + 2]);
    }
    let norm = 0.5f64.powi(n as i32);
    let value = sum.value().scale(norm);
    let scale = scale.value() * norm;
    if !value.is_finite() || !scale.is_finite() {
        return Err(Error::NonFinite);
    }
    Ok(PermResult { value, scale, positions, blocks: threads })
}

/// Subsets per hafnian thread (Gray-code order, added with Neumaier summation in the thread); 16 measured faster
/// than 64 on the GB10 (more threads in flight).
pub const HAF_GROUP: u64 = 16;
/// Hafnian threads per launch: at most 16384, and at most HAF_WORKSPACE_BYTES of workspace (allocated per launch).
const HAF_LAUNCH: u64 = 16384;
const HAF_WORKSPACE_BYTES: u64 = 256 << 20;
/// The largest number of pairs on the GPU (the workspace grows as (2n)^2 per thread).
pub const HAF_MAX_PAIRS: usize = 32;

/// complex entries of a hafnian thread's workspace for n pairs (the layout in gpu/haf.cu)
fn haf_workspace(n: usize) -> usize {
    let d = 2 * n;
    2 * d * d + (d + 1) * (d + 2) / 2 + 3 * d + 3 * n + 2
}

/// haf(A), or the loop hafnian with `loops` on the diagonal, of a symmetric m x m matrix on the GPU: the subset terms
/// of hafnian.rs bit for bit, regrouped in threads of HAF_GROUP subsets
pub fn hafnian(a: &[C64], m: usize, loops: Option<&[C64]>, cancel: &dyn Fn() -> bool) -> Result<HafResult, Error> {
    if a.len() != m * m || loops.is_some_and(|l| l.len() != m) {
        return Err(Error::BadShape("matrix and loop weights"));
    }
    if a.iter().chain(loops.unwrap_or(&[])).any(|z| !z.is_finite()) {
        return Err(Error::NonFinite);
    }
    if m % 2 == 1 {
        return match loops {
            None => Ok(HafResult { value: C64::ZERO, scale: 0.0, subsets: 0 }),
            Some(l) => {
                // as the CPU kernel: one more vertex with loop weight 1 and no edges
                let mp = m + 1;
                let mut ap = vec![C64::ZERO; mp * mp];
                for i in 0..m {
                    for j in 0..m {
                        ap[i * mp + j] = a[i * m + j];
                    }
                }
                let mut lp = l.to_vec();
                lp.push(C64::ONE);
                hafnian(&ap, mp, Some(&lp), cancel)
            }
        };
    }
    let n = m / 2;
    if n == 0 {
        return Ok(HafResult { value: C64::ONE, scale: 1.0, subsets: 0 });
    }
    if n > HAF_MAX_PAIRS {
        return Err(Error::BadShape("at most 32 pairs on the GPU"));
    }
    let total = 1u64 << n;
    let threads = total.div_ceil(HAF_GROUP);
    let fa: Vec<f64> = a.iter().flat_map(|z| [z.re, z.im]).collect();
    let fl: Option<Vec<f64>> = loops.map(|l| l.iter().flat_map(|z| [z.re, z.im]).collect());
    let stride = 2 * haf_workspace(n) as u64;
    let dev = driver::get()?;
    let mut partial = Vec::with_capacity(3 * threads as usize);
    let mut first = 0;
    while first < threads {
        if cancel() {
            return Err(Error::Cancelled);
        }
        let count = HAF_LAUNCH.min(HAF_WORKSPACE_BYTES / (8 * stride)).max(1).min(threads - first);
        partial.extend(dev.hafnian_terms(&fa, m, fl.as_deref(), n, HAF_GROUP, first, count, stride)?);
        first += count;
    }
    let mut sum = CSum::new(threads);
    let mut scale = Sum::new(threads);
    for t in 0..threads as usize {
        sum.add(C64::new(partial[3 * t], partial[3 * t + 1]));
        scale.add(partial[3 * t + 2]);
    }
    let value = sum.value();
    if !value.is_finite() {
        return Err(Error::NonFinite);
    }
    Ok(HafResult { value, scale: scale.value(), subsets: total })
}

/// The FFI boundary to the CUDA driver API (libcuda.so.1), opened once per process.
#[allow(unsafe_code)]
mod driver {
    use super::{Error, OnceLock, HAF_PTX, PTX};
    use libloading::{Library, Symbol};
    use std::ffi::{c_char, c_int, c_uint, c_void, CString};
    use std::sync::Mutex;

    type CuResult = c_int;
    type CuDevice = c_int;
    type CuContext = *mut c_void;
    type CuModule = *mut c_void;
    type CuFunction = *mut c_void;
    type CuDevPtr = u64;

    pub struct Device {
        pub name: String,
        // the library must outlive every function pointer taken from it
        _lib: Library,
        ctx: CuContext,
        func: CuFunction,
        haf: CuFunction,
        api: Api,
        lock: Mutex<()>,
    }

    // the context and function handles are used only under `lock`, with the context made current first
    unsafe impl Send for Device {}
    unsafe impl Sync for Device {}

    struct Api {
        ctx_set_current: unsafe extern "C" fn(CuContext) -> CuResult,
        mem_alloc: unsafe extern "C" fn(*mut CuDevPtr, usize) -> CuResult,
        mem_free: unsafe extern "C" fn(CuDevPtr) -> CuResult,
        htod: unsafe extern "C" fn(CuDevPtr, *const c_void, usize) -> CuResult,
        dtoh: unsafe extern "C" fn(*mut c_void, CuDevPtr, usize) -> CuResult,
        launch: unsafe extern "C" fn(CuFunction, c_uint, c_uint, c_uint, c_uint, c_uint, c_uint, c_uint, *mut c_void,
                                     *mut *mut c_void, *mut *mut c_void) -> CuResult,
        sync: unsafe extern "C" fn() -> CuResult,
    }

    static DEVICE: OnceLock<Result<Device, String>> = OnceLock::new();

    pub fn get() -> Result<&'static Device, Error> {
        DEVICE.get_or_init(open).as_ref().map_err(|e| Error::Unavailable(e.clone()))
    }

    fn check(r: CuResult, what: &str) -> Result<(), String> {
        if r == 0 { Ok(()) } else { Err(format!("{what}: CUDA error {r}")) }
    }

    fn open() -> Result<Device, String> {
        unsafe {
            let lib = Library::new("libcuda.so.1").map_err(|e| format!("libcuda.so.1: {e}"))?;
            macro_rules! sym {
                ($name:expr, $t:ty) => {{
                    let s: Symbol<$t> = lib.get($name).map_err(|e| format!("{}: {e}", String::from_utf8_lossy($name)))?;
                    *s
                }};
            }
            let init = sym!(b"cuInit", unsafe extern "C" fn(c_uint) -> CuResult);
            let count = sym!(b"cuDeviceGetCount", unsafe extern "C" fn(*mut c_int) -> CuResult);
            let device_get = sym!(b"cuDeviceGet", unsafe extern "C" fn(*mut CuDevice, c_int) -> CuResult);
            let device_name = sym!(b"cuDeviceGetName", unsafe extern "C" fn(*mut c_char, c_int, CuDevice) -> CuResult);
            let retain = sym!(b"cuDevicePrimaryCtxRetain", unsafe extern "C" fn(*mut CuContext, CuDevice) -> CuResult);
            let ctx_set_current = sym!(b"cuCtxSetCurrent", unsafe extern "C" fn(CuContext) -> CuResult);
            let module_load = sym!(b"cuModuleLoadData", unsafe extern "C" fn(*mut CuModule, *const c_void) -> CuResult);
            let get_function = sym!(b"cuModuleGetFunction", unsafe extern "C" fn(*mut CuFunction, CuModule, *const c_char) -> CuResult);
            let api = Api {
                ctx_set_current,
                mem_alloc: sym!(b"cuMemAlloc_v2", unsafe extern "C" fn(*mut CuDevPtr, usize) -> CuResult),
                mem_free: sym!(b"cuMemFree_v2", unsafe extern "C" fn(CuDevPtr) -> CuResult),
                htod: sym!(b"cuMemcpyHtoD_v2", unsafe extern "C" fn(CuDevPtr, *const c_void, usize) -> CuResult),
                dtoh: sym!(b"cuMemcpyDtoH_v2", unsafe extern "C" fn(*mut c_void, CuDevPtr, usize) -> CuResult),
                launch: sym!(b"cuLaunchKernel", unsafe extern "C" fn(CuFunction, c_uint, c_uint, c_uint, c_uint, c_uint, c_uint, c_uint,
                                                                       *mut c_void, *mut *mut c_void, *mut *mut c_void) -> CuResult),
                sync: sym!(b"cuCtxSynchronize", unsafe extern "C" fn() -> CuResult),
            };
            check(init(0), "cuInit")?;
            let mut n: c_int = 0;
            check(count(&mut n), "cuDeviceGetCount")?;
            if n == 0 {
                return Err("no CUDA device".into());
            }
            let mut dev: CuDevice = 0;
            check(device_get(&mut dev, 0), "cuDeviceGet")?;
            let mut buf = [0 as c_char; 128];
            check(device_name(buf.as_mut_ptr(), 128, dev), "cuDeviceGetName")?;
            let name = std::ffi::CStr::from_ptr(buf.as_ptr()).to_string_lossy().into_owned();
            let mut ctx: CuContext = std::ptr::null_mut();
            check(retain(&mut ctx, dev), "cuDevicePrimaryCtxRetain")?;
            check(ctx_set_current(ctx), "cuCtxSetCurrent")?;
            let ptx = CString::new(PTX).map_err(|_| "PTX contains a NUL byte".to_string())?;
            let mut module: CuModule = std::ptr::null_mut();
            check(module_load(&mut module, ptx.as_ptr() as *const c_void), "cuModuleLoadData")?;
            let fname = CString::new("glynn").expect("name");
            let mut func: CuFunction = std::ptr::null_mut();
            check(get_function(&mut func, module, fname.as_ptr()), "cuModuleGetFunction")?;
            let hptx = CString::new(HAF_PTX).map_err(|_| "PTX contains a NUL byte".to_string())?;
            let mut hmodule: CuModule = std::ptr::null_mut();
            check(module_load(&mut hmodule, hptx.as_ptr() as *const c_void), "cuModuleLoadData (hafnian)")?;
            let hname = CString::new("hafnian_terms").expect("name");
            let mut haf: CuFunction = std::ptr::null_mut();
            check(get_function(&mut haf, hmodule, hname.as_ptr()), "cuModuleGetFunction (hafnian)")?;
            Ok(Device { name, _lib: lib, ctx, func, haf, api, lock: Mutex::new(()) })
        }
    }

    impl Device {
        /// one launch of hafnian_terms over threads first .. first + threads - 1 (3 doubles each), with a workspace of
        /// `stride` doubles per thread
        #[allow(clippy::too_many_arguments)]
        pub fn hafnian_terms(&self, a: &[f64], m: usize, loops: Option<&[f64]>, n: usize, group: u64, first: u64, threads: u64,
                             stride: u64) -> Result<Vec<f64>, Error> {
            let _guard = self.lock.lock().map_err(|_| Error::Unavailable("device lock poisoned".into()))?;
            let fail = |e: String| Error::Unavailable(e);
            let mut out = vec![0.0f64; 3 * threads as usize];
            unsafe {
                check((self.api.ctx_set_current)(self.ctx), "cuCtxSetCurrent").map_err(fail)?;
                let sizes = [a.len() * 8, loops.map_or(0, |l| l.len() * 8), (stride * threads) as usize * 8, out.len() * 8];
                let mut ptrs = [0u64; 4];
                let mut taken = 0;
                let mut alloc_err = None;
                for (k, &bytes) in sizes.iter().enumerate() {
                    if bytes == 0 {
                        continue;
                    }
                    if let Err(e) = check((self.api.mem_alloc)(&mut ptrs[k], bytes), "cuMemAlloc") {
                        alloc_err = Some(e);
                        break;
                    }
                    taken = k + 1;
                }
                let free_all = |ptrs: &[u64; 4]| {
                    for &p in ptrs.iter().take(taken) {
                        if p != 0 {
                            (self.api.mem_free)(p);
                        }
                    }
                };
                if let Some(e) = alloc_err {
                    free_all(&ptrs);
                    return Err(fail(e));
                }
                let mut run = || -> Result<(), String> {
                    check((self.api.htod)(ptrs[0], a.as_ptr() as *const c_void, sizes[0]), "cuMemcpyHtoD")?;
                    if let Some(l) = loops {
                        check((self.api.htod)(ptrs[1], l.as_ptr() as *const c_void, sizes[1]), "cuMemcpyHtoD")?;
                    }
                    let (mut pa, mut mm, mut pl, mut nn) = (ptrs[0], m as c_int, ptrs[1], n as c_int);
                    let (mut grp, mut fst, mut thr, mut pw, mut st, mut po) = (group, first, threads, ptrs[2], stride, ptrs[3]);
                    let mut args: [*mut c_void; 10] = [
                        &mut pa as *mut u64 as *mut c_void,
                        &mut mm as *mut c_int as *mut c_void,
                        &mut pl as *mut u64 as *mut c_void,
                        &mut nn as *mut c_int as *mut c_void,
                        &mut grp as *mut u64 as *mut c_void,
                        &mut fst as *mut u64 as *mut c_void,
                        &mut thr as *mut u64 as *mut c_void,
                        &mut pw as *mut u64 as *mut c_void,
                        &mut st as *mut u64 as *mut c_void,
                        &mut po as *mut u64 as *mut c_void,
                    ];
                    let block: u64 = 64;
                    let grid = threads.div_ceil(block) as c_uint;
                    check((self.api.launch)(self.haf, grid, 1, 1, block as c_uint, 1, 1, 0, std::ptr::null_mut(),
                                            args.as_mut_ptr(), std::ptr::null_mut()), "cuLaunchKernel (hafnian)")?;
                    check((self.api.sync)(), "cuCtxSynchronize")?;
                    check((self.api.dtoh)(out.as_mut_ptr() as *mut c_void, ptrs[3], sizes[3]), "cuMemcpyDtoH")
                };
                let r = run();
                free_all(&ptrs);
                r.map_err(fail)?;
            }
            Ok(out)
        }

        /// one launch of the glynn kernel over threads first .. first + threads - 1: 3 doubles per thread
        pub fn glynn(&self, u: &[f64], n: usize, positions: u64, first: u64, threads: u64) -> Result<Vec<f64>, Error> {
            let _guard = self.lock.lock().map_err(|_| Error::Unavailable("device lock poisoned".into()))?;
            let fail = |e: String| Error::Unavailable(e);
            let mut out = vec![0.0f64; 3 * threads as usize];
            unsafe {
                check((self.api.ctx_set_current)(self.ctx), "cuCtxSetCurrent").map_err(fail)?;
                let (ubytes, obytes) = (u.len() * 8, out.len() * 8);
                let mut du: u64 = 0;
                let mut dout: u64 = 0;
                check((self.api.mem_alloc)(&mut du, ubytes), "cuMemAlloc").map_err(fail)?;
                if let Err(e) = check((self.api.mem_alloc)(&mut dout, obytes), "cuMemAlloc") {
                    (self.api.mem_free)(du);
                    return Err(fail(e));
                }
                let mut run = || -> Result<(), String> {
                    check((self.api.htod)(du, u.as_ptr() as *const c_void, ubytes), "cuMemcpyHtoD")?;
                    let mut nn: c_int = n as c_int;
                    let mut pos: u64 = positions;
                    let mut fst: u64 = first;
                    let mut thr: u64 = threads;
                    let mut pu = du;
                    let mut po = dout;
                    let mut args: [*mut c_void; 6] = [
                        &mut pu as *mut u64 as *mut c_void,
                        &mut nn as *mut c_int as *mut c_void,
                        &mut pos as *mut u64 as *mut c_void,
                        &mut fst as *mut u64 as *mut c_void,
                        &mut thr as *mut u64 as *mut c_void,
                        &mut po as *mut u64 as *mut c_void,
                    ];
                    let block: u64 = 128;
                    let grid = threads.div_ceil(block) as c_uint;
                    check((self.api.launch)(self.func, grid, 1, 1, block as c_uint, 1, 1, 0, std::ptr::null_mut(),
                                            args.as_mut_ptr(), std::ptr::null_mut()), "cuLaunchKernel")?;
                    check((self.api.sync)(), "cuCtxSynchronize")?;
                    check((self.api.dtoh)(out.as_mut_ptr() as *mut c_void, dout, obytes), "cuMemcpyDtoH")
                };
                let r = run();
                (self.api.mem_free)(du);
                (self.api.mem_free)(dout);
                r.map_err(fail)?;
            }
            Ok(out)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matrix(n: usize, seed: u64) -> Vec<C64> {
        let mut s = crate::rng::Stream::new(seed, 0);
        (0..n * n).map(|_| C64::new(s.uniform() - 0.5, s.uniform() - 0.5)).collect()
    }

    /// on a machine with a GPU: equal to the CPU kernel within a few ulps of the scale, bit-identical on repetition;
    /// without one: Unavailable, never a wrong value
    #[test]
    fn against_the_cpu_kernel() {
        match device() {
            Err(Error::Unavailable(why)) => eprintln!("no GPU ({why}): skipped"),
            Err(e) => panic!("{e:?}"),
            Ok(name) => {
                for n in [1usize, 2, 5, 13, 14, 18] {
                    let a = matrix(n, n as u64);
                    let cpu = crate::perm::permanent(&a, n, 1, &|| false).unwrap();
                    let gpu = permanent(&a, n).unwrap();
                    let err = (gpu.value.re - cpu.value.re).abs().max((gpu.value.im - cpu.value.im).abs());
                    assert!(err <= 8.0 * f64::EPSILON * cpu.scale, "{name} n = {n}: {err} vs scale {}", cpu.scale);
                    assert!((gpu.scale - cpu.scale).abs() <= 8.0 * f64::EPSILON * cpu.scale);
                    let again = permanent(&a, n).unwrap();
                    assert_eq!(gpu.value, again.value);
                }
            }
        }
    }

    fn symmetric(m: usize, seed: u64) -> Vec<C64> {
        let a = matrix(m, seed);
        let mut s = vec![C64::ZERO; m * m];
        for i in 0..m {
            for j in 0..m {
                s[i * m + j] = a[i * m + j].add(a[j * m + i]).scale(0.5);
            }
        }
        s
    }

    /// the hafnian and loop hafnian on the GPU against the CPU kernel: within a few ulps of the scale, odd sizes
    /// included, and bit-identical on repetition
    #[test]
    fn hafnian_against_the_cpu_kernel() {
        let Ok(_) = device() else { return };
        for m in [2usize, 4, 7, 10, 16, 22] {
            let a = symmetric(m, 50 + m as u64);
            let l = matrix(1, 70 + m as u64).into_iter().cycle().take(m).enumerate().map(|(i, z)| z.scale(1.0 + i as f64 * 0.1)).collect::<Vec<_>>();
            for loops in [None, Some(l.as_slice())] {
                let cpu = crate::hafnian::hafnian(&a, m, loops, 1, &|| false).unwrap();
                let gpu = hafnian(&a, m, loops, &|| false).unwrap();
                let tol = 8.0 * f64::EPSILON * cpu.scale.max(1e-300);
                assert!((gpu.value.re - cpu.value.re).abs() <= tol && (gpu.value.im - cpu.value.im).abs() <= tol,
                        "m = {m} loops {}: {:?} vs {:?} (scale {})", loops.is_some(), gpu.value, cpu.value, cpu.scale);
                assert_eq!(gpu.subsets, cpu.subsets);
                assert_eq!(gpu.value, hafnian(&a, m, loops, &|| false).unwrap().value);
            }
        }
    }

    /// spike G1: time against the CPU kernel (1 and `threads` threads) at n = 20..28, and the difference in units of
    /// eps * S. cargo test --release --lib gpu::tests::spike_g1 -- --ignored --nocapture
    #[test]
    #[ignore]
    fn spike_g1() {
        let Ok(name) = device() else { return };
        let threads = std::thread::available_parallelism().map(|x| x.get()).unwrap_or(1);
        eprintln!("device {name}; CPU threads {threads}");
        for n in [16usize, 20, 22, 24, 26, 28] {
            let a = matrix(n, 100 + n as u64);
            let _ = permanent(&a, n).unwrap(); // the first launch JIT-compiles
            let t = std::time::Instant::now();
            let gpu = permanent(&a, n).unwrap();
            let tg = t.elapsed().as_secs_f64();
            let t = std::time::Instant::now();
            let cpu_n = crate::perm::permanent(&a, n, threads, &|| false).unwrap();
            let tn = t.elapsed().as_secs_f64();
            let t1 = if n <= 24 {
                let t = std::time::Instant::now();
                let _ = crate::perm::permanent(&a, n, 1, &|| false).unwrap();
                t.elapsed().as_secs_f64()
            } else {
                f64::NAN
            };
            let err = (gpu.value.re - cpu_n.value.re).abs().max((gpu.value.im - cpu_n.value.im).abs());
            eprintln!("n = {n}: GPU {tg:.4} s, CPU x{threads} {tn:.4} s, CPU x1 {t1:.4} s; |GPU - CPU| = {:.2} eps S",
                      err / (f64::EPSILON * cpu_n.scale));
        }
    }

    /// spike G1 for the hafnian: cargo test --release --lib gpu::tests::spike_g1_hafnian -- --ignored --nocapture
    #[test]
    #[ignore]
    fn spike_g1_hafnian() {
        let Ok(name) = device() else { return };
        let threads = std::thread::available_parallelism().map(|x| x.get()).unwrap_or(1);
        eprintln!("device {name}; CPU threads {threads}");
        for m in [16usize, 24, 28, 32, 36] {
            let a = symmetric(m, 200 + m as u64);
            let _ = hafnian(&a, m, None, &|| false).unwrap();
            let t = std::time::Instant::now();
            let gpu = hafnian(&a, m, None, &|| false).unwrap();
            let tg = t.elapsed().as_secs_f64();
            let t = std::time::Instant::now();
            let cpu = crate::hafnian::hafnian(&a, m, None, threads, &|| false).unwrap();
            let tc = t.elapsed().as_secs_f64();
            let err = (gpu.value.re - cpu.value.re).abs().max((gpu.value.im - cpu.value.im).abs());
            eprintln!("m = {m}: GPU {tg:.3} s, CPU x{threads} {tc:.3} s; |GPU - CPU| = {:.2} eps S", err / (f64::EPSILON * cpu.scale));
        }
    }
}
