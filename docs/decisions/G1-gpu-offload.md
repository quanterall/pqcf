# G1: GPU offload of the permanent

**Decision: Glynn's permanent gets a CUDA kernel behind the existing NIF contract, selected by `device => cpu | gpu |
auto` (default `cpu`).**

- The kernel is CUDA C (`native/pq_nif/gpu/perm.cu`), compiled once to PTX with `-fmad=false`. The PTX is committed
  and embedded in the crate.
- The CUDA driver library (`libcuda.so.1`) is opened at run time, so building the crate needs nothing from CUDA, and
  a machine without a GPU behaves as before.
- Each GPU thread computes exactly the positions between two refresh points of the CPU kernel, so every term is
  bit-identical to the CPU's.

Phase 4, item 5 (FS 15 item 5). Released after 1.0.0-rc.1, purely additive under the freeze: a new option and a new
`device` entry in the kernels' info map.

## Why PTX and the driver API

- **No CUDA build dependency:** the crate builds on CI and the reference machines without a CUDA toolkit. The driver
  JIT-compiles compute_80 PTX for any newer device; on the GB10 that is the Blackwell `sm_121` target, with driver
  580 and CUDA 13.0. `gpu/build_ptx.sh` regenerates the PTX when the source changes.
- **No CUDA runtime library and no binding crate:** the driver API needs about ten functions, opened with
  `libloading`, which is ISC-licensed and already in the dependency tree through Rustler.
- **Unavailable is a value:** without the library, a device, or a loadable module, `device => gpu` returns
  `{error, {unsupported_input, #{device => gpu, reason => ...}}}`, never a silent CPU run. `auto` falls back to the
  CPU.

## Determinism

The CPU kernel (perm.rs) cuts the 2^(n−1) Gray-code positions into blocks of 2^16 and recomputes the column sums w
directly every 4096 positions. GPU thread t takes positions [4096 t, 4096 (t + 1)):

- it starts with w computed directly, as the CPU does at that position;
- it moves the digits with the same rule;
- it updates w with the same operations, so every term is the CPU kernel's term, bit for bit (`-fmad=false`: no
  fused multiply-add, AN 2.2).

The thread adds its terms with Neumaier summation, and the host adds the threads' sums in thread order. Only that
grouping differs from the CPU's blocks of 2^16. It does not depend on the launch configuration, so GPU results are
reproducible bit for bit (D1 on the GPU path).

Against the CPU they differ only by the regrouped rounding. The tests allow 8 ε S; the measured difference is
0.00 ε S (rounded) at every n from 16 to 28.

Long enumerations run as launches of 2^16 threads (2^28 positions), with the cancel token checked between launches.

## Scope

- Plain permanents, `permanent_repeated` with every multiplicity 0 or 1, and `bs_probability` for collision-free
  patterns run on the GPU. Those cases reduce to a submatrix.
- Repeated rows or columns run on the CPU, and the info map says which device ran.
- `pq_api:probabilities` with the boson backend passes `device` through.
- The Clifford–Clifford sampler's minor permanents stay on the CPU: they are small and numerous, so launch overhead
  would dominate.

## Spike G1: speed on the GB10 (shared with other work; best of one run after a warm-up launch)

| n | GPU | CPU, 20 threads | CPU, 1 thread |
| --- | --- | --- | --- |
| 16 | 0.024 s | 0.0014 s | 0.0014 s |
| 20 | 0.034 s | 0.0071 s | 0.030 s |
| 22 | 0.037 s | 0.016 s | 0.13 s |
| 24 | 0.037 s | 0.056 s | 0.59 s |
| 26 | 0.079 s | 0.26 s | — |
| 28 | 0.30 s | 1.01 s | — |

- **Fixed cost:** about 25 ms per call (context, allocation and synchronisation).
- **Crossover:** the GPU passes all 20 CPU threads at n = 24, which is where `auto` switches; at n = 28 it is 3.3×
  faster.
- **Peak speed:** the GB10's FP64 rate limits it. A double-single or integer-assisted kernel would trade the
  bit-identity for speed. Not done: the bit-identical terms are what make the GPU path testable against the CPU
  kernel.

## The hafnian (follow-up, 5 Nov 2026)

`gpu/haf.cu`: each subset term of the power-trace formula is computed independently, so a thread takes a group of
16 subsets in Gray-code order. The thread computes, operation for operation like `hafnian.rs`:

- E = (A X)_Z;
- the Householder reduction to Hessenberg form;
- La Budde's characteristic polynomial;
- Newton's identities;
- the loop terms;
- the exponential's power series.

With `-fmad=false` and IEEE division and square root, every term is bit-identical to the CPU's. The thread adds
its terms with Neumaier summation, and the host adds the threads in order, so the result is reproducible bit for bit
and equal to the CPU within the regrouped rounding (0.00 ε S measured, m = 16–40). The matrices of a thread, up to
2n × 2n, live in a global-memory workspace; launches are capped at 16,384 threads and 256 MB of workspace.

- **Covered:** `hafnian/2`, `hafnian_repeated/3` and `loop_hafnian/4` with `device => gpu`, odd sizes included.
  Repeated indices are expanded on the GPU unless Kan's formula is cheaper, in which case they run on the CPU and
  the info map says so.
- **Speed: no gain on the GB10.** A smaller group (16 instead of 64) and bigger launches helped, but against 20 CPU
  threads, on a machine shared with other work:

  | m | GPU | CPU, 20 threads |
  | --- | --- | --- |
  | 28 | 0.16 s | 0.09 s |
  | 32 | 0.27 s | 0.24 s |
  | 36 | 1.07 s | 0.45–1.5 s |
  | 40 | 6.2 s | 5.2 s |

  The work is O(d³) complex FP64 per subset, with its matrices in global memory, which is the GB10's weak spot
  (low FP64 rate, memory shared with the CPU). So `auto` keeps the hafnian on the CPU. The GPU path is exact and
  tested, ready for a GPU with high FP64 throughput.
- **Tests:**
  - GT-H01 through the GPU path: 52 cases of the golden file up to m = 28, 42 of them on the GPU (the others on
    Kan's formula), all within 1e-10 of the scale;
  - Rust: hafnians and loop hafnians at m = 2–22 against the CPU kernel, bit-identical on repetition;
  - `auto` choosing the CPU, and the refusal without a GPU;
  - mutation check: one Householder reflection factor 2 → 1 in the CUDA source, caught by GT-H01.

## Open

- **Torontonian on the GPU:** its Cholesky walk carries state from subset to subset, a third kernel shape. Not
  done: on the GB10 it would meet the same FP64 limit as the hafnian.
- **`auto` for the hafnian:** stays on the CPU until a GPU with high FP64 throughput is measured, such as the
  dedicated machine's, if it has one.
- **Fock gates on the GPU:** these are memory-bound, and the GB10's memory is shared with the CPU, so a gain is
  unlikely there. To be measured on the dedicated machine.

## Tests (`test/phase4_gpu_tests.erl`, Rust `gpu::tests`)

- **With a GPU** (`docker run --rm --gpus all pqcf-dev`):
  - GT-P01, every permanent and probability of the golden file against thewalrus within 1e-10 of the scale;
  - n = 2–22 against the CPU kernel within 8 ε S, the scale too, and bit-identical on repetition;
  - `auto` choosing the CPU at n = 10 and the GPU at n = 24;
  - repeated multiplicities on the CPU, multiplicities of 0 and 1 reduced to the GPU;
  - boson probabilities through the facade equal on both devices;
  - a set cancel token returning `cancelled`;
  - mutation check: the Glynn sign tested on bit 1 instead of bit 0 in the CUDA source, caught by three tests.
- **Without one** (CI, `docker run --rm pqcf-dev`): `device => gpu` refused with the reason, and `auto` running on
  the CPU.
