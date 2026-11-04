// Glynn's permanent on the GPU (Phase 4, FS 15 item 5; docs/decisions/G1-gpu-offload.md).
//
// The same enumeration as perm.rs with every multiplicity 1: rows 1..n-1 are the Gray digits (fastest first), row 0
// is fixed (k = 0, coefficient 2). Thread t computes positions [t L, (t+1) L), L = 4096 = REFRESH: it starts from w
// computed directly, exactly where the CPU kernel refreshes, and moves the digits like it, so every term is the CPU
// kernel's term bit for bit (built with -fmad=false: no fused multiply-add, AN 2.2). The thread's terms are added
// in order with Neumaier summation; the host adds the threads' sums in thread order. Only that grouping differs
// from the CPU's blocks of 2^16, and it is fixed, so results do not depend on the launch.
//
// u: n x n row-major, (re, im) interleaved; out: 3 doubles per thread of the launch (sum re, sum im, sum of |term|).
// A long enumeration is cut into launches of consecutive threads (`first`), so the host can check cancellation.

struct c64 { double re, im; };

__device__ __forceinline__ c64 cmul(c64 a, c64 b) {
    c64 r; r.re = a.re * b.re - a.im * b.im; r.im = a.re * b.im + a.im * b.re; return r;
}

struct neumaier { double s, c; };

__device__ __forceinline__ void nadd(neumaier &n, double x) {
    double t = n.s + x;
    if (fabs(n.s) >= fabs(x)) n.c += (n.s - t) + x; else n.c += (x - t) + n.s;
    n.s = t;
}

#define MAXN 40
#define L 4096ULL

extern "C" __global__ void glynn(const double *u, int n, unsigned long long positions, unsigned long long first,
                                 unsigned long long threads, double *out) {
    // threads first .. first + threads - 1 of the enumeration; out is indexed from first
    unsigned long long local = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (local >= threads) return;
    unsigned long long t = first + local;
    unsigned long long start = t * L;
    unsigned long long end = start + L < positions ? start + L : positions;
    int d = n - 1;                       // moving digits: rows 1..n-1
    int g[MAXN], dir[MAXN];
    // digits_at(start): digit k is bit k of start, reflected when the bits above are odd
    for (int k = 0; k < d; k++) {
        unsigned long long a = (start >> k) & 1ULL, high = start >> (k + 1);
        if ((high & 1ULL) == 0) { g[k] = (int)a; dir[k] = 1; } else { g[k] = 1 - (int)a; dir[k] = -1; }
    }
    c64 w[MAXN];
    // w_direct: rows in significance order (1..n-1, then 0 with k = 0)
    for (int j = 0; j < n; j++) {
        c64 acc; acc.re = 0.0; acc.im = 0.0;
        for (int k = 0; k < d; k++) {
            double f = 1.0 - 2.0 * (double)g[k];
            int i = k + 1;
            acc.re = acc.re + u[2 * (i * n + j)] * f;
            acc.im = acc.im + u[2 * (i * n + j) + 1] * f;
        }
        acc.re = acc.re + u[2 * j] * 1.0;
        acc.im = acc.im + u[2 * j + 1] * 1.0;
        w[j] = acc;
    }
    int parity = 0;
    for (int k = 0; k < d; k++) parity += g[k];
    neumaier sr = {0.0, 0.0}, si = {0.0, 0.0}, sa = {0.0, 0.0};
    for (unsigned long long p = start; p < end; p++) {
        if (p > start) {
            int k = 0;
            for (;;) {
                int next = g[k] + dir[k];
                if (next >= 0 && next < 2) break;
                dir[k] = -dir[k];
                k++;
            }
            g[k] = g[k] + dir[k];
            parity += dir[k];
            double f = -2.0 * (double)dir[k];
            int i = k + 1;
            for (int j = 0; j < n; j++) {
                w[j].re = w[j].re + u[2 * (i * n + j)] * f;
                w[j].im = w[j].im + u[2 * (i * n + j) + 1] * f;
            }
        }
        c64 prod; prod.re = 1.0; prod.im = 0.0;
        for (int j = 0; j < n; j++) prod = cmul(prod, w[j]);
        double coef = (parity & 1) ? -2.0 : 2.0;
        double tr = prod.re * coef, ti = prod.im * coef;
        nadd(sr, tr); nadd(si, ti); nadd(sa, sqrt(tr * tr + ti * ti));
    }
    out[3 * local] = sr.s + sr.c;
    out[3 * local + 1] = si.s + si.c;
    out[3 * local + 2] = sa.s + sa.c;
}
