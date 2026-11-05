// The hafnian and loop hafnian on the GPU (Phase 4 follow-up; docs/decisions/G1-gpu-offload.md).
//
// The same subset terms as hafnian.rs (the power-trace formula of Bjorklund, Gupt and Quesada), operation for
// operation: E = (A X)_Z, Householder reduction to Hessenberg form, La Budde's characteristic polynomial, Newton's
// identities, the exponential's power series (built with -fmad=false and IEEE division and square root, so each
// term equals the CPU kernel's bit for bit). Thread t adds the terms of subsets [G t, G (t+1)) in Gray-code order
// with Neumaier summation; the host adds the threads' sums in thread order.
//
// a: m x m row-major (re, im); loops: m weights or null; work: per thread `stride` doubles of scratch;
// out: 3 doubles per thread of the launch (sum re, sum im, sum of |term|).

struct c64 { double re, im; };

__device__ __forceinline__ c64 mk(double r, double i) { c64 z; z.re = r; z.im = i; return z; }
__device__ __forceinline__ c64 cadd(c64 a, c64 b) { return mk(a.re + b.re, a.im + b.im); }
__device__ __forceinline__ c64 csub(c64 a, c64 b) { return mk(a.re - b.re, a.im - b.im); }
__device__ __forceinline__ c64 cmul(c64 a, c64 b) { return mk(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re); }
__device__ __forceinline__ c64 cscale(c64 a, double s) { return mk(a.re * s, a.im * s); }
__device__ __forceinline__ c64 cconj(c64 a) { return mk(a.re, -a.im); }
__device__ __forceinline__ c64 cdivr(c64 a, double x) { return mk(a.re / x, a.im / x); }
__device__ __forceinline__ double cabs(c64 a) { return sqrt(a.re * a.re + a.im * a.im); }

struct neumaier { double s, c; };
__device__ __forceinline__ void nadd(neumaier &n, double x) {
    double t = n.s + x;
    if (fabs(n.s) >= fabs(x)) n.c += (n.s - t) + x; else n.c += (x - t) + n.s;
    n.s = t;
}

#define MAXP 32

// the workspace layout for D = 2n (complex entries): E D*D | P (D+1)(D+2)/2 | v D | s n+1 | g n | h n+1 | w D |
// nw D | H D*D (the Hessenberg copy of E; E itself is still needed by the loop terms)
__device__ c64 subset_term(const c64 *a, int m, const c64 *loops, int n, unsigned long long gray, c64 *ws) {
    int idx[2 * MAXP];
    int d = 0;
    for (int i = 0; i < n; i++)
        if ((gray >> i) & 1ULL) { idx[d++] = 2 * i; idx[d++] = 2 * i + 1; }
    int D = 2 * n;
    c64 *e = ws;
    c64 *P = e + D * D;
    c64 *v = P + (D + 1) * (D + 2) / 2;
    c64 *s = v + D;
    c64 *g = s + (n + 1);
    c64 *h = g + n;
    c64 *w = h + (n + 1);
    c64 *nw = w + D;
    for (int r = 0; r < d; r++)
        for (int c = 0; c < d; c++) e[r * d + c] = a[idx[r] * m + (idx[c] ^ 1)];
    // power traces tr(E^k), k = 1..n (all zero for the empty subset)
    if (d == 0) {
        for (int k = 1; k <= n; k++) s[k] = mk(0.0, 0.0);
    } else {
        c64 *hm = nw + D;
        for (int x = 0; x < d * d; x++) hm[x] = e[x];
        for (int k = 0; k + 2 < d; k++) {
            int len = d - k - 1;
            double norm = 0.0;
            for (int t = 0; t < len; t++) { c64 z = hm[(k + 1 + t) * d + k]; norm = norm + (z.re * z.re + z.im * z.im); }
            norm = sqrt(norm);
            if (norm == 0.0) continue;
            c64 x0 = hm[(k + 1) * d + k];
            double a0 = cabs(x0);
            c64 phase = (a0 == 0.0) ? mk(1.0, 0.0) : cdivr(x0, a0);
            for (int t = 0; t < len; t++) v[t] = hm[(k + 1 + t) * d + k];
            v[0] = cadd(x0, cscale(phase, norm));
            double vn = 0.0;
            for (int t = 0; t < len; t++) vn = vn + (v[t].re * v[t].re + v[t].im * v[t].im);
            vn = sqrt(vn);
            if (vn == 0.0) continue;
            for (int t = 0; t < len; t++) v[t] = cdivr(v[t], vn);
            for (int j = 0; j < d; j++) {
                c64 dot = mk(0.0, 0.0);
                for (int t = 0; t < len; t++) dot = cadd(dot, cmul(cconj(v[t]), hm[(k + 1 + t) * d + j]));
                for (int t = 0; t < len; t++) {
                    int id = (k + 1 + t) * d + j;
                    hm[id] = csub(hm[id], cscale(cmul(v[t], dot), 2.0));
                }
            }
            for (int i = 0; i < d; i++) {
                c64 dot = mk(0.0, 0.0);
                for (int t = 0; t < len; t++) dot = cadd(dot, cmul(hm[i * d + k + 1 + t], v[t]));
                for (int t = 0; t < len; t++) {
                    int id = i * d + k + 1 + t;
                    hm[id] = csub(hm[id], cscale(cmul(dot, cconj(v[t])), 2.0));
                }
            }
        }
        // La Budde: p_i in P at offset i(i+1)/2, coefficients of x^0 .. x^i
        P[0] = mk(1.0, 0.0);
        for (int i = 1; i <= d; i++) {
            c64 *nx = P + i * (i + 1) / 2;
            c64 *pp = P + (i - 1) * i / 2;
            for (int t = 0; t <= i; t++) nx[t] = mk(0.0, 0.0);
            c64 hii = hm[(i - 1) * d + (i - 1)];
            for (int t = 0; t < i; t++) {
                c64 q = pp[t];
                nx[t + 1] = cadd(nx[t + 1], q);
                nx[t] = csub(nx[t], cmul(hii, q));
            }
            c64 prod = mk(1.0, 0.0);
            for (int mm = 1; mm < i; mm++) {
                prod = cmul(prod, hm[(i - mm) * d + (i - mm - 1)]);
                c64 coef = cmul(hm[(i - mm - 1) * d + (i - 1)], prod);
                c64 *pq = P + (i - mm - 1) * (i - mm) / 2;
                for (int t = 0; t <= i - mm - 1; t++) nx[t] = csub(nx[t], cmul(coef, pq[t]));
            }
        }
        c64 *pd = P + d * (d + 1) / 2;   // c[k] = pd[d - k]
        s[0] = mk(0.0, 0.0);
        for (int k = 1; k <= n; k++) {
            c64 acc = mk(0.0, 0.0);
            if (k <= d) acc = cscale(pd[d - k], -(double)k);
            int lim = k < d + 1 ? k : d + 1;
            for (int i = 1; i < lim; i++) acc = csub(acc, cmul(pd[d - i], s[k - i]));
            s[k] = acc;
        }
    }
    for (int k = 1; k <= n; k++) g[k - 1] = cdivr(s[k], 2.0 * (double)k);
    if (loops) {
        for (int r = 0; r < d; r++) { v[r] = loops[idx[r]]; w[r] = v[r]; }
        for (int k = 0; k < n; k++) {
            c64 acc = mk(0.0, 0.0);
            for (int r = 0; r < d; r++) acc = cadd(acc, cmul(v[r], w[r ^ 1]));
            g[k] = cadd(g[k], cscale(acc, 0.5));
            for (int r = 0; r < d; r++) {
                c64 sm = mk(0.0, 0.0);
                for (int c = 0; c < d; c++) sm = cadd(sm, cmul(e[r * d + c], w[c]));
                nw[r] = sm;
            }
            for (int r = 0; r < d; r++) w[r] = nw[r];
        }
    }
    h[0] = mk(1.0, 0.0);
    for (int mm = 1; mm <= n; mm++) {
        c64 acc = mk(0.0, 0.0);
        for (int k = 1; k <= mm; k++) acc = cadd(acc, cmul(cscale(g[k - 1], (double)k), h[mm - k]));
        h[mm] = cdivr(acc, (double)mm);
    }
    return h[n];
}

extern "C" __global__ void hafnian_terms(const double *a, int m, const double *loops, int n, unsigned long long group,
                                         unsigned long long first, unsigned long long threads, double *work,
                                         unsigned long long stride, double *out) {
    unsigned long long local = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (local >= threads) return;
    unsigned long long t = first + local;
    unsigned long long total = 1ULL << n;
    unsigned long long start = t * group;
    unsigned long long end = start + group < total ? start + group : total;
    c64 *ws = (c64 *)(work + local * stride);
    neumaier sr = {0.0, 0.0}, si = {0.0, 0.0}, sa = {0.0, 0.0};
    for (unsigned long long p = start; p < end; p++) {
        unsigned long long gray = p ^ (p >> 1);
        int size = __popcll(gray);
        c64 term = subset_term((const c64 *)a, m, (const c64 *)loops, n, gray, ws);
        if ((n - size) % 2 == 1) term = cscale(term, -1.0);
        nadd(sr, term.re); nadd(si, term.im); nadd(sa, cabs(term));
    }
    out[3 * local] = sr.s + sr.c;
    out[3 * local + 1] = si.s + si.c;
    out[3 * local + 2] = sa.s + sa.c;
}
