# S3 and S4: the hafnian family and the torontonian

**Decisions** (closed 30 October 2026): both algorithms of AN 4 as specified; the hafnian formula is written out
below (AN 4.1 left it to this spike), and no prefix caching is needed for the threshold sampler.

## S3: hafnian and loop hafnian (`native/pq_nif/src/hafnian.rs`)

For a symmetric M x M matrix, M = 2n, with the fixed pairing (2i, 2i+1) and X swapping the members of each pair
(Bjorklund, Gupt and Quesada 2019):

    haf(A)  = sum over subsets Z of the n pairs  (-1)^(n-|Z|) [eta^n] exp( sum_(k=1..n) g_k eta^k )
    g_k     = tr(E^k) / (2k),           E = (A X) restricted to the indices of Z
    lhaf(A, v):  g_k += (1/2) v_Z^T X E^(k-1) v_Z

`tr(E^k)` from the characteristic polynomial (Householder reduction to Hessenberg form, La Budde's recurrence,
Newton's identities), so each subset costs O(|Z|^3); the coefficient of `eta^n` from the power series of the
exponential. Odd M: the hafnian is 0; the loop hafnian pads a vertex with loop weight 1. Repeated indices expand
the matrix. Subsets in Gray-code order of their index, blocks of 2^12, block sums in block order.

Exit test: equal to the matching recursion for M <= 12 within 1e-12 of S (unit test, with and without loops);
equal to thewalrus 0.22.0 within 1e-10 of S for M <= 40, with repeated indices (`native/pq_nif/tests/gt_h01.rs`,
fixtures `test/fixtures/gt_h01.json`). **Result: worst error 1.2e-15 of S** (`spikes/s3_s4/s3_hafnian.txt`).

| M | hafnian, ours | thewalrus | error / S | S / value |
| --- | --- | --- | --- | --- |
| 20 | 0.003 s | 0.007 s | 7e-17 | 1.5e4 |
| 30 | 0.24 s | 0.25 s | 2e-16 | 1.2e7 |
| 36 | 3.0 s | 2.3 s | 3e-17 | 4.4e7 |
| 40 | 15.4 s | 11.6 s | 8e-17 | 7.2e8 |

One core of a shared aarch64 machine for ours; thewalrus as installed in the oracle container. The last column is
why tolerances are relative to S (amendment A-1): at M = 40 the terms cancel by nine orders of magnitude, so an
error of 1e-16 of S is about 1e-7 of the value.

## S4: torontonian (`native/pq_nif/src/toronto.rs`)

AN 4.2 as written: `tor(O) = sum_Z (-1)^(c-|Z|) / sqrt(det(I - O_Z))`, each `I - O_Z` a principal submatrix of
`Q^-1`, its Cholesky factor grown by two bordered rows per included mode in a depth-first search, so a term is
`1 / prod(L_ii)`. The first min(6, c) decisions give at most 64 subtrees, summed in order. A matrix that is not of
the form `I - Q^-1` (not positive definite) returns `{non_physical_state, not_positive_definite}`.

Exit test: equal to thewalrus within 1e-10 of S for c <= 16; c = 20 in under 2 s; the threshold chain rule under
60 s per sample at m = 20 or else prefix caching. **Result**: worst error 4e-14 of S up to c = 20 (random mixed
states from thewalrus, `tests/gt_tor.rs`, `test/fixtures/gt_tor.json`); P(B) = tor / sqrt(det Q) equals
thewalrus's `threshold_detection_prob` within 4e-14; FS A.10's check P(click) = 1 - 1/cosh r holds to 1e-14;
the definition (determinants by LU) is matched for c <= 10.

| c | ours | thewalrus |
| --- | --- | --- |
| 16 | 0.019 s | 0.081 s |
| 18 | 0.088 s | 0.35 s |
| 20 | 0.40 s | 1.5 s |

The threshold chain rule at m = 20 evaluates one torontonian per mode on the clicks so far; with every mode
clicked that sums to about twice the c = 20 cost, about 1 s per sample, far below 60 s: **no prefix caching**.
