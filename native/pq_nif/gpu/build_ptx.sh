#!/bin/sh
# Regenerates gpu/perm.ptx and gpu/haf.ptx from their .cu sources (needs nvcc; the PTX is committed, so building the crate does not).
# compute_80 PTX is JIT-compiled by the driver for any newer device; -fmad=false keeps the CPU kernel's
# operations (AN 2.2: no fused multiply-add), and the divisions and square roots IEEE-rounded.
set -e
cd "$(dirname "$0")"
for k in perm haf; do
    nvcc -ptx -arch=compute_80 -fmad=false -prec-div=true -prec-sqrt=true -ftz=false -o $k.ptx $k.cu
done
