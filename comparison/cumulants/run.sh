#!/bin/sh
# QuantumCumulants.jl's results for PQCF's cumulant tests (test/pq_cumulant_tests.erl): the single-atom laser and the
# driven Jaynes-Cummings model at orders 2-4, written to test/fixtures/qc_cumulants.json. QuantumCumulants 0.3.7 and
# its dependencies are pinned in the image (julia/Dockerfile; the versions are printed at build time).
#   sh comparison/cumulants/run.sh            (about 2 minutes once the image is built; the build about 5)
set -e
cd "$(dirname "$0")"
docker build -q -t pqcf-cumulants julia >/dev/null
L=$(docker run --rm -v "$PWD":/c pqcf-cumulants julia /c/laser.jl 2>/dev/null | tail -1)
J=$(docker run --rm -v "$PWD":/c pqcf-cumulants julia /c/jc_orders.jl 2>/dev/null | tail -1)
V=$(docker run --rm pqcf-cumulants julia -e 'using QuantumCumulants; print(pkgversion(QuantumCumulants))')
printf '{"quantumcumulants": "%s", "laser": %s, "jc_orders": %s}\n' "$V" "$L" "$J" > ../../test/fixtures/qc_cumulants.json
cat ../../test/fixtures/qc_cumulants.json
