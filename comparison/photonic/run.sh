#!/bin/sh
# The photonic head-to-head (comparison/photonic/README.md): PQCF against The Walrus, Piquasso, Perceval, DeepQuantum
# and Strawberry Fields on the same inputs, the same machine and the same core, one program at a time. Benchmarks
# the committed HEAD. Steps:
#   1. the inputs from fixed seeds (gen_inputs.py)
#   2. PQCF (photonic_pqcf.erl), then every peer (peers.py), each pinned to core CORE with one thread
#   3. the answers checked against each other and the tables written (compare.py)
# Writes comparison/photonic/results/<host>-<date>.json and RESULTS.md, appends to history.csv; exits 1 if any
# answer disagrees.
#   PQ_BENCH_HOST=<label> sh comparison/photonic/run.sh [core]
set -e
cd "$(dirname "$0")/../.."
CORE=${1:-8}
HOST=${PQ_BENCH_HOST:-$(hostname -s)}
DATE=$(date -u +%Y-%m-%dT%H%MZ)
COMMIT=$(git rev-parse --short HEAD)
CPU=$(lscpu | sed -n 's/^Model name: *//p' | sort -u | paste -sd+)
docker build -q -t pqcf-dev . >/dev/null
docker build -q -t pqcf-photonic-peers comparison/photonic/peers >/dev/null
docker build -q -t pqcf-photonic-sf comparison/photonic/sf >/dev/null
WORK=$(mktemp -d "${TMPDIR:-/tmp}/pqcf-photonic-XXXXXX")
mkdir -p "$WORK/src" "$WORK/out"
git archive HEAD | tar -x -C "$WORK/src"
ONE="--cpuset-cpus=$CORE -e OMP_NUM_THREADS=1 -e NUMBA_NUM_THREADS=1 -e MKL_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1"
peer() {
    docker run --rm $ONE -v "$WORK/src/comparison/photonic":/p -v "$WORK/out":/out "$1" python /p/peers.py "$2" /out
}
docker run --rm -v "$WORK/src/comparison/photonic":/p -v "$WORK/out":/out pqcf-photonic-peers \
    sh -c 'python /p/gen_inputs.py > /out/inputs.json'
docker run --rm $ONE -v "$WORK/src":/src -v "$WORK/out":/out -v pqcf-bench-build:/src/_build \
    -v pqcf-bench-target:/src/native/pq_nif/target -w /src pqcf-dev sh -c \
    "rebar3 compile >/dev/null && mkdir -p /tmp/b && erlc -I include -o /tmp/b comparison/photonic/photonic_pqcf.erl &&
     erl -noshell -pa _build/default/lib/*/ebin -pa /tmp/b -eval 'application:load(pq), photonic_pqcf:main(\"/out\"), halt().'"
echo "PQCF done"
for T in thewalrus piquasso perceval deepquantum; do peer pqcf-photonic-peers $T; done
peer pqcf-photonic-sf strawberryfields
printf '{"host": "%s", "date": "%s", "commit": "%s", "cpu": "%s", "core": %s}\n' "$HOST" "$DATE" "$COMMIT" "$CPU" "$CORE" \
    > "$WORK/out/meta.json"
mkdir -p comparison/photonic/results
CHECK=0
docker run --rm -v "$WORK/out":/out -v "$PWD/comparison/photonic":/c --user "$(id -u):$(id -g)" pqcf-photonic-peers \
    python /c/compare.py /out "/c/results/$HOST-$DATE.json" /c/RESULTS.md /c/history.csv /out/meta.json || CHECK=1
docker run --rm -v "$WORK":/w pqcf-photonic-peers chown -R "$(id -u):$(id -g)" /w
rm -rf "${WORK:?}"
exit $CHECK
