#!/bin/sh
# G0 of the genome research (docs/research/genome/README.md): the classical baselines on the genomic instances,
# one core, one thread. Writes research/genome/results/g0-<host>-<date>.json.
#   sh research/genome/run.sh [core]
set -e
cd "$(dirname "$0")"
CORE=${1:-9}
HOST=${PQ_BENCH_HOST:-$(hostname -s)}
DATE=$(date -u +%Y-%m-%dT%H%MZ)
docker build -q -t pqcf-genome docker >/dev/null
mkdir -p results
OUT=$(mktemp -d "${TMPDIR:-/tmp}/pqcf-genome-XXXXXX")
CHECK=0
docker run --rm --cpuset-cpus="$CORE" -e OMP_NUM_THREADS=1 -e NUMBA_NUM_THREADS=1 --user "$(id -u):$(id -g)" \
    -v "$PWD":/g -v "$OUT":/out -w /g pqcf-genome python g0.py /out || CHECK=1
cp "$OUT/g0.json" "results/g0-$HOST-$DATE.json"
rm -rf "${OUT:?}"
echo "written results/g0-$HOST-$DATE.json"
exit $CHECK
