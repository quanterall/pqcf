#!/bin/sh
# G4 of the genome research (README.md): genomes as matrix product states, the bond dimension they need.
# Writes research/genome/results/g4-<host>-<date>.json and RESULTS_G4.md.
#   sh research/genome/data/fetch.sh && sh research/genome/g4_run.sh ["5-8"]
set -e
cd "$(dirname "$0")"
CORES=${1:-"5-8"}
HOST=${PQ_BENCH_HOST:-$(hostname -s)}
DATE=$(date -u +%Y-%m-%dT%H%MZ)
docker build -q -t pqcf-genome docker >/dev/null
WORK=$(mktemp -d "${TMPDIR:-/tmp}/pqcf-g4-XXXXXX")
G="--user=$(id -u):$(id -g)"
CHECK=0
docker run --rm $G --cpuset-cpus="$CORES" -e OMP_NUM_THREADS=4 -e OPENBLAS_NUM_THREADS=4 -v "$PWD":/g -v "$WORK":/out \
    -w /g pqcf-genome python g4.py /out || CHECK=1
mkdir -p results
cp "$WORK/g4.json" "results/g4-$HOST-$DATE.json"
docker run --rm $G -v "$PWD":/g -w /g pqcf-genome python report_g4.py "results/g4-$HOST-$DATE.json" > RESULTS_G4.md
rm -rf "${WORK:?}"
echo "written results/g4-$HOST-$DATE.json"
exit $CHECK
