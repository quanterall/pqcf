#!/bin/sh
# G2 of the genome research (README.md): coherent-Ising-machine solvers against simulated annealing and parallel
# tempering, each tuned with the same budget, by TTS* on one core each. Writes research/genome/results/g2-*.json.
#   sh research/genome/g2_run.sh ["15 16 17 18"]
set -e
cd "$(dirname "$0")"
CORES=${1:-"15 16 17 18"}
HOST=${PQ_BENCH_HOST:-$(hostname -s)}
DATE=$(date -u +%Y-%m-%dT%H%MZ)
docker build -q -t pqcf-genome docker >/dev/null
WORK=$(mktemp -d "${TMPDIR:-/tmp}/pqcf-g2-XXXXXX")
RUN="docker run --rm --user $(id -u):$(id -g) -e OMP_NUM_THREADS=1 -e NUMBA_NUM_THREADS=1 -v $PWD:/g -v $WORK:/out -w /g"
$RUN --cpuset-cpus=5-8,15-18 pqcf-genome python g2.py exact /out
set -- $CORES
for S in sa pt simcim ahc; do
    $RUN --cpuset-cpus="$1" pqcf-genome python g2.py run "$S" /out > "$WORK/$S.log" 2>&1 &
    shift
done
wait
cat "$WORK"/*.log
mkdir -p results
CHECK=0
$RUN pqcf-genome python g2.py merge /out "/g/results/g2-$HOST-$DATE.json" || CHECK=1
rm -rf "${WORK:?}"
echo "written results/g2-$HOST-$DATE.json"
exit $CHECK
