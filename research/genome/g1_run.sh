#!/bin/sh
# G1 of the genome research (README.md): Gaussian boson sampling simulated by PQCF, against Oh et al.'s
# quantum-inspired sampler and uniform sampling, for dense subgraphs and cliques. PQCF samples in parallel, one
# container per core, each with its own random streams. Writes research/genome/results/g1-<host>-<date>.json.
#   sh research/genome/g1_run.sh ["15 16 17 18 5 6 7 8"]
set -e
cd "$(dirname "$0")"
CORES=${1:-"15 16 17 18 5 6 7 8"}
HOST=${PQ_BENCH_HOST:-$(hostname -s)}
DATE=$(date -u +%Y-%m-%dT%H%MZ)
N=$(echo $CORES | wc -w)
docker build -q -t pqcf-genome docker >/dev/null
WORK=$(mktemp -d "${TMPDIR:-/tmp}/pqcf-g1-XXXXXX")
G=--user=$(id -u):$(id -g)
docker run --rm $G -v "$PWD":/g -v "$WORK":/out -w /g pqcf-genome python g1.py specs /out "$N"
# compile once: compiling in every container rewrites the NIF library under the others (SIGBUS)
docker run --rm -v "$PWD/../..":/src -v pqcf-genome-build:/src/_build -v pqcf-genome-target:/src/native/pq_nif/target \
    -w /src pqcf-dev rebar3 compile >/dev/null
I=0
for C in $CORES; do
    docker run --rm --cpuset-cpus="$C" -v "$PWD/../..":/src -v "$WORK/chunk-$I":/out \
        -v pqcf-genome-build:/src/_build -v pqcf-genome-target:/src/native/pq_nif/target -w /src pqcf-dev sh -c \
        "mkdir -p /tmp/b && erlc -I include -o /tmp/b research/genome/g1_pqcf.erl &&
         erl -noshell -pa _build/default/lib/*/ebin -pa /tmp/b -eval 'g1_pqcf:main(\"/out\"), halt().'" &
    I=$((I + 1))
done
wait
docker run --rm -v "$WORK":/w pqcf-genome chown -R "$(id -u):$(id -g)" /w
CHECK=0
docker run --rm $G --cpuset-cpus="$(echo $CORES | cut -d' ' -f1)" -v "$PWD":/g -v "$WORK":/out -w /g pqcf-genome \
    python g1.py analyse /out || CHECK=1
mkdir -p results
cp "$WORK/g1.json" "results/g1-$HOST-$DATE.json"
docker run --rm $G --cpuset-cpus="$(echo $CORES | cut -d' ' -f1)" -v "$PWD":/g -w /g pqcf-genome \
    python g1_equal.py "results/g1-$HOST-$DATE.json"
docker run --rm $G -v "$PWD":/g -w /g pqcf-genome python report_g1.py "results/g1-$HOST-$DATE.json" > RESULTS_G1.md
rm -rf "${WORK:?}"
echo "written results/g1-$HOST-$DATE.json"
exit $CHECK
