#!/bin/sh
# G3 of the genome research (README.md): the GBS graph kernel against classical features, for telling which
# organism a window of DNA comes from. PQCF checks the GBS features on two window graphs first.
# Writes research/genome/results/g3-<host>-<date>.json and RESULTS_G3.md.
#   sh research/genome/data/fetch.sh && sh research/genome/g3_run.sh [core]
set -e
cd "$(dirname "$0")"
CORE=${1:-15}
HOST=${PQ_BENCH_HOST:-$(hostname -s)}
DATE=$(date -u +%Y-%m-%dT%H%MZ)
docker build -q -t pqcf-genome docker >/dev/null
WORK=$(mktemp -d "${TMPDIR:-/tmp}/pqcf-g3-XXXXXX")
G="--user=$(id -u):$(id -g)"
docker run --rm $G -v "$PWD":/g -v "$WORK":/out -w /g pqcf-genome python g3.py specs /out
docker run --rm -v "$PWD/../..":/src -v pqcf-genome-build:/src/_build -v pqcf-genome-target:/src/native/pq_nif/target \
    -w /src pqcf-dev rebar3 compile >/dev/null
docker run --rm --cpuset-cpus="$CORE" -v "$PWD/../..":/src -v "$WORK/chunk-0":/out -v pqcf-genome-build:/src/_build \
    -v pqcf-genome-target:/src/native/pq_nif/target -w /src pqcf-dev sh -c \
    "mkdir -p /tmp/b && erlc -I include -o /tmp/b research/genome/g1_pqcf.erl &&
     erl -noshell -pa _build/default/lib/*/ebin -pa /tmp/b -eval 'g1_pqcf:main(\"/out\"), halt().'"
docker run --rm -v "$WORK":/w pqcf-genome chown -R "$(id -u):$(id -g)" /w
CHECK=0
docker run --rm $G --cpuset-cpus="$CORE" -e OMP_NUM_THREADS=1 -v "$PWD":/g -v "$WORK":/out -w /g pqcf-genome \
    python g3.py analyse /out || CHECK=1
mkdir -p results
cp "$WORK/g3.json" "results/g3-$HOST-$DATE.json"
docker run --rm $G -v "$PWD":/g -w /g pqcf-genome python report_g3.py "results/g3-$HOST-$DATE.json" > RESULTS_G3.md
rm -rf "${WORK:?}"
echo "written results/g3-$HOST-$DATE.json"
exit $CHECK
