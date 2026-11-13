#!/bin/sh
# The worked example on public genomes (README.md): assembly of SARS-CoV-2 and an E. coli region, phasing a
# simulated individual on human chromosome 21, dense modules of E. coli's DNA replication and repair network with
# PQCF's simulated GBS device, and the genomes' quantum encoding (from G4). Writes example.json and EXAMPLE.md.
#   sh research/genome/data/fetch.sh && sh research/genome/example/run.sh ["15 16 17 18 5 6 7 8"]
set -e
cd "$(dirname "$0")"
CORES=${1:-"15 16 17 18 5 6 7 8"}
N=$(echo $CORES | wc -w)
docker build -q -t pqcf-genome ../docker >/dev/null
WORK=$(mktemp -d "${TMPDIR:-/tmp}/pqcf-example-XXXXXX")
G="--user=$(id -u):$(id -g) -e NUMBA_CACHE_DIR=/tmp/numba"
docker run --rm $G --cpuset-cpus="$(echo $CORES | tr ' ' ',')" -v "$PWD/..":/g -v "$WORK":/out -w /g/example \
    pqcf-genome python example.py prepare /out "$N"
# PQCF: compile once, then one sampler per core
docker run --rm -v "$PWD/../../..":/src -v pqcf-genome-build:/src/_build -v pqcf-genome-target:/src/native/pq_nif/target \
    -w /src pqcf-dev rebar3 compile >/dev/null
I=0
for C in $CORES; do
    docker run --rm --cpuset-cpus="$C" -v "$PWD/../../..":/src -v "$WORK/chunk-$I":/out -v pqcf-genome-build:/src/_build \
        -v pqcf-genome-target:/src/native/pq_nif/target -w /src pqcf-dev sh -c \
        "mkdir -p /tmp/b && erlc -I include -o /tmp/b research/genome/g1_pqcf.erl &&
         erl -noshell -pa _build/default/lib/*/ebin -pa /tmp/b -eval 'g1_pqcf:main(\"/out\"), halt().'" &
    I=$((I + 1))
done
wait
docker run --rm -v "$WORK":/w pqcf-genome chown -R "$(id -u):$(id -g)" /w
CHECK=0
docker run --rm $G -v "$PWD/..":/g -v "$WORK":/out -w /g/example pqcf-genome python example.py analyse /out || CHECK=1
cp "$WORK/example.json" example.json
docker run --rm $G -v "$PWD/..":/g -w /g/example pqcf-genome python report_example.py example.json > EXAMPLE.md
rm -rf "${WORK:?}"
echo "written example.json and EXAMPLE.md"
exit $CHECK
