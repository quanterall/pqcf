#!/bin/sh
# The quantum optics head-to-head (comparison/quantum-optics/README.md): PQCF against QuTiP on the tasks of tasks.json,
# on the same machine and the same core, one program at a time. Benchmarks the committed HEAD. Writes
# comparison/quantum-optics/results/<host>-<date>.json and RESULTS.md, appends to history.csv; exits 1 if any answer
# is wrong.
#   PQ_BENCH_HOST=<label> sh comparison/quantum-optics/run.sh [core]
set -e
cd "$(dirname "$0")/../.."
CORE=${1:-8}
HOST=${PQ_BENCH_HOST:-$(hostname -s)}
DATE=$(date -u +%Y-%m-%dT%H%MZ)
COMMIT=$(git rev-parse --short HEAD)
CPU=$(lscpu | sed -n 's/^Model name: *//p' | sort -u | paste -sd+)
docker build -q -t pqcf-dev . >/dev/null
docker build -q -t pqcf-qutip comparison/quantum-optics/qutip >/dev/null
WORK=$(mktemp -d "${TMPDIR:-/tmp}/pqcf-qo-XXXXXX")
mkdir -p "$WORK/src" "$WORK/out"
git archive HEAD | tar -x -C "$WORK/src"
ONE="--cpuset-cpus=$CORE -e OMP_NUM_THREADS=1 -e MKL_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1"
docker run --rm $ONE -v "$WORK/src":/src -v "$WORK/out":/out -v pqcf-bench-build:/src/_build \
    -v pqcf-bench-target:/src/native/pq_nif/target -w /src pqcf-dev sh -c \
    "rebar3 compile >/dev/null && mkdir -p /tmp/b && erlc -o /tmp/b comparison/quantum-optics/qo_pqcf.erl &&
     erl -noshell -pa _build/default/lib/*/ebin -pa /tmp/b -eval
       'application:load(pq), qo_pqcf:main(\"comparison/quantum-optics/tasks.json\", \"/out\"), halt().'"
echo "PQCF done"
docker run --rm $ONE -v "$WORK/src/comparison/quantum-optics":/p -v "$WORK/out":/out pqcf-qutip \
    python /p/qutip_side.py /p/tasks.json /out || echo "QuTiP failed: its results are missing from the tables"
printf '{"host": "%s", "date": "%s", "commit": "%s", "cpu": "%s", "core": %s}\n' "$HOST" "$DATE" "$COMMIT" "$CPU" "$CORE" \
    > "$WORK/out/meta.json"
mkdir -p comparison/quantum-optics/results
CHECK=0
docker run --rm -v "$WORK/out":/out -v "$PWD/comparison/quantum-optics":/c --user "$(id -u):$(id -g)" pqcf-qutip \
    python /c/compare.py /out "/c/results/$HOST-$DATE.json" /c/RESULTS.md /c/history.csv /out/meta.json \
    /c/tasks.json || CHECK=1
docker run --rm -v "$WORK":/w pqcf-qutip chown -R "$(id -u):$(id -g)" /w
rm -rf "${WORK:?}"
exit $CHECK
