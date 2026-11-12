#!/bin/sh
# PQCF against FORM + Forcer on FORM's own Forcer checks (comparison/forcer/README.md). Benchmarks the committed
# HEAD. Steps:
#   1. PQCF reduces every case (forcer_suite:pqcf/1)
#   2. Forcer evaluates PQCF's masters and each target, and FORM runs its own checks, timed
#   3. the answers are compared: PQCF composed with Forcer's master values against Forcer's results, and the
#      MS-bar sums against FORM's assertions (forcer_suite:check/1)
# Writes comparison/forcer/results/<host>-<date>.json; exits 1 if anything differs.
#   docker build -t pqcf-form -f comparison/form/Dockerfile . && sh comparison/forcer/run.sh [cores]
set -e
cd "$(dirname "$0")/../.."
CORES=${1:-0-15}
HOST=${PQ_BENCH_HOST:-$(hostname -s)}
DATE=$(date -u +%Y-%m-%dT%H%MZ)
COMMIT=$(git rev-parse --short HEAD)
docker build -q -t pqcf-forcer comparison/forcer >/dev/null
WORK=$(mktemp -d "${TMPDIR:-/tmp}/pqcf-forcer-XXXXXX")
mkdir -p "$WORK/src" "$WORK/out"
git archive HEAD | tar -x -C "$WORK/src"
dev() {
    docker run --rm --cpuset-cpus="$CORES" -v "$WORK/src":/src -v "$WORK/out":/out -v pqcf-bench-build:/src/_build \
        -v pqcf-bench-target:/src/native/pq_nif/target -w /src pqcf-dev sh -c \
        "rebar3 compile >/dev/null && mkdir -p /tmp/b && erlc -o /tmp/b comparison/forcer/*.erl && erl -noshell -pa _build/default/lib/*/ebin -pa /tmp/b -eval '$1, halt().'"
}
dev 'forcer_suite:pqcf("/out")'
# Forcer: PQCF's masters and the targets (exact, in Forcer's basis), then FORM's checks timed (one core)
CORE1=$(echo "$CORES" | cut -d- -f1 | cut -d, -f1)
docker run --rm --cpuset-cpus="$CORE1" -v "$WORK/out":/out -w /out pqcf-forcer sh -c '
    cp /opt/forcer/forcer/Masters.prc /opt/formchecks/forcer.frm /out/
    for f in *_masters.frm *_target.frm; do form -q "$f" > "${f%.frm}.out" 2>&1; done
    cd /opt/formchecks && echo "[" > /out/forcer.json && SEP=""
    for T in Forcer_example Forcer_1 Forcer_2; do
        S0=$(python3 -c "import time; print(time.time())")
        form -q -D TEST=$T forcer.frm > /out/$T.log 2>&1
        SECS=$(python3 -c "import time; print(round(time.time() - $S0, 3))")
        printf "%s{\"check\": \"%s\", \"form_forcer_seconds\": %s}\n" "$SEP" "$T" "$SECS" >> /out/forcer.json; SEP=","
        echo "FORM + Forcer, $T: $SECS s"
    done
    echo "]" >> /out/forcer.json'
CHECK=0
dev 'forcer_suite:check("/out")' || CHECK=1
docker run --rm -v "$WORK":/w pqcf-dev chown -R "$(id -u):$(id -g)" /w
mkdir -p comparison/forcer/results
# PQCF's symbolic results (reductions, expansions, MS-bar results) of the latest run, readable
rm -rf comparison/forcer/results/symbolic && cp -r "$WORK/out/symbolic" comparison/forcer/results/symbolic
OUT=comparison/forcer/results/$HOST-$DATE.json
python3 - "$WORK/out" "$OUT" "$HOST" "$DATE" "$COMMIT" "$CORES" "$(cut -d' ' -f1-3 /proc/loadavg)" <<'PY'
import json, sys
work, out, host, date, commit, cores, load = sys.argv[1:]
rec = {"host": host, "date": date, "commit": commit, "cores": cores, "load_average": load,
       "pqcf": json.load(open(f"{work}/pqcf.json")), "form_forcer": json.load(open(f"{work}/forcer.json")),
       "check": json.load(open(f"{work}/check.json"))}
json.dump(rec, open(out, "w"), indent=1)
print("written", out)
PY
rm -rf "${WORK:?}"
exit $CHECK
