#!/bin/sh
# The reduction benchmark suite (decision record R1, "A benchmark suite"): bench/ibp_suite.erl's cases, the ones FIRE
# ships as examples also run by FIRE 6.5.2, both on the same core back to back; FIRE's answers and PQCF's compared
# exactly. Benchmarks the committed HEAD (git archive), so that a result names its commit. Writes
# bench/results/ibp-<host>-<date>.json, prints it against the previous run on this host, and exits 1 if an answer
# differs from FIRE's. Needs docker and the images pqcf-dev (Dockerfile) and pqcf-fire2 (comparison/fire/Dockerfile).
#   bench/ibp.sh [core]            (default core 0; the machine's other load is recorded, not controlled)
# To run it regularly (here weekly, Sunday 03:00): crontab -e, then
#   0 3 * * 0  cd $HOME/pqcf && bench/ibp.sh 0 >> bench/results/ibp.log 2>&1
set -e
cd "$(dirname "$0")/.."
CORE=${1:-0}
HOST=${PQ_BENCH_HOST:-$(hostname -s)}
DATE=$(date -u +%Y-%m-%dT%H%MZ)
COMMIT=$(git rev-parse --short HEAD)
WORK=$(mktemp -d "${TMPDIR:-/tmp}/pqcf-ibp-XXXXXX")
mkdir -p "$WORK/src" "$WORK/out"
git archive HEAD | tar -x -C "$WORK/src"
LOAD=$(cut -d' ' -f1-3 /proc/loadavg)
DEV="docker run --rm -v $WORK/src:/src -v $WORK/out:/out -v pqcf-bench-build:/src/_build -v pqcf-bench-target:/src/native/pq_nif/target -w /src pqcf-dev"
ERL='erl -noshell -pa _build/default/lib/*/ebin -pa /tmp'
COMPILE='rebar3 compile >/dev/null && erlc -o /tmp bench/ibp_suite.erl comparison/fire/fire_check.erl'
# PQCF on one core, then on all of them; FIRE (one thread) on the same core; then the check
$DEV --cpuset-cpus="$CORE" sh -c "$COMPILE && $ERL +S 1 -eval 'ibp_suite:fire_cases(\"/out\"), ibp_suite:pqcf(\"/out\", \"one\"), halt().'"
$DEV sh -c "$COMPILE && $ERL -eval 'ibp_suite:pqcf(\"/out\", \"all\"), halt().'"
docker run --rm --cpuset-cpus="$CORE" -v "$WORK/out":/out -v "$PWD/bench":/b:ro pqcf-fire2 sh /b/fire_run.sh
CHECK=0
$DEV sh -c "$COMPILE && $ERL -eval 'ibp_suite:check(\"/out\"), halt().'" || CHECK=1
# the files the containers wrote are root's: hand them back before removing the work directory
docker run --rm -v "$WORK":/w pqcf-dev chown -R "$(id -u):$(id -g)" /w
OUT=bench/results/ibp-$HOST-$DATE.json
python3 - "$WORK/out" "$OUT" "$HOST" "$DATE" "$COMMIT" "$CORE" "$LOAD" "$(lscpu | sed -n 's/^Model name: *//p' | head -1)" "$(nproc)" <<'PY'
import glob, json, os, sys
work, out, host, date, commit, core, load, cpu, nproc = sys.argv[1:]
one = {r["case"]: r for r in json.load(open(f"{work}/pqcf-one.json"))}
allc = {r["case"]: r for r in json.load(open(f"{work}/pqcf-all.json"))}
fire = {r["case"]: r for r in json.load(open(f"{work}/fire.json"))}
check = {r["case"]: r for r in json.load(open(f"{work}/check.json"))} if os.path.exists(f"{work}/check.json") else {}
rows = []
for case, r in one.items():
    row = dict(r, seconds_one_core=r.pop("seconds"), seconds_all_cores=allc[case]["seconds"])
    row.pop("seconds", None)
    if case in fire:
        row["fire_seconds"] = fire[case]["fire_seconds"]
        row["equal_to_fire"] = check.get(case, {}).get("equal", False)
    rows.append(row)
rec = {"host": host, "cpu": cpu, "cores": int(nproc), "core": int(core), "load_average": load, "date": date,
       "commit": commit, "fire": "6.5.2 (FLINT), one thread", "cases": rows}
# the previous run on this host
prev = sorted(p for p in glob.glob(os.path.join(os.path.dirname(out), f"ibp-{host}-*.json")) if p != out)
json.dump(rec, open(out, "w"), indent=1)
before = {r["case"]: r for r in json.load(open(prev[-1]))["cases"]} if prev else {}
print(f"{date} {commit} on {host}, core {core}, load {load}" + (f"; against {os.path.basename(prev[-1])}" if prev else ""))
print(f"{'case':16} {'1 core':>9} {'all':>8} {'FIRE':>8} {'ratio':>6} {'adds':>5}  check   change")
for r in rows:
    f = r.get("fire_seconds")
    ratio = f"{r['seconds_one_core'] / f:6.2f}" if f else "     -"
    eq = {True: "equal", False: "DIFFER", None: "-"}[r.get("equal_to_fire")]
    note = ""
    b = before.get(r["case"])
    if b:
        t = r["seconds_one_core"] / max(b["seconds_one_core"], 1e-9)
        note = f"time x{t:.2f}" + (" SLOWER" if t > 1.25 else "") + \
               ("" if b["additions"] == r["additions"] else f", additions {b['additions']} -> {r['additions']}")
    print(f"{r['case']:16} {r['seconds_one_core']:8.3f}s {r['seconds_all_cores']:7.3f}s "
          f"{(f'{f:7.3f}s' if f else '       -')} {ratio} {r['additions']:5}  {eq:7} {note}")
print("written", out)
PY
rm -rf "$WORK"
exit $CHECK
