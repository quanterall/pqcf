#!/bin/sh
# NFR-01 baseline (spike S1): our single-thread permanent against thewalrus on the same machine, both pinned to
# one core. Writes bench/results/<host>-<date>.json and prints the ratios. Needs docker and the two images:
#   docker build -t pqcf-dev . && docker build -t pqcf-oracle test/oracle && bench/run.sh [core]
set -e
cd "$(dirname "$0")/.."
CORE=${1:-0}
HOST=${PQ_BENCH_HOST:-$(hostname -s)}   # a label for the machine (PQ_BENCH_HOST), else its host name
DATE=$(date -u +%Y-%m-%dT%H%MZ)
OUT=bench/results/$HOST-$DATE.json
OURS=$(docker run --rm --cpuset-cpus="$CORE" pqcf-dev sh -c 'cd native/pq_nif && cargo run -q --release --example bench_permanent')
THEIRS=$(docker run --rm --cpuset-cpus="$CORE" -e NUMBA_NUM_THREADS=1 -e OMP_NUM_THREADS=1 -v "$PWD/bench":/b pqcf-oracle python /b/thewalrus_permanent.py)
GOV=$(cat /sys/devices/system/cpu/cpu$CORE/cpufreq/scaling_governor 2>/dev/null || echo unknown)
CPU=$(lscpu | sed -n 's/^Model name: *//p' | head -1)
python3 - "$OUT" "$HOST" "$DATE" "$CORE" "$GOV" "$CPU" "$(uname -m)" "$(uname -r)" "$(nproc)" "$OURS" "$THEIRS" <<'PY'
import json, sys
out, host, date, core, gov, cpu, arch, kernel, nproc, ours, theirs = sys.argv[1:]
ours, theirs = json.loads(ours), json.loads(theirs)
rows = []
for o, t in zip(ours, theirs["results"]):
    best = min(t["bbfg"], t["glynn"], t["ryser"])
    ok = abs(complex(*o["value"]) - complex(*t["value"])) <= 1e-8 * abs(complex(*t["value"]))
    rows.append({"n": o["n"], "ours": o["median_seconds"], "thewalrus_best": best, "thewalrus": t,
                 "ratio": o["median_seconds"] / best, "values_agree": ok})
rec = {"hardware": {"host": host, "cpu": cpu, "arch": arch, "kernel": kernel, "cores": int(nproc),
                    "pinned_core": int(core), "governor": gov},
       "date": date, "thewalrus": theirs["thewalrus"], "nfr01_limit": 1.5, "results": rows}
json.dump(rec, open(out, "w"), indent=1)
print(f"{cpu} ({arch}), core {core}, governor {gov}")
for r in rows:
    print(f"n = {r['n']}: ours {r['ours']*1e3:.2f} ms, thewalrus best {r['thewalrus_best']*1e3:.2f} ms, "
          f"ratio {r['ratio']:.2f} ({'within' if r['ratio'] <= 1.5 else 'OVER'} 1.5), values agree: {r['values_agree']}")
print("written", out)
PY
