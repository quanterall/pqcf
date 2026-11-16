#!/bin/sh
# The four-point function D0 (pq_oneloop:d0/11, plan "theoretical physics instruments", F3): references for
# test/fixtures/gt_d0.json.
#  - three Euclidean configurations (all p^2, s, t < 0: F > 0 on the whole simplex): mpmath's direct 3D integral of
#    1/F^2 at 20 digits;
#  - four configurations above threshold: LoopTools 2.16 and QCDLoop 2.1.0 (double precision), recorded as they are
#    (QCDLoop returns nan for two of them).
#   docker build -t pqcf-oracle test/oracle && docker build -t pqcf-looptools comparison/oneloop/looptools
#   docker build -t pqcf-qcdloop comparison/oneloop/qcdloop && sh test/oracle/gt_d0.sh
set -e
cd "$(dirname "$0")/../.."
W=$(mktemp -d "${TMPDIR:-/tmp}/pqcf-d0-XXXXXX")
cat > "$W/eu.py" <<'PY'
import json, sys
import mpmath as mp
mp.mp.dps = 22
cfgs = [[-1, -2, -1, -3, -2, -1, 1, 2, 1, 3], [-1, -1, -1, -1, -3, -2, 1, 1, 1, 1], [-2, -1, -3, -1, -1, -4, 2, 1, 3, 1]]
out = []
for c in cfgs:
    p1, p2, p3, p4, s12, s23, m1, m2, m3, m4 = [mp.mpf(v) for v in c]
    m = [m1, m2, m3, m4]
    p = {(0, 1): p1, (1, 2): p2, (2, 3): p3, (0, 3): p4, (0, 2): s12, (1, 3): s23}
    S = [[(m[i] + m[j] - (0 if i == j else p[(min(i, j), max(i, j))]))/2 for j in range(4)] for i in range(4)]
    F = lambda x: sum(S[i][j]*x[i]*x[j] for i in range(4) for j in range(4))
    v = mp.quad(lambda x: mp.quad(lambda y: mp.quad(lambda z: 1/F([1-x, x-y, y-z, z])**2, [0, y]), [0, x]), [0, 1])
    out.append({"args": [str(a) for a in c], "mpmath": mp.nstr(v, 20)})
    print(c, mp.nstr(v, 20), flush=True)
json.dump(out, open(sys.argv[1], "w"))
PY
docker run --rm -v "$W":/w pqcf-oracle python /w/eu.py /w/eu.json
printf 'D 1 1 1 1 9 -2 1 1 1 1\nD 2 3 1 4 20 -5 1 2 1 3\nD 0.5 0.5 0.5 0.5 12 -3 1 1 1 1\nD 1 2 3 4 30 -10 2 1 2 1\n' > "$W/d.in"
docker run --rm -i pqcf-looptools lt_eval < "$W/d.in" 2>/dev/null | grep -E '^-?[0-9n]' > "$W/lt.out"
docker run --rm -i pqcf-qcdloop ql_eval < "$W/d.in" 2>/dev/null | grep -E '^-?[0-9n]' > "$W/ql.out"
python3 - "$W" <<'PY'
import json, sys
w = sys.argv[1]
eu = json.load(open(f"{w}/eu.json"))
args = [["1", "1", "1", "1", "9", "-2", "1", "1", "1", "1"], ["2", "3", "1", "4", "20", "-5", "1", "2", "1", "3"],
        ["1/2", "1/2", "1/2", "1/2", "12", "-3", "1", "1", "1", "1"], ["1", "2", "3", "4", "30", "-10", "2", "1", "2", "1"]]
lt = [l.split() for l in open(f"{w}/lt.out")]
ql = [l.split() for l in open(f"{w}/ql.out")]
phys = [{"args": a, "looptools": l, "qcdloop": q} for a, l, q in zip(args, lt, ql)]
json.dump({"euclidean": eu, "physical": phys}, open("test/fixtures/gt_d0.json", "w"), indent=0)
print(len(eu), "Euclidean,", len(phys), "physical")
PY
rm -rf "${W:?}"
