#!/bin/sh
# The three-point function C0 (pq_oneloop:c0/7, plan "theoretical physics instruments", F3): 40 configurations drawn
# (p1^2, p2^2, p3^2; m1^2, m2^2, m3^2), exact fractions fixed by a seed, over space-like, below-threshold and
# above-threshold kinematics, lambda(p1^2, p2^2, p3^2) of both signs, massless lines (infrared-divergent and
# lambda = 0 configurations left out, collinear-divergent ones dropped after the draw). References:
#  - mpmath (test/oracle/gt_c0.py): contour deformation for all, the direct integral at 30 digits where D > 0;
#  - QCDLoop 2.1.0 (comparison/oneloop/qcdloop) and LoopTools 2.16 (comparison/oneloop/looptools), recorded as they
#    are: in our aarch64 builds each of them is wrong in some lambda < 0 configurations (README there).
#   docker build -t pqcf-qcdloop comparison/oneloop/qcdloop && docker build -t pqcf-looptools comparison/oneloop/looptools
#   sh test/oracle/gt_c0.sh
set -e
cd "$(dirname "$0")/../.."
W=$(mktemp -d "${TMPDIR:-/tmp}/pqcf-c0-XXXXXX")
python3 - "$W" <<'PY'
import json, random, sys
from fractions import Fraction as F
w = sys.argv[1]
rnd = random.Random(20261116)
masses = [F(0), F(1, 4), F(1), F(2), F(3)]
moms = [F(-3), F(-1), F(-1, 2), F(0), F(1, 2), F(1), F(3), F(5), F(9), F(20)]
def lam(x, y, z): return x*x + y*y + z*z - 2*(x*y + y*z + x*z)
def ir(p1, p2, p3, m1, m2, m3):
    return any(m == 0 and pa == ma and pb == mb
               for m, pa, ma, pb, mb in [(m1, p1, m2, p3, m3), (m2, p1, m1, p2, m3), (m3, p3, m1, p2, m2)])
cfgs = []
seen = set()
while len(cfgs) < 40:
    c = tuple([rnd.choice(moms) for _ in range(3)] + [rnd.choice(masses) for _ in range(3)])
    if c in seen or lam(*c[:3]) == 0 or ir(*c) or c[:3] == (0, 0, 0):
        continue
    seen.add(c)
    cfgs.append(c)
# collinear-divergent ones (p^2 = 0 between two massless propagators) are left out after the draw
def coll(p1, p2, p3, m1, m2, m3):
    return any(p == 0 and mi == 0 and mj == 0 for p, mi, mj in [(p1, m1, m2), (p2, m2, m3), (p3, m1, m3)])
cfgs = [c for c in cfgs if not coll(*c)]
json.dump([[str(x) for x in c] for c in cfgs], open(f"{w}/cfg.json", "w"))
with open(f"{w}/c.in", "w") as f:
    for c in cfgs:
        f.write("C " + " ".join(repr(float(x)) for x in c) + "\n")
PY
docker run --rm -i pqcf-qcdloop ql_eval < "$W/c.in" 2>/dev/null | grep -E '^-?[0-9n]' > "$W/ql.out"
timeout 60 docker run --rm -i pqcf-looptools lt_eval < "$W/c.in" 2>/dev/null | grep -E '^-?[0-9n]' > "$W/lt.out" || true
docker run --rm -v "$W":/w -v "$PWD/test/oracle":/o:ro pqcf-oracle python /o/gt_c0.py /w/cfg.json /w/mp.json
python3 - "$W" <<'PY'
import json, sys
w = sys.argv[1]
mp_ = json.load(open(f"{w}/mp.json"))
ql = [l.split() for l in open(f"{w}/ql.out")]
lt = [l.split() for l in open(f"{w}/lt.out")]
assert len(ql) == len(mp_)
out, bad_ql, bad_lt = [], 0, 0
for k, rec in enumerate(mp_):
    ref = complex(float(rec["deformed"][0]), float(rec["deformed"][1]))
    def agree(v):
        try:
            z = complex(float(v[0]), float(v[1]))
        except (ValueError, IndexError):
            return False
        return abs(z - ref) <= 1e-9*max(1, abs(ref))
    rec["qcdloop"] = ql[k]
    rec["qcdloop_agrees"] = agree(ql[k])
    rec["looptools"] = lt[k] if k < len(lt) else None
    rec["looptools_agrees"] = agree(lt[k]) if k < len(lt) else None
    bad_ql += not rec["qcdloop_agrees"]
    bad_lt += rec["looptools_agrees"] is False or rec["looptools_agrees"] is None
    out.append(rec)
json.dump(out, open("test/fixtures/gt_c0.json", "w"), indent=0)
print(len(out), "configurations;", sum(r["stable"] for r in out), "stable under the deformation;",
      sum("direct" in r for r in out), "with the direct integral; QCDLoop disagrees on", bad_ql,
      "; LoopTools disagrees or gives nothing on", bad_lt)
PY
rm -rf "${W:?}"
