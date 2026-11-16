#!/bin/sh
# QCD's renormalisation group (pq_qcd, plan "theoretical physics instruments", QCD1): two references for
# test/fixtures/gt_qcd.json.
#  - mpmath (test/oracle/gt_qcd.py): the beta and gamma_m coefficients for nf = 3..6 to 40 digits, alpha_s and the
#    MS-bar mass by an ODE solver at 45 digits, the decoupling summed directly;
#  - CRunDec 3.1.2 (comparison/qcd/crundec): AsmMSrunexact (its exact running of alpha_s and m, to about its step
#    tolerance 1e-10), DecAsDownMS and DecAsUpMS.
#   docker build -t pqcf-oracle test/oracle && docker build -t pqcf-crundec comparison/qcd/crundec
#   sh test/oracle/gt_qcd.sh
set -e
cd "$(dirname "$0")/../.."
W=$(mktemp -d "${TMPDIR:-/tmp}/pqcf-qcd-XXXXXX")
cat > "$W/cases.json" <<'J'
{"running": [
 {"as0": "59/500", "mu0": "911876/10000", "mu1": "10", "nf": 5, "loops": 5, "m0": "209/50"},
 {"as0": "59/500", "mu0": "911876/10000", "mu1": "1000", "nf": 5, "loops": 4, "m0": "209/50"},
 {"as0": "2253/10000", "mu0": "209/50", "mu1": "911876/10000", "nf": 5, "loops": 5, "m0": "209/50"},
 {"as0": "2253/10000", "mu0": "209/50", "mu1": "911876/10000", "nf": 5, "loops": 1, "m0": "209/50"},
 {"as0": "2253/10000", "mu0": "209/50", "mu1": "911876/10000", "nf": 5, "loops": 2, "m0": "209/50"},
 {"as0": "2253/10000", "mu0": "209/50", "mu1": "911876/10000", "nf": 5, "loops": 3, "m0": "209/50"},
 {"as0": "38/100", "mu0": "3", "mu1": "2", "nf": 4, "loops": 5, "m0": "1"},
 {"as0": "1/4", "mu0": "3", "mu1": "127/100", "nf": 3, "loops": 5, "m0": "93/1000"},
 {"as0": "1/10", "mu0": "200", "mu1": "3000", "nf": 6, "loops": 5, "m0": "163"},
 {"as0": "1/10", "mu0": "200", "mu1": "100", "nf": 6, "loops": 3, "m0": "163"}],
 "decoupling": [
 {"direction": "down", "as": "59/500", "mh": "209/50", "mu": "209/50", "nl": 4, "loops": 4},
 {"direction": "down", "as": "59/500", "mh": "163", "mu": "911876/10000", "nl": 5, "loops": 4},
 {"direction": "down", "as": "59/500", "mh": "163", "mu": "326", "nl": 5, "loops": 3},
 {"direction": "down", "as": "21/100", "mh": "127/100", "mu": "3", "nl": 3, "loops": 4},
 {"direction": "up", "as": "2253/10000", "mh": "209/50", "mu": "209/50", "nl": 4, "loops": 4},
 {"direction": "up", "as": "2253/10000", "mh": "209/50", "mu": "8", "nl": 4, "loops": 2},
 {"direction": "up", "as": "3/10", "mh": "127/100", "mu": "2", "nl": 3, "loops": 4},
 {"direction": "up", "as": "1/10", "mh": "163", "mu": "100", "nl": 5, "loops": 3}]}
J
docker run --rm -v "$W":/w -v "$PWD/test/oracle":/o:ro pqcf-oracle python /o/gt_qcd.py /w/cases.json /w/mp.json
python3 - "$W" <<'PY'
import json, sys
from fractions import Fraction as F
w = sys.argv[1]
c = json.load(open(f"{w}/cases.json"))
f = lambda x: repr(float(F(x)))
with open(f"{w}/crd.in", "w") as o:
    for r in c["running"]:
        o.write(f"X {f(r['m0'])} {f(r['as0'])} {f(r['mu0'])} {f(r['mu1'])} {r['nf']} {r['loops']}\n")
    for d in c["decoupling"]:
        k = "D" if d["direction"] == "down" else "U"
        o.write(f"{k} {f(d['as'])} {f(d['mh'])} {f(d['mu'])} {d['nl']} {d['loops']}\n")
PY
docker run --rm -i pqcf-crundec crd_eval < "$W/crd.in" > "$W/crd.out"
python3 - "$W" <<'PY'
import json, sys
w = sys.argv[1]
mp_ = json.load(open(f"{w}/mp.json"))
lines = [l.split() for l in open(f"{w}/crd.out")]
nr = len(mp_["running"])
for r, l in zip(mp_["running"], lines[:nr]):
    r["crundec"] = {"mass": l[0], "alphas": l[1]}
    assert abs(float(r["alphas"]) - float(l[1])) < 1e-8 * float(l[1]), (r, l)
    assert abs(float(r["mass"]) - float(l[0])) < 1e-8 * float(l[0]), (r, l)
for d, l in zip(mp_["decoupling"], lines[nr:]):
    d["crundec"] = l[0]
    assert abs(float(d["value"]) - float(l[0])) < 1e-14, (d, l)
json.dump(mp_, open("test/fixtures/gt_qcd.json", "w"), indent=0)
print(len(mp_["coefficients"]), "coefficient sets,", nr, "running cases,", len(mp_["decoupling"]),
      "decoupling cases; mpmath and CRunDec agree (1e-8 running, 1e-14 decoupling)")
PY
rm -rf "${W:?}"
