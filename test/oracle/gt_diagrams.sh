#!/bin/sh
# Feynman diagrams of QED and QCD by FORM 5's generator (diagrams_, Kaneko's grcc), the reference for pq_diagrams
# (plan "theoretical physics instruments", F4). For each case FORM's diagrams are terms whose coefficient is
# sign / S (sign = (-1)^(closed fermion loops)); the script writes test/fixtures/gt_diagrams.json with, per case, the
# sorted list of those coefficients (the multiset pq_diagrams must reproduce) and the number of diagrams.
#   docker build -t pqcf-form -f comparison/form/Dockerfile . && sh test/oracle/gt_diagrams.sh
set -e
cd "$(dirname "$0")/../.."
W=$(mktemp -d "${TMPDIR:-/tmp}/pqcf-diagrams-XXXXXX")
python3 - "$W" <<'PY'
import json, sys
w = sys.argv[1]
# name, model, incoming, outgoing, loops, options
cases = []
for L in (1, 2, 3):
    cases += [("qed_pho_se", "QED", ["pho"], ["pho"], L, "1pi"), ("qed_ele_se", "QED", ["ele"], ["ele"], L, "1pi")]
for L in (1, 2):
    cases += [("qed_vertex", "QED", ["ele", "pho"], ["ele"], L, "1pi")]
for L in (2, 3, 4):
    cases += [("qed_vacuum", "QED", [], [], L, "none")]
for L in (1, 2, 3):
    cases += [("qcd_glu_se", "QCD", ["glu"], ["glu"], L, "1pi"), ("qcd_gho_se", "QCD", ["gho"], ["gho"], L, "1pi"),
              ("qcd_qua_se", "QCD", ["qua"], ["qua"], L, "1pi")]
for L in (1, 2):
    cases += [("qcd_ggg", "QCD", ["glu", "glu", "glu"], [], L, "1pi"), ("qcd_qqg", "QCD", ["qua", "glu"], ["qua"], L, "1pi")]
json.dump(cases, open(f"{w}/cases.json", "w"))
with open(f"{w}/d.frm", "w") as f:
    f.write("""Off Statistics;
Vectors p1,...,p4,k1,...,k20;
Set kk:k1,...,k20;
Symbols e,g;
Model QED;
  Particle ele,ELE,-2;
  Particle pho,+3;
  Vertex ELE,ele,pho:e;
EndModel;
Model QCD;
  Particle qua,QUA,-2;
  Particle gho,GHO,-1;
  Particle glu,+3;
  Vertex QUA,qua,glu:g;
  Vertex GHO,gho,glu:g;
  Vertex glu,glu,glu:g;
  Vertex glu,glu,glu,glu:g^2;
EndModel;
""")
    for k, (name, model, inc, out, L, opt) in enumerate(cases):
        n = len(inc) + len(out)
        moms = ",".join(f"p{i}" for i in range(1, n + 1))
        o = "`ONEPI_'" if opt == "1pi" else "0"
        f.write(f"L F{k} = diagrams_({model},{{{','.join(inc)}}},{{{','.join(out)}}},{{{moms}}},kk,{L},{o});\n")
    f.write(".sort\nid topo_(?x) = 1;\nid node_(?x) = 1;\n* each diagram is now its coefficient times e^n g^m: keep them apart\n")
    f.write("Format nospaces;\nFormat 255;\n.sort\n")
    # count and coefficients: print every term of every expression separately, before the products collapse
    f.write(".end\n")
PY
# FORM collapses equal terms once node_ and topo_ are gone, so the coefficients are read from the full output
sed -i 's/^id topo_(?x) = 1;$//; s/^id node_(?x) = 1;$//' "$W/d.frm"
sed -i 's/^\.end$/Print +s;\n.end/' "$W/d.frm"
docker run --rm -v "$W":/w -w /w pqcf-form sh -c 'form -q d.frm' > "$W/out.txt"
python3 - "$W" <<'PY'
import json, re, sys
from fractions import Fraction
w = sys.argv[1]
cases = json.load(open(f"{w}/cases.json"))
text = open(f"{w}/out.txt").read()
out = []
for k, (name, model, inc, out_, L, opt) in enumerate(cases):
    m = re.search(rf"\n\s*F{k}\s*=(.*?);", text, re.S)
    body = re.sub(r"\s+", "", m.group(1)) if m else ""
    coefs = []
    if body and body != "0":
        # split into terms at top-level signs
        terms, depth, cur = [], 0, ""
        for ch in body:
            if ch in "+-" and depth == 0 and cur:
                terms.append(cur); cur = ch
            else:
                cur += ch
                depth += (ch == "(") - (ch == ")")
        terms.append(cur)
        for t in terms:
            mm = re.match(r"([+-]?)(\d+(?:/\d+)?\*)?", t)
            sgn = -1 if mm.group(1) == "-" else 1
            c = Fraction(mm.group(2)[:-1]) if mm.group(2) else Fraction(1)
            coefs.append(sgn * c)
    coefs.sort()
    out.append({"name": name, "model": model, "in": inc, "out": out_, "loops": L, "options": opt,
                "count": len(coefs), "coefficients": [str(c) for c in coefs]})
json.dump(out, open("test/fixtures/gt_diagrams.json", "w"), indent=0)
for c in out:
    print(c["name"], c["loops"], c["count"])
PY
rm -rf "${W:?}"
