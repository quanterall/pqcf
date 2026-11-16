#!/bin/sh
# d-dimensional Dirac traces by FORM (tracen), the reference for pq_dirac (plan "theoretical physics instruments", F2).
# 60 fully contracted traces of up to 12 gamma matrices (momenta p1..p4, repeated indices), fixed by a seed; FORM's
# results written to test/fixtures/gt_dirac.json with the dot products as dot_pI_pJ and the dimension as d.
#   docker build -t pqcf-dev . && docker build -t pqcf-form -f comparison/form/Dockerfile . && sh test/oracle/gt_dirac.sh
set -e
cd "$(dirname "$0")/../.."
W=$(mktemp -d "${TMPDIR:-/tmp}/pqcf-dirac-XXXXXX")
python3 - "$W" <<'PY'
import json, random, sys
w = sys.argv[1]
rnd = random.Random(20261116)
cases = []
while len(cases) < 60:
    n_idx = rnd.randint(0, 5)
    n_vec = rnd.randrange(0, 13 - 2 * n_idx)
    if (2 * n_idx + n_vec) % 2 or 2 * n_idx + n_vec < 2:
        continue
    items = [f"i{k}" for k in range(1, n_idx + 1)] * 2 + [f"p{rnd.randint(1, 4)}" for _ in range(n_vec)]
    rnd.shuffle(items)
    cases.append(items)
json.dump(cases, open(f"{w}/cases.json", "w"))
with open(f"{w}/t.frm", "w") as f:
    f.write("Dimension d;\nVectors p1,p2,p3,p4;\nIndices i1,i2,i3,i4,i5;\n")
    for k, c in enumerate(cases):
        f.write(f"Local F{k} = g_(1,{','.join(c)});\n")
    f.write("tracen,1;\nFormat nospaces;\nFormat 255;\nPrint;\n.end\n")
PY
docker run --rm -v "$W":/w -w /w pqcf-form sh -c 'form -q t.frm' > "$W/out.txt"
python3 - "$W" <<'PY'
import json, re, sys
w = sys.argv[1]
cases = json.load(open(f"{w}/cases.json"))
text = open(f"{w}/out.txt").read()
exprs = {}
for m in re.finditer(r"F(\d+)=\s*(.*?);", text, re.S):
    e = re.sub(r"\s+", "", m.group(2))
    e = re.sub(r"(p\d)\.(p\d)", lambda x: "dot_" + "_".join(sorted([x.group(1), x.group(2)])), e)
    exprs[int(m.group(1))] = e
out = [{"gammas": c, "form": exprs[k]} for k, c in enumerate(cases)]
json.dump(out, open("test/fixtures/gt_dirac.json", "w"), indent=0)
print(len(out), "traces; longest", max(len(c) for c in cases))
PY
rm -rf "${W:?}"
