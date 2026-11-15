#!/bin/sh
# The three-loop QCD beta function (comparison/beta/README.md), every step from the diagrams, into results/:
#   1. FORM + Forcer: the gluon and ghost bare self-energies at one, two and three loops (se.frm);
#   2. FORM writes the one- and two-loop integrands (se.frm -D DUMP); PQCF reduces them (beta_pqcf:reduce/3);
#   3. Forcer gives the values of PQCF's master integrals;
#   4. PQCF renormalises both and extracts the beta function (renorm.erl, beta_run.erl).
# Uses the working tree. Images: pqcf-dev, pqcf-form, pqcf-forcer, built here. FORM's colour package (color.h, van
# Ritbergen, Schellekens and Vermaseren) is fetched at a fixed commit and checked by its SHA-256.
# Afterwards `git diff comparison/beta/results` shows any change; the script prints whether the beta functions are
# the known ones and exits 1 if not.
#   sh comparison/beta/run.sh [cores]        # default 0-7; about 15 minutes on 8 cores, most of it the 3-loop gluon
set -e
cd "$(dirname "$0")/../.."
CORES=${1:-0-7}
docker build -q -t pqcf-form -f comparison/form/Dockerfile . >/dev/null && docker build -q -t pqcf-forcer comparison/forcer >/dev/null
NW=$(docker run --rm --cpuset-cpus="$CORES" pqcf-forcer nproc)
COLOR_REF=47f4739b732e6aa30d213fbff0aac817aede2f11
COLOR_SHA=616e225880fda51c5e1d37b8730a520b1994f05d70713e7e8b2fca2e2ab3fd00
docker build -q -t pqcf-dev . >/dev/null
WORK=$(mktemp -d "${TMPDIR:-/tmp}/pqcf-beta-XXXXXX")
cp comparison/beta/se.frm "$WORK/"
U="-u $(id -u):$(id -g)"
FORM="docker run --rm $U --cpuset-cpus=$CORES -v $WORK:/w -w /w pqcf-forcer"
PQ="docker run --rm $U --cpuset-cpus=$CORES -v $PWD:/src -v $WORK:/w -w /src pqcf-dev"
$FORM sh -c "wget -q -O color.h https://raw.githubusercontent.com/form-dev/form-packages/$COLOR_REF/color/color.h &&
             echo '$COLOR_SHA  color.h' | sha256sum -c - >/dev/null && cp /opt/forcer/forcer/Masters.prc ."
for L in 1 2 3; do
    for F in glu gho; do
        echo "FORM + Forcer: $F, $L loops"
        $FORM sh -c "tform -w$NW -D LOOPS=$L -D FIELD=$F -D LANDAU=1 se.frm > res_${F}_$L.log 2>&1"
    done
done
for L in 1 2; do
    for F in glu gho; do
        echo "PQCF: $F, $L loops"
        $FORM sh -c "form -D LOOPS=$L -D FIELD=$F -D LANDAU=1 -D DUMP=dump_${F}_$L.txt se.frm > dump_${F}_$L.log 2>&1"
    done
done
$PQ sh -c 'rebar3 compile >/dev/null && mkdir -p /tmp/b && erlc -o /tmp/b comparison/forcer/*.erl comparison/beta/*.erl &&
           erl -noshell -pa _build/default/lib/*/ebin -pa /tmp/b -eval "
             [beta_pqcf:reduce(\"/w\", F, L) || L <- [1, 2], F <- [\"glu\", \"gho\"]], halt()." >/dev/null'
$FORM sh -c 'for f in masters_glu_1 masters_gho_1 masters_glu_2 masters_gho_2; do form -q $f.frm > $f.out 2>&1; done'
$PQ sh -c 'mkdir -p /tmp/b && erlc -o /tmp/b comparison/forcer/*.erl comparison/beta/*.erl &&
           erl -noshell -pa _build/default/lib/*/ebin -pa /tmp/b -eval "
             beta_run:write_results(\"/w\", \"/src/comparison/beta/results\"), halt()."'
rm -rf "${WORK:?}"
# the known beta functions (Tarasov, Vladimirov, Zharkov 1980; Larin, Vermaseren, hep-ph/9302208)
CHECK=0
cat > "${TMPDIR:-/tmp}/pqcf-beta-known.txt" <<'EOF'
beta_0 = -4/3*I2R*nf + 11/3*cA
beta_1 = -20/3*I2R*cA*nf + -4*I2R*cR*nf + 34/3*cA^2
beta_2 = -205/9*I2R*cA*cR*nf + -1415/27*I2R*cA^2*nf + 2*I2R*cR^2*nf + 158/27*I2R^2*cA*nf^2 + 44/9*I2R^2*cR*nf^2 + 2857/54*cA^3
EOF
grep '^beta_' comparison/beta/results/renormalisation_form_forcer_3loops.txt | cmp -s - "${TMPDIR:-/tmp}/pqcf-beta-known.txt" ||
    { echo "the three-loop beta functions differ from the known ones"; CHECK=1; }
rm -f "${TMPDIR:-/tmp}/pqcf-beta-known.txt"
grep -q "Finite: true" comparison/beta/results/renormalisation_form_forcer_3loops.txt || CHECK=1
grep -q "Finite: true" comparison/beta/results/renormalisation_pqcf_2loops.txt || CHECK=1
git --no-pager diff --stat comparison/beta/results
[ $CHECK = 0 ] && echo "beta functions: as known" || echo "beta functions: CHECK FAILED"
exit $CHECK
