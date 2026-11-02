#!/bin/sh
# The standing comparison of pq_symbolic with FORM, by hand (the weekly workflow does the same on GitHub's runners):
# builds the dev and the FORM images, runs pq_form's tests with FORM, then the reference tasks; rewrites
# comparison/form/RESULTS.md and appends to comparison/form/history.csv.
#   COMPARISON_MACHINE=<label> [CPUSET=8] sh comparison/form/run.sh
set -e
cd "$(dirname "$0")/../.."
docker build -q -t pqcf-dev . >/dev/null
docker build -q -t pqcf-form -f comparison/form/Dockerfile . >/dev/null
docker run --rm pqcf-form rebar3 eunit --module=pq_form_tests
docker run --rm ${CPUSET:+--cpuset-cpus=$CPUSET} -e COMPARISON_MACHINE="${COMPARISON_MACHINE:-local}" \
  -e PQCF_COMMIT="$(git rev-parse --short HEAD)" \
  -e COMPARISON_CPU="$(lscpu | sed -n 's/^Model name: *//p' | sort -u | paste -sd+)" -e HOME=/tmp --user "$(id -u):$(id -g)" \
  -v "$PWD/comparison/form":/out pqcf-form \
  sh -c 'erlc -o /tmp /out/compare.erl && erl +S 1 -pa _build/default/lib/*/ebin -pa /tmp -noshell -eval "compare:main(\"/out\")"'
