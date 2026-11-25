#!/bin/sh
# PQCF results certified by Lean 4: a lean-check node (Lean 4 + Mathlib, image lean-check) and a PQCF node (image
# pqcf-dev) on an internal Docker network (no outside access), joined by Erlang distribution. PQCF writes the Lean
# statement and proof for each result (pq_lean), lean-check returns the verdict.
#   git clone https://github.com/quanterall/erlang-lean4 ../erlang-lean4
#   docker build -t pqcf-dev . && (cd ../erlang-lean4 && docker build -t lean-check docker)
#   LEAN_CHECK_SRC=../erlang-lean4 sh examples/lean/run.sh                         # the three examples
#   LEAN_CHECK_SRC=../erlang-lean4 DEMO=certify_published sh examples/lean/run.sh  # published results (T7, T8)
set -e
cd "$(dirname "$0")/../.."
LC=$(cd "${LEAN_CHECK_SRC:?path to the erlang-lean4 checkout}" && pwd)
DEMO=${DEMO:-pq_lean_demo}
ARGS=${ARGS:-}
NET=pqcf-lean-$$
COOKIE=pqcf-lean-$$
docker network create --internal "$NET" >/dev/null
trap 'docker rm -f "lc-$$" >/dev/null 2>&1; docker network rm "$NET" >/dev/null 2>&1' EXIT
docker run -d --rm --name "lc-$$" --hostname lc --network "$NET" -v "$LC":/src:ro lean-check sh -c "
  mkdir /w && cd /src && tar --exclude=./_build -cf - . | (cd /w && tar xf -) && cd /w && rebar3 compile >/dev/null &&
  exec erl -sname lean_check -setcookie $COOKIE -noshell -pa _build/default/lib/*/ebin \
    -eval 'application:ensure_all_started(lean_check)'" >/dev/null
docker run --rm --hostname pq --network "$NET" -v "$PWD":/src:ro pqcf-dev sh -c "
  mkdir /w && cd /src && tar --exclude=./_build -cf - . | (cd /w && tar xf -) && cd /w && rebar3 compile >/dev/null &&
  erlc -o /tmp examples/lean/pq_lean_demo.erl examples/lean/certify_published.erl comparison/form/compare.erl &&
  erl -sname pqcf -setcookie $COOKIE -noshell -pa _build/default/lib/*/ebin -pa /tmp -eval '
    Node = lean_check@lc,
    Wait = fun W(0) -> error(lean_check_not_ready);
               W(N) -> case catch erpc:call(Node, lean_check, status, []) of
                           #{free := F} when F > 0 -> ok;
                           _ -> timer:sleep(1000), W(N - 1)
                       end end,
    Wait(900),
    $DEMO:main(Node$ARGS), halt().'"
