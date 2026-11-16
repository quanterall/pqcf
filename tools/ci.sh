#!/bin/sh
# The suite in groups, for CI's separate jobs; `docker run --rm pqcf-dev` still runs all of it.
#   sh tools/ci.sh static        Dialyzer, xref, clippy
#   sh tools/ci.sh fast          every EUnit module not in another group: unit and property tests, a few minutes
#   sh tools/ci.sh integration   nodes, executors, the service, the NFRs; the Python client; the guide's examples
#   sh tools/ci.sh scientific    the validations against published results and other programs; cargo test
#   sh tools/ci.sh list          the EUnit modules of each group
# A test module belongs to `fast` unless it is named below, so a new module is never left out.
set -eu

INTEGRATION="phase1_jobs_tests phase2_distributed_tests pq_rexec_tests pq_task_tests pq_service_tests nfr_tests
  guide_examples_tests phase3_api_tests phase3_requirements_tests"
SCIENTIFIC="stat_samplers_tests pq_holo_tests pq_ode_tests pq_qo_exact_tests pq_lindblad_tests pq_rpoly_tests
  pq_ibp_tests pq_recon_tests pq_tensor_tests pq_geometry_tests pq_series_tests pq_mp_tests gt_h02_tests
  gt_p01_tests phase4_noise_tests phase4_feed_forward_tests phase4_dm_blocks_tests sample_boson_tests
  pq_cumulant_tests pq_dicke_tests pq_slh_tests pq_lie_tests pq_perturb_tests pq_spectra_tests pq_dirac_tests
  pq_oneloop_tests"

all_modules() {
    for f in test/*_tests.erl; do basename "$f" .erl; done
}

# the modules of a group, as rebar3 wants them: a,b,c
modules() {
    case $1 in
        integration) echo $INTEGRATION ;;
        scientific) echo $SCIENTIFIC ;;
        fast) all_modules | grep -vxF -f /dev/fd/3 3<<EOF
$(printf '%s\n' $INTEGRATION $SCIENTIFIC)
EOF
        ;;
    esac | tr -s ' \n' ',' | sed 's/,$//'
}

# every named module must exist: a renamed file would otherwise drop out of its group silently
for m in $INTEGRATION $SCIENTIFIC; do
    [ -f "test/$m.erl" ] || { echo "tools/ci.sh: no test/$m.erl" >&2; exit 1; }
done

case ${1:-} in
    static)
        rebar3 dialyzer && rebar3 xref
        (cd native/pq_nif && cargo clippy --release --all-targets -- -D warnings) ;;
    fast)
        rebar3 eunit --module="$(modules fast)" ;;
    integration)
        rebar3 eunit --module="$(modules integration)"
        rebar3 compile
        PQCF_QUIET=1 python3 -W error::ResourceWarning -m unittest discover -s python/tests ;;
    scientific)
        rebar3 eunit --module="$(modules scientific)"
        (cd native/pq_nif && cargo test --release -q) ;;
    list)
        for g in fast integration scientific; do echo "$g: $(modules $g)"; done ;;
    *)
        echo "usage: sh tools/ci.sh static|fast|integration|scientific|list" >&2; exit 2 ;;
esac
