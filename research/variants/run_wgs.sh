#!/bin/sh
# The whole-genome benchmark (README.md): every chromosome of GIAB HG002, sliced and distributed, then checked on
# chromosome 20 (1 worker against 8, sliced against unsliced) and evaluated against the truth set (chr1-22).
# Needs prepare.sh and prepare_wgs.sh. Writes research/variants/results/wgs-<host>-<date>.json and RESULTS_WGS.md.
#   sh research/variants/run_wgs.sh [cpuset] [threads] [call workers]
set -e
cd "$(dirname "$0")"
CPUS=${1:-5-9,15-19}
THREADS=${2:-10}
SLOTS=${3:-8}
HOST=${PQ_BENCH_HOST:-$(hostname -s)}
DATE=$(date -u +%Y-%m-%dT%H%MZ)
docker build -q -t pqcf-variants docker >/dev/null
REF=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
REGIONS=$(awk '$1 ~ /^chr([0-9]+|X|Y)$/ {printf "[\"%s\",%s],", $1, $2}' data/$REF.fai)
LEN20=$(awk '$1 == "chr20" {print $2}' data/$REF.fai)
# WORK=<dir> SKIP_ALIGN=true: resume in a kept work directory, with its chunks' alignments
WORK=${WORK:-$(mktemp -d "${TMPDIR:-/tmp}/pqcf-wgs-XXXXXX")}
SKIP_ALIGN=${SKIP_ALIGN:-false}
RUN="docker run --rm --user $(id -u):$(id -g) --cpuset-cpus=$CPUS -v $PWD:/v -v $WORK:/work -w /v pqcf-variants"
CHUNKS=$(ls data/reads.all.chunk*.R1.fastq.gz | sed 's|.*chunk\([0-9]*\).*|\1|' | sort | while read i; do
    printf '["/v/data/reads.all.chunk%s.R1.fastq.gz","/v/data/reads.all.chunk%s.R2.fastq.gz"],' $i $i; done)
spec() { # name, call_slots, regions json, slice, margin, skip_align
    printf '{"exec":"local","align_slots":1,"call_slots":%s,"threads":%s,"ref":"/v/data/%s","index":"/v/data/%s.mmi",
"sample":"HG002","chunks":[%s],"regions":[%s],"slice":%s,"margin":%s,"work":"/work","out":"/work/%s.vcf.gz",
"skip_align":%s}' "$2" "$THREADS" "$REF" "$REF" "${CHUNKS%,}" "$3" "$4" "$5" "$1" "$6" > "$WORK/$1.json"
}
erl_run() {
    $RUN sh -c "mkdir -p /tmp/b && erlc -o /tmp/b pq_genome.erl && erl -noshell -pa /tmp/b -eval 'pq_genome:main(\"/work/$1.json\"), halt().'"
}
spec wgs "$SLOTS" "${REGIONS%,}" 1000000 1000 "$SKIP_ALIGN";       erl_run wgs
spec chr20_sliced1 1 "[\"chr20\",$LEN20]" 1000000 1000 true;       erl_run chr20_sliced1
spec chr20_unsliced 1 "[\"chr20\",$LEN20]" "$LEN20" 0 true;        erl_run chr20_unsliced
CHECK=0
$RUN python3 evaluate_wgs.py /work "$REF" "results/wgs-$HOST-$DATE.json" "$HOST" "$DATE" "$SLOTS" "$THREADS" || CHECK=1
$RUN python3 report_wgs.py "results/wgs-$HOST-$DATE.json" > RESULTS_WGS.md
if [ -n "$KEEP" ]; then echo "kept $WORK"; else
    docker run --rm -v "$WORK":/w pqcf-variants sh -c 'rm -rf /w/*'
    rm -rf "${WORK:?}"
fi
echo "written results/wgs-$HOST-$DATE.json"
exit $CHECK
