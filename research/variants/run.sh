#!/bin/sh
# The variant benchmark (README.md) on one chromosome of GIAB HG002: the sliced pipeline, its determinism and
# equivalence checks, and its accuracy against the truth set. Needs prepare.sh first.
# Writes research/variants/results/<host>-<date>.json and RESULTS.md.
#   sh research/variants/run.sh [chrom] [cpuset] [threads]
set -e
cd "$(dirname "$0")"
CHROM=${1:-chr20}
CPUS=${2:-5-9,15-19}
THREADS=${3:-10}
HOST=${PQ_BENCH_HOST:-$(hostname -s)}
DATE=$(date -u +%Y-%m-%dT%H%MZ)
docker build -q -t pqcf-variants docker >/dev/null
REF=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
LEN=$(awk -v c="$CHROM" '$1 == c {print $2}' data/$REF.fai)
WORK=$(mktemp -d "${TMPDIR:-/tmp}/pqcf-variants-XXXXXX")
RUN="docker run --rm --user $(id -u):$(id -g) --cpuset-cpus=$CPUS -v $PWD:/v -v $WORK:/work -w /v pqcf-variants"
CHUNKS=$(ls data/reads.$CHROM.chunk*.R1.fastq | sed 's|.*chunk\([0-9]*\).*|\1|' | sort | while read i; do
    printf '["/v/data/reads.%s.chunk%s.R1.fastq","/v/data/reads.%s.chunk%s.R2.fastq"],' $CHROM $i $CHROM $i; done)
spec() { # name, call_slots, slice, margin, skip_align
    printf '{"exec":"local","align_slots":1,"call_slots":%s,"threads":%s,"ref":"/v/data/%s","index":"/v/data/%s.mmi",
"sample":"HG002","chunks":[%s],"chrom":"%s","length":%s,"slice":%s,"margin":%s,"work":"/work",
"out":"/work/%s.vcf.gz","skip_align":%s}' "$2" "$THREADS" "$REF" "$REF" "${CHUNKS%,}" "$CHROM" "$LEN" "$3" "$4" "$1" "$5" \
        > "$WORK/$1.json"
}
erl_run() {
    $RUN sh -c "mkdir -p /tmp/b && erlc -o /tmp/b pq_genome.erl && erl -noshell -pa /tmp/b -eval 'pq_genome:main(\"/work/$1.json\"), halt().'"
}
spec sliced8 8 1000000 1000 false;   erl_run sliced8
spec sliced1 1 1000000 1000 true;    erl_run sliced1
spec unsliced 1 "$LEN" 0 true;       erl_run unsliced
CHECK=0
$RUN python3 evaluate.py /work "$CHROM" "$REF" "results/$HOST-$DATE.json" "$HOST" "$DATE" || CHECK=1
$RUN python3 report.py "results/$HOST-$DATE.json" > RESULTS.md
if [ -n "$KEEP" ]; then echo "kept $WORK"; else
    docker run --rm -v "$WORK":/w pqcf-variants sh -c 'rm -rf /w/*'
    rm -rf "${WORK:?}"
fi
echo "written results/$HOST-$DATE.json"
exit $CHECK
