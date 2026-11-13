#!/bin/sh
# Whole-genome preparation (README.md), after prepare.sh (reference, indexes, truth set): the HG002 NovaSeq 35x
# whole-genome reads (46 GB BAM from the DeepVariant team's public case-study data) streamed back into read pairs
# and cut into numbered gzipped chunks of PAIRS read pairs, without uncompressed intermediates.
#   sh research/variants/prepare_wgs.sh [pairs per chunk] [cpuset]
set -e
cd "$(dirname "$0")/data"
PAIRS=${1:-6000000}
CPUS=${2:-5-9,15-19}
BAM=HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam
[ -s reads.all.chunk000.R1.fastq.gz ] && { echo "chunks exist"; exit 0; }
mkdir -p collate-tmp
docker run --rm --user "$(id -u):$(id -g)" --cpuset-cpus="$CPUS" -v "$PWD":/d -v "$PWD/..":/s:ro -w /d pqcf-variants \
  bash /s/split_fastq.sh "$BAM" /d/collate-tmp "$PAIRS"
rmdir collate-tmp 2>/dev/null || true
N=$(ls reads.all.chunk*.R1.fastq.gz 2>/dev/null | wc -l)
echo "chunks $N"
[ "$N" -gt 0 ]
