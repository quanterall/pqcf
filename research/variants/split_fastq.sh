#!/bin/bash
# Used by prepare_wgs.sh inside the pqcf-variants image: BAM -> pairs -> numbered gzipped chunks of PAIRS read pairs
# (reads.all.chunkNNN.R1/R2.fastq.gz in the current directory). The two splitters read named pipes and are waited
# for by PID, so the last chunk is complete when this returns.
#   split_fastq.sh BAM TMPDIR PAIRS
set -euo pipefail
BAM=$1 TMP=$2 PAIRS=$3
mkfifo "$TMP/r1" "$TMP/r2"
split -d -a 3 -l $((PAIRS * 4)) --filter='pigz -p 2 > $FILE.R1.fastq.gz' - reads.all.chunk < "$TMP/r1" & P1=$!
split -d -a 3 -l $((PAIRS * 4)) --filter='pigz -p 2 > $FILE.R2.fastq.gz' - reads.all.chunk < "$TMP/r2" & P2=$!
samtools collate -@ 4 -O -u "$BAM" "$TMP/c" \
  | samtools fastq -@ 2 -F 0x900 -n -0 /dev/null -s /dev/null -1 "$TMP/r1" -2 "$TMP/r2" -
wait "$P1"; wait "$P2"
rm -f "$TMP/r1" "$TMP/r2"
