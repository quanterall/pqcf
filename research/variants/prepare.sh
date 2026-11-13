#!/bin/sh
# One-off preparation of the benchmark (README.md), in research/variants/data (not in git):
#   the GRCh38 no-alt analysis set, checked by NCBI's MD5, indexed for samtools, minimap2 (short reads) and RTG;
#   the GIAB HG002 v4.2.1 truth set and its confident regions, cut to the chromosome;
#   the HG002 reads of the chromosome turned back into read pairs and split into numbered chunks, so the alignment
#   is measured too.
#   sh research/variants/prepare.sh [chrom] [chunks]
set -e
cd "$(dirname "$0")/data"
CHROM=${1:-chr20}
CHUNKS=${2:-8}
REF=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
grep "${REF}.gz" md5checksums.txt | sed 's| \./| |' | md5sum -c -
RUN="docker run --rm --user $(id -u):$(id -g) -v $PWD:/d -w /d pqcf-variants"
[ -s "$REF" ] || $RUN sh -c "pigz -dc $REF.gz > $REF"
[ -s "$REF.fai" ] || $RUN samtools faidx "$REF"
[ -s "$REF.mmi" ] || $RUN minimap2 -x sr -t 8 -d "$REF.mmi" "$REF"
[ -d "$REF.sdf" ] || $RUN rtg format -o "$REF.sdf" "$REF"
T=HG002_GRCh38_1_22_v4.2.1_benchmark
[ -s "truth.$CHROM.vcf.gz" ] || $RUN sh -c "bcftools view -r $CHROM -Oz -o truth.$CHROM.vcf.gz $T.vcf.gz && tabix -p vcf truth.$CHROM.vcf.gz"
[ -s "truth.$CHROM.bed" ] || awk -v c="$CHROM" '$1 == c' ${T}_noinconsistent.bed > "truth.$CHROM.bed"
BAM=HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.$CHROM.bam
if [ ! -s "reads.$CHROM.chunk$(printf %03d $((CHUNKS - 1))).R2.fastq" ]; then
    # read pairs, both mates present, no secondary or supplementary records
    $RUN sh -c "samtools collate -@ 8 -O -u $BAM /tmp/collate | samtools fastq -@ 4 -F 0x900 -1 /tmp/r1.fq -2 /tmp/r2.fq -0 /dev/null -s /dev/null -n - &&
        P=\$(( \$(wc -l < /tmp/r1.fq) / 4 )) && L=\$(( (P + $CHUNKS - 1) / $CHUNKS * 4 )) &&
        split -d -a 3 -l \$L /tmp/r1.fq reads.$CHROM.chunk && for f in reads.$CHROM.chunk[0-9][0-9][0-9]; do mv \$f \$f.R1.fastq; done &&
        split -d -a 3 -l \$L /tmp/r2.fq reads.$CHROM.chunk && for f in reads.$CHROM.chunk[0-9][0-9][0-9]; do mv \$f \$f.R2.fastq; done &&
        echo pairs \$P"
fi
ls -la
