#!/bin/sh
# The large public genomes used by G3, G4 and the worked example, fetched from NCBI (not kept in git) and checked.
#   sh research/genome/data/fetch.sh
set -e
cd "$(dirname "$0")"
for id in NC_000913.3 NC_000021.9; do
    [ -s "$id.fasta" ] || curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$id&rettype=fasta&retmode=text" -o "$id.fasta"
done
sha256sum -c SHA256SUMS
