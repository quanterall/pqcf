# Data sources

| file | what | source | licence |
| --- | --- | --- | --- |
| `phix174.fasta` | Bacteriophage φX174, complete genome, 5,386 bp | NCBI RefSeq [NC_001422.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1), fetched 12 Nov 2026 with E-utilities `efetch` | public domain (NCBI) |
| `NC_045512.2.fasta` | SARS-CoV-2 isolate Wuhan-Hu-1, complete genome, 29,903 bp | NCBI RefSeq [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2), fetched 13 Nov 2026 | public domain (NCBI) |
| `NC_000913.3.fasta` (not in git: `fetch.sh`) | Escherichia coli K-12 MG1655, complete genome, 4,641,652 bp | NCBI RefSeq [NC_000913.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3), fetched 13 Nov 2026, SHA-256 in `SHA256SUMS` | public domain (NCBI) |
| `NC_000021.9.fasta` (not in git: `fetch.sh`) | Homo sapiens chromosome 21, GRCh38.p14, 46,709,983 bp (40,088,619 known bases) | NCBI RefSeq [NC_000021.9](https://www.ncbi.nlm.nih.gov/nuccore/NC_000021.9), fetched 13 Nov 2026, SHA-256 in `SHA256SUMS` | public domain (NCBI) |
| `string_dna_repair.tsv` | Protein–protein associations among 45 human DNA-repair genes, combined score ≥ 0.7 | [STRING](https://string-db.org) 12.5, API `network`, species 9606, `required_score=700`, fetched 12 Nov 2026 | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) |

The 45 genes and the complexes they form (the known modules used as a reference for dense-subgraph searches):
- MRN complex: MRE11, RAD50, NBN.
- BRCA1–BARD1; BRCA2–PALB2; RAD51 and its paralogs RAD51C, XRCC2.
- Fanconi anaemia core complex: FANCA, FANCB, FANCC, FANCE, FANCF, FANCG, FANCL, FANCM; ID2 complex: FANCD2, FANCI.
- Mismatch repair: MLH1–PMS2 (MutLα), MSH2–MSH6 (MutSα), MSH2–MSH3 (MutSβ).
- Nucleotide excision repair: XPA, XPC–RAD23B, ERCC1–ERCC4 (XPF), ERCC2 (XPD), ERCC3 (XPB).
- Non-homologous end joining: XRCC5–XRCC6 (Ku), PRKDC (DNA-PKcs), LIG4–XRCC4, NHEJ1.
- Checkpoint signalling: ATM, ATR, CHEK1, CHEK2, TP53, MDC1, H2AX.

STRING: D. Szklarczyk et al., *The STRING database in 2025*, Nucleic Acids Res. 53, D730 (2025),
[doi:10.1093/nar/gkae1113](https://doi.org/10.1093/nar/gkae1113).
