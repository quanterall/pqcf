# Finding the variants in a human genome, sliced and distributed

Research branch `research/variants`, developed in parallel with the core platform. The goal is speed: find every
variant of a sequenced person against the reference genome, by cutting the work into numbered pieces, running them
in parallel on as many machines as are available, and gathering one deterministic result. The algorithms are
classical, the established ones; what this branch adds is the slicing, the distribution and the checks.

**Latest run:** [RESULTS.md](RESULTS.md). **Every run:** [`results/`](results).

## The pipeline

1. **Prepare once** (`prepare.sh`):
   - the reference, GRCh38's no-alt analysis set from NCBI, checked by NCBI's MD5 and indexed for samtools, minimap2
     and RTG;
   - the truth set, GIAB HG002 v4.2.1 with its confident regions (NIST);
   - the reads: the DeepVariant team's public HG002 NovaSeq 35× reads of one chromosome, turned back into raw read
     pairs so that the alignment is measured too, and split into numbered chunks.
2. **Scatter A, alignment.** Each chunk of read pairs is one task: minimap2 in short-read mode against the whole
   genome, then sorted. All chunks carry the same read group sample, so the calling treats them as one person.
3. **Scatter B, calling.** The chromosome is cut into numbered slices, 1 Mb each. Each slice is one task:
   - it reads all chunks' alignments in the slice plus a 1-kb margin on either side;
   - it calls variants with bcftools (`mpileup | call -mv | norm`);
   - it keeps the variants whose position lies in the slice's own range, so a variant in a margin is counted once,
     by the slice that owns it.
4. **Gather.** The slices' calls are concatenated in slice order into one VCF.

**One task interface, several executors** (`pq_genome.erl`, the pattern of PQCF's `pq_rexec`):
- A task is plain data, and its result is a pure function of it.
- `map(Tasks, Executor, Slots)` returns the results in task order, whichever worker ran what.
- Executors:
  - `local`: slots on this machine;
  - `{nodes, Nodes}`: one owner's Erlang nodes, called by `erpc` in a star, with no mesh between the workers;
  - the FVLLMP federation: the same tasks over its HTTPS gate (below).
- A task whose worker fails is retried elsewhere.

## The checks

Every run checks that:
- **8 workers and 1 worker give byte-identical calls,** so how the work is spread never changes the result;
- **the sliced calls equal one unsliced call over the whole chromosome,** so slicing with margins loses or doubles
  nothing;
- **the calls agree with the GIAB truth set:** precision, recall and F1 for SNPs and indels, by RTG vcfeval in the
  truth set's confident regions.

## Running it

```sh
docker build -t pqcf-dev . && docker build -t pqcf-variants research/variants/docker
sh research/variants/prepare.sh chr20 8          # downloads are listed in prepare.sh's header
PQ_BENCH_HOST=dgx-spark sh research/variants/run.sh chr20 5-9,15-19 10
```

## Results on chromosome 20

13 November 2026, [RESULTS.md](RESULTS.md). One DGX Spark, 10 cores. **All 4 checks hold.**

**Data:** GIAB HG002, the DeepVariant team's NovaSeq 35× reads of chromosome 20 (64 million bases), realigned from
scratch:
- 8 chunks of about 1.4 million read pairs;
- 65 slices of 1 Mb;
- 126,066 variant records.

**Time:**

| step | seconds |
| --- | --- |
| alignment, 8 chunks one after another (10 threads each) | 302 |
| calling, 65 slices, 8 workers | **29** |
| calling, 65 slices, 1 worker | 218 |
| calling, one unsliced run | 205 |
| gathering | 0.5 |

- Slicing makes the calling **7.5 times faster with 8 workers,** nearly linear.
- The slices are even: median 3.2 s, slowest 7.7 s (the deep regions around the centromere).
- Alignment is 91% of the total time, so it is what more machines must divide.

**Accuracy against the GIAB v4.2.1 truth set** (confident regions, RTG vcfeval):

| variants | precision | recall | F1 |
| --- | --- | --- | --- |
| SNPs | 99.14% | 99.28% | 99.21% |
| indels | 97.24% | 96.37% | 96.80% |

That is bcftools' accuracy, a fast but basic caller. Deep-learning callers such as DeepVariant are more accurate,
especially for indels. Any caller can be swapped in per slice, and the slicing spreads its extra cost.

**A finding: what "sliced = unsliced" means with bcftools.**
- Sliced and unsliced calls are **identical at every one of the 124,789 sites with depth ≤ 255.**
- They differ at 1,402 records, about 1%, and every one is at a site deeper than 255 reads.
- The cause is bcftools' error model. At such sites it keeps a random 255 of the bases, from one random stream for
  the whole run, so the draw depends on where the run started.
- None of these sites lies in the truth set's confident regions. They cluster around the centromere, where repeats
  pile up reads.
- Raising mpileup's depth caps does not change it: the subsampling is inside the error model.
- **The worker count never changes anything:** 8 workers and 1 give byte-identical calls, because every slice starts
  at the same place.

**What this projects to.** These are estimates, not measured; chromosome 20 is about 2.1% of the genome.

| setting | alignment | calling |
| --- | --- | --- |
| whole genome on this one machine | about 4 hours | about 25 minutes |
| whole genome on N such machines | divided by N | divided by N |

Data transfer is extra. A known saving: every alignment task now loads the 6 GB index again. A worker that keeps it
in memory across its chunks would avoid that.

## Towards the whole genome and the federation

- **The whole genome** is the same pipeline with more tasks:
  - about 3,100 one-megabase slices instead of 64;
  - the reads of a 30× genome, about 100 GB, in more chunks.
- **On many machines,** both scatters divide by the number of machines. What does not divide is moving data:
  - each alignment task needs the 6 GB minimap2 index and its reads;
  - each calling task needs the alignments of its slice.
- **On the FVLLMP federation,** where machines reach each other only through each node's HTTPS gate:
  - **Bring the work to the data.** Each site aligns and calls the genomes it holds, and only the result (a VCF of a
    few hundred megabytes) travels. That also keeps a person's reads at the site that holds them.
  - **When one genome must be spread over sites,** the alignment chunks travel once with their reads. The calling
    slices then need the alignments of their region, which is about 1/3,000 of the alignments per slice.
- **What it needs from FVLLMP:** a gate route that accepts a task (plain data, as here) and returns its result, plus
  the credit unit for a task. As for PQCF's other distributed work, that is for the FVLLMP team to provide; this
  branch does not change FVLLMP.
