# Genome research, G0: data, instances and classical baselines

Step G0 of the plan in [docs/research/genome/README.md](../../docs/research/genome/README.md). It does four things:
- turns public genomic data into optimisation problems;
- finds their exact optima;
- measures the classical heuristics that any photonic method (G1: Gaussian boson sampling; G2: a coherent Ising
  machine) will have to beat;
- measures what translating each problem into the binary quadratic form (QUBO) costs, the form annealers, Ising
  machines and QAOA need.

That is the benchmarking Maurizio & Mazzola ask for ([arXiv:2507.04111](https://arxiv.org/abs/2507.04111)): exact
optima, tuned classical heuristics, and the optimal time to solution (TTS\*) as the problem grows.

**Latest run:** [RESULTS.md](RESULTS.md). **Every run:** [`results/`](results).

## Data

From [`data/`](data/SOURCES.md):
- **φX174 genome** (NCBI NC_001422.1);
- **STRING 12.5 network** of 45 human DNA-repair genes (CC BY 4.0).

## The three problem families

`instances.py` builds each problem in two forms: native and QUBO.

1. **Haplotype phasing.**
   - *Input:* reads from simulated diploid haplotypes over n heterozygous sites (coverage 10, 2–6 sites per read,
     5% sequencing errors).
   - *Model:* an Ising model, the Max-Cut relaxation of HapCUT
     ([doi:10.1093/bioinformatics/btn298](https://doi.org/10.1093/bioinformatics/btn298)).
   - Natively an Ising model, so its QUBO costs nothing extra.
   - Scored against the true haplotypes: switch errors, and the minimum-error-correction count.
2. **Genome assembly** (the layout step of overlap–layout–consensus).
   - *Input:* reads of 100 bp from a segment of φX174, coverage 4, shuffled.
   - *Model:* their overlap graph, and the Hamiltonian path of largest total overlap.
   - Native: CP-SAT's circuit constraint.
   - QUBO: Boev et al.'s two encodings ([arXiv:2004.06719](https://arxiv.org/abs/2004.06719)), positions (n²
     variables) and edges (one per overlap), each with a reward for overlap length.
   - Every optimum must be the reads' true order.
3. **Dense modules in a gene network.**
   - *Model:* the densest k-subgraph of the STRING network (links with score ≥ 0.9).
   - Native: CP-SAT with a cardinality constraint.
   - QUBO: a penalty for the size, which is the form G1's comparison needs.

## The solvers

All in `solvers.py`, on one core with one thread:
- **Exact:** OR-Tools CP-SAT, on the native form and, separately, on the QUBO. A QUBO's penalties must keep the
  native optimum.
- **Simulated annealing (SA):** single-spin Metropolis with dwave-neal's default temperature range.
- **Parallel tempering (PT):** 16 replicas.

## How TTS\* is measured

- **Run lengths:** each heuristic runs at 10 to 10,000 sweeps, 40 runs per length (SA) or 20 (PT).
- **Success probability p(t):** the fraction of runs at length t that reach the exact optimum (and, for assembly,
  decode to the true read order).
- **TTS:** TTS(t) = t_run · ln 0.01 / ln(1 − p), the expected time to succeed with 99% probability.
- **TTS\*:** the minimum of TTS(t) over the run lengths.

Two caveats:
- When the minimum sits at the shortest (10) or longest (10,000) run length, TTS\* is only an upper bound. That
  happens for 16-site phasing (p = 1 at 10 sweeps) and for 256-site phasing (p = 0.05 at 10,000 sweeps).
- With 20–40 runs per length, p is coarse.

## Results of the first run

12 November 2026, one Cortex-X925 core. **All 19 checks hold.**

**1. Haplotype phasing gets genuinely hard, and is the best target for the coherent Ising machine (G2).**

| sites | exact (CP-SAT) | SA TTS\* | PT TTS\* |
| --- | --- | --- | --- |
| 16 | 0.002 s | 0.00002 s | 0.00007 s |
| 64 | 0.054 s | 0.0025 s | 0.0022 s |
| 128 | 1.0 s | 0.096 s | 0.063 s |
| 256 | 44 s | 14 s (upper bound) | 54 s |

- Exact solving grows about 20,000-fold from 16 to 256 sites, and the heuristics grow similarly.
- The optimum phases the reads with 0–2 switch errors against the truth.
- This is the kind of problem the review names as a candidate: native Ising, with no translation overhead, hard at
  a few hundred variables.

**2. Assembly is easy natively and hard as a QUBO**, which is the review's warning, now measured on real data.

| reads | native exact | QUBO: CP-SAT | QUBO: SA / PT |
| --- | --- | --- | --- |
| 6 | about 1 ms | reaches the optimum | reach it |
| 8 | about 1 ms | misses it on the n² encoding (64 variables) within 60 s | reach it |
| 16 | about 1 ms | misses it on both encodings | n² encoding: SA 20 s, PT 17 s; edge encoding: only PT, 0.6 s |
| 24 | about 1 ms | misses it on both encodings | none reach it |

- The native exact solver finds the true read order at every size.
- At these sizes, assembly as a QUBO is a poor target for annealers or Ising machines.

**3. The densest modules of the DNA-repair network are cliques of hub genes.**
- Up to k = 10 the densest module is a complete clique (10 genes, all 45 links): BRCA1, BARD1, RAD51, MRE11,
  RAD50, ATM, ATR, CHEK2, TP53 and H2AX. At k = 12 it has 64 of 66 links, proven optimal in 250 s.
- **The optimum is not unique.** Several cliques of the same size exist, and two runs returned different ones for
  k = 3, 5 and 8.
- **It mixes complexes.** It is the hub of the DNA-damage response, not one known complex. Unweighted density
  favours hubs. To find complexes, G1 should also use weighted or degree-normalised density.
- **The QUBO's penalty landscape defeats SA and PT at k = 10.** Neither reaches the optimum, while the native exact
  solver needs 2.4 s. GBS samples dense subgraphs directly from the graph, with no penalty term, so G1 compares it
  with both forms.

## Running it

```sh
PQ_BENCH_HOST=dgx-spark sh research/genome/run.sh [core]   # about 45 minutes on one core
python research/genome/report.py research/genome/results/<run>.json > research/genome/RESULTS.md
```

## Next: G1

G1 adds the graph → GBS encoding to PQCF:
- the Takagi decomposition of the adjacency matrix, its mean photon number, and the interferometer;
- sampling with PQCF's threshold and photon-counting samplers;
- GBS-seeded search for dense modules and cliques on this network.

It is compared at equal time with the exact optimum, SA/PT on the QUBO, uniform random seeding, and the
quantum-inspired classical sampler of Oh et al. ([arXiv:2302.00536](https://arxiv.org/abs/2302.00536)).
