# Genome research: G0 (data, instances, classical baselines) and G1 (Gaussian boson sampling)

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

## G1: Gaussian boson sampling for dense modules and cliques

**Latest run:** [RESULTS_G1.md](RESULTS_G1.md), from `g1_run.sh` (PQCF sampling on 8 cores, then the analysis).

**The question** (from the plan): does sampling from a Gaussian boson sampler (GBS) find dense subgraphs and
cliques faster than classical samplers? Faster than uniform random sampling, and faster than the quantum-inspired
classical sampler of Oh et al. ([arXiv:2302.00536](https://arxiv.org/abs/2302.00536))? And does the gap grow with
size?

### What was built

- **A graph encoded as a GBS device** (`gbs_graph.py`):
  - the Takagi decomposition A = U diag(d) Uᵀ of the adjacency matrix;
  - squeezing tanh rᵢ = c dᵢ on every mode, then the interferometer U.
  - A pattern of k photons on k distinct modes then has probability ∝ |haf(A_S)|², the number of perfect
    matchings of that k-vertex subgraph, squared (Arrazola & Bromley, [arXiv:1803.10730](https://arxiv.org/abs/1803.10730)).
- **The device simulated exactly by PQCF** (`g1_pqcf.erl`): squeezing, the interferometer and photon counting
  through PQCF's Gaussian backend and chain-rule sampler. Samples are post-selected on exactly k photons, one per
  mode at most.
- **An exact low-photon regime.** PQCF caps the photon number step by step, which is exact only if almost no
  probability lies above the cap.
  - The post-selected distribution over k-subsets does not depend on the squeezing scale c.
  - So the device is simulated at the scale where the probability above the cap is below 10⁻⁵, computed exactly
    from the squeezing parameters. That means a mean photon number of 0.4–1.0, and 400–10,500 shots per kept sample.
  - At the scale an experiment would use (mean photon number k), the simulation never finished: a 45-mode run with
    up to 32 photons was stopped after 9 hours.
- **Oh et al.'s sampler**, exactly as the paper constructs it for a nonnegative matrix. It draws k/2 edges with
  probability ∝ weight and keeps the draw if they are disjoint, so a set appears with probability ∝ haf(A_S),
  not its square.
- **The same search for every sampler:**
  - random search for the densest k-subgraph;
  - shrink-and-search for the maximum clique (Banchi et al., [arXiv:1902.00462](https://arxiv.org/abs/1902.00462)).
- **Instances:**
  - a 30-vertex graph with a planted dense 10-vertex subgraph of low-degree vertices, like Arrazola & Bromley's;
  - the STRING DNA-repair network (45 genes).

### Checks: 28 of 28 hold

- **PQCF's GBS samples follow |haf(A_S)|² exactly** on a 10-vertex graph (chi-square p = 0.09 at k = 4, 0.25 at
  k = 6). The encoding, PQCF's squeezing and interferometer conventions, and the low-photon regime are all right.
- **Oh et al.'s sampler follows haf(A_S).**
- **Every optimum** is proven by CP-SAT; the maximum cliques come from exhaustive clique enumeration.
- **The probability above PQCF's cap** is below 10⁻⁴ in every chunk.

### Results

**1. Per sample, GBS is the best sampler for dense subgraphs, and its lead grows with k.** Samples needed to reach
the densest k-subgraph:

| graph | k | GBS | Oh et al. | uniform |
| --- | --- | --- | --- | --- |
| STRING | 4 | 2 | 5 | 62 |
| STRING | 6 | 20 | 80 | 2,458 |
| planted | 6 | 10 | 78 | 3,290 |
| planted | 8 | 50 | 1,040 | 221,086 |

- GBS needs 2 to 21 times fewer samples than Oh et al.'s sampler. That is in line with Oh et al.'s own estimate
  of "about 10 times more samples", and the gap widens with k.
- Uniform sampling, the baseline most GBS papers use, is 13 to 4,400 times worse than GBS. It overstates GBS's
  advantage.

**2. For maximum cliques there is no advantage on these graphs.** After shrink-and-search, every sampler (uniform
too) finds the maximum clique (size 9 planted, 11 STRING) in every block within 10–200 samples, mostly within 50. These graphs are too small and
easy to separate the samplers.

**3. In time on a classical computer, Oh et al.'s sampler wins by about 800–3,000 times.** A GBS sample costs PQCF
5 ms to 3 s of CPU to simulate, while Oh et al.'s costs 1–5 µs. For example, the planted 8-subgraph takes 8.6 s
of simulated GBS against 2.9 ms of Oh et al.'s sampler.

**4. The break-even for real hardware.** To beat Oh et al.'s sampler on one CPU core, a GBS device would have to
deliver 17,000–450,000 post-selected samples per second. Current GBS experiments (Jiuzhang, Borealis) have fixed
interferometers, so no programmable device that accepts an arbitrary graph exists today to try it.

**5. Biology.** The densest modules of the DNA-repair network are cliques of hub genes (BRCA1, BARD1, RAD51, the
MRE11–RAD50–NBN complex, ATM, ATR, H2AX, MSH6 at k = 10). Several equally dense modules exist, so the one reported
can change between runs (as in G0). Unweighted density favours hubs over individual complexes.

### Answer to G1's go/no-go question

- **Per sample: yes.** GBS beats both classical samplers, by a factor that grows with k.
- **At equal time on classical hardware: no, by 800–3,000 times.** A classical simulation of GBS (PQCF's or
  any exact one) is not a practical tool for this task.
- **The advantage exists only for a physical, programmable GBS device**, at the sample rates above.

**Our recommendation:**
- keep G1 as a hardware-readiness result, rerunnable whenever a programmable device appears;
- move on to G2, the coherent Ising machine on haplotype phasing. There a photonic heuristic competes natively on
  an Ising problem, and G0 showed the instances get genuinely hard.

### Running it

```sh
PQ_BENCH_HOST=dgx-spark sh research/genome/g1_run.sh ["15 16 17 18 5 6 7 8"]   # about 30 minutes on 8 cores
```

