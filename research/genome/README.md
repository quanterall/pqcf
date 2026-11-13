# Genome research: G0 (baselines), G1 (Gaussian boson sampling), G2 (coherent Ising machines)

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

**1. Haplotype phasing gets slow for single-worker CP-SAT** (*corrected by G2*: it is not hard. CP-SAT with 8 workers
proves 256 sites in 0.09 s, and a banded dynamic programme in 0.011 s; see G2).

| sites | exact (CP-SAT) | SA TTS\* | PT TTS\* |
| --- | --- | --- | --- |
| 16 | 0.002 s | 0.00002 s | 0.00007 s |
| 64 | 0.054 s | 0.0025 s | 0.0022 s |
| 128 | 1.0 s | 0.096 s | 0.063 s |
| 256 | 44 s | 14 s (upper bound) | 54 s |

- Exact solving on one CP-SAT worker grows about 20,000-fold from 16 to 256 sites, and the heuristics grow
  similarly.
- The optimum phases the reads with 0–2 switch errors against the truth.
- G0 read this as the kind of problem the review names as a candidate. G2 showed it is not: the growth was single-
  worker CP-SAT's, and short-read phasing is solved exactly in linear time.

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


## G2: coherent-Ising-machine solvers for haplotype phasing and assembly

**Latest run:** [RESULTS_G2.md](RESULTS_G2.md), from `g2_run.sh`: exact optima on 8 cores, then the four solvers in
parallel, one core each.

**The question:** do solvers derived from the coherent Ising machine (CIM), an optical Ising solver, find the
optimum of genomic Ising problems faster than tuned classical heuristics, and does the gap grow with size?

### What was built

- **Two CIM solvers** (`cim.py`, Numba, one core, sparse couplings):
  - **SimCIM**, as Tiunov, Ulanov & Lvovsky ([arXiv:1901.08927](https://arxiv.org/abs/1901.08927)) define it. It is
    the simulated CIM Boev et al. used for genome assembly.
  - **AHC**, the CIM with amplitude-heterogeneity correction of Leleu, Yamamoto & McMahon
    ([arXiv:1810.12565](https://arxiv.org/abs/1810.12565), Eqs. 1–2). Its error variables destabilise local
    minima. We keep its target amplitude constant; the paper also modulates it.
- **The classical heuristics from G0:** simulated annealing (SA) and parallel tempering (PT).
- **The same tuning budget for every solver:** 8 parameter settings, chosen by geometric-mean TTS\* on separate
  tuning instances for each problem family.
- **Exact optima:**
  - CP-SAT with 8 workers, run before the solvers, alone;
  - for short reads also a **banded dynamic programme** (`solvers.exact_banded`), whose cost is linear in the
    number of sites.
- **Three problem families:**
  1. haplotype phasing from **short reads** (2–6 sites per read), 16–256 sites;
  2. haplotype phasing from **long reads** (10–40 sites per read), 64–512 sites;
  3. Boev et al.'s assembly QUBO for φX174 (one variable per overlap edge), 6–16 reads.

  The phasing families have 3 instances per size, each with a proven optimum.

### Checks: 46 of 46 hold

Every optimum is proven by CP-SAT, and the dynamic programme gives the same optimum on every short-read instance.

### Two corrections to G0, found on the way

- **G0's exact times were CP-SAT on one worker.** On one worker, CP-SAT spent 44 s on a 256-site instance it proves
  in 0.09 s with 8. On another it stopped after 300 s at a solution far from the optimum (−3,272 against −3,486).
  An earlier G2 run had used that unproven value as its target; the check caught it, and the exact phase now uses
  8 workers.
- **Short-read phasing is not a hard problem.** A read couples only the sites it covers, so the couplings form a
  band 5 sites wide.
  - A dynamic programme over that band solves it exactly in linear time: 0.011 s at 256 sites, 0.18 s at 4,096.
  - The standard phasing tools rely on the same kind of structure. WhatsHap solves the weighted
    minimum-error-correction problem exactly, with a cost exponential only in the read coverage (Patterson et al.,
    [J. Comput. Biol. 22, 498 (2015)](https://doi.org/10.1089/cmb.2014.0157)).
  - G0's statement that phasing "gets genuinely hard at a few hundred sites" was an artefact of a single-worker
    exact solver.

### Results

**1. Short reads: every heuristic loses to the exact dynamic programme.**

| method | time at 256 sites |
| --- | --- |
| dynamic programme, 1 core, exact | 0.011 s |
| CP-SAT, 8 workers, exact | 0.09 s |
| AHC (CIM), best heuristic | 0.07–0.32 s (TTS\*) |
| SA | about 6.6 s on two instances, not reached on the third |
| PT and SimCIM | not reached |

Among the heuristics, AHC scales best: its TTS\* grows about 14 times per doubling of the sites, against 88 for SA.
On this family that only ranks the heuristics, since none of them is needed.

**2. Long reads, where exact solving is costly, are won by a classical heuristic.** Median TTS\* at 512 sites, and its
growth per doubling over 128–512 sites:

| method | time at 512 sites | growth per doubling |
| --- | --- | --- |
| PT | 0.0085 s | 5.3 |
| AHC (CIM) | 0.037 s | 9.9 |
| SA | 0.072 s | 32 |
| SimCIM | 0.46 s | 14 |
| CP-SAT, 8 workers, exact | 23.5 s | 9.3 |

- Every heuristic is far faster than proving the optimum, but the best is parallel tempering. It is about 4 times
  faster than AHC at 512 sites and has the flattest growth.
- Many PT times sit at the shortest run length tried, so they are upper bounds and PT is, if anything, faster still.
- SimCIM, the solver used in the assembly paper, is the slowest heuristic here.

**3. Assembly QUBO: no consistent winner, and none is needed.**
- AHC is fastest at 6–12 reads.
- At 16 reads only SA (0.033 s) and PT (0.0035 s) reach the optimum; neither CIM solver does.
- The native exact solver needs under 1 ms at every size (G0).

### Answer to G2's go/no-go question

**No.** On no genomic instance class here does a CIM-derived solver beat the best classical method:
- **Short-read phasing:** solved exactly in linear time by dynamic programming.
- **Long-read phasing:** won by parallel tempering.
- **Assembly:** a QUBO only makes a native one-millisecond problem harder.

Of the two, AHC is a strong Ising heuristic: it beats SA and SimCIM on both phasing families. SimCIM is not.

**Limits of this test, stated plainly:**
- AHC runs without its target-amplitude modulation.
- Tuning used small instances and only 8 settings per solver.
- Everything ran on one CPU core; a physical CIM's speed is not modelled.
- Phasing was solved in its Ising (Max-Cut) form, not as the minimum-error-correction problem that phasing tools
  actually solve.

None of these would change the main conclusions:
- the short-read problem is easy;
- the best solver for the long-read Ising problem is classical.

### Running it

```sh
PQ_BENCH_HOST=dgx-spark sh research/genome/g2_run.sh ["15 16 17 18"]   # about 40 minutes
```
