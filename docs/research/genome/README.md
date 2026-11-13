# Can PQCF be used for genome research?

Research branch `research/genome`, developed in parallel with the core platform. This document is a summary of
five papers, of where a photonic simulator like PQCF fits, and of what we would build. Started 12 November 2026.

## The short answer

**Yes, but not the way most quantum-genomics papers propose.** Those papers mostly propose:
- Grover-style search for read alignment;
- quantum "images" of DNA;
- loading whole genomes into qubits.

None of these gives a speedup once the cost of loading the data is counted, and the field's own critical review
says so (Maurizio & Mazzola, [arXiv:2507.04111](https://arxiv.org/abs/2507.04111)). PQCF is a photonic
(continuous-variable) simulator, not a qubit machine, so they are not native to it either.

**Where PQCF does fit** are the problems that review names as the realistic candidates: hard combinatorial problems
on graphs and Ising models, of modest size, where the classical heuristics struggle. Photonics has two tools for
exactly these, and PQCF already simulates the physics behind both:

1. **Gaussian boson sampling (GBS) on graphs.** A graph's adjacency matrix can be programmed into a GBS device,
   which then preferentially samples dense subgraphs and cliques, and gives features for comparing graphs.
   - The genomic uses are dense modules in gene co-expression and protein-interaction networks, sets of unrelated
     individuals for association studies (GWAS), and comparing assembly or pangenome graphs.
   - PQCF has the GBS samplers, hafnians and torontonians. What's missing is the graph encoding and the search
     around the samples.
2. **Coherent Ising machines (CIM).** These are optical Ising solvers, and a simulated CIM is the
   "quantum-inspired" solver one of the papers used for genome assembly.
   - The genomic uses are haplotype phasing (a Max-Cut problem), overlap-graph assembly, pangenome-guided
     assembly and gene-panel selection.
   - PQCF's Gaussian backend with homodyne measurement and feed-forward is the physics of a measurement-feedback
     CIM. A fast mean-field CIM solver and the genomics-to-Ising builders would be new.

**The most valuable thing we can contribute may be the method, not a solver.** The review's central demand is
honest benchmarking: time-to-solution scaling against the best classical methods, on real instances. We now run
exactly that kind of head-to-head every week (`comparison/photonic`). A public, reproducible "quantum genomics
benchmark" built the same way would be new.

**What a simulator cannot do.** PQCF runs on classical hardware, so any speed it shows is classical speed. Its
roles here are to:
- test photonic proposals exactly at small scale, with cross-checked answers;
- run the quantum-inspired classical algorithms that come out of them (GBS-inspired samplers, simulated CIM);
- prepare problem encodings for photonic hardware (GBS devices, CIMs) when it is available.

## Status after G0–G2 (13 November 2026)

We measured the two photonic routes of this plan on genomic data, each against the best classical methods, with
every answer checked against a proven optimum. **Neither gives a practical advantage on these problems.**

1. **Gaussian boson sampling (G1).** Simulated exactly by PQCF, it needs 2–21 times fewer samples than the best
   classical sampler to find dense gene modules, a lead growing with size.
   - But simulating it costs 800–3,000 times more time than that classical sampler.
   - An advantage would need a programmable GBS device delivering 17,000–450,000 post-selected samples per second,
     and none exists.
2. **Coherent Ising machines (G2).** The CIM-derived solvers lose on every genomic problem we tried:
   - short-read phasing is solved exactly, in linear time, by dynamic programming;
   - long-read phasing is won by classical parallel tempering;
   - assembly is easier in its native form than as a QUBO.
3. **The data loading and problem translation that the critical review warns about are real.** Assembly as a QUBO
   turns a 1 ms exact problem into one no solver finishes at 24 reads. And the problems that look hard often have
   structure that classical exact methods exploit.

**What PQCF contributed:** an exact, cross-checked simulation of the photonic devices (GBS), and a benchmark
discipline that can be rerun when hardware appears. **What remains open:** G3 (the GBS graph kernel) and G4 (an MPS
engine). Given G1 and G2, we would only continue them with a genomics partner and a concrete dataset where
classical methods are known to struggle.

## The five papers

### 1. Maurizio & Mazzola, *Quantum computing for genomics: conceptual challenges and practical perspectives* (2025), [arXiv:2507.04111](https://arxiv.org/abs/2507.04111)

The critical review, and the frame for everything below.

- **Grover search for alignment has no practical speedup.**
  - Without quantum RAM, loading the genome into superposition costs at least N operations, which erases the
    √N advantage. Whether quantum RAM can ever be built is open.
  - Even with it, matching a 100-bp read against the human genome (3·10⁹ bp) takes classical tools under a
    minute. A quadratic speedup then needs each Grover step in under 1 ms. At a realistic 10 kHz logical gate
    rate, that is a circuit of depth 10.
- **Pairwise alignment** (Needleman–Wunsch, Smith–Waterman, O(NM), and BLAST, nearly linear) is not a quantum
  target.
- **Optimization is the realistic place to look**, with caveats:
  - Recasting a problem as a QUBO costs a polynomial overhead. Assembly as a Hamiltonian path needs N² variables,
    and embedding on sparse hardware adds another quadratic factor.
  - Strong classical heuristics exist, and some NP-hard problems have good polynomial-time approximations.
- **Candidates the review names:** Max-Cut (haplotype phasing), maximum independent set (unrelated individuals
  for GWAS, gene selection, non-repetitive promoters), knapsack (gene-panel design), the quadratic assignment
  problem, and low-autocorrelation binary sequences (guide-RNA design).
- **Quantum-enhanced Markov chain Monte Carlo** (Layden et al., [arXiv:2203.12497](https://arxiv.org/abs/2203.12497))
  is named for inference tasks such as phylogenetics and population genetics.
- **Method:** any claimed advantage must come from the scaling of the *optimal* time to solution, TTS\*(N),
  against the best classical solver. Earlier claims for annealers were refuted that way.
- **Quantum machine learning** has no shown advantage on real data, and suffers from trainability (barren
  plateaus) and data loading.

### 2. Boev et al., *Genome assembly using quantum and quantum-inspired annealing* (2020/21), [arXiv:2004.06719](https://arxiv.org/abs/2004.06719)

Overlap–layout–consensus (OLC) assembly, with the layout step as a Hamiltonian path written as a QUBO.

- **Two formulations:** N² variables in general; one variable per edge, M of them, for acyclic overlap graphs.
- **Two solvers:** a D-Wave 2000Q annealer, and **SimCIM**, a simulated coherent Ising machine (Tiunov, Ulanov,
  Lvovsky, [arXiv:1901.08927](https://arxiv.org/abs/1901.08927)), on a laptop GPU.
- **Results:**
  - φX174 (5,386 bp; 50 simulated error-free reads of 600 bp; a 248-vertex overlap graph): solved after
    splitting the graph by hand into 3 parts and using D-Wave's hybrid mode.
  - Synthetic sequences: D-Wave reconstructs up to 7 nucleotides (25 graph nodes); SimCIM goes further.
  - Time to solution for φX174: D-Wave about 8.5 ms of CPU plus 0.5 ms of QPU, SimCIM 0.26 ms (mean).
- **Our reading:** a proof of principle with no sequencing errors and no repeats. The simulated photonic machine
  did better than the quantum hardware. Its encodings and instances are a ready-made test set for PQCF.

### 3. Creevey, Jing, Hollenberg, *Implementation of a quantum sequence alignment algorithm* (2025), [arXiv:2506.22775](https://arxiv.org/abs/2506.22775)

Hollenberg's quantum sequence alignment from 2000 ([arXiv:quant-ph/0002076](https://arxiv.org/abs/quant-ph/0002076)),
implemented.

- **How it works:** a Grover-type search (BHHT, for an unknown number of solutions) over a superposition of all
  subsequences, for the smallest Hamming distance to the target. The data registers are prepared approximately
  by a genetic algorithm (GASP).
- **Scale:** simulated on 3–8-qubit data registers, with 2 to 32 sequences.
- **Result:** the accuracy stays around 80% even when the data is loaded with low fidelity, with large spread.
- **Our reading:**
  - It measures matches by Hamming distance only, with no gaps.
  - It does not address the cost of loading the data that paper 1 raises.
  - It is a qubit algorithm, not native to PQCF.

### 4. Creevey, Hassan, McCafferty, Hollenberg, Strelchuk, *Scalable quantum state preparation for encoding genomic data with matrix product states* (2025/26), [arXiv:2508.06184](https://arxiv.org/abs/2508.06184)

From the Wellcome Sanger Institute, Melbourne and Oxford.

- **Encoding:** a genome as Σᵢ |i⟩|bᵢ⟩, position and base. The state is written as a matrix product state (MPS),
  which is compiled into layers of nearest-neighbour two-qubit gates.
- **φX174:** fits in 15 qubits, and needs bond dimension 98 for an error below 10⁻⁶. Halving the bond dimension
  gives about 20% error.
- **Cost:** gate counts grow about as 2ⁿ, i.e. linearly in the genome's length. The authors note this means the
  genome has almost no redundancy to compress.
- **Comparison:** better than Qiskit's general state preparation.
- **What it points to next:** pangenome-graph "tangle" resolution with QAOA (Cudby et al.,
  [arXiv:2508.08200](https://arxiv.org/abs/2508.08200); Cudby & Strelchuk on closest-string problems,
  [arXiv:2510.15529](https://arxiv.org/abs/2510.15529)).
- **Our reading:**
  - The MPS method is solid, general tensor-network engineering.
  - The data-loading cost stays linear in the data, so it does not answer paper 1's objection.
  - MPS techniques are also used to simulate lossy photonic circuits. An MPS engine in PQCF would serve both.

### 5. Kösoglu-Kind et al., *A biological sequence comparison algorithm using quantum computers* (2023), [arXiv:2303.13608](https://arxiv.org/abs/2303.13608)

From IBM, CERN and FOM.

- **Method:** nucleotides as rotation angles in the FRQI quantum-image encoding. Two sequences are compared by
  measuring one "strip" qubit, with similarity = 1 − 2·P(1).
- **Scale:** run on a 5-qubit IBM machine with four-nucleotide sequences. An example gives similarity 0.246.
- **Our reading:** the paper claims lower time and memory than Needleman–Wunsch. It does not count the
  encoding's multi-controlled gates, whose cost grows with the sequence, nor the many shots a precise estimate
  needs. It is a toy demonstration, and not a direction for PQCF.

## What follows for PQCF: the ground rules

1. **No claims of speed from loading data into quantum states.** That holds for Grover, FRQI, quantum RAM and
   amplitude encoding (papers 1, 4).
2. **Target small but hard instances.** These are combinatorial problems where the best classical heuristics
   struggle at tens to hundreds of variables, so that scaling can be measured within reach.
3. **Compare against the best classical method, by TTS\*(N).** That means exact solvers where possible, plus
   simulated annealing, parallel tempering and problem-specific tools (HapCUT for phasing, real assemblers for
   assembly). Never only a uniform random baseline.
4. **Prefer native formulations.** Use graph problems that GBS encodes directly, and Ising problems without heavy
   QUBO overhead. Avoid N²-variable encodings where a sparse one exists.
5. **Cross-check every answer**, as in the head-to-head suite, and publish the instances, the code and the results.

## Where PQCF fits: tasks, methods, what exists, what to build

| genomics task | photonic method | PQCF today | to build | honest outlook |
| --- | --- | --- | --- | --- |
| Dense modules in gene co-expression or protein-interaction networks (densest k-subgraph, max-weight clique) | GBS seeded search (Arrazola & Bromley, [arXiv:1803.10730](https://arxiv.org/abs/1803.10730)) | GBS photon-counting and threshold samplers, hafnian, torontonian, loss | Graph → GBS encoding (Takagi decomposition, mean-photon tuning), sample → subgraph post-processing, seeded random search and simulated annealing | Oh et al. found GBS's advantage "not significant in general", with a classical quantum-inspired sampler as good ([arXiv:2302.00536](https://arxiv.org/abs/2302.00536)). Worth testing on real networks, against that sampler. |
| Unrelated individuals for GWAS; non-redundant gene sets (maximum independent set) | Max clique on the complement graph via GBS (as in docking, Banchi et al., [arXiv:1902.00462](https://arxiv.org/abs/1902.00462)) | as above | as above, plus independent-set ↔ clique conversion and real relatedness graphs | Same caveat. Relatedness graphs are sparse and large, so large instances will be hard. |
| Comparing assembly or pangenome graphs; classifying sequences by their graphs | GBS graph kernel (Schuld et al., [arXiv:1905.12646](https://arxiv.org/abs/1905.12646); isomorphism, Bradler et al., [arXiv:1810.10644](https://arxiv.org/abs/1810.10644)) | samplers; exact probabilities for small graphs | Feature vectors from samples or exact orbit probabilities; a kernel SVM pipeline | Plausible as a research result; the classical graph kernels it must beat are strong. |
| Haplotype phasing (Max-Cut on the fragment graph; HapCUT, [doi:10.1093/bioinformatics/btn298](https://doi.org/10.1093/bioinformatics/btn298)) | Coherent Ising machine (McMahon et al., [doi:10.1126/science.aah5178](https://doi.org/10.1126/science.aah5178); 100,000 spins, Honjo et al., [doi:10.1126/sciadv.abh0952](https://doi.org/10.1126/sciadv.abh0952)) | Gaussian backend, homodyne, feed-forward, loss (a small exact quantum CIM model) | Fast mean-field CIM solver in Rust (as SimCIM, [arXiv:1901.08927](https://arxiv.org/abs/1901.08927)); fragment-graph builder from reads; HapCUT-style baseline | Max-Cut is the standard Ising benchmark, and phasing instances are natural Ising problems. A useful test of CIM-style heuristics against HapCUT. |
| OLC genome assembly (Hamiltonian path) | CIM / QUBO (paper 2) | as above | Overlap-graph builder, both QUBO encodings, φX174 instances | Paper 2's instances are reproducible. The N² encoding is a known handicap (paper 1). |
| Pangenome-guided assembly; closest string | Binary optimisation (Cudby et al., [arXiv:2508.08200](https://arxiv.org/abs/2508.08200)) | — | Their formulation as Ising; solve with the CIM and classical solvers | Real data, and close to the Sanger group's work. A good external comparison. |
| Gene-panel design (knapsack), guide-RNA design (low-autocorrelation sequences) | CIM / QUBO | — | Builders | Speculative in the review itself; low priority. |
| Read alignment, sequence comparison (papers 3, 5) | Grover, FRQI (qubit) | — (qubit algorithms) | Only as small cross-checks, on Quanterall's qubit simulators or a future MPS engine | No practical advantage (paper 1). Not a goal. |
| Encoding genomes into quantum states (paper 4) | MPS → circuits | — | MPS engine (also useful for lossy photonic simulation); reproduce φX174 at bond dimension 98 | An engineering result. Its value for PQCF is the MPS engine itself. |
| Phylogenetics and population-genetics likelihoods | Exact symbolic computation (not quantum) | Symbolic engine, finite-field reconstruction | Likelihood of small trees as exact polynomials in the branch parameters; phylogenetic invariants | A PQCF-unique capability: exact answers to check approximate methods. Classical, exploratory. |

## Proposed plan for the branch

Each step ends with published instances, code and cross-checked results, like the comparisons. Each also has a
go/no-go question.

- **G0. Data and baselines.** **Done** (12 November 2026): [research/genome](../../../research/genome/README.md),
  results in [RESULTS.md](../../../research/genome/RESULTS.md).
  - *Haplotype phasing* looked hard at a few hundred sites (exact 44 s, SA 14 s at 256 sites); G2 corrected this:
    the 44 s were a single-worker exact solver, and short-read phasing is solved exactly in linear time.
  - *Assembly* is trivial natively (1 ms) but hard as a QUBO (no solver reaches the optimum at 24 reads).
  - *The densest modules* of the DNA-repair network are cliques of hub genes, and SA/PT on the QUBO fail at k = 10.

  What G0 was planned to contain:
  - **Public data:** φX174 (NCBI NC_001422.1), simulated reads with errors, and a small public gene co-expression
    or protein-interaction network.
  - **Builders:** graph and Ising/QUBO builders from these.
  - **Classical baselines:** an exact solver for small sizes, simulated annealing, parallel tempering, and Oh et
    al.'s quantum-inspired sampler.
  - **Measurement:** a TTS\*(N) harness on the pattern of `comparison/photonic`.
- **G1. GBS for dense subgraphs and cliques.** **Done** (13 November 2026):
  [research/genome](../../../research/genome/README.md#g1-gaussian-boson-sampling-for-dense-modules-and-cliques),
  results in [RESULTS_G1.md](../../../research/genome/RESULTS_G1.md).
  - Per sample, GBS (simulated exactly by PQCF) needs 2–21 times fewer samples than Oh et al.'s classical
    sampler to find the densest subgraph, a lead growing with k.
  - At equal time on a classical computer, Oh et al.'s sampler wins by about 800–3,000 times.
  - So an advantage exists only for a programmable GBS device (break-even 17,000–450,000 post-selected samples
    per second). For cliques, all samplers do equally well on these graphs.

  What G1 was planned to contain:
  - Graph → GBS encoding in PQCF, and seeded search.
  - Reproduce Arrazola & Bromley's 30-vertex result, then Oh et al.'s comparison.
  - Then apply it to gene-network modules.
  - *Go/no-go:* at equal time, does GBS seeding beat both the uniform and the quantum-inspired classical sampler,
    and does the gap grow with N?
- **G2. Coherent Ising machine for phasing and assembly.** **Done** (13 November 2026):
  [research/genome](../../../research/genome/README.md#g2-coherent-ising-machine-solvers-for-haplotype-phasing-and-assembly),
  results in [RESULTS_G2.md](../../../research/genome/RESULTS_G2.md). **No advantage on any instance class.**
  - Short-read phasing is solved exactly in linear time by a banded dynamic programme (0.011 s at 256 sites);
    G0's "hard" was a single-worker exact solver.
  - Long-read phasing is won by classical parallel tempering, about 4 times faster than the best CIM-derived
    solver (AHC) at 512 sites.
  - The assembly QUBO only makes a 1 ms native problem harder.
  - AHC is a strong Ising heuristic, ahead of SA and SimCIM; SimCIM, used in the assembly paper, is the weakest.

  What G2 was planned to contain:
  - A mean-field CIM solver at scale, and the exact Gaussian CIM model at small N to validate it.
  - Haplotype phasing from simulated reads against a HapCUT-style method.
  - Paper 2's φX174 and synthetic assembly instances.
  - *Go/no-go:* is TTS\*(N) better than tuned simulated annealing or parallel tempering on any instance class?
- **G3. The GBS graph kernel** for assembly and pangenome graphs, against standard graph kernels.
- **G4 (optional). An MPS engine.** Reproduce paper 4's φX174 encoding, and use the engine for lossy photonic
  simulation.

## Risks, stated plainly

- **The likely outcome is a negative or a narrow result.** Classical quantum-inspired algorithms match GBS on
  graph problems in general (Oh et al.). SimCIM, a classical simulation of a photonic machine, beat the hardware it
  imitates. PQCF's results would be honest measurements either way, and that is still worth publishing.
- **Genomics expertise.** Real instances (relatedness graphs, read fragments, pangenome graphs) need a genomics
  partner to choose meaningful data and judge results. The Sanger group behind papers 4 and the pangenome work is
  the natural contact.
- **Scale.** Exact GBS simulation is limited to tens of modes. That is enough for the method, but not for large
  graphs, where the quantum-inspired samplers would take over.

## References

- A. Maurizio, G. Mazzola, [arXiv:2507.04111](https://arxiv.org/abs/2507.04111)
- A.S. Boev et al., [arXiv:2004.06719](https://arxiv.org/abs/2004.06719)
- F.M. Creevey, M. Jing, L.C.L. Hollenberg, [arXiv:2506.22775](https://arxiv.org/abs/2506.22775)
- F.M. Creevey, H.T. Hassan, J. McCafferty, L.C.L. Hollenberg, S. Strelchuk, [arXiv:2508.06184](https://arxiv.org/abs/2508.06184)
- B. Kösoglu-Kind et al., [arXiv:2303.13608](https://arxiv.org/abs/2303.13608)
- L.C.L. Hollenberg, quantum sequence alignment, [arXiv:quant-ph/0002076](https://arxiv.org/abs/quant-ph/0002076)
- J.M. Arrazola, T.R. Bromley, GBS for dense subgraphs, [arXiv:1803.10730](https://arxiv.org/abs/1803.10730)
- M. Schuld, K. Brádler, R. Israel et al., GBS graph kernel, [arXiv:1905.12646](https://arxiv.org/abs/1905.12646)
- K. Brádler, S. Friedland, J. Izaac et al., graph isomorphism and GBS, [arXiv:1810.10644](https://arxiv.org/abs/1810.10644)
- L. Banchi, M. Fingerhuth, T. Babej et al., molecular docking with GBS, [arXiv:1902.00462](https://arxiv.org/abs/1902.00462)
- S. Jahangiri, J.M. Arrazola, N. Quesada et al., point processes with GBS, [arXiv:1906.11972](https://arxiv.org/abs/1906.11972)
- C. Oh, B. Fefferman et al., quantum-inspired classical algorithm for GBS graph problems, [arXiv:2302.00536](https://arxiv.org/abs/2302.00536)
- E.S. Tiunov, A.E. Ulanov, A.I. Lvovsky, simulating the coherent Ising machine (SimCIM), [arXiv:1901.08927](https://arxiv.org/abs/1901.08927)
- P.L. McMahon et al., a 100-spin coherent Ising machine, [Science 354, 614 (2016)](https://doi.org/10.1126/science.aah5178)
- T. Honjo et al., a 100,000-spin coherent Ising machine, [Sci. Adv. 7, eabh0952 (2021)](https://doi.org/10.1126/sciadv.abh0952)
- D. Layden, G. Mazzola et al., quantum-enhanced Markov chain Monte Carlo, [arXiv:2203.12497](https://arxiv.org/abs/2203.12497)
- J. Cudby et al., pangenome-guided assembly via binary optimisation, [arXiv:2508.08200](https://arxiv.org/abs/2508.08200)
- J. Cudby, S. Strelchuk, closest string problems, [arXiv:2510.15529](https://arxiv.org/abs/2510.15529)
- V. Bansal, V. Bafna, HapCUT, [Bioinformatics 24, i153 (2008)](https://doi.org/10.1093/bioinformatics/btn298)
