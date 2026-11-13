# Worked example: a virus, a bacterium and a human chromosome

This example runs the genome research (G0–G4, [../README.md](../README.md)) on public genomes, end to end. It shows
what each method does with real data, and where quantum and photonic methods stand against classical ones.
**Results:** [EXAMPLE.md](EXAMPLE.md), written by `run.sh`.

```sh
sh research/genome/data/fetch.sh            # E. coli and human chromosome 21 from NCBI, checked by SHA-256
sh research/genome/example/run.sh           # about 10 minutes on 8 cores
```

## The genomes

All from NCBI RefSeq (public domain); see [../data/SOURCES.md](../data/SOURCES.md).

| genome | accession | size | used for |
| --- | --- | --- | --- |
| SARS-CoV-2, isolate Wuhan-Hu-1 (virus) | NC_045512.2 | 29,903 bases | assembly, quantum encoding |
| *Escherichia coli* K-12 MG1655 (bacterium) | NC_000913.3 | 4.6 million bases | assembly of a region, its protein network, quantum encoding |
| *Homo sapiens* chromosome 21, GRCh38 | NC_000021.9 | 46.7 million bases (40.1 million known) | phasing, quantum encoding |

E. coli's protein network is from STRING 12.5 (CC BY 4.0): 45 DNA replication and repair proteins.

## The steps

### 1. Assembling a genome from reads

Sequencing machines read short pieces of a genome, and assembly puts them back in order. Here:
- **The reads:** simulated, error-free reads of 300 bases at 6-fold coverage, from the whole SARS-CoV-2 genome and
  from a 50,000-base region of E. coli (the lac operon).
- **The overlap graph:** found through an index of read beginnings.
- **The layout:** the order of reads that maximises the total overlap (a Hamiltonian path), found by an exact solver
  on that graph, then joined into one sequence.

**What it shows:**
- The exact solver finds the true order in a fraction of a second, and the assembled sequence equals the reference.
- The same problem written for a quantum annealer would need hundreds of thousands to a million binary variables
  (the n² encoding of Boev et al., [arXiv:2004.06719](https://arxiv.org/abs/2004.06719)). G0 found no solver that
  finishes that form beyond 24 reads.
- Real assemblers also handle sequencing errors and repeats; error-free reads keep the example about the layout
  step itself.

### 2. Phasing a human chromosome region

Each person has two copies of each chromosome. Phasing decides which variants sit together on the same copy. Here:
- **The individual:** simulated on a real 1-million-base region of chromosome 21, with a heterozygous variant about
  every 1,500 bases (typical for humans) in a random phase.
- **The reads:** paired short reads (2 × 150 bases) and long reads (15,000 bases, 1% errors), each from one copy.

**What it shows:**
- **Short reads cannot phase a chromosome:** they rarely cover two variants, so the region falls apart into many
  small phase blocks.
- **Long reads join the whole region into one block.** The resulting Ising problem is solved exactly in well under a
  second, with no switch errors.
- The heuristics tuned in G2, parallel tempering and the CIM-derived AHC, do not reach that optimum here.

### 3. Dense modules in E. coli's DNA replication and repair network

The densest k-gene modules of the protein network, from three samplers each followed by the same search:
- a Gaussian boson sampler (GBS), simulated exactly by PQCF;
- Oh et al.'s quantum-inspired classical sampler ([arXiv:2302.00536](https://arxiv.org/abs/2302.00536));
- uniform random sampling.

Each found module is checked against the exact optimum (G1's method).

**What it shows:**
- **The modules are biologically meaningful:**
  - at k = 6, the replication machinery (DnaA, DnaC, DnaG, DnaE, DnaN, DnaX);
  - at k = 8, the RecF/RecO/RecR–RuvAB recombination-repair pathway with RecA, RecQ and RecJ.
- **Per sample, GBS finds them most reliably,** ahead of Oh et al.'s sampler and far ahead of uniform sampling.
- **But each simulated GBS sample costs PQCF milliseconds to a quarter of a second of CPU,** against microseconds for
  Oh et al.'s sampler (G1). An advantage would need a programmable GBS device.

### 4. The genomes as quantum states

From G4: the virus, the bacterium and the human chromosome encoded as quantum states (Creevey et al.,
[arXiv:2508.06184](https://arxiv.org/abs/2508.06184)). The table gives the bond dimension a matrix product state
needs, and with it the size of the circuit that loads the genome.

**What it shows:** every genome needs exactly what a random sequence of the same length needs. Loading a genome into
a quantum computer costs as much as the genome is long, which is the data-loading objection of Maurizio & Mazzola
([arXiv:2507.04111](https://arxiv.org/abs/2507.04111)).

## In one paragraph

On real genomes, the classical exact methods solve these problems in a fraction of a second:
- assembly of a virus and a bacterial region;
- phasing of a human chromosome region from long reads;
- the densest modules of a bacterial protein network.

The quantum and photonic routes either need hardware that does not exist (a programmable boson sampler fast
enough) or lose to classical methods. Underneath every quantum route sits the cost of loading the genome itself,
which this example measures to be as large as the genome. What PQCF provides here is an exact simulation of the
photonic devices and measurements that can be checked. Today's answer to "can PQCF be used for genome research?"
is therefore: as a precise testbed for photonic ideas, not as a faster genome tool.
