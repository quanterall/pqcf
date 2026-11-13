# G4: genomes as matrix product states, latest run

From [`results/g4-dgx-spark-2026-11-13T1129Z.json`](results/g4-dgx-spark-2026-11-13T1129Z.json), written by `g4_run.sh`. The method is in [README.md](README.md). The encoding is Creevey et al.'s: |ψ⟩ = L^(-1/2) Σ_i |i⟩|b_i⟩.

**Checks: 2 of 2 hold.**

## Small genomes in full: the whole MPS, its exact fidelity

| sequence | bases | qubits | amplitudes | error ≤ 1e-2: largest bond, MPS parameters | error ≤ 1e-4 | error ≤ 1e-6 |
| --- | --- | --- | --- | --- | --- | --- |
| phiX174 | 5,386 | 15 | 32,768 | 96, 38,956 | 98, 39,930 | 98, 39,930 |
| random, 5386 bases | 5,386 | 15 | 32,768 | 97, 39,224 | 98, 39,930 | 98, 39,930 |
| periodic (a 100-base unit repeated), 5386 bases | 5,386 | 15 | 32,768 | 26, 7,822 | 26, 8,140 | 26, 8,220 |
| SARS-CoV-2 | 29,903 | 17 | 131,072 | 231, 211,440 | 234, 215,542 | 234, 215,578 |
| random, 29903 bases | 29,903 | 17 | 131,072 | 231, 211,440 | 234, 215,542 | 234, 215,578 |
| periodic (a 100-base unit repeated), 29903 bases | 29,903 | 17 | 131,072 | 26, 10,890 | 26, 11,424 | 26, 11,604 |

## Large genomes at the middle cut: the Schmidt values needed

*Most possible*: the number of blocks of the sequence on the left of the cut (each block one row of the state's matrix; the register's remaining rows are padding).

| sequence | bases | qubits | most possible at the cut | needed for error ≤ 1e-2 | 1e-4 | 1e-6 | entanglement (bits) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| E. coli K-12 | 4,641,652 | 25 | 2,267 | 2,118 | 2,265 | 2,267 | 8.94 |
| random, 4,641,652 bases | 4,641,652 | 25 | 2,267 | 2,123 | 2,266 | 2,267 | 8.97 |
| periodic (a 100-base unit repeated), 4,641,652 bases | 4,641,652 | 25 | 2,267 | 25 | 26 | 26 | 4.07 |
| human chromosome 21 (N removed) | 40,088,619 | 28 | 9,788 | 7,957 | 9,695 | 9,787 | 10.17 |

## Every check

| check | holds | detail |
| --- | --- | --- |
| phiX174: bond dimension 98 for an error below 1e-6, as Creevey et al. report | yes | max bond 98, infidelity 3.3e-15 |
| the block-overlap spectrum equals the SVD of the state (SARS-CoV-2, middle cut) | yes | largest difference 1.1e-16 |
