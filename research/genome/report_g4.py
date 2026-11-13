"""Tables of a G4 run as Markdown.
    python report_g4.py results/g4-<host>-<date>.json > RESULTS_G4.md"""
import json
import os
import sys

path = sys.argv[1]
r = json.load(open(path))
ok = sum(c["ok"] for c in r["checks"])
out = ["# G4: genomes as matrix product states, latest run", "",
       f"From [`results/{os.path.basename(path)}`](results/{os.path.basename(path)}), written by `g4_run.sh`. The "
       "method is in [README.md](README.md). The encoding is Creevey et al.'s: |ψ⟩ = L^(-1/2) Σ_i |i⟩|b_i⟩.", "",
       f"**Checks: {ok} of {len(r['checks'])} hold.**", "",
       "## Small genomes in full: the whole MPS, its exact fidelity", "",
       "| sequence | bases | qubits | amplitudes | error ≤ 1e-2: largest bond, MPS parameters | error ≤ 1e-4 | "
       "error ≤ 1e-6 |", "| --- | --- | --- | --- | --- | --- | --- |"]
for x in r["small"]:
    cells = [f"{f['max_bond']}, {f['parameters']:,}" for f in x["full"]]
    out.append(f"| {x['sequence']} | {x['bases']:,} | {x['qubits']} | {x['amplitudes']:,} | " + " | ".join(cells) + " |")
out += ["", "## Large genomes at the middle cut: the Schmidt values needed", "",
        "*Most possible*: the number of blocks of the sequence on the left of the cut (each block one row of the "
        "state's matrix; the register's remaining rows are padding).", "",
        "| sequence | bases | qubits | most possible at the cut | needed for error ≤ 1e-2 | 1e-4 | 1e-6 | "
        "entanglement (bits) |", "| --- | --- | --- | --- | --- | --- | --- | --- |"]
for x in r["large"]:
    nd = x["needed"]
    block = 2 ** (x["qubits"] - 2 - x["cut_after_qubits"])
    most = min(x["max_possible"], -(-x["bases"] // block))
    out.append(f"| {x['sequence']} | {x['bases']:,} | {x['qubits']} | {most:,} | {nd['0.01']:,} | "
               f"{nd['0.0001']:,} | {nd['1e-06']:,} | {x['entropy_bits']:.2f} |")
out += ["", "## Every check", "", "| check | holds | detail |", "| --- | --- | --- |"]
out += [f"| {c['what']} | {'yes' if c['ok'] else '**NO**'} | {c['detail']} |" for c in r["checks"]]
print("\n".join(out))
