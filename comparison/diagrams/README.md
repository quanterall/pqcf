# Feynman diagrams and Feynman rules: PQCF against FORM, zero-dimensional field theory and Peskin and Schroeder

`pq_diagrams` generates every Feynman graph of a theory for given external legs and loops, with its weight
sign/S (S the symmetry factor, the sign (-1) per closed fermion loop), with QGRAF's filters (connected, one-particle
irreducible, no tadpoles, no snails). `pq_lagrangian` turns a Lagrangian into Feynman rules: the vertices in
momentum space, canonical up to the symmetries of the metric and the structure constants, and the propagators.

## The checks

| check | against | agreement |
| --- | --- | --- |
| QED photon and electron self-energies (1-3 loops), vertex (1-2), vacuum graphs (2-4); QCD gluon, ghost and quark self-energies (1-3 loops), three-gluon and quark-gluon vertices (1-2): 24 cases, up to 324 diagrams | FORM 5's diagram generator (`diagrams_`, Kaneko's grcc), `test/oracle/gt_diagrams.sh` | the same number of diagrams and the same multiset of weights sign/S in every case |
| phi^4 vacuum graphs to 6 loops, two- and four-point graphs (connected and one-particle irreducible); phi^3 vacuum graphs to 5 loops (eight vertices); a two-field theory | zero-dimensional field theory: ln Z, the cumulants of Z(J), Sigma = 1 - 1/G, G_c4/G^4, computed in the test in exact rationals (Cvitanovic, Lautrup and Pearson, [doi:10.1103/PhysRevD.18.1939](https://doi.org/10.1103/PhysRevD.18.1939)) | exact |
| QED vacuum graphs with N fermion flavours, up to 8 vertices, N = 1, 2, 3 | ln <(1 - e A)^N>, the zero-dimensional generating function | exact (the fermion-loop signs) |
| the orderly generator | Wick's pairings counted one by one (the first method, kept as a check) | graph by graph |
| phi^4, phi^3, Yukawa, QED, scalar QED, QCD (quark-gluon, three-gluon, four-gluon, ghost-gluon vertices) | Peskin and Schroeder, *An Introduction to Quantum Field Theory* (1995) | exact, term by term |
| scalar, Dirac and gauge-field (R_xi) propagators | -Gamma_2^-1 of their own quadratic Lagrangians (Gamma_2 from the vertex code), evaluated with the Minkowski metric and Dirac matrices at a momentum | exact |

Speed (one core): the three-loop gluon self-energy (324 diagrams) in 0.1 s; phi^3 vacuum graphs at five loops (71
graphs from eight vertices) in 1.8 s, at six loops (388 graphs) in 100 s.

## Reproduce

```sh
git clone https://github.com/quanterall/pqcf.git && cd pqcf
docker build -t pqcf-dev .
docker build -t pqcf-form -f comparison/form/Dockerfile .   # FORM 5, about 2 minutes
sh test/oracle/gt_diagrams.sh                                # FORM's diagrams, seconds; rewrites test/fixtures/gt_diagrams.json
git diff --stat test/fixtures/gt_diagrams.json               # no change expected
docker run --rm pqcf-dev rebar3 eunit --module=pq_diagrams_tests,pq_lagrangian_tests
```

Expected: the script lists the 24 cases with their numbers of diagrams (e.g. `qcd_glu_se 3 324`); the tests end with
`23 tests, 0 failures` (about 15 seconds).
