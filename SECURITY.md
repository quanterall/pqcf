# Security policy

## Reporting a vulnerability

Please report vulnerabilities privately, through GitHub: the repository's **Security** tab, **Report a
vulnerability** (<https://github.com/quanterall/pqcf/security/advisories/new>). Do not open a public issue for them.

Please include what is affected (module, version or commit), how to reproduce it, and what an attacker gains. We
acknowledge a report within five working days, agree a disclosure date with you, and credit you in the advisory
unless you prefer otherwise.

## Supported versions

| version | security fixes |
| --- | --- |
| 1.4.x, and `main` | yes |
| earlier | no: please upgrade (the API is frozen and additive until 2.0, `docs/decisions/API-freeze.md`) |

## What is in scope

PQCF runs as an Erlang/OTP application with numerical kernels in Rust, loaded as NIFs. The parts that take input
from outside a trusted node are the most important:

- the JSON service (`pq_service`) and the Python client;
- the JSON codec (`pq_json`) and the expression parser (`pq_sym_parse`);
- the NIFs (`native/pq_nif`): a crash there takes down the Erlang node;
- distributed execution: the node handshake (`pq_node_registry`) and the task executors (`pq_rexec`, `pq_task`),
  where a node runs only the task modules it allows.

Distributed Erlang trusts every connected node completely: PQCF's nodes are meant for one owner's machines on a
trusted network. A threat model for federation peers (nodes of different owners) is to come, before federation is
released.

PQCF implements no cryptographic algorithms; it is not post-quantum cryptography.
