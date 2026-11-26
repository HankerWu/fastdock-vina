# Transfer of a different bound ligand's receptor contact footprint

The query starts and searches were already frozen. For each experimental query, the reference is the different ligand deposited in the receiving structure, in that receptor's coordinates. The query crystal pose is used solely to judge RMSD. These known-site, ligand-informed comparisons require an independent bound ligand in the same pocket. [Predefined score and controls](PROTOCOL.md) · [Marcou and Rognan's interaction fingerprint study](https://pubmed.ncbi.nlm.nih.gov/17238265/)

Top-1 is a whole-ligand receptor-frame symmetry-aware RMSD ≤ 2 Å. Each count is over one ligand/receptor direction at one seed. All scores within a row operate on identical coordinates. The element column is the existing bidirectional closest-same-element template with the same strength and identical reference data.

| Frozen source | Conditions | Native | Element | Contact | Attainable |
|---|---:|---:|---:|---:|---:|
| Five human targets, original fixed cap | 20 | 7 | 7 | 7 | 9 |
| Same five, cohort-calibrated pools | 20 | 10 | 11 | 12 | 14 |
| Same five, separately calibrated pools and new seeds | 20 | 11 | 13 | 13 | 14 |
| BPTF dry | 18 | 0 | 1 | 1 | 1 |
| BRDT dry | 6 | 1 | 2 | 3 | 4 |
| BCL6, 12k / 120k | 4 / 4 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| ESR1, PDE4D, PPARG, 3k | 18 | 5 | 6 | 6 | 7 |
| Same three, 12k | 12 | 5 | 5 | 5 | 6 |

The four extra hits across the two calibrated five-target rows originate from FKBP1A (two ligands and repeated seeds). BRDT covers two ligands and three seeds, BPTF six query/receiver directions. The 3k external-row rescue is PPARG. The nine rows reuse chemistry. They should not be combined to estimate population accuracy. A new template can change within-threshold and outside-threshold RMSDs in either direction even when a binary success count stays equal. The dry pools used here are each one emitted set, unlike the earlier two-run water unions.

From the repository root, run the following to regenerate the per-pose scores and all selections. Each input archive is already part of the repository. The separate C++/Python check compares the real CLI output for an existing element template and the new contact template, the unmodified source pool, the native score and energy window, and the independent numeric evaluation at serialized coordinates.

```bash
python benchmarks/frozen_core/contacts/evaluate.py /tmp/frozen-contact-rerun
python benchmarks/frozen_core/contacts/check_cli.py --help
```
