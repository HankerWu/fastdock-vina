# Archived search-budget and score-aggregation replay (2026-11-26)

This is a retrospective analysis of serialized scores and RMSDs, **not** a new
docking run or a wall-time-equivalent experiment. Source snapshot: 301db5c.
Run `python3 benchmarks/search_budget/replay.py` from the repository checkout.
RMSD threshold is 2 Å; top-1 is the lowest reported score across included runs
(ties by seed). The oracle is any emitted mode, not every generated pose.

Each of three held-out-seed folds starts with the other two seeds. The adaptive
policy spends a third run on floor(N/2) directions in each panel, chosen by the
largest absolute difference between the two reported best scores, ties by
target/query/receiver. The random comparator selects the same number uniformly;
its expectation is computed analytically. No RMSD is supplied to the policy.

| Panel | Targets | Directions | One top1/any | Two top1/any | Three top1/any | Adaptive top1/any | Uniform expected top1/any |
|---|---:|---:|---:|---:|---:|---:|---:|
| cubic human | 6 | 11 | 14/17 | 14/22 | 5/8 | 14/24 | 14.4545/22.9091 |
| torsion bridge | 4 | 8 | 2/4 | 3/7 | 1/3 | 3/9 | 3/8 |

One/two/adaptive/random denominators are 33 and 24; three denominators are
11 and 8. Three-seed aggregation yields 6/19 correct top1 and 11/19 oracle.
The adaptive policy improves oracle coverage but not top-1 relative to random;
it is not supported as a default policy. Cases/seeds are dependent; no
population-level significance claim is made.

For the affinity experiment, use the earlier interpolation-only archived
`benchmarks/cubic_grid/affinity/cells.csv` (not the subsequent wall-fix dataset).
Davis et al. same-assay Kd labels in panel.csv; one structure and two seeds per
compound. Pair credit is one for a correct order, one-half for a tie, for Kd
ratios >= 3. See https://doi.org/10.1038/nbt.1990.

| Target | Pairs | Native min/mean/max | Cubic min/mean/max | Size baseline |
|---|---:|---:|---:|---:|
| SYK | 68 | 52/43/38 | 52/47/47 | 39.5 |
| FAK | 47 | 9/9/8 | 18/19/21 | 37.5 |

The minimum over random optimizer starts is not a thermodynamic ensemble
free energy. FAK remains worse than size-only. Neither policy is promoted.
New independent targets and endpoint-qualified assays are necessary before
claiming improved generalization.
