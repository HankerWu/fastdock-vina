# Common-field scoring of equal-budget dry and hydrated pose pools

The confirmation runs already contain two rigid receptor hypotheses for human BPTF and BRDT BD1. We asked whether using one search in each field would preserve near-crystal poses under an equal two-search budget, and whether the common energy function could select those poses. All calculations use existing searches at 120,000 evaluations per trajectory with three trajectories. In total, 72 source searches emitted 637 candidates. We evaluated every coordinate under both a common dry and a common hydrated receptor, giving 1,274 new fixed-coordinate score calculations.

For seeds 95701, 95713, 95717, the three ordered folds are (1,2), (2,3), (3,1) in that sequence. Compare dry(i)+dry(j), dry(i)+polar(j), and polar(i)+polar(j). Each candidate pool can emit nine poses. We preserve all emitted candidates, up to 18 per fold. The RMSD oracle uses the minimum across both source pools. For rank one, we rerank the whole union using the original prepared ligand topology, identical explicit receptor and zero unbound reference (within query), the original 24 Å search box and the emitted coordinates. There is no reoptimization. Exact ties at the reported 0.001 kcal/mol precision use first source, then original rank. None of the tied best scores spans the 2 Å hit threshold. RMSDs come from the hash-matched source outputs and the independently audited protein-frame alignments and equivalent-atom mappings.

| Target and task | Two searches | Folds | Oracle hits | Dry common score hits | Polar common score hits |
|---|---|---:|---:|---:|---:|
| BPTF cross | dry + dry | 18 | 2 | 0 | 0 |
| BPTF cross | dry + polar | 18 | 3 | 0 | 0 |
| BPTF cross | polar + polar | 18 | 6 | 0 | 0 |
| BRDT cross | dry + dry | 6 | 5 | 1 | 5 |
| BRDT cross | dry + polar | 6 | 5 | 2 | 5 |
| BRDT cross | polar + polar | 6 | 6 | 3 | 6 |
| BPTF cognate diagnostic | dry + dry | 6 | 5 | 3 | 0 |
| BPTF cognate diagnostic | dry + polar | 6 | 6 | 3 | 0 |
| BPTF cognate diagnostic | polar + polar | 6 | 3 | 0 | 0 |
| BRDT cognate diagnostic | dry + dry | 6 | 3 | 0 | 3 |
| BRDT cognate diagnostic | dry + polar | 6 | 3 | 0 | 2 |
| BRDT cognate diagnostic | polar + polar | 6 | 3 | 0 | 3 |

On BRDT, ranking the dry-search pool with explicit polar sites recovers five of six, at fixed candidates, compared with one under dry scoring. Switching search fields adds at most one more top hit under polar scoring. On BPTF, hydrated search puts near-crystal candidates into six of 18 all-hydrated cross folds, but neither common scorer selects one. Adding a hydrated search to a dry search yields three oracle hits and zero top hits. Overall, mixed search provides no consistent cross-target rank-one advantage. Results are conditional on these fixed protonation states and recurrent, nondisplaceable waters. The shared-seed cyclic folds reuse searches. The 18 BPTF cross folds correspond to six query/receiver directions, the six BRDT folds to two directions. They are not statistically independent experiments.

Actual mean global evaluations per cross fold (dry+dry, dry+polar, polar+polar) were 720149, 720143 and 720137 for BPTF; 720193, 720201 and 720208 for BRDT. The sum of the archived source run times per fold averaged 2.04, 2.12 and 2.20 s on BPTF and 3.52, 3.26 and 3.01 s on BRDT. Timings include source map generation and refinement on the recorded host and exclude this offline rescoring. Additional explicit-refinement counts and candidate counts are in summary.csv.

## Replay

The script extracts the existing bptf_water, pose and water_transfer archives directly. A fresh Python 3.12 environment uses only the standard library, plus the branch's Vina executable. In energies.csv the common score is reported with the CLI's three decimal places. It is comparable across poses for the same ligand and field. The original reported score can use a run-specific unbound reference.

```bash
B=benchmarks/geometry_comparison/pooled_water
W=$(mktemp -d)
python "$B/compare.py" /absolute/path/to/vina "$W" --jobs 4
python "$B/audit.py" "$W"
cmp "$B/energies.csv" "$W/energies.csv"
cmp "$B/folds.csv" "$W/folds.csv"
cmp "$B/summary.csv" "$W/summary.csv"
```
