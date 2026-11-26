# Mixed Vina/Vinardo search (experimental, not validated)

Hypothesis: at fixed nominal MC evaluation budget, allocating even (1-based) trajectories to Vinardo improves basin coverage across held-out protein families. Odd trajectories retain the identical Vina RNG stream and objective. All saved poses are rescored and refined with the native Vina objective before selection. This is not a combined binding free-energy model; reported scores are Vina scores. No new fitted constants. Vinardo uses its standard pair weights and the same cutoff and macrocycle closure weight as the native search. The additional grid is built per docking; report runtime and realized evaluations rather than claiming exact wall-clock or evaluation parity.

```sh
python tools/build.py --output-dir build/mixed-vinardo --jobs 4
python benchmarks/mixed_vinardo/check.py build/mixed-vinardo/vina
build/mixed-vinardo/vina --receptor R.pdbqt --ligand L.pdbqt --config box.txt --pose_selection late --mixed_vinardo_search --exhaustiveness 8 --max_evals 15000 --num_modes 20 --refine_candidates 60 --energy_range 100 --seed 20261126 --search_counters --cpu 1 --out mixed.pdbqt
# matched control: omit --mixed_vinardo_search and change output path only
```

## External paired diagnostic (not a test of this implementation)
Source: https://github.com/Qasim-Hussain-Code/vina_gnina_pose_benchmark_pipeline/blob/313851a0941f70da32f190efb66790a0c8cbb74a/results/scored/run_scores.tsv . Pair by dataset, arm, case key and seed; retain two successful runs with recorded RMSD. Count top1_success_2a and top5_success_2a (physical validity plus symmetry-corrected RMSD <= 2 A). Top5 counts: crossdock 191 pairs, 40 both/9 Vina only/9 Vinardo only/133 neither; PoseBusters predicted box 289, 64/16/11/198; PoseBusters reference box 293, 171/31/20/71; Astex reference box 76, 52/5/9/10. The union uses two full searches and oracle reference labels; it is only evidence of complementarity. Cases are not independent protein families. No affinity conclusion follows.

## Preregistered acceptance
Use new families, clustered by protein family and scaffold. Freeze receptor/ligand preparation, non-query pocket definition, seeds (at least three), mode/refinement caps and nominal evaluation budget. Compare native Late, mixed, and separate higher-budget native control; report actual hunt/second/explicit evaluations, wall time, top1 and top5 symmetry-corrected heavy-atom RMSD without ligand superposition, physical validity (PoseBusters), and paired target-level bootstrap. Report all failures. Independently validate atom correspondence and output native scores; do not change defaults before actual held-out results. For affinity use only assay-compatible equilibrium Ki or Kd, never mix with IC50, keep receptor/conditions fixed, report per-target rank and pairwise sign by experimental uncertainty; no tuning on the evaluation targets. Pose improvement need not improve ranking.

References: https://doi.org/10.1002/jcc.21334 ; https://doi.org/10.1371/journal.pone.0155183 .

Status: source implementation and regression harness; not yet compiled or newly docked in the editing environment.
