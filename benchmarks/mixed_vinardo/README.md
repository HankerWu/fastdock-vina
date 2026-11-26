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

## Frozen first-pass external panels

PoseBusters: deterministic FNV-1a selection on case identifiers (seed prefix 20261126:, sorted ascending, one per PDB, first 16). Exclude 7WL4 (BRD4, historical development family) before running; no outcome-based replacement. The remaining 15 identifiers are embedded in paired_panel.py, including FABP3, JAK2 JH2, bacterial FabI, KDNase, PPARdelta, CMP-sialic acid transporter and N-myristoyltransferase. Some are metabolites/cofactors rather than druglike ligands; do not call these 15 independent protein families without sequence verification. The official data archive (MD5 pinned) supplies receptor and reference ligand. Generated ligand starts, 25 A reference-centered box, three seeds, four trajectories and cap 3000 per trajectory, 10 final modes, 40 refinement candidates; explicitly a low-budget known-site redock pilot. Compare native/mixed with identical preparation and settings. Strict receptor preparation can fail on cofactors/modified residues; record failures, never quietly drop them. All per-case logs, ligand poses, physical checks and independent no-alignment graph mapping RMSD checks saved in CI artifacts. Report per-pair actual gradient counts and wall time. Not a blind/cross-docking or fully held-out family benchmark.

Affinity: pinned OpenFF benchmark fd88824f9114244f95a14b485e6d6c96c1de716d. MCL1 25 Ki entries (Friberg et al., 10.1021/jm301448p, Table 3); PTP1B 22 Ki entries (Wilson et al., 10.1021/jm0702478, Tables 1-3); TYK2 13 Ki entries (10.1016/j.ejmech.2013.03.070 Table 4 and 10.1021/jm400266t Table 2). Exclude IC50 and other endpoint types; group by target/DOI/table/endpoint, not pooled. One fixed benchmark-aligned anchor per target defines the pocket; no query poses used as starting coordinates. Strict comparable-assay adjudication (pH, substrate, construct, conditions, censored measurements and errors) is not yet complete: these are preliminary rank diagnostics, not validated binding affinity claims. Require all planned pair failures reported. Repeated seed pairs do not constitute independent measurements. Receptor coordinate sources and all raw rows retained. Do not use benchmark-generated ligand conformers as crystal RMSD truth.

Workflow is scoped to this research branch. Results are experimental; do not merge/enable by default before review.
