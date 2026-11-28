# First-divergence and locked-proposal BFGS replay

The historical commits e7031e9 and 2326430 differ in one source header. We traced both, without changing the search inputs, across the largest reciprocal improvement, loss and closest tie. The six resulting native PDBQT files match the historical archives byte-for-byte. All results here concern three purposefully selected development cells. Cartesian comparisons use ligand graph symmetry and the frozen pocket Cα alignment in the receptor frame. A model's output coordinates are rounded to millångström precision for geometry, and its full conformation is recorded at double round-trip precision. Crystal coordinates enter the analysis only.

## Historical failure retention

The native proposal counts total 7,130 across the six runs; each proposal was independently evaluated by both historical solvers. All 7,130 same-solver controls reproduce the original conformation, energy, model geometry, evaluation count and failure information exactly. Full-precision source states and mutated proposals agree at the first state split on all nine paired trajectories. Six split in the hunt optimizer; three split in the second, save-stage optimizer. First divergence occurs well before either run exhausts its budget. Step indices are zero-based.

| Query → receptor, seed | First state split in trajectories 1 / 2 / 3 | Published Top-1 old → repaired, Å | Published oracle old → repaired, Å |
|---|---|---:|---:|
| FABP4 2NNQ → 3FR2, 85237 | Hunt 6 / hunt 0 / hunt 0 | 2.814 → 15.954 | 2.814 → 12.919 |
| FKBP1A 1FKG → 1FKH, 85237 | Hunt 28 / hunt 0 / save 0 | 5.509 → 2.449 | 3.942 → 2.449 |
| EPAS1 5TBM → 4GHI, 85231 | Save 1 / save 8 / hunt 0 | 0.661 → 0.671 | 0.661 → 0.671 |

* **FABP4.** On trajectory 1, step 6, the same source and proposal return 10.295 Å (79 evaluations; accepted and sent to the save stage) and 15.308 Å (11 evaluations; rejected). In the old driver’s 606-proposal library, the best hunt is 2.768 Å under either solver. On the repaired driver’s 2,691 proposals it is 5.102 Å with the former solver and 10.835 Å with the repaired one. Even the cross-solve cannot recover the old library’s best. The repaired policy's per-trajectory best hunts on old-length prefixes are 12.227, 10.835 and 11.957 Å. Further proposals lower only the first to 11.633 Å. For the released Top-1 difference of 13.140 Å, published pool deterioration accounts arithmetically for 10.104 Å and additional regret for 3.035 Å. These are exact case-level endpoint differences; a percentage attributable solely to optimizer vs MC cannot be identified because a changed local return feeds subsequent proposals and their acceptance.
* **FKBP1A.** The repaired driver generates four hunts within 2 Å; three are accepted, all on trajectory 2 after the common 234-proposal prefix. The best saved state in that trajectory is 2.066 Å; the final trajectory pool's best is 2.120 Å. The published best at 2.449 Å comes from trajectory 2 step 255, at 2.479 Å after save. The old run's best hunt is 3.666 Å. The four near hunts occurred only on the longer path, and the saved-state threshold and final ranking further changed the observed endpoint.
* **EPAS1.** Although the released Top-1 changes by only 0.010 Å, the old driver visits ≤2 Å on two trajectories (259 of 465 hunts). The repaired driver reaches ≤2 Å only on trajectory 1 (202 of 1,582 hunts); its other two trajectories' best hunts are 9.715 and 12.579 Å. The longest consecutive near-native runs and each acceptance, save, insertion and eviction denominator are in the attrition table. The old EPAS1 run shows seven modes in its human-readable log but writes six under the 100 kcal/mol energy window; the seventh is excluded from published metrics.

The historical corrected solver on the adverse FABP4 library exhausts all 10 trials on 2,677 hunts, and stops on a non-descent direction on the other 14. On 1,282 of the 2,691 hunts even the last attempted angular increment exceeds π. Similar extreme angles also occur in the improvement and near-tie controls. That observation motivated the separate candidate experiment and cannot determine the sign of final pose accuracy.

## Conditional test of the principal-angle step

The addendum applies stable ed0d42d and candidate 0bb827e to the same 7,130 historical proposals. Both use the optimizer object configured by the production MC, including the candidate's new angle limit. They return identical input energies and initial directional derivatives on every proposal. Across all trials the number of ten-trial failures falls from 5,079 to 2,868, and the number whose last rejected angular displacement exceeds π falls from 2,210 to zero. All exhausted directions have negative directional derivatives and neither positive-curvature solver terminates on a non-descent direction. The last angular displacement of any exhausted candidate call is at most π/512. These observations establish fewer line-search failures from common starting proposals. The paths use an optimization objective and a finite search budget.

* On the 2,691 repaired-driver FABP4 proposals, the median paired grid energy change is −13.354 search-energy units. Best hunt improves from 6.824 to 3.917 Å, with no returned hunt within 2 Å. Evals rise from 32,653 to 72,732. On the 606 old-driver proposals the best hunt improves from 2.768 to 1.789 Å and the near count rises from zero to one.
* Across all six libraries, 4,559 proposals have lower energy by more than 1e−7 and 1,186 have higher energy. Median paired energy change is nonpositive in each library. Both geometric signs occur in each. The 7,130 hunts require 144,229 evaluations with stable and 196,297 with the candidate. A fixed total search budget can therefore change the number of proposals. Each library is conditional on a driver and all observations within a trajectory are dependent.

## Evidence and reproduction

* [Protocol](PROTOCOL.md) and [candidate addendum](CANDIDATE_PROTOCOL.md) were registered before their respective runs.
* [First divergences](results/first_divergence.csv), [all stage counts](results/path_attrition.csv), [released pose lineage](results/published_lineage.csv), [historical cross-solves](results/locked_proposals.csv) and [stable/candidate cross-solves](results/candidate_locked.csv) are reproducible from the raw event log.
* [Historical costs and hashes](raw/costs.json), [candidate replay costs and hashes](raw/candidate_costs.json) and [archive manifest](raw/files.json) preserve all inputs via the existing four-version panel and record per-cell output hashes. The six instrumented native runs took 2.548 CPU s (2.632 sequential wall s, peak RSS 47,460 KiB), twelve historical shadow runs 6.932 CPU s (7.112 wall s, 52,908 KiB), and twelve stable/candidate shadow runs 4.083 CPU s (4.220 wall s, 54,772 KiB). These times include grid initialization and diagnostic recording, with input preparation reused. Benchmark engine timing should use uninstrumented binaries.
* The raw ASCII event form is tab-delimited. T starts a CPU=1 trajectory. S gives the MC source and its energy. P gives the proposal (the energy slot is stale until the optimizer runs). H is the hunt return, A is an accepted Metropolis transition, D is the second optimizer return, and M enumerates the pool immediately after that save attempt. I/O link merged inputs to explicit-refinement outputs by index; F records the final mode-table ordering, and the actual PDBQT applies its energy window. Pose dumps include position, quaternion components and torsions. A detailed return also gives evaluation count, ten-trial exhaustion count, non-descent count, first exhausted unscaled max angular increment or −1, and chemical-heavy Cartesian coordinates in original PDBQT record order. In the candidate addendum the final four columns are starting energy, starting directional slope, terminating failure slope and actual last rejected angle. Missing failure data are NaN. The three-case probes use no hydrogens or receptor atoms in the RMSD. For historical experiments the angular scalar always starts at one. A reported old non-descent count is the number of early returns on that guard; the pre-fix default can attempt non-descent directions without this guard.
* Native run and shadow reproducibility requires the compiler, Boost setup and fixed PDBQT/metadata specified by the upstream four-version benchmark. To regenerate all data from the repository root, use a fresh output path for each command sequence.

```sh
D=benchmarks/bfgs_source_replay
bash "$D/build.sh" /tmp/bfgs-trace-build "$BOOST_PREFIX"
python "$D/run_replay.py" /tmp/bfgs-work /tmp/bfgs-trace-build/e7031e9/vina /tmp/bfgs-trace-build/2326430/vina
python "$D/analyze.py" /tmp/bfgs-work /tmp/bfgs-tables
python "$D/candidate_replay.py" /tmp/bfgs-work /tmp/bfgs-tables --binaries /tmp/bfgs-trace-build/ed0d42d/vina /tmp/bfgs-trace-build/0bb827e/vina
diff -qr "$D/results" /tmp/bfgs-tables
```

For an offline audit, restore the lossless bundle, which checks the SHA-256 of every shard, the joined stream, member files, independently executed historical same-solver returns, and every instrumented trace. Then recalculate both tables directly from the extracted Cartesian coordinates and energies.

```sh
D=benchmarks/bfgs_source_replay
python "$D/archive.py" restore "$D/raw" /tmp/bfgs-offline
python "$D/analyze.py" /tmp/bfgs-offline /tmp/bfgs-offline-tables
python "$D/candidate_replay.py" /tmp/bfgs-offline /tmp/bfgs-offline-tables
diff -qr "$D/results" /tmp/bfgs-offline-tables
```

An [independent verifier](verify.py) uses a NumPy SVD pocket fit, an independent event parser and direct PDBQT coordinates. [Verification results](verification.json) cover all nine first state divergences, all 18 native trajectories, 51 published ranks, 7,130 same-solver controls and 7,130 stable–candidate proposal pairs. The largest absolute discrepancy on checked aggregates is 4.45 × 10⁻¹⁵ or less. Run `python "$D/verify.py" /tmp/bfgs-offline "$D/results"` after restoring the bundle.
