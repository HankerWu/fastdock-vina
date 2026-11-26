# Alternating receptor repulsion on human thrombin

We applied the same `--soft_receptor_search` tested on ESR1, PPARG and PDE4D to the twelve compounds and three structures in the parent benchmark. Only even-numbered search trajectories use half-strength receptor repulsion. All final rankings in the native-score column use full-strength Late Vina. All first-panel hard-search results are reused from the parent fixtures. The soft searches and confirmation searches add 216 fresh docking runs. The two seed panels together contain 288 independent process runs (12 compounds × three receivers × four seeds × two arms), and repeat the same twelve experimental ligand identities. Searches use the same 50000-evaluation-per-trajectory setting, four trajectories, full prepared ligand/receptor files and boxes as the parent benchmark. The protein pocket does not move in these searches. Build SHA-256 for the modified binary used here is c26d269353976b0cc1909d7f9398363ec26e6c83750b7cacd1207ddfbde167a4.

Two additional seeds (88091, 302123) were frozen after the first two seeds (83519, 194027) showed better pose and affinity outcomes. The two arms of the confirmation used the same modified binary, same seeds and same search settings, with the soft-receptor flag as the sole intervention. The confirmation searches also recorded trajectory and explicit-refinement counts. Repeated seeds sample uncertainty in the search of the same compounds. A single thrombin scaffold provides no independent new-scaffold validation.

| Native final ranking, hard → soft search | First seeds | Confirmation seeds | All four seeds |
|---|---:|---:|---:|
| Top-1 ≤2 Å, crystal-labeled complexes | 14→22 / 30 | 11→17 / 30 | 25→39 / 60 |
| Top-1 ≤2 Å, different receptor/ligand deposits | 11→19 / 26 | 9→15 / 26 | 20→34 / 52 |
| Best returned ≤2 Å, different deposits | 16→24 / 26 | 15→20 / 26 | 31→44 / 52 |
| Correct unique ≥threefold ITC Kd pairs | 11→19 / 37 | 18→16 / 37 | 17→17 / 37 |
| Same-charge subset of unique Kd pairs | 1→14 / 26 | 10→9 / 26 | 8→12 / 26 |

On distinct deposits the first seeds gain ten hits and lose two. The second seeds gain eight and lose two. All five crystallographic ligands have a net gain across the two seed sets. Their counts are 1a 5→10/12, 1b 2→5/8, 3a 2→5/12, 4 5→6/12, and 5 6→8/8. The four-seed unique Kd scores on each receiver are 15→17/37 for 3D49, 21→13/37 for 2ZFF and 10→23/37 for 2ZC9. The 17/37 pooled score applies to the signs of medians of all 12 matched gaps per ligand pair. Omitting any one compound changes the all-seed accuracy difference to a value between −6.90 and +15.38 percentage points. The native new-seed affinity result fails to reproduce the first-panel rise. In the parent benchmark the ligand-only heavy atom count scores 31/37 on precisely these experimental pairs.

As an additional fixed-pool postprocessing control, applying the parent single-water steric exclusion to the four-seed soft-search results reduces the distinct-deposit pose count from 34 to 25/52, while the number of correct unique Kd pairs moves from 17 to 20/37. The first and confirmation seeds separately yield 19→24 and 16→16/37 on this rescoring contrast. These numbers characterize Vina candidates and manual prior information about one water. They do not give a water displacement free energy. See `expected/soft_summary.csv` for every arm, receiver and seed panel.

## Reproduction

Start in the fresh work directory with unzipped `raw/`, `prepared/`, `results/` and `series.csv` following the parent README. Let `T` be the full path to the parent benchmark directory and `VINA` the full path to the newly built binary supporting `--soft_receptor_search`. The following commands audit all coordinates and compare every search arm from archived outputs.

```bash
S="$T/soft_search"
unzip -q "$S/first_soft.zip"
unzip -q "$S/confirmation.zip"
sha256sum -c "$S/files.sha256" > soft-fixture-checksums.log
python "$S/evaluate_soft_search.py" results first_soft confirmation verified
diff -u "$S/expected/soft_summary.csv" verified/soft_summary.csv
```

To regenerate the additional 216 raw runs, use a new output directory and rerun the evaluator with those freshly generated outputs.

```bash
python "$S/soft_challenge.py" "$VINA" "$PWD" rerun/soft --jobs 3
python "$S/confirm_search.py" "$VINA" "$PWD" rerun/confirm --jobs 3
python "$S/evaluate_soft_search.py" results rerun/soft/results rerun/confirm recomputed
```

The evaluator reads each raw PDBQT, reconstructs crystal correspondence against named CCD bonds, and uses independent SciPy rotations and the canonical SMILES serial map from the parent audit. `expected/` contains pose and pair results, matched gains/losses and leave-one-ligand-out sensitivity. `files.sha256` covers all 432 additional raw logs and output models.
