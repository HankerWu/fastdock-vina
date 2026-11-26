# Fixed-coordinate polar transfer cycle

The primary experiment tests receptor/ligand Coulomb attraction together with the cost of changing the dielectric boundary. Each source PDBQT contains docked coordinates, merged-H Meeko charges and explicit donor hydrogens. At those coordinates the model evaluates

$$\Delta E_{\rm polar} = E_{\rm Coul}^{R,L}/4 + G_{\rm OBC2}(RL)-G_{\rm OBC2}(R)-G_{\rm OBC2}(L).$$

Here the numerator of the first term is the vacuum receptor/ligand charge interaction. OpenMM computes OBC-II polar energies with internal dielectric 4 and external dielectric 78.5, no salt, no pair cutoff and zero nonpolar surface-area term. All three legs use the same charge and cavity conventions. The receptor coordinates are fixed and the free ligand is held in each candidate's own internal geometry. Radii and screening factors follow OpenMM mbondi2/OBC2, retaining the benchmark's united-atom PDBQT representation. This approximation omits nonpolar hydrogen volumes and water structure, permits no dielectric or conformational relaxation and uses charges distinct from the all-atom force fields ordinarily paired with GB. Receptor Zn is removed from both OBC legs and the new Coulomb term. The native contribution still scores it. This avoids a spurious +2 point-charge contact but omits metal-induced polarization, charge transfer and its dielectric volume. [Onufriev, Bashford & Case](https://doi.org/10.1002/prot.20033) · [OpenMM GBSAOBCForce](https://docs.openmm.org/latest/api-python/generated/openmm.openmm.GBSAOBCForce.html)

For each pool, compute the normalized native Vina score with a common unbound reference equal to the intra energy of the raw-native winner. Add the cycle in kcal/mol after this normalization. Pose selection and same-assay ligand comparisons use the resulting score. Calibration used zero affinity or RMSD labels. Separate controls add cross-Coulomb only, GB difference only, cycle at internal dielectric 2 or 8, cycle before the Vina transformation, or subtract the engine's reported uncapped native N/O intermolecular H-bond component. Vina's empirical terms can overlap the new energy. Subtraction is a component diagnostic because per-atom caps are nonlinear. Another control scores the original native winner with the cycle to separate selection from fixed-coordinate offsets. With retained scaffold geometry, the new composite affinity contains the geometry term; its affinity result has a different definition from the previous round's native-at-selected report. Both are saved separately.

| Fixed pool | Native → full cycle top-1 | Native → full cycle ≥3-fold affinity ordering |
|---|---:|---:|
| hCA II two zinc-bound tertiary-N input arrangements, scaffold-generated, 8 valid candidate pools | 1 → 1 / 20, heavy-atom RMSD ≤2 Å and correct N arrangement | 9 → 2 / 28 |
| BRD4 native | 1 → 1 / 24 | 10 → 12 / 30 |
| BRD4 scaffold-generated | 7 → 3 / 24 | 9 → 9 / 30 |
| T4L native | 28 → 28 / 36 | 50 → 46 / 76 |
| T4L scaffold-generated | 28 → 28 / 36 | 49 → 47 / 76 |
| CA I/XII zinc-point, 36 Å | 9 → 4 / 10 | 2 → 2 / 4, CA I only |
| Same small-ligand panel, 20 Å | 5 → 2 / 10 | 3 → 0 / 4, CA I only |

The 2,067 poses include 538 long hCA II poses, 422 BRD4, 588 T4L, and 262 + 257 short-ligand poses. Repeated receptor/seed observations and alternate scoring of the same coordinates are correlated. Distinct pair counts are reported separately. The short-ligand labels use the original CA I DSF comparison; CA XII has no affinity pair. Timing and temperature of the experimental measurements remain panel-specific.

On BRD4 native candidates, Coulomb alone produces 27/30 correct, including 18/18 comparisons between different formal charge states, while losing both cells' possible near-native first pose (native top-1 1/24 to 0/24; candidate-positive cells 2/24). On hCA II, GB alone gives 20/28. Adding its cancelling Coulomb partner yields 2/28. For the strongest hCA II ligand 46 (6WQ8, charge −1), median Coulomb and polarization are +24.96 and −15.53 kcal/mol, with net polar-transfer median +9.16. For ligand 41 (6WQ5, +1), the corresponding values are −12.24, +21.19, +8.36. These medians describe individual components; a sum of medians need not equal the median of a sum. The relative change between per-ligand residuals is large enough to overwhelm several native separations. The two-state hCA II pool has only one ≥3-fold same-charge pair, repeated four times. Its native and cycle credit are both 0/4. The remaining six pairs contribute 9/24 → 2/24. There is no evidence for promoting this endpoint to a general affinity or pose score.

A separate blockwise float64 implementation evaluates the HCT descreening integrals, OBC2 Born radii, pair and self energies on three poses per panel. The 538 hCA II native INTER, INTRA and torsion factors agree with the previous independently reconstructed fixed-topology scores within 5×10⁻⁹. See `results/independent.csv` and `results/charge_pairs.csv` for audit values.

## Recompute

The study uses the previous anchor evidence files. Arguments for `run.py` are panel, pool directory, preparation directory, native binary, and output CSV path. Install NumPy, SciPy and OpenMM, using one OpenMM CPU thread per worker. For hCA II use the double-state scaffold directory with its `inputs` and `candidates.csv`. BRD4/T4L use native and scaffold pools with `modes.csv`. For the small ligands, point to `runs/` under each holdout experiment.

```bash
python tools/polar_endpoint.py receptor.pdbqt prepared_ligand.pdbqt native_poses.pdbqt
python benchmarks/polar_endpoint/run.py hcaii "$state_runs" "$ca2_inputs/prepared" build/vina "$results/hcaii.csv"
python benchmarks/polar_endpoint/run.py brd4 "$brd4_runs" "$brd4_inputs/prepared" build/vina "$results/brd4.csv"
python benchmarks/polar_endpoint/run.py t4l "$t4l_runs" "$t4l_inputs/prepared" build/vina "$results/t4l.csv"
python benchmarks/polar_endpoint/run.py ca_holdout36 "$hold36/runs" "$hold36/prepared" build/vina "$results/ca_holdout36.csv"
python benchmarks/polar_endpoint/run.py ca_holdout20 "$hold20/runs" "$hold20/prepared" build/vina "$results/ca_holdout20.csv"
python benchmarks/polar_endpoint/evaluate.py "$results" "$state_runs/rank_energies.csv"
python benchmarks/polar_endpoint/charge_audit.py "$results" --hcaii "$ca2_inputs/prepared" --brd4 "$brd4_inputs/prepared" --t4l "$t4l_inputs/prepared" --ca_holdout36 "$hold36/prepared" --ca_holdout20 "$hold20/prepared"
```
