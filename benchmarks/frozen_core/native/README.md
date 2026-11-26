# Independent native-search confirmation of the rotor annotation

The chemistry, seven receptors and exact Ki/Kd measurements come from the [five-target affinity panel](../../geometry_comparison/affinity/README.md). Seeds 96431 and 96457 are new. Native search and native refinement used two trajectories, 12,000 evaluations per trajectory and one CPU. The correction was defined before generating these poses. Each of 140 jobs was run twice, with and without the annotation. After the single additional remark was removed every output was byte-identical. All 140 search trajectory histories and summary counters agreed.

The single-pose readout is the native first pose. Aggregate the two seed-specific values by their median, then count correct orders only among same-target, same-receptor, same-assay pairs at least threefold apart. Shared ligands, seeds and receptor copies introduce dependence. Only the Native column uses the original native torsion divisor. Raw removes that divisor. Uniform adds RT ln(3) per eligible rotor to raw, and contact adds it per rotor having a qualifying contact on each side. Upper values are better in this table.

| Target | Eligible comparisons | Native | Raw | Uniform | Contact |
|---|---:|---:|---:|---:|---:|
| MAPK14 | 58 | 29 | 29 | 36 | 35 |
| ALK | 19 | 15 | 10 | 15 | 15 |
| RET | 19 | 13 | 9 | 13 | 13 |
| TNK2 | 19 | 9 | 8 | 13 | 13 |
| BACE1 | 72 | 23 | 25 | 25 | 28 |
| Total | 187 | 89 | 81 | 102 | 104 |

Contact reversed 15 native mistakes, with no correct native decisions lost. The contact-specific increment over the uniform penalty was two pairs. This is a descriptive paired result. These are new seeds on previously tested chemicals. Factor Xa in the [external series](../affinity/) is a separate transfer test.

The independent verifier constructs the original PDBQT branch cuts and uses a SciPy KD-tree for contacts. C++ cuts inferred rotatable bonds in the chemical graph and uses explicit receptor distances. Across all 1,190 modes their counts agreed exactly. The largest score discrepancy was 0.000968 after serialization to PDBQT.

## Reproduction

From the repository root, prepare the existing affinity inputs, then either replay the deposited outputs or run all 280 CLI invocations. The run writes annotated PDBQT files, paired logs, and a hash of each unannotated pool.

```bash
B=benchmarks/frozen_core/native
A=benchmarks/geometry_comparison/affinity
mkdir -p native-inputs native-work
tar -xJf "$A/inputs.tar.xz" -C native-inputs
cp "$A/panel.csv" native-inputs/
tar -xJf "$B/outputs.tar.xz" -C native-work
python "$B/confirm.py" native-inputs native-work
cmp "$B/cells.csv" native-work/cells.csv
cmp "$B/pairs.csv" native-work/pairs.csv
python "$B/check_fixed.py" build/vina
python "$B/confirm.py" native-inputs native-rerun --run --binary build/vina --jobs 4
```
