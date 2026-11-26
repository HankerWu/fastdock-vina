# Transfer of a frozen contact ranker to two additional human targets

The prespecified contact ranker recovered a near-crystal roscovitine pose against the independently determined human CDK2–staurosporine receptor in both new seed runs. Native Vina and its nearest-element template selected poses over 6.5 Å from the reference in those runs. The reverse CDK2 direction and both directions for PARP1 supplied no pose within 2 Å at this search budget. The data support one directional cross-compound recovery with two seed repeats.

## Design

Two human targets, four deposited complexes, both reciprocal directions, one independent experimental ligand in the receiving pocket and two new seeds (98207 and 98213) were specified before docking. The published two-channel ranker and all thresholds in the [frozen protocol](../contacts/PROTOCOL.md) were used without recalibration. These are additional receptor/ligand structures for this contact method. Other repository studies previously included PARP1 or CDK2. [PARP1–olaparib, 7KK4](https://www.rcsb.org/structure/7KK4) · [PARP1–talazoparib, 4UND](https://www.rcsb.org/structure/4UND) · [CDK2–roscovitine, 2A4L](https://www.rcsb.org/structure/2A4L) · [CDK2–staurosporine, 1AQ1](https://www.rcsb.org/structure/1AQ1)

All four first-choice structures passed the existing single-chain dry preparation protocol. Both CDK2 asymmetric units have no omitted nonwater heavy atoms. For PARP1, the nearest such atom is 11.66 Å from olaparib and 20.41 Å from talazoparib. No missing protein heavy atoms were recorded in the prepared receptors. Starting query conformations were regenerated from CCD chemical graphs with the fixed ETKDGv3/MMFF recipe. Each receiver used its experimental ligand to center the 24 Å box and to provide the template for conditional ranking. The query ligand's experimental coordinates enter only the subsequent RMSD calculation. Meeko's default prepared protonation, dry rigid receptors and incomplete sampling affect these outcomes.

Each CLI execution used the branch's native selection and optimizer defaults, four single-threaded trajectories, max_evals 60,000 per trajectory, nine output modes and energy_range 100. Native, element-template, combined general/polar contact, general-channel alone and overlap-only scores were evaluated on exactly the same nine printed poses per cell. The two primary comparisons use the same fixed prior strength 6 on the pre-divisor scale. Native order breaks ties. The footprint is an empirical prior provided by another compound with a known ligand binding site.

## Results

| Query → receiver | Seed | Native RMSD (Å) | Element RMSD (Å) | Contact RMSD (Å) | Contact-selected native rank | Best in pool (Å) |
|---|---:|---:|---:|---:|---:|---:|
| CDK2 staurosporine → roscovitine structure | 98207 | 4.129 | 4.129 | 4.129 | 1 | 4.129 |
| CDK2 staurosporine → roscovitine structure | 98213 | 5.970 | 5.970 | 5.970 | 1 | 5.970 |
| CDK2 roscovitine → staurosporine structure | 98207 | 6.518 | 6.518 | 1.927 | 2 | 1.927 |
| CDK2 roscovitine → staurosporine structure | 98213 | 6.550 | 6.550 | 1.793 | 3 | 1.793 |
| PARP1 talazoparib → olaparib structure | 98207 | 2.184 | 2.184 | 2.184 | 1 | 2.184 |
| PARP1 talazoparib → olaparib structure | 98213 | 2.179 | 2.179 | 2.179 | 1 | 2.179 |
| PARP1 olaparib → talazoparib structure | 98207 | 2.906 | 3.109 | 3.109 | 2 | 2.906 |
| PARP1 olaparib → talazoparib structure | 98213 | 2.854 | 2.854 | 2.854 | 1 | 2.854 |

Top-one hits are 0/8 native, 0/8 element, 2/8 contact, 0/8 general alone and 2/8 contact overlap alone. Two pools contain an attainable hit. The stricter 2 Å endpoint counts six pool-coverage failures. In the one worsening PARP1 run a template arm selected 3.109 Å, while the native 2.906 Å pose was best available.

In the recovered CDK2 direction the wrong native poses had *higher* general overlap (0.822 and 0.830) than the two near-crystal poses (0.717 and 0.723). Their N/O-channel overlaps were only 0.138 and 0.144, versus 0.482 for each near-crystal pose. The combined score selected a pose that sacrificed 0.381 or 0.351 native energy units. The total winning margins were approximately 0.099 and 0.113 corrected score units. This identifies the contribution of the N/O proximity channel in these pools. Neither donor–acceptor compatibility nor hydrogen bond direction is encoded in this channel. No affinity inference follows from this pose-only comparison. [Ros­covitine structural study](https://doi.org/10.1111/j.1432-1033.1997.0518a.x) · [Staurosporine structural study](https://doi.org/10.1038/nsb1097-796)

## Independent verification and replay

Protein-frame alignment used 51 shared pocket Cα atoms for CDK2 (0.570 Å Cα RMSD) and 50 for PARP1 (1.071 Å). Native-output pose RMSD uses each complete symmetry-equivalent, stereo-compatible CCD atom mapping, with the receptor-fixed transform and no subsequent ligand superposition. A separate implementation rebuilt output molecules from Meeko, downloaded RCSB ModelServer ligand-instance SDFs independently of CCD atom names, reread the original protein deposits, used SciPy's vector alignment, and compared unmoved ligand coordinates using RDKit CalcRMS. All 72 independent measurements agree to 9.55×10⁻⁷ Å. Cross-docking RMSD measures disagreement with the ligand pose observed in a different receptor structure after pocket registration. Experimental uncertainties and true induced pose changes contribute.

```sh
B=benchmarks/frozen_core/transfer
W=$(mktemp -d)
(cd "$B" && sha256sum -c SHA256SUMS)
tar -xJf "$B/inputs.tar.xz" -C "$W"
tar -xJf "$B/outputs.tar.xz" -C "$W"
python "$B/score.py" "$W" --out "$W"
python "$B/audit_pose.py" "$W" --out "$W"
cmp "$W/cells.csv" "$B/cells.csv"
cmp "$W/modes.csv" "$B/modes.csv"
# Optional reproduction of the eight original docking calculations
python "$B/run.py" dock "$W" --binary /absolute/path/to/vina
```

The input archive contains raw PDBs and CCD records, the four independent receiver/query ligand-instance SDFs, fully prepared inputs, original source and preparation metadata, and the original binary fingerprint. The output archive contains all generated coordinates and CLI logs. The scripts use Python, NumPy, SciPy, RDKit, Meeko, Gemmi, OpenMM and PDBFixer from the original geometry panel. The fixed pose statistics select a winner from each common returned pool and have no influence on how that pool was generated.
