# Restrained internal relaxation of the reported top pose

CHEK1 produced experimental and predicted N···N distances close to the PoseBusters intramolecular-distance threshold. We tested whether relaxing covalent geometry could restore physically acceptable poses without losing their placement. Each existing first-ranked Meeko pose was reconstructed with explicit hydrogens, retaining its graph, stereochemistry and receptor frame. Optimization uses MMFF94s with a harmonic positional restraint on each heavy atom at the reported coordinates, zero flat-bottom radius and a fixed spring coefficient of 10 kcal mol⁻¹ Å⁻². All ligand atoms may move. No receptor potential or crystal coordinates enter this operation. The comparison arm removes the harmonic springs. The score from the source file applies only to the original pose. [RDKit force-field interface](https://www.rdkit.org/docs/source/rdkit.ForceField.rdForceField.html)

The diagnostic retained the 32 paired native and prepared-start top poses of the 250,000-cap cross-docking experiment. An external transfer check used one additional panel of 10 distinct human targets, three native top poses per target. We excluded BRD9, which also appears in the cross-docking set. The same ligand/control chemistry and asymmetric-unit preparation limitations as the halogen panel apply. Thirty repetitions represent 10 target–ligand cases. The conditions were chosen after the CHEK1 observation. Each crystal evaluation is receptor-frame, graph-symmetry-corrected RMSD. The output includes PoseBusters full docking checks against the conditioned protein. Every relaxation converged.

| Panel / operation | Top pose ≤2 Å | PoseBusters valid | Both | Median heavy-atom shift |
|---|---:|---:|---:|---:|
| Cross, native coordinates | 4/32 | 28/32 | 1/32 | 0 Å |
| Cross, spring 10 | 4/32 | 27/32 | 1/32 | 0.126 Å |
| Cross, unrestrained | 4/32 | 9/32 | 1/32 | 0.730 Å |
| Additional systems, native coordinates | 4/30 | 30/30 | 4/30 | 0 Å |
| Additional systems, spring 10 | 4/30 | 30/30 | 4/30 | 0.098 Å |
| Additional systems, unrestrained | 4/30 | 12/30 | 2/30 | 0.864 Å |

Three near-crystal CHEK1 77A top poses still triggered the intraligand clash check after both relaxations. The weak 1,5 contact can correspond to an internal hydrogen bond. The deposited CHEK1 A42 reference also triggers the check. An empirically chosen cutoff can remove true geometries. The current experiment does not establish a ligand-only repair for this case. [Every before-and-after coordinate metric and failed check](results.csv)

A single operation on any supported native Meeko pool is available with `python tools/relax_pose_geometry.py native.pdbqt relaxed_top.sdf`. The default operation writes one top geometry. A full replay needs the two frozen benchmark packs and their separate Python environments. From the repository root, set POSE and HAL to two empty directories. The commands below use PY3 for the structure-preparation Python and PB for the Python with PoseBusters 0.6.3.

```bash
P=benchmarks/torsion_proposals
H=benchmarks/assembly_halogen
mkdir -p "$POSE" "$HAL/results"
for a in "$P"/fixtures/*.tar.xz; do tar -xJf "$a" -C "$POSE"; done
cp "$P"/results/modes-*-250000.csv "$POSE"/
"$PY3" "$P/export_frames.py" . "$POSE"
"$PB" "$P/check_validity.py" "$POSE" --deep
tar -xJf "$H/fixtures.tar.xz" -C "$HAL"
cp "$H/results/cells.csv" "$HAL/results/"
"$PB" benchmarks/torsion_search/geometry/check.py . "$POSE" "$HAL" geometry.csv
```
