# Same-assay binding affinity on SYK and focal adhesion kinase

The original Davis et al. competition-binding panel provides 19 exact Kd entries for unmodified human SYK and 18 for FAK (PTK2). The fixed name/chemical-identity dictionary from the two earlier kinase benchmarks covers 14 and 12, respectively. All 26 entries have assigned graph stereochemistry. They represent 16 distinct compounds across both targets and a distinct typed Bemis–Murcko scaffold per compound within each target. Measured Kd and molecular inputs are in [panel.csv](panel.csv), and missing identities appear in [exclusions.csv](exclusions.csv). Table 1 confirms each assay and Table 4 supplies its Kd in nM. [Original paper](https://doi.org/10.1038/nbt.1990)

## Paired results

Both arms used the same 26 Å box, ligand PDBQT, rigid receptor, seeds 76123 and 76231, two trajectories, cap of 12,000 evaluations per trajectory, nine requested modes, native explicit Late refinement, and energy range 100. Each ligand's arithmetic mean top-pose score over the two seeds defines its primary ordering. Count a measured pair if the weaker Kd is at least three times the stronger Kd. Equal predicted scores earn half credit. Each pair counts once in the consensus. The controls use the preceding release and a descending heavy-atom-count ranking.

| Target (structure) | Compounds | ≥3-fold pairs | Native | Cubic | Larger compound stronger | Spearman native → cubic | Changed pairs, gains / losses |
|---|---:|---:|---:|---:|---:|---|---:|
| SYK (6ZCY) | 14 | 68 | 43 | 47 | 39.5 | 0.367 → 0.402 | 7 / 3 |
| FAK (6YT6) | 12 | 47 | 9 | 19 | 37.5 | −0.601 → −0.175 | 15 / 5 |

Both seeds exhibit the same direction on these receptor preparations. Across 10,000 whole-compound bootstrap draws, the descriptive 95% intervals for the increase in pair accuracy are −11.8 to +24.3 percentage points for SYK and −22.9 to +63.8 for FAK. Each target contributes one experimental receptor conformation. Assay-specific ionization populations and bound waters were unavailable. The native and cubic runs used 1,252,604 and 1,252,700 measured trajectory evaluations and produced 451 and 458 modes. End-to-end elapsed times from overlapping CPU jobs are diagnostic. This comparison controls the evaluation cap. [All rows](summary.csv) · [Individual pairs](pairs.csv) · [Scores, timing and evaluation counts](cells.csv)

An independent script extracted assay identifiers and labels directly from both original Excel files, read every archived model pool, and recovered all four consensus counts. Every serialized first pose was rescored with its original topology, original ligand-specific unbound reference and the native explicit scoring function. The maximum absolute discrepancy was 0.007 kcal/mol, and all four correct-pair counts persisted. All 52 native pools were byte-identical to outputs from revision ce9a89f. This study used engine revision 23b24f6c53892e3edef8fbfe726cd2b34a080a89. [Raw-output checks](independent.json) · [Native replay](native_parity.json)

## Experimental coordinates and repair sensitivity

Input protein coordinates and the separate box-defining ligand SDF come from the PoseBusters authors' collection (archive MD5 f004ac7c4e68317b5348497d2bb6bee6). The compounds were prepared from their molecular graphs with the existing independent ETKDG/MMFF workflow. Its one chemical form per molecule and one dry protein structure per target limit comparison to measured binding equilibria. [PoseBusters data](https://zenodo.org/records/8278563) · [Method](https://doi.org/10.1039/D3SC04185A)

FAK binds PKE in chain B. The other chain's closest heavy atom lies 21.82 Å from PKE, and sulfate lies 38.38 Å away. All retained protein heavy atoms match the input coordinates. The nearest observed Cα bounding a residue-number gap is 11.88 Å away. [Structure 6YT6](https://www.rcsb.org/structure/6YT6)

SYK binds QF8 in chain A. The prepared model contains 48 added heavy-atom positions, the nearest 6.09 Å from QF8; the nearest Cα bounding a numbering gap lies 6.06 Å away. The previous explicit protein model contained a 1.68 Å contact between modeled Lys386 NZ and an experimental backbone oxygen. Adjusting the dihedrals of Lys386 and Ser379 changes three reconstructed heavy-atom coordinates. All retained observed heavy atoms have zero displacement. After correction the recorded check finds zero contacts below 2 Å between reconstructed heavy atoms and non-adjacent residues. The discarded coordinates gave 49 → 44 correct pairs on the same 68 SYK pairs, versus 43 → 47 in the corrected model. This change in sign makes the sensitivity to structural preparation part of the interpretation. The crystal pose only defines the search box and preparation-distance audit. Protein-only clash calculations selected modeled rotamers. [Structure 6ZCY](https://www.rcsb.org/structure/6ZCY) · [Coordinate comparison](coordinate_audit.json) · [All preparation details](structure_audit.json)

## Replay

Run commands from the repository root, substituting this benchmark's directory for B. The archive includes the original assay tables, chemical identity records, original coordinate pairs, prepared PDBQT and all accepted raw pools/logs. The reference archive additionally includes the discarded-repair protein, 56 raw pools and byte-parity verification. Use the archived prepared inputs to replay docking directly.

```bash
B=benchmarks/cubic_grid/affinity
W=$(mktemp -d)
(cd "$B" && sha256sum -c evidence.sha256)
tar -xJf "$B/evidence.tar.xz" -C "$W"
cp "$B"/panel.csv "$W"
tar -xJf "$B/native-and-cubic-output.tar.xz" -C "$W"
python "$B/analyze.py" "$W" build/vina
python "$B/verify.py" "$W" "$W/sources/table1.xls" "$W/sources/table4.xls"
python "$B/coordinate_audit.py" "$W"

F=$(mktemp -d)
tar -xJf "$B/evidence.tar.xz" -C "$F"
cp "$B"/panel.csv "$F"
python "$B/run.py" "$F" build/vina
python "$B/analyze.py" "$F" build/vina
python "$B/check_native.py" "$F" /path/to/ce9a89f/vina
```

For regeneration of the prepared ligand and protein inputs use Python 3.12, RDKit 2026.09.2, Meeko 0.8.0, PDBFixer 1.12.0, OpenMM 8.6.1 and xlrd 2.0.2. Source coordinates are in the two original PoseBusters folders under sources. Preparation logs record an intermediate SYK repair superseded by the final protein-based rotamer selection.

```bash
D=$(mktemp -d)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OPENMM_CPU_THREADS=1
python "$B/prepare_panel.py" . "$W/sources" "$D"
python "$B/prepare_receptor.py" . "$W/sources" "$D"
python "$B/coordinate_audit.py" "$D"
```
