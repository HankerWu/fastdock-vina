# Human thrombin S1 hydration and affinity panel

This benchmark evaluates one falsifiable transfer hypothesis on human alpha-thrombin (F2, UniProt P00734). If the main omitted effect were steric competition with a recurrent S1 water, penalizing access to its space should improve the selected pose and the ranking of ligands in the same calorimetry assay. The water is localized in a crystal obtained without a modeled active-site inhibitor. The alternative and native rankings use precisely the same candidate coordinates. The experiment is retrospective and the known S1 chemistry informed selection of this system.

## Experimental reference

Baum et al., *Journal of Molecular Biology* **390** (2009), 56–69, Figure 1 and Tables 1 and 3. The authors fitted a one-site isotherm to ITC heats. They report ΔG° derived from the fitted Kd. We transcribed Table 1 in kJ/mol and reconstructed Kd (1 M standard state, T = 298.15 K); the paper's one-decimal rounding carries through. All 12 assay compounds are included. All are from the same D-Phe–Pro scaffold with systematic S1 substituents. The ITC conditions were 25 °C, 50 mM Tris-HCl pH 7.8, 100 mM NaCl, 2.5% DMSO and 0.1% PEG 8000. The authors report at least triplicate measurements, with compound 5 measured in duplicate. The 26 inhibitor Ki values were measured by a separate substrate competition assay at pH 7.4, 154 mM NaCl and 5% DMSO. This panel scores only the ITC endpoint. Primary data and conditions are in the [article](https://doi.org/10.1016/j.jmb.2009.04.051).

Five ligands have a crystallographic reference (5, H, 2ZFF; 1a, m-F, 2ZDV; 3a, m-Me, 2ZF0; 1b, m-Cl, 2ZC9; 4, p-amidinium, 2ZDA). Three receiving structures are 3D49, 2ZFF and 2ZC9. The 3D49 model has another benzamidine molecule; the minimum distance from its **A** conformer's heavy atoms to the S1 water used here is 13.06 Å. The experimental structures of the five benchmark ligands were withheld from ligand preparation, selection and docking. Two receivers come from benchmark ligands; they were included to compare open/occupied-pocket conformations, with cognate receptor/ligand cells reported separately. The binding-site landmarks and the water-site hypothesis use prior structural knowledge. Source URLs and file hashes are in `sources.json`. wwPDB coordinate files are provided under the [wwPDB public-domain policy](https://www.wwpdb.org/about/usage-policies).

## Search and score

All input ligands were rebuilt from the isomeric graph with an ammonium D-Phe terminus and, for compound 4, protonated amidinium. Eight ETKDGv3 starts and MMFF94 optimization selected the lowest-energy converged conformer per ligand, seed 408151. Only RDKit-generated coordinates entered PDBQT preparation. Two stereo centers were verified as R and S for every graph. This common fixed protomer choice has not been validated by an independent titration. The PDBQT ligands and dry protein coordinates were prepared using Meeko 0.8.0. Receptors keep chains H and L, remove crystallization/solvent components, and select altloc blank/A. Incomplete residues omitted by Meeko are recorded in the preparation logs. Minimum distance from any omitted observed atom to the five aligned crystal ligands is 9.21 Å in 3D49, 21.93 Å in 2ZFF and 21.89 Å in 2ZC9.

A 25 Å cube is centered on the mean of 11 H-chain Cα landmarks (57, 60A, 60D, 99, 174, 189, 195, 215, 216, 227 and 228). The same landmarks, and the protein alone, are used to align each crystal to its receiving structure. Each receiver/ligand pair has two independent Late Vina searches (seeds 83519 and 194027), single CPU, exhaustiveness 4, max_evals 50000, 12 output modes and 16 refinement candidates. There are 72 process logs and 864 returned poses. The binary was built from a checkout with tree f04fa353be2aff5f8dd32acc1e2fd537cf2a4a26, identical to the GitHub baseline commit 79c793afb3f3ffe257a07845917e3d2ba762474e. Its binary SHA-256 was 7e26c1abd31ccbcdba8ce88a62e1c9cf5763d7bffe0019e2d24ef96319b4d39e.

The source water is the 3D49 water nearest the midpoint between the backbone carbonyl oxygens of H-chain Trp215 and Phe227, H:2006. Its distances from those oxygens are 3.13 and 3.09 Å. It is mapped to each receiver by the same protein alignment. Let d be the shortest distance between its position and any ligand heavy atom. Let t = clip((d−2.0)/0.6, 0, 1) and A = t²(3−2t). Add −RT ln[(1−p)+pA] to the printed native Vina affinity, with R in kcal/mol/K, T = 298.15 K, and p = (1+0.5)/(1+1) = 0.75. The last number is a one-source transfer of the prior benchmark’s pseudocount convention; there is no measured 75% solution occupancy. The implemented correction has a maximum 0.8214 kcal/mol. This tests only the single-site **exclusion component**, using its previous distance window. It is a one-source mechanistic challenge and does not satisfy the historical multi-donor full hydration pipeline's two-source site filter. Results cannot separate solvent displacement free energy, halogen interactions and ligand desolvation.

The native/modified ranking uses the lowest score, with model order breaking ties. Pose RMSD is calculated in the receptor frame, allowing element- and connectivity-preserving ligand symmetry, with no ligand superposition. Within each receptor and seed the affinity gap is the weaker ligand's selected energy minus the stronger one's selected energy. For each distinct ligand pair the reported prediction is the median of its six gaps (or the two gaps for a single receiver), with half credit for ties. We include only experimental gaps of at least RT ln 3. The 37 pairs share ligands.

## Results

| Metric | Native | Exclusion |
|---|---:|---:|
| Top pose ≤2 Å, all five labeled ligands × three receivers × two seeds | 14/30 | 12/30 |
| Top pose ≤2 Å, distinct ligand and receiver deposits only | 11/26 | 9/26 |
| Best candidate ≤2 Å, all / distinct-deposit cells | 19/30, 16/26 | 19/30, 16/26 |
| Correct unique threefold-separated Kd pairs, all receivers | 11/37 | 10/37 |
| Correct unique Kd pairs, 3D49 | 12/37 | 16/37 |
| Correct unique Kd pairs, 2ZFF | 18/37 | 21/37 |
| Correct unique Kd pairs, 2ZC9 | 10/37 | 7/37 |

The closest crystal ligand atom to the transferred water lies 0.26 Å away for m-Cl, 0.65 Å for m-Me and 0.63 Å for m-F; the nearest reported water in each bound structure is 4.95–5.17 Å away. The unsubstituted and amidinium structures retain waters 0.70 and 1.01 Å from the transferred site. In the unsubstituted case one ligand atom is also just 2.21 Å from the transferred site; a rigid exclusion field assigns cost even when a shifted water was observed. The alignments' Cα RMSD spans 0.09–0.19 Å. The paper also reports differences in S1 ligand mobility, desolvation and binding enthalpy among these substituents. The observed replacement and the failures of the ranking jointly constrain the interpretation.

Simple ligand-only controls emphasize series composition. More heavy atoms predicts 31/37 with half credit on ties, greater molecular mass 29/37, and greater RDKit Crippen logP 26/37. Removing the 11 pairs involving the additional positive charge of compound 4 leaves 26 pairs among 11 monocations; logP scores 26/26 and native Vina 1/26. Extending native Vina by a full octanol-transfer term (−RT ln(10) × calculated logP, no fit; a deliberately crude reference) gives 14/37. Calculated logP on these fixed charged graphs is a descriptor, not a measured distribution coefficient, and the association within this small scaffold is not a protein-specific affinity prediction. A positive or negative change on these labels should always be compared to the ligand-only controls.

## Reproduction

Python 3.12.14; dependencies used here were RDKit 2026.9.2, Meeko 0.8.0, NumPy 2.3.5, SciPy 1.17.0 and Gemmi 0.7.5. The archived experimental inputs, prepared fixtures and results can also be audited without network access. From the repository root (a fresh work directory is recommended):

```bash
T="$PWD/benchmarks/thrombin_water"
VINA="$PWD/build/vina"
workdir=$(mktemp -d)
cd "$workdir"
unzip -q "$T/raw.zip"
unzip -q "$T/prepared.zip"
unzip -q "$T/results.zip"
sha256sum -c "$T/files.sha256" > checksums.log
cp "$T/expected/series.csv" series.csv
python "$T/analyze_panel.py"
python "$T/controls.py"
python "$T/independent_audit.py"
```

To redo all 72 searches with the exact prepared fixtures, run `python "$T/run_panel.py" "$VINA" --prepared --jobs 3` and rerun the three analysis commands. Omitting `--prepared` rebuilds all inputs from the six archived PDB files and the formulae and stereochemistry in the script. This command regenerates the local `series.csv` and overwrites existing work-directory results. Use `expected/` to compare numerical outputs. Cross-platform and cross-version differences in conformation generation, floating point operations or compilation can change estimates. The independent analysis rereads all models, overlays proteins using SciPy Rotation.align_vectors, reconstructs atom correspondence from CCD bond tables and Meeko's canonical SMILES-index map (the primary analysis uses coordinate-derived bonds and the original input atom index), and calculates pair eligibility from Kd ratios directly. It reproduced the results above with maximum per-candidate RMSD disagreement below 4 × 10⁻¹⁵ Å on the 360 labeled poses and maximum corrected-score disagreement below 8 × 10⁻¹⁵ kcal/mol on all 864 poses.
