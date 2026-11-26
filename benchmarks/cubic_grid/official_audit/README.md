# Independent check of official Vina force fields and side-chain flexibility

The verifier reads ligand chemistry from the PDBQT SMILES and atom maps, and reference chemistry/coordinates from the original CCD and PDB. It obtains pocket frames by an SVD fit. It identifies ligand atoms from their residue name even when a flexible protein residue reuses ligand atom serials. The reference reader is the independent deposition reader in `human_panel/verify.py`. Neither the acid-state nor the flexibility run script supplies calculations to this check.

| Experiment | Search outputs | Modes compared | Maximum absolute RMSD difference |
|---|---:|---:|---:|
| Three acids, two receptor structures each, two force fields, two chemical states, three seeds | 72 | 637 | < 3 × 10⁻¹⁴ Å |
| Four ligands, two receptors per target, rigid/flexible, three seeds | 48 | 425 | < 3 × 10⁻¹⁴ Å |

The heavy-atom coordinates of the two prepared chemical states matched exactly for all three acids, including all 36 heavy atoms of T4B, 23 of 8CA and 28 of 2T1. PDBQT summed charges for the respective anions were −0.998, −1.001 and −1.000 in elementary-charge units. Each anion also loses one torsion associated with the OH hydrogen. Tests of carboxylate oxygen equivalence explicitly exchange the two terminal oxygens after all graph mappings. Applying the same exchange to a protonated acid supplies a sensitivity analysis for uncertain crystal proton assignments. In the neutral chemical graph its labeled carbonyl and hydroxyl oxygens are chemically distinct. Oxygen exchange changed no 2 Å cell classification at the top 1, 5 or all returned modes. The native Vina HA → A first-position hits are 6 → 5 of 18. AD4 gives 4 → 6 of 18.

The verifier checked equality of every prepared receptor atom's identity, position, type and charge when the rigid and flexible input portions are reunited. Across 24 paired cells rigid and flexible each recovered 5 first-ranked poses. Side-chain RMSDs omit Cα and hydrogen and use the deposited frame. In top-ranked flex outputs, the median side-chain RMSD to the starting receiver was 0.446 Å for Thr60 and 1.988 Å for Gln98; the median against the corresponding query-bound side-chain conformation following pocket alignment was 1.006 and 2.328 Å. Distortions from a single modeled side-chain starting geometry are not solution-state free energies. The per-model values are given in `sidechain_rmsd.csv`.

To reproduce, extract the raw structures and prepared inputs from `human_panel` to `fixtures`. Extract each of the two official experiment archives into its own output folder, including input state PDBQT and `modes.csv`. Use the Python environment specified for the human panel.

```bash
python verify_official.py "$repo/benchmarks/cubic_grid/human_panel" "$fixtures" "$acid_output" "$flex_output" "$audit_output"
```
