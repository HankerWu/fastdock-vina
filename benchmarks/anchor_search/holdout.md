# Ligand-independent site transfer to human carbonic anhydrases I and XII

The primary transfer test was specified before docking these ligands. Two thiadiazole-family drugs, acetazolamide (AZM) and methazolamide (MZM), differ from the longer phenylsulfonamide series used on CA II. CA I tested the two drugs on both experimental receptor conformations (1AZM and 1BZM). CA XII tested AZM on 1JD0. Each receptor–ligand combination used two matched seeds, 314159 and 271828. The proteins share the three-histidine catalytic zinc motif. This comparison covers related human isoforms, one short-ligand series and one CA I affinity pair.

[CA I crystal structures, including both ligands](https://doi.org/10.1006/jmbi.1994.1655) · [CA XII crystal structures](https://doi.org/10.1073/pnas.161301298)

## Protein, ligand and anchor preparation

* One rigid protein chain A and its experimental zinc were retained. PDBFixer added missing heavy atoms in observed residues with seed 2026, including the Asp4 side chain in 1AZM. Unresolved residues, waters, additives and other chains were omitted. The closest protein heavy atom of the second chain in the 1JD0 dimer is 12.34 Å from the chain-A experimental ligand. Meeko used HID for His94/His96 and HIE for His119, with the experimentally coordinating histidine nitrogen free of hydrogen.
* The ligand was generated once per drug from the CCD canonical isomeric SMILES, using the bound sulfonamide monoanion, ETKDG seed 2026 and MMFF optimization. The methazolamide imine stereochemistry was retained. The query donor was identified by the sulfonamide chemical pattern and mapped to its PDBQT serial. Independent ligand-instance SDF coordinates confirm the stereochemistry. Building the search inputs uses no coordinates from the experimental ligand.
* The current three-coordinate zinc geometry check in `tools/select_zinc_pose.py` defines the vacant direction from the protein donors. The directional point is 2 Å from Zn on this ray. The radial shell has the same Zn origin and a 2 Å target radius. Both arms leave native Late ranking and refinement in place. The search box is centered 5 Å from Zn along the protein-derived ray. Across the five primary ligand/receptor pairs all reference atoms are at least 12.05 Å inside the 36 Å box. CA I's shared receptor conformations differ by 0.149 Å over matched Cα positions.
* Each arm used one thread, exhaustiveness 4, 50,000 nominal evaluations per trajectory, up to 9 output modes and an energy range of 100. Ten native outputs from the previous binary reproduced byte-for-byte with the new binary. Nineteen input/manifest files from the primary preparation reproduced byte-for-byte in a separate directory. Actual primary hunt plus second-stage grid evaluations were 2,000,809 (native), 2,000,676 (directional point) and 2,000,730 (radial shell).
* Pose evaluation uses heavy-atom RMSD in the receptor frame, minimizing over graph-equivalent atom assignments. Cross-docks use common Cα atoms for receptor superposition. An independent reconstruction from RCSB ligand-instance SDF, Gemmi, SciPy, Meeko and RDKit CalcRMS agreed with every primary RMSD to 0.00000109 Å. The largest pose-classification sensitivity occurs near the 2 Å threshold, e.g. 1AZM–AZM at seed 314159 gives 1.996709 Å with the point.

## Experimental strength ordering

Bai et al., Table 1, reports independent triplicate measurements on the same commercial human CA I preparation. At 57 °C in 20 mM Tris, 100 mM NaCl, pH 6.1 and 1% DMSO, their isothermal analysis of DSF yielded a Kd of 7.2 ± 0.8 µM for acetazolamide and 1.2 ± 0.1 µM for methazolamide. This sixfold pair meets the study threshold. The separately measured esterase Ki values in the same paper were 2.1 ± 0.2 and 1.9 ± 0.2 µM, which do not meet it. The rank test compares the two drugs within identical receptor and seed. Scores from different isoforms are excluded. One unique pair is repeated four times. The docking model uses the specified bound-state charges and does not include a solution microstate or thermal-unfolding model. [Isothermal DSF and enzymatic inhibition assays](https://doi.org/10.1038/s41598-018-37072-x)

### Independent pH-linked affinity control at 25 °C

A separate experiment by Morkūnaitė et al. measured the same two ligands by ITC against CA I in phosphate at pH 7.0, 25 °C. Table 1 gives observed association constants of 1.19 × 10⁶ M⁻¹ for AZM and 1.14 × 10⁷ M⁻¹ for MZM. Their reciprocals are 840 and 87.7 nM, a 9.58-fold difference corresponding to 1.339 kcal/mol. Table 3 gives intrinsic association constants of 4.2 × 10⁶ and 3.3 × 10⁷ M⁻¹ for the deprotonated ligands binding to the Zn–water receptor state (7.86-fold). These measured and derived quantities give an independent assay for the same single chemical pair. Table 1 reports an uncertainty of ±4 kJ/mol for the individual observed free energies. [Original ITC and protonation data, Tables 1–3](https://doi.org/10.3109/14756366.2014.908291)

Table 2 lists sulfonamide pKa values of 7.3 for AZM and 7.1 for MZM. Its footnote attributes both tabulated values to Matulis and Todd (2004), the original authors' reference 25. The paper additionally shows an independent example MZM potentiometric titration at 24 °C in Fig. 6C and uses the tabulated protonation parameters in its 25 °C analysis. The pKa values are reported to a single decimal place; a value of 0.167 kcal/mol should be interpreted as an estimate of about 0.17 kcal/mol. Under the paper's single-acid model, the free anion fraction is `f = 1/(1 + 10^(pKa - pH))`. At pH 7.0 these fractions are 0.334 and 0.443. At 298.15 K the associated penalties `−RT ln(f)` are 0.650 and 0.483 kcal/mol. The additive difference on the AZM-minus-MZM score is +0.167 kcal/mol. Treating each tabulated pKa as rounded to the nearest 0.1 alone allows a differential correction of roughly +0.084 to +0.251 kcal/mol; this rounding interval is not an experimental confidence interval. For example, the wide-box shell has score gaps of −0.150 and −0.098 on 1AZM, so their nominal reversals depend on the precision of this correction. The common receptor protonation factor cancels in this same-target pair. This correction assumes the bound-state Vina score approximates an intrinsic binding term; Vina's empirical training does not guarantee that interpretation.

Adding that fixed offset to the existing 25 °C comparison leaves the directional point's ranking at 2/4 in the 36 Å box and 3/4 at 20 Å. Native Late changes from 1 to 1 and 2 to 4, and radial shell from 2 to 4 and 2 to 4. Several raw errors are smaller than the offset, so these changes can be driven by a small physical term on a single pair. The majority of the observed binding difference remains after the free-ligand protonation correction. No extrapolation to the Bai 57 °C condition was used. `holdout_linkage*.csv` retains the raw and adjusted per-case gaps.

## Results

| Primary 36 Å box | Native Late | Point-guided hunt | Radial-shell hunt |
|---|---:|---:|---:|
| CA I top 1 / 8 | 4 | 7 | 4 |
| CA I any returned pose / 8 | 4 | 8 | 8 |
| CA XII top 1 / 2 | 0 | 2 | 1 |
| CA XII any returned pose / 2 | 1 | 2 | 1 |
| CA I Kd ordering / 4 | 1 | 2 | 2 |

Top-1 hits over both isoforms were 3, 5 and 3 at 1.5 Å, and 4, 9 and 7 at 2.5 Å, in native, point and shell order. Directional guidance generated more crystal-like poses than radial-only guidance in the primary broad box. Some poses still had inaccurate coordination even when their whole-ligand RMSD met the 2 Å cutoff. The scores reverse the experimental affinity order on both 1AZM-receptor trials in all three arms.

## Pocket-size sensitivity

The primary 36 Å cubes worked for direct comparison with the long CA II ligands, but the short thiadiazoles leave over 12 Å to the nearest face. After inspecting the primary experiment, all three arms were repeated on 20 Å cubes with the same center, two seeds, scores and search budget. The smallest reference-to-face margin was 4.05 Å. Results at this smaller size are a separate sensitivity test.

| Additional 20 Å box | Native Late | Point-guided hunt | Radial-shell hunt |
|---|---:|---:|---:|
| CA I top 1 / 8 | 6 | 5 | 3 |
| CA I any returned pose / 8 | 8 | 8 | 8 |
| CA XII top 1 / 2 | 0 | 0 | 2 |
| CA XII any returned pose / 2 | 1 | 2 | 2 |
| CA I Kd ordering / 4 | 2 | 3 | 2 |

In CA XII, both point-guided runs now returned a candidate within 0.99 Å of the crystal pose while ranking poses at 3.094 and 3.268 Å first. Radial-shell hunt returned first poses at 1.563 and 1.565 Å. Tightening the box alone also recovered all eight native CA I pools. Accurate localization is a major contributor to the primary wide-box effect, and the directional point shows no first-pose benefit across both box sizes. All 257 additional poses passed the independent check with maximum difference 0.00000093 Å. The two search volumes are correlated tests on the same compounds and targets.

## Fixed-candidate zinc geometry check

As a post-hoc test of why selection failed, all original output pools were re-ranked using the existing `tools/select_zinc_pose.py` hard, continuous and radial rules. All geometry cutoffs and score parameters in that selector predate these ligands. The hard rule requires a compatible monodentate direction and Zn–N distance, and otherwise abstains. Counting each abstention as a missed whole-ligand pose, the hard rule on the 36 Å native/point/shell pools gave 4/10, 10/10 and 6/10 (5, 0 and 0 abstentions). At 20 Å it gave 8/10, 6/10 and 5/10 (1, 0 and 0 abstentions). The continuous directional rule gave 4, 8 and 8 hits at 36 Å, and 8, 7 and 7 at 20 Å. Enforcing plausible metal geometry does not determine the correct placement of the rest of the ligand. The per-mode selected model, geometry rule and RMSD are recorded in `holdout_rerank*.csv`.

## Reproduction

Dependencies recorded for these tests include Python 3.12, RDKit 2026.9.2, Meeko 0.8.0, PDBFixer 1.12.0, OpenMM 8.6.1, Gemmi 0.7.5, NumPy 2.3.5 and SciPy 1.17.0. The preparation checks SHA-256 for three RCSB PDBs and two CCD files. The independently downloaded instance-SDF hashes are printed by the verification script.

```bash
work=/tmp/fastdock-ca-holdout
python benchmarks/anchor_search/holdout.py prepare "$work"
for arm in native anchor shell; do
  python benchmarks/anchor_search/holdout.py "$arm" "$work" --binary build/vina
done
python benchmarks/anchor_search/holdout.py evaluate "$work"
python benchmarks/anchor_search/holdout_verify.py "$work"
python benchmarks/anchor_search/holdout_check_linkage.py
```

The per-mode RMSDs, scores and Zn geometry are in `holdout_modes.csv`, experimental-condition summaries in `holdout_cells.csv` and the within-CA-I rank decisions in `holdout_pairs.csv`. For the smaller-box sensitivity use a separate directory and add `--box-size 20` at preparation.
