# Screening intraligand electrostatics on fixed native pools

## Frozen design

Use the native Late pools from the previous cutoff experiment at seeds 271828 and 314159. Keep all shared-receptor T4L and hCAII cells and all BRD4 cells. BRD4 includes extra cognate receptors for pose success only. Use the published ligands and their original prepared microstates. Compare model selection at symmetry-corrected RMSD ≤2 Å and same-assay order for pairs separated by at least threefold, within receptor and seed. Repeated appearances of one compound pair are dependent.

A screened internal potential might reduce the bias toward compact conformations of charged ligands. A separate possibility is that the energy published as strain mixes the penalty for distorting a pose within its local minimum with the choice of the conformer against which other molecules are compared. Test both by retaining the repository's 0.1 replacement weight, Vina flexibility normalization, Meeko atom mapping and hydrogen-only minimization. The primary new arm changes the dielectric from 1 to a fixed 4 in both bound energies and the common eight-conformer reference. The basin sensitivity changes the reference to the free MMFF94s minimum reached by minimizing the very same pose, keeping dielectric 4. This local difference estimates relaxation cost and omits costs of occupying a rare basin. A native-interaction-only control drops the native intraligand difference without replacing it. The historical gas-phase arm and original native score complete the five arms. There is no fit, search, pooling of experimental concentrations across targets or choice of parameters from BRD4 outcomes. Scores are model units. A constant dielectric does not represent explicit solvent.

The test is retrospective and every target was examined in earlier cutoff experiments. BRD4 is a transfer evaluation of the frozen strain hypothesis, not a blinded external benchmark. Refinement operates on private copies for ligand internal energy evaluation only. All geometric outcomes use the originally saved receptor-frame coordinates. Exact same-pool tests cannot recover missing candidates.

## Matched results

Each row gives the number of top-ranked poses within 2 Å and the credit for correctly ordered, at least threefold-spaced same-assay pairs. Ties receive 0.5. Pools contain 304 T4L, 179 hCAII and 249 BRD4 poses. The three panels cover 18 distinct compounds. Reranking leaves all coordinates fixed.

| Score | T4L top1 / Kd | hCAII top1 / Ki | BRD4 top1 / Kd |
|---|---|---|---|
| Native | 28/36 · 50/76 | 0/20 · 14/28 | 3/28 · 10/30 |
| Inter only | 28/36 · 51/76 | 0/20 · 12/28 | 3/28 · 10/30 |
| Vacuum, global reference | 28/36 · 28/76 | 0/20 · 16/28 | 3/28 · 6/30 |
| Screened, global reference | 28/36 · 28/76 | 0/20 · 13/28 | 3/28 · 7/30 |
| Screened, local reference | 26/36 · 30/76 | 0/20 · 15/28 | 3/28 · 11/30 |
| Exploratory torsion-only addition | 28/36 · 45/76 | 0/20 · 10/28 | 3/28 · 7/30 |

The global reference is the minimum among eight independently embedded conformers for each molecule. The screened global and local trials keep identical frozen-pose energy at dielectric 4. The reference subtraction changes across-compound ordering. Values below a finite-sampling global reference are retained. The number of independent ligand pairs is 19, 7 and 5, respectively. BRD4 outcomes include four additional cognate cells beyond the 24 shared-receptor cells used for rank comparisons. These three ligands sets and their assays remain small and reuse the same ligands.

In a separate hindsight comparison we added an independently native-refined crystal-initialized pose to each matched hCAII cognate decoy pool (seeds 42 and 1729, nine native decoys per pool). Same-ligand reference and normalization constants cancel. Explicit graph scoring reproduced all ten historical native gaps within 0.000002 raw energy units. All ten native selections favored decoys. Each of the three MMFF variants favored the crystal-like pose for compound 34/6WQ4 in seed 42 alone. Other seed/compound combinations retained a lower-scoring decoy. This diagnoses scoring; the crystal-derived poses are never candidates in the main 84 cells.

The 22 T4L affinity-order reversals under screened global strain all involve the two strongest ligands. Six reversals involve isobutylbenzene and 16 involve n-butylbenzene. Recomputing the seven MMFF contributions from their own minimum-energy conformer to the two 186L native n-butylbenzene top poses found 10.61 and 11.91 kcal/mol total strain. Intramolecular van der Waals accounts for 7.13 and 7.94 kcal/mol, while the respective torsional terms contribute 0.63 and 0.97. These poses have 0.87 Å crystal RMSD. Applying exactly the same heavy-coordinate procedure to the deposited 186L ligand gave 4.12 kcal/mol. Crystal coordinates have experimental and model uncertainty. This comparison illustrates the sensitivity of all-atom strain to fine geometry of docked poses refined under another potential. It does not establish a gas-phase force-field energy as the real binding cost.

Following this diagnosis, we tested adding only the MMFF torsional contribution at the same fixed weight to the original complete native score. Each torsional reference comes from the same lowest *total*-energy vacuum conformer. This is a post-hoc sensitivity on the same pools; it reduces rank concordance on all three targets. The results support keeping these terms outside default scoring. Related original studies discuss bound-state strain and its dependence on coordinates and reference conformers. [Perola and Charifson, 2004](https://doi.org/10.1021/jm030563w) · [Gu et al., 2021](https://doi.org/10.1021/acs.jcim.1c00368)

## Reproduction

Extract the preceding experiment's FastDock-Vina-Pro-results.zip and provide the extracted directory as the first argument. Use the environment documented in the cutoff experiment. The scripts recompute energy from original coordinate files and verify the archived native pool hashes. Column names eps1 and eps4 indicate dielectric constants 1 and 4. Kd is in µM for T4L/BRD4 and Ki in nM for hCAII.

```bash
python benchmarks/strain_transfer/evaluate.py EVIDENCE build/mixed/vina results --workers 2
python benchmarks/strain_transfer/crystal_diagnostic.py EVIDENCE build/mixed/vina results/crystal_contrasts.csv
python benchmarks/strain_transfer/alkylbenzene_diagnostic.py EVIDENCE results/alkylbenzene_terms.csv
python benchmarks/strain_transfer/torsion_diagnostic.py EVIDENCE results
```
