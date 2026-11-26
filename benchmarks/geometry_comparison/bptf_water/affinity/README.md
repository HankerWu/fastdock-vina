# BPTF water-site controls and experimental affinity ranking

Six ligands were docked independently into two known BPTF binding sites, 7F5D and 7RWO. Source molecular graphs and independently generated ETKDGv3/MMFF conformers were used. None of the graphs matches either water donor ligand (7K6R/GSK4027, 7M2E/compound 19) or either receiver's crystallographic ligand. The four BROMOscan analogs are related to the donor scaffolds. All runs use the same analytic geometry refinement, Late ranking, a 24 Å cube, three trajectories and one CPU. Three matched seeds per budget give one median top score per ligand and receptor. Protein and water models are fixed throughout each run.

| Matched Kd series | Compounds and concentration |
|---|---|
| Zahid et al., Figure 5B, BROMOscan | BZ1 6.3 nM, 21 9.3 nM, 22 70 nM, 24 200 nM |
| Ycas et al., Table 1, **unlabeled BPTF** SPR | TP-238 200 ± 10 nM, bromosporine 9,000 ± 1,000 nM |

The sources and fingerprint of each input are recorded in identities.json. Ycas gives separate values for 5FW-labeled protein and other assay formats. The four BROMOscan values were in the same reported panel. For compounds 21, 22 and 24 the figure reports one experiment with two technical replicates. Each assay contributes only its within-assay comparisons. No IC50-to-Kd conversion or cross-assay ranking was used. [Zahid et al.](https://doi.org/10.1021/acs.jmedchem.1c01294) · [Ycas et al.](https://doi.org/10.1039/D0OB00506A)

**Identity discrepancy for 24.** ChEMBL record CHEMBL4863631 contains a nitrile, whereas the original synthesis for 24 specifies prop-2-yn-1-yl and the paper identifies an alkyne as the click handle. The calculations use the published terminal-alkyne structure. The manuscript's observed HRMS value also appears transposed relative to its calculated value. Every summary therefore includes the BROMOscan series with compound 24 omitted. BZ1 has a primary aminoethyl substituent. Its related water-donor ligand 19 has a dimethylaminoethyl substituent. The chemically plausible cation at one aliphatic amine was fixed for TP-238 and the four BROMOscan compounds. Bromosporine retained its uncharged input form. No microstate populations or differential uptake free energies were fitted. [Paper and experimental section](https://pmc.ncbi.nlm.nih.gov/articles/PMC9012132/) · [Deposited donor ligand YOV](https://www.rcsb.org/ligand/YOV)

## Initial tests

The table shows the number of correctly ordered pairs with experimental Kd differing at least threefold. Each BROMOscan row has four eligible pairs. Omitting 24 gives two eligible pairs with the same numerator, since 24 was wrongly ranked against both stronger molecules in all initial conditions. TP-238 and bromosporine contribute one separate SPR pair. Every arm ranked this 45-fold pair incorrectly in both receptors at both budgets.

| Receptor | Per-trajectory cap | Dry | Steric sites | Polar sites |
|---|---:|---:|---:|---:|
| 7F5D | 48,000 | 1/4 | 0/4 | 2/4 |
| 7RWO | 48,000 | 2/4 | 2/4 | 2/4 |
| 7F5D | 120,000 | 2/4 | 0/4 | 1/4 |
| 7RWO | 120,000 | 0/4 | 2/4 | 2/4 |

Initial polar sites are the recurrent donor coordinates defined by the parent experiment. The steric arm retains sites as nondirectional excluded volume. Coordinates and site-selection rules are held fixed before comparing scores. The receptor- and seed-dependent rankings provide no reproducible gain in Kd ordering. Repeated runs on these six ligands do not provide evidence of generalization to additional targets.

## Fresh-seed confirmation

With receiver-water support and registration defined in the parent experiment, the six compounds were docked again for all four arms at cap 120,000 with three new seeds (95701, 95713, 95717). Supported retains eligible polar sites at their mean donor coordinates. Registered moves the eligible sites to their observed receiver coordinates. The 144 runs freeze the same ligands and use identical search budgets in every arm.

| Endpoint | Receptor | Eligible pairs | Dry | Polar | Supported | Registered |
|---|---|---:|---:|---:|---:|---:|
| BROMOscan | 7F5D | 4 | 2 | 2 | 2 | 2 |
| BROMOscan | 7RWO | 4 | 2 | 2 | 2 | 2 |
| Unlabeled SPR | 7F5D | 1 | 0 | 0 | 0 | 0 |
| Unlabeled SPR | 7RWO | 1 | 1 | 0 | 0 | 1 |

In 7RWO, registration restored the dry arm's ordering of TP-238 and bromosporine. The weak-minus-strong margins were +0.425, −2.472, −2.472 and +0.101 kcal/mol. In the BROMOscan sensitivity analysis excluding 24, the counts for 7F5D became 2, 0, 2 and 2 of two, while 7RWO gave 2 of two in every arm. The fresh seeds test search variability; the molecules and receptor structures remain the same.

## Reproduction

The evidence archives include raw Vina CLI logs, emitted poses, prepared ligands and all receptor arms used. From the repository root, first run the commands below for the original budget and seed sets. Pass the corresponding empty work directory when repeating the follow-up.

```bash
B=benchmarks/geometry_comparison/bptf_water/affinity
W=$(mktemp -d)
tar -xJf "$B/evidence.tar.xz" -C "$W"
python "$B/experiment.py" "$W"
python "$B/audit.py" "$W"
# To repeat the docking itself instead of replaying the frozen scores:
python "$B/experiment.py" "$W" --binary /absolute/path/to/vina --receptors "$W/receptors"
```

Ligand generation can be repeated by adding --prepare; the method reuses the kinase-state affinity ligand preparer. The preparation used RDKit 2026.09.2 and Meeko 0.8.0. Source chemical graphs (before protonation) and postulated microstate rules are shared across all arms. The displayed unit is nM for measured affinity and kcal/mol for predicted score margins. Positive weak-minus-strong scores in pairs.csv indicate the correct order.

For fresh-seed confirmation, extract confirmation/evidence.tar.xz in an empty directory and add --confirmation to experiment.py. The audit command uses the corresponding cells.csv and summary.csv written by that analysis.
