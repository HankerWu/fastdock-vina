# Paired affinity ranking under geometric BFGS preconditioning

Armijo and analytic geometry use the same Vina potential, search box, input and random seed. Full Late docking changes only the final explicit-receptor optimizer. Every search trajectory counter, hunt-chain record, starting point, prescreen count and number of inputs to explicit refinement matched across arms. These checks cover 280 paired comparisons and 560 full searches over both budgets. The geometry mode used 98,893 explicit score evaluations in total, versus 62,384 for Armijo. These counters exclude grid search and map construction.

Each receptor and assay has a separate ordering. For each compound, take the median rank-one Vina energy across the two fixed seeds at that budget. Report concordance over exact measured concentrations differing by at least threefold, with lower energy predicting lower Ki or Kd. A score tie receives half credit. CSVs also give all unequal experimental pairs, tenfold pairs and within-stratum Spearman correlations. Seed repeats and shared ligands/receptors make pair counts dependent.

| Target and assay | Receptor | Ligands | Eligible pairs | 12k Armijo / geometry | 60k Armijo / geometry |
|---|---|---:|---:|---:|---:|
| MAPK14, Davis Kd | 3GCP | 9 | 29 | 15 / 12 | 13 / 13 |
| MAPK14, Davis Kd | 1KV2 | 9 | 29 | 19 / 19 | 6 / 6 |
| ALK, Davis Kd | 7BTT | 8 | 19 | 17 / 16 | 15 / 15 |
| RET, Davis Kd | 7DUA | 8 | 19 | 8 / 8 | 14 / 14 |
| TNK2, Davis Kd | 7KP6 | 8 | 19 | 12 / 12 | 13 / 12 |
| BACE1, CHEMBL2015090 Ki | 4DJW | 6 | 14 | 2 / 2 | 2 / 2 |
| BACE1, CHEMBL2015090 Ki | 4GID | 6 | 14 | 4 / 4 | 12 / 11 |
| BACE1, CHEMBL2155063 Ki | 4DJW | 8 | 22 | 2 / 1 | 2 / 2 |
| BACE1, CHEMBL2155063 Ki | 4GID | 8 | 22 | 14 / 13 | 9 / 9 |

At 12k, geometry recovered one of the pairs missed by Armijo and lost seven. At 60k, it lost two pairs. The budgets used different seeds. Both high-budget flips involved small score margins. For bosutinib (29 nM) and sunitinib (8,900 nM) on TNK2, the predicted weak-minus-strong margin moved from +0.0085 to −0.0065 kcal/mol. For BACE1 on 4DJW, 57 nM CHEMBL2151138 scored worse than 9,500 nM CHEMBL2151136 with both optimizers at both budgets. The high-budget medians were −8.4555 and −10.0885 kcal/mol under Armijo, and −8.5025 and −10.1730 under geometry. Relative affinity accuracy remained strongly dependent on the receptor and assay series.

## Sources, crystal structures and chemistry

- Davis et al. supplementary Table 1 identifies unmutated constructs and Table 4 gives measured Kd in nM. Empty assay cells and uncertain stereochemistry were excluded. For MAPK14, exact concentrations were available for ten agents. SB-203580 has unspecified sulfoxide stereochemistry and was excluded. ALK, RET and TNK2 use the first eight eligible numeric PubChem IDs per target from the recorded identity dictionary. These three targets were selected as an extension after running the initial MAPK14/BACE1 panel. The choices used WT human structure quality and label availability. No selections used docking predictions. [Primary experiment and tables](https://www.nature.com/articles/nbt.1990)
- Cumming et al. and Stamford et al. describe two separate human BACE1 inhibitor series. The exact nM Ki rows come from the frozen ChEMBL records for assays CHEMBL2015090 and CHEMBL2155063, identified by document, endpoint, assay and species. The first assay left six eligible ligands after filtering stereochemistry. The second used the first eight eligible numeric ChEMBL IDs. The first assay describes a 30-minute preincubation, quenched peptide substrate and a 90-minute FRET read. Neither time-series assays nor measurements across studies were pooled. [Cumming et al.](https://doi.org/10.1016/j.bmcl.2012.02.013) · [Stamford et al.](https://doi.org/10.1021/ml3001165)
- In the extension, the nearest omitted nonwater heavy atom to the crystallographic ligand in the asymmetric unit was 13.21 Å for RET and 10.36 Å for TNK2. ALK chain A had only its own ligand and water. BRAF 7K0V was excluded before docking because an omitted chain contacted the ligand at 6.54 Å. No reconstructed heavy atom lay in the three accepted pockets. ALK had three truncated, incomplete distal side chains more than 18 Å from the bound ligand. The crystallographic pocket center determines the box. The model omits waters, symmetry mates and nonlocal domains. It samples rigid deposited receptors. [ALK 7BTT](https://www.rcsb.org/structure/7BTT) · [RET 7DUA](https://www.rcsb.org/structure/7DUA) · [TNK2 7KP6](https://www.rcsb.org/structure/7KP6)
- MAPK14 uses the repository's existing ring-ensemble receptor fixtures. Both BACE1 chain A receptors had complete recorded side-chain coordinates. Distances to the closest omitted nonwater atom were 17.07 Å in 4DJW and 9.96 Å in 4GID. Each BACE1 ligand used one fixed cationic cyclic guanidinium representation. Standard Meeko protein typing was shared between arms. Alternative ligand microstates, catalytic aspartate protonation, linked proton uptake and binding-site water populations were not enumerated. The 47 assay-compound rows cover 32 distinct compounds and 27 distinct Bemis–Murcko scaffolds. Some kinase compounds appear across targets. Ligand starts use ETKDGv3 plus MMFF from the deposited chemical identities, without crystallographic pose coordinates. [4DJW](https://www.rcsb.org/structure/4DJW) · [4GID](https://www.rcsb.org/structure/4GID)

## Reproduction

With Python 3.12, RDKit 2026.09.2, Meeko 0.8.0, OpenMM 8.6.1, PDBFixer 1.12, xlrd 2.0.1, SciPy 1.17.0 and NumPy 2.3.5, from the repository root:

```bash
B=benchmarks/geometry_comparison/affinity
W=$(mktemp -d)
H=$(mktemp -d)
V=/absolute/path/to/analytic-build/vina
(cd "$B" && sha256sum -c SHA256SUMS)
tar -xJf "$B/inputs.tar.xz" -C "$W"
cp "$B/panel.csv" "$W/"
# Optional reconstruction of ligands, receptors, panel and exclusions from raw inputs.
OPENMM_CPU_THREADS=1 python "$B/prepare_affinity.py" "$W"
cp "$W/panel.csv" "$H/"
ln -s "$W/prepared" "$H/prepared"
python "$B/run.py" "$W" "$V" --cap 12000 --seeds 89513 89527
python "$B/run.py" "$H" "$V" --cap 60000 --seeds 89609 89611
python "$B/summarize.py" "$W" --seeds 89513 89527
python "$B/summarize.py" "$H" --seeds 89609 89611
python "$B/audit.py" "$W" "$W/raw" 89513 89527
python "$B/audit.py" "$H" "$W/raw" 89609 89611
```

The release includes all 560 raw PDBQT outputs and CLI logs, predicted scores, pair decisions, included and excluded labels, compiled-run inputs and frozen source data. For offline verification, replace the two docking commands with extraction of the corresponding cap's outputs.tar.xz into its work directory, then run both summarizers and audits. The audit reads the source spreadsheet and two ChEMBL records directly and recomputes every reported concordance count from the raw PDBQT energies.
