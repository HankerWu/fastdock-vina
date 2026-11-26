# Tricubic grid search on human protein–ligand complexes

Changing the rigid search-map interpolant at 0.375 Å reduced first-ranked pose recovery from 32 to 30 of 66 paired cells. Five first-ranked poses crossed the 2 Å threshold in the favorable direction and seven crossed in the unfavorable direction. Recovery in the complete nine-mode pool stayed at 38 of 66. A favorable first FKBP12 affinity ranking disappeared on three fresh seeds.

## Pose experiment

Each query started from a single ETKDGv3 embedding of its chemical graph, fixed seed 43103, followed by MMFF relaxation and Meeko preparation. Crystal coordinates served to locate the receptor pocket and evaluate predictions. Every target contributed two different co-crystallized ligands. Eleven receptor preparations supported 11 cognate and 11 reciprocal cross-ligand cells. Three seeds (14591, 30553, 75431) gave 66 paired searches per arm. Repeated seeds estimate search variability within a structural case. They supply no additional independent target or scaffold.

| Target | PDB and ligand | PDB and ligand | Native → cubic top-1 | Gain / loss |
|---|---|---|---:|---:|
| EPAS1, HIF-2α | 4GHI / 0X3 | 5TBM / 79A | 9 → 8 of 12 | 1 / 2 |
| FKBP1A, FKBP12 | 1FKG / SB3 | 1FKH / SBX | 0 → 2 of 12 | 2 / 0 |
| RBP4 | 1RBP / RTL | 3FMZ / 2T1 | 4 → 2 of 12 | 0 / 2 |
| FABP4 | 2NNQ / T4B | 3FR2 / 8CA | 6 → 6 of 12 | 2 / 2 |
| DYRK1A | 3ANR / HRM | 7AKL / RK5 | 12 → 12 of 12 | 0 / 0 |
| JAK2 | 6VS3 / R6V | 6WTN / RXT (query only) | 1 → 0 of 6 | 0 / 1 |

Both EPAS1 structures use an engineered R247E variant. The experiment evaluated each deposited chain A, omitting waters and preserving the ligand protonation encoded by the chemical component dictionary. No excluded nonwater atom from the deposited asymmetric unit was within 8 Å of a studied ligand. This check covers the asymmetric unit. It does not enumerate space-group symmetry contacts. In 6WTN, the deposited Lys883 lacks its side chain beyond Cβ, with observed atoms reaching 6.80 Å from ruxolitinib. Attempts to rebuild its side chain produced a spurious covalent bond during PDBQT preparation. The ligand had complete coordinates and remained eligible for cross-docking against 6VS3. The other 11 receptor inputs passed preparation. The closest rebuilt atom in these 11 receptors was 8.58 Å from its ligand.

All searches used one CPU, a 24 Å cubic box centered on the ligand belonging to the receiver, three trajectories capped at 24,000 evaluations each, Late explicit refinement, 15 requested refinement candidates and at most nine output modes. The common cap is a search budget, not a wall-time match. Per-run outputs include logs and elapsed seconds, with concurrently scheduled processes. The cubic arm adds only `--cubic_search_grid`. The default-arm outputs were reproduced byte for byte in a separate 66-run replay with the experimental binary.

Symmetry-aware heavy-atom RMSD was measured in the receptor frame. For cross-docking, corresponding Cα positions in the two known pockets determined the transform. Only atom-name and sequence correspondence, chemical graph, chirality and receptor alignment entered structural evaluation. Cross-docking measures agreement with transfer of a known crystal pose to a different receptor. The actual bound pose of that query in the alternate conformation is unobserved. `verify.py` independently read the raw PDB and CCD with Gemmi, mapped a fresh RDKit graph to the PDBQT SMILES indices, and computed protein transforms with SciPy. Agreement on all 1,188 modes was within 1.78 × 10⁻¹⁴ Å. Each output contained nine modes.

| Set | Paired cells | Top-1 ≤2 Å | Top-5 ≤2 Å | Any returned mode ≤2 Å |
|---|---:|---:|---:|---:|
| Cognate | 33 | 18 → 16 | 20 → 19 | 21 → 19 |
| Cross-ligand | 33 | 14 → 14 | 17 → 17 | 17 → 19 |
| All | 66 | 32 → 30 | 37 → 36 | 38 → 38 |

Two FKBP12 first-ranked successes coexist with a worsening per-target median from 4.80 to 6.19 Å. On RBP4 the median improved from 5.61 to 2.45 Å while first-ranked hits fell from four to two. These targets illustrate the sensitivity of binary success to the error distribution around the threshold. The six targets cover only 12 molecules. Within each target a particular receptor shape and ligand series can still dominate the results.

## Same-assay FKBP12 ranking

Holt and colleagues measured apparent inhibition constants in recombinant human FKBP12 with a chymotrypsin-coupled PPIase assay. The standard assay used 35 mM HEPES at pH 7.8, 10 °C, 25 µM succinyl-Ala-Phe-Pro-Phe-4-nitroanilide and 0.5% ethanol. This benchmark includes the complete eight-compound Table II, using the structures in the synthesis section. The experimental values, reported errors and explicit stereochemical graphs appear in `affinity/ligands.csv`. The 8/9 graphs equal the deposited SB3/SBX descriptors, and 10/11 are the corresponding carbinol epimers. All eight compounds share a pipecolate scaffold. [Holt et al., JACS 1993, pp. 9927 and 9933–9936](https://doi.org/10.1021/ja00075a008)

The affinity protocol used the same binary settings and both receptor structures for all molecules. For each molecule we first took the median top score across seeds within each receiver, then averaged the two medians. This prespecified consensus treats receptor contributions equally and lacks a thermodynamic population calibration. If the minimum possible ratio of the weaker to the stronger measurement was at least three, each distinct chemical pair contributed one comparison. The 300–600 nM range for compound 10 was kept as an interval. The 15 qualifying pairs all compare one of 6, 8 and 9 against a member of the other five. They share compounds and therefore are statistically dependent.

| FKBP12 seeds | Paired cell count | Correct orders, native → cubic | Changed orders, gain / loss | Spearman ρ on seven point values |
|---|---:|---:|---:|---:|
| Primary three | 48 per arm | 4 → 11 of 15 | 7 / 0 | 0.071 → 0.714 |
| Three additional seeds (92461, 148721, 211693) | 48 per arm | 9 → 9 of 15 | 1 / 1 | 0.214 → 0.536 |

Both algorithms got the 9 versus 11 matched epimer order right in both cohorts. For 8 versus 10, cubic alone succeeded in the primary cohort and both succeeded on the additional seeds. Compound 6 introduced multiple methoxy substituents and remained a difficult strong binder. The new-seed replication reuses molecules and receptors and tests only sampling sensitivity. `audit_affinity.py` checked all eight molecular formulae and reconstructed all pair orders directly from serialized scores.

## Reproduction

The study measured the candidate in commit 911ebb2 against the baseline from ce9a89f. After the explicit box derivative was corrected, all 132 structural searches and 192 FKBP12 searches were rerun with identical inputs. All top scores, first-ranked pose RMSDs, per-cell top-five and oracle RMSDs, and all affinity rankings reproduced. Minor differences occurred in secondary modes. The complete before/after SHA-256 inventory and every changed output are in `wall_replay.csv` and `wall-replay-differences.tar.xz`. The archives contain the downloaded PDB/CCD, prepared molecules, mode pools, logs and timings. `sources.csv` records retrieval URLs and SHA-256 hashes. Each receiver's `preparation.json` contains the box center, ligand atom maps, observed references, alignment coordinates and the details of any repaired atoms. For 6WTN the archived log and partial chain document the failed receptor preparation.

The commands below assume this directory is `benchmarks/cubic_grid/human_panel`. Set `repo` to the current checkout and `vina_bin` to its independently built Vina executable. Use Python 3.12 with RDKit 2026.09.2, Meeko 0.8.0, Gemmi 0.7.5, NumPy and SciPy. Preparing receptors from raw coordinates also requires OpenMM 8.6.1 and PDBFixer.

```bash
here="$repo/benchmarks/cubic_grid/human_panel"
workdir="$(mktemp -d)"
for archive in depositions prepared pose-pools affinity-pools affinity_confirmation-pools; do
    tar -xJf "$here/$archive.tar.xz" -C "$workdir"
done
cp "$here"/modes_*.csv "$workdir"
for cohort in affinity affinity_confirmation; do
    cp "$here/$cohort"/*.csv "$workdir/$cohort/"
done
python "$here/verify.py" "$workdir"
python "$here/audit_affinity.py" "$workdir"
python "$here/summarize.py" "$workdir"
OPENMM_CPU_THREADS=1 python "$here/panel.py" prepare --repo "$repo" --work "$workdir" --jobs 3
python "$here/panel.py" native --repo "$repo" --work "$workdir" --binary "$vina_bin"
python "$here/panel.py" cubic --repo "$repo" --work "$workdir" --binary "$vina_bin"
python "$here/affinity.py" --work "$workdir" --binary "$vina_bin"
python "$here/affinity.py" --work "$workdir" --binary "$vina_bin" --confirmation
```
