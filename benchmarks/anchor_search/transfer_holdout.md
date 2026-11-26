# Ligand-excluded scaffold transfer on human BRD4 and the T4L cavity

The 13 query compounds started from the independently prepared ETKDG/MMFF geometries of the preceding matched cutoff experiment. For each query the donor was a different compound with a published co-crystal on the same target. Searches used two previously used seeds and the three BRD4 or two T4L shared receptors. Donor selection used the graph-only maximum common substructure procedure in `transfer_mapping.py`. The Meeko SMILES atom map supplies serials; the donor crystal position comes from the earlier receptor-framework alignment. For graph symmetries the first stored donor crystal graph map is fixed before model evaluation. Similar ligands can bind in different orientations; forcing a single correspondence is a model choice. These are leave-one-ligand-out tests under access to other co-crystals, including analogues within the same chemical series.

The new search and the paired native Late control used four trajectories, 50,000 maximum evaluations per trajectory, nine output modes, one CPU and the original search boxes. The final energies are native. The query crystal is used only to evaluate output and perform the separately saved `mapping_audit.csv` check after docking. A receptor may itself be the query's stripped holo receptor. Every native output SHA-256 matched the corresponding original artifact.

| Shared target panel | Top-1 ≤2 Å | Candidate pool ≤2 Å | Median top-1 RMSD, Å | Same-assay ≥3-fold Kd pairs | Native → transfer grid evaluations |
|---|---|---|---|---|---|
| BRD4, 4 compounds, 3 states, 2 seeds | 1 → 7 / 24 | 2 → 12 / 24 | 6.738 → 6.734 | 10 → 9 / 30 | 4,802,803 → 4,802,491 |
| T4L, 9 compounds, 2 states, 2 seeds | 28 → 28 / 36 | 36 → 35 / 36 | 1.180 → 1.012 | 50 → 49 / 76 | 7,202,403 → 7,203,758 |

The BRD4 gains all occurred on two congeners that transfer 29 of 36 or 37 heavy atoms. TG-101209 reached top 1 in four of six cells, up from one. TG-101348 reached three of six, up from zero. Every transfer pool for these two ligands contained a near-experimental pose. For the 7-atom donor match on SB-202190, the only native pool hit was lost. BI-2536's 13-atom match left zero near-native poses. In the 3MXF receptor frame the later reference audit found symmetry-adjusted donor-to-query scaffold RMSDs of 0.96 Å for either congener, 7.35 Å for SB-202190 and 5.36 Å for BI-2536. The last two have respective scaffold centroid errors of 6.85 and 4.72 Å.

In T4L, two benzofuran native near misses (2.020 and 2.007 Å) became hits at 0.615 and 0.594 Å. Two o-xylene hits became misses (1.973 → 2.044 and 1.058 → 2.809 Å). The o-xylene match covers seven of eight heavy atoms of each ligand. Despite a 0.33 Å centroid error in the 1L90 frame, the matched atoms as a group lie 2.73 Å RMSD from their true positions after allowed whole-ligand symmetries. Substructure coverage by itself cannot guarantee a transferable orientation. Affinity credit uses five distinct BRD4 and 19 distinct T4L molecule pairs, repeated across states and seeds.

An independent calculation used RCSB ligand-instance SDFs, receptor alignments from SciPy, and Meeko/RDKit pose reconstruction with receptor-frame symmetry-aware CalcRMS. Over both arms it checked all 422 BRD4 and 588 T4L poses; maximum absolute discrepancies were 0.00000073 and 0.00000078 Å. Input structures, preparation and assay sources are documented in the cutoff experiment.

## Reproduction

Use the `brd4/` and `t4l/` directories from the previous evidence archive, containing `prepared/` with the manifest and `raw/` with each original PDB and RCSB ligand-instance SDF. From the repository root, substitute `brd4` or `t4l` for PANEL and the corresponding absolute filesystem paths for DATA, OUT and VINA.

```bash
python benchmarks/anchor_search/transfer_holdout.py native PANEL DATA/prepared OUT --binary VINA
python benchmarks/anchor_search/transfer_holdout.py transfer PANEL DATA/prepared OUT --binary VINA
python benchmarks/anchor_search/transfer_holdout.py evaluate PANEL DATA/prepared OUT
python benchmarks/anchor_search/transfer_holdout_verify.py PANEL DATA OUT
```

## Conditional final ranking diagnostic

A further matched calculation reused the final coordinates of the transferred pools. We first kept their released native ranks, then ranked by printed native INTER + INTRA plus the same mean point restraint used in search, and finally ranked by restraint alone. For displacement d in Å, each point contributes 4d² / (sqrt(1+d²) + 1). The means are over transferred atoms in the published mapping. These comparisons use no fitted weight. The composite is a conditional within-ligand pose ranking, with no affinity interpretation.

| Panel on transferred poses | Released top-1 | Composite top-1 | Restraint-only top-1 | Pools with any hit |
|---|---|---|---|---|
| BRD4 | 7/24 | 10/24 | 11/24 | 12/24 |
| T4L | 28/36 | 22/36 | 22/36 | 35/36 |

Reapplying the restraint in BRD4 rescued four cases and lost one. TG-101209 in 3MXF went from 6.875 to 1.185 Å and from 3.905 to 1.142 Å in the two seeds. TG-101348 in 4OGJ at seed 314159 went from 1.841 to 2.180 Å. T4L gained three and lost nine, with all four benzofuran cases changing from 0.59–1.01 to 2.97–3.09 Å. Keeping the donor's single orientation in the final score can therefore undo a successful search. All 493 new guided poses were separately re-scored with the prepared topology and no anchor: maximum absolute deviations from saved native INTER and INTRA were 0.04654 and 0.00749 energy units at the saved coordinate precision. Native score terms have been released from the search prior.

```bash
python benchmarks/anchor_search/transfer_holdout_rerank.py OUT RESULTS
```
