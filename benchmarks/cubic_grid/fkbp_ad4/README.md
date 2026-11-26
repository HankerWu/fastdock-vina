# Official Vina and AD4 applied to a matched FKBP12 kinetic series

The native Vina and AD4 modes of the same official AutoDock Vina binary correctly ordered 7 and 6 of 15 unique, threefold-separated FKBP12 apparent inhibition constants. Using AD4 to select among the nine poses from each native Vina search gave 3/15. The eight compounds in Holt et al. Table II share a pipecolate ester scaffold. Both receptor conformations are human FKBP12. The three seeds (131071, 196613, 262147) were fixed before evaluation and differ from the two prior cubic-grid cohorts.

| Search and reported score | Correct pairs | Spearman ρ, seven point measurements | Correct matched epimer orders |
|---|---:|---:|---:|
| Official Vina search and Vina score | 7/15 | 0.286 | 1/2 |
| Official AD4 search and AD4 score | 6/15 | −0.179 | 2/2 |
| Official Vina poses, lowest independent AD4 score | 3/15 | −0.071 | 1/2 |

Changing from native Vina to AD4 gained two pair orders and lost three. AD4 recovered both 8 > 10 and 9 > 11, where “>” denotes stronger measured inhibition. Their predicted score differences favoring the stronger compound were +0.421 and +0.193 kcal/mol. Vina gave −0.060 and +0.195. The two epimer pairs have identical TORSDOF within each pair. On this series, explicit scoring terms were insufficient to yield a consistent benefit. The crossover evaluates existing Vina candidates at their written coordinates, with no AD4 minimization.

## Physical and experimental controls

Holt et al. measured same-assay apparent Ki via recombinant human FKBP12 and a chymotrypsin-coupled PPIase assay at 10 °C, pH 7.8. Compound 10 retained its published 300–600 nM interval. A pair qualified when the lowest possible weaker/stronger ratio was at least three. A prediction was the mean of the two receptor-specific medians over three seed-level top scores. The 15 pairs share eight compounds and provide one scaffold-level observation for transfer. Absolute values across scoring functions were not compared. [Holt et al., JACS 1993](https://doi.org/10.1021/ja00075a008) · [1FKG](https://www.rcsb.org/structure/1FKG) · [1FKH](https://www.rcsb.org/structure/1FKH)

Both arms used the same prepared 3D ligand and receptor PDBQT inputs, grid center rounded to 0.001 Å, 24 Å box, 0.375 Å spacing, one CPU, three trajectories, 24,000 nominal evaluation cap per trajectory, nine output modes, and identical seeds for each condition. This was a matched-cap comparison. AD4 used AutoGrid 4.2.6 (Debian 4.2.6-10) with the compiled default atom parameters, 0.5 Å smoothing, standard distance-dependent dielectric (-42), four ligand atom-type maps (A, C, N, OA), an electrostatic map and a charge-dependent desolvation map. All prepared compounds were neutral, and both receptors and ligands had nonzero partial charges. Vina computed its own empirical maps and used its default explicit refinement. AD4 used its own map scoring and torsion treatment. Median docking wall times per run under concurrent execution were 1.02 s for Vina and 1.87 s for AD4; the reusable AD4 maps were generated separately. [Official Vina AD4 protocol](https://autodock-vina.readthedocs.io/en/latest/docking_basic.html) · [Vina 1.2.0 methods paper](https://doi.org/10.1021/acs.jcim.1c00203)

All 63 C++ source and header files compiled from the official source checkout were verified against official Git object IDs at commit `bb64e5e26b520569a86360d69248069b5e665b3f`. The historical version string was supplied by the build helper and still contains 3c65c0b. Binary checksums, map checksums and exact input identities are in `results/metadata.json`. Replaying the first seed using the FastDock fork's normal AD4 mode yielded byte-identical complete PDBQT files in all 16 ligand–receptor cells. [Official comparison commit](https://github.com/forlilab/AutoDock-Vina/tree/bb64e5e26b520569a86360d69248069b5e665b3f)

## Coordinate-based independent audit

There were 96 docking runs (48 per arm) and 864 written poses. SciPy interpolated the six AutoGrid maps from the written coordinates of every pose independently of the search binary. For rigid FKBP12 the bound AD4 intramolecular energy cancels against its chosen unbound reference. The total therefore consists of each atom's typed affinity map, the electrostatic map multiplied by its signed charge, the desolvation map multiplied by its absolute charge, and 0.2983 × ligand TORSDOF. Each map value undergoes the same positive-energy saturation at 1000 as in upstream Vina. The typed map already contains charge-independent desolvation. All 432 AD4 scores reproduced within 0.00285 kcal/mol from coordinates and energies printed to three decimals. Model counts and first-mode scores also agree with all 96 raw logs/pools, with explicit checks for the engine and seed.

An exploratory reaggregation of the same AD4 poses without the ligand-specific TORSDOF term gave 10/15. This constant cannot change the pose order within a molecule or the matched epimer orders; it can change the cross-molecule comparison because the ligands span 10–13 TORSDOF. The result is sensitive to this crude flexibility correction. Removing a fitted physical term based on this single series would need evaluation on unseen scaffolds and targets.

## Reproduction

The existing human panel provides identical prepared inputs and measured labels. The additional archives contain AD4 maps and grid-generation logs, all 96 output pools and logs, and score audits. The published paper and prepared ligand/receptor inputs are not duplicated in these archives. For a replay from the repository root, set `official_vina` and `autogrid` to the two binaries and use Python with NumPy and SciPy.

```bash
src=benchmarks/cubic_grid/human_panel
here=benchmarks/cubic_grid/fkbp_ad4
workdir=$(mktemp -d)
tar -xJf "$src/prepared.tar.xz" -C "$workdir"
tar -xJf "$src/affinity-pools.tar.xz" -C "$workdir" --wildcards --no-wildcards-match-slash 'affinity/*.pdbqt'
python "$here/run_official.py" "$workdir" "$workdir/affinity" "$workdir/official" "$autogrid" "$official_vina" --jobs 4
python "$here/analyze.py" "$workdir/official" "$src/affinity/ligands.csv" "$workdir/computed"
python "$here/check_branch.py" "$workdir" "$workdir/affinity" "$workdir/official" build/vina "$workdir/replay.json"
```

For an offline recomputation using archived poses and maps, extract the three archives into one working directory, copy `results/cells.csv` there, then pass that directory to `analyze.py` as its first argument.
