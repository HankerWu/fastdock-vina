# Acyclic amide torsions in ligand docking

Five three-tailed hCA II inhibitors each encode a tertiary carboxamide C–N as a movable PDBQT branch. All five deposited O–C–N–substituent torsions are within 2.1° of 0° or 180°. The unrestrained Vina potential permits substantially twisted poses. In the new continuous-search pools, 505/699 candidates and 16/20 native winners deviate by over 30° from the nearest planar well. The crystal angles for the same distinguished substituent in 6WQ4/5/7/8/9 are +2.01°, −179.85°, +1.83°, −179.44° and −179.52°. Geometry starts in the opposite well from the crystal in three of the five. No ligand crystal coordinate enters either the generation of alternate inputs or the ranking formula. The pre-existing independent same-target scaffold supplies the search anchor.

The chemically motivated test adds B sin²(ω) to the native score on fixed poses. A π system with maximum orbital overlap in the planar geometry motivates a twofold torsion profile. B is an effective test amplitude, in kcal/mol, chosen independently of the affinity labels. Experimental rotation barriers in simple dimethylamides lie around 15–19 kcal/mol and depend on the substituents, solvent and temperature. A kinetic barrier from those molecules cannot establish the equilibrium torsional potential of an individual solvated protein-bound amide. We compare 0, 5 and 15. Both planar orientations receive equal energy. Ligand internal relaxation, nitrogen pyramidalization, non-identical free rotamer populations and amide conjugation already captured implicitly by Vina remain sources of error. The protocol leaves the default score intact. [Dynamic NMR measurements](https://doi.org/10.1246/nikkashi.1975.579) · [Spin-saturation transfer measurement](https://doi.org/10.1021/ed074p978) · [hCA II experiments](https://doi.org/10.1021/acs.jmedchem.0c00733)

## Results

| hCA II panel and endpoint | Native | B=5 | B=15 |
|---|---:|---:|---:|
| Historical two-state candidates, joint top / 20 | 0 | 6 | 6 |
| Historical candidates, cross-docking only / 16 | 0 | 4 | 4 |
| Historical same-assay Ki edges / 7 | 3 | 3 | 2 |
| Fresh equal-budget continuous search, joint top / 20 | 2 | 6 | 8 |
| Fresh continuous search, cross-docking only / 16 | 2 | 4 | 5 |
| Fresh same-assay Ki edges / 7 | 3 | 6 | 5 |

A joint hit requires native heavy-atom RMSD ≤2 Å in the registered receptor frame and the deposited local amine geometry. Each 20-cell panel contains five ligands, two receivers and two seeds. The historical and new panels overlap in chemistry, receptor structures, priors and some stochastic paths. They are sensitivity checks on the same series. With the 15-unit penalty in the new panel, 6WQ8 contributes four additional hits; 6WQ4 and 6WQ5 each contribute one. In the continuous pools 6WQ7 and 6WQ9 still contain zero joint hits. The new and historical candidate ceilings are 12/20 and 10/20.

A displaced-control penalty of the same strength on the neighboring O–C–CH₂–N⁺ torsion gave 4/20 and 2/20 hits in the new and historical panels. Rewarding perpendicular amides through 15 cos²(ω) gave 0/20 and 1/20. The historical affinity signs, including the better signs from the adjacent control, do not support adopting the amide term as a general affinity correction. Full cell and control outputs are under [results](results).

Each new hCA II comparison has four independent Vina processes per receptor/ligand/seed: two solution-declared geometries of the protonated amine, each with two one-trajectory starts and at most 9 returned modes per start. The native input, its nearest planar rotamer and the opposite planar rotamer have identical atom identities and covalent bonds. For the latter two, the C–N branch is removed, retaining the descendant branches, coordinates outside the rotated side, and the graph-declared stereochemistry. Aligning the chosen O–C–N–X torsion to 0 or 180° preserves the prepared nitrogen geometry. It does not flatten every substituent on a pyramidal nitrogen. In the two mixed arms each amine geometry gets both amides, with the seed-to-amide assignment swapped between arms. The union within each arm has equal search budget. Every returned candidate is transformed back into the original PDBQT representation using only its Cartesian coordinates, then rescored on the same explicit receptor, ligand graph, torsion factor, intramolecular pairs and native baseline unbound reference for that ligand. The open rotor gave 2/20 native top and a 12/20 ceiling. Equal-budget two-well unions gave 1/20 top and an 11/20 ceiling in either assignment, without a torsion penalty. Locking out twisted intermediates during search lost useful candidates.

BRD4 4OGI is a structural control with a secondary amide that was already fixed in the input tree. The continuous/input label for this ligand denotes repeated use of that original fixed rotamer. None of three receivers × two seeds had a whole-ligand ≤2 Å candidate in the original or alternative pools. FXa tests all 28 members of the existing same-assay Ki series on both 1FJS and 2P16, including the previously repaired CHEMBL369480 geometry. Three ligands have the eligible movable amide. The 25 common controls reuse the same new native pools across all arms. The 3 affected ligands have experimental Ki of 300, 370 and 417 nM, generating no pair with a threefold separation among themselves. The 66 ranking edges touching these weak ligands went from 62/66 to 66/66 with a 5- or 15-unit penalty in the original pools. Across all 295 ≥3-fold edges that is 231 to 235. The 4-point gain depends on affected-to-control comparisons. It does not establish the relative order of the three affected compounds. [FXa same-assay series](https://doi.org/10.1016/S0960-894X(00)00139-6)

## Running the tools

The endpoint takes one prepared Meeko ligand and one matching native Vina pool. It checks the chemical records, torsion tree, atom identities and charges. Its native energies come from the printed VINA RESULT records, so use it on one comparable pool at a time. To pool across preparations, use the common-topology rescoring protocol above. The penalty sums over movable, acyclic, neutral, distinguishable carboxamides. For each nitrogen the selected carbon substituent is identified by graph symmetry rank; equivalent arms are skipped. Primary amides, carbamates, ureas, lactams and amides with a symmetric pair of N substituents are outside the present protocol. No automatic change is made to input protonation states.

```sh
python tools/rank_amide_pool.py PREPARED.pdbqt POOL.pdbqt --barrier 15
python tools/amide_rotamers.py PREPARED.pdbqt
python tools/amide_rotamers.py PREPARED.pdbqt --bond C_SERIAL N_SERIAL --lock --planar --keep-orientation --output near.pdbqt
python tools/amide_rotamers.py PREPARED.pdbqt --bond C_SERIAL N_SERIAL --lock --planar --output far.pdbqt
```

Use the source branch's standard Python preparation dependencies (NumPy, RDKit, Meeko and Gemmi). Compile the CLI with the branch build script. The public fixture includes the inputs for this study. See [fixture contents and reproduction](fixtures.md).

```sh
python benchmarks/amide_search/run.py "$BIN" "$BASE" "$WORK" --jobs 4
python benchmarks/amide_search/evaluate.py "$BIN" "$BASE" "$WORK" "$RESULTS" --jobs 4
python benchmarks/amide_search/archive.py "$BASE" "$RESULTS"
python benchmarks/amide_search/controls.py "$BASE" "$WORK" "$RESULTS"
```

The 548 successful searches produced 4,665 models. An independent check recovers the whole-ligand structural RMSDs from experimental coordinates and residue registrations to within 1.72×10⁻⁷ Å. The ordinary autoboxed fixed-pose scoring route independently reproduces all 4,665 explicit energies within its 0.001 kcal/mol print precision. The independently recomputed top-model selections, candidate ceilings and unique affinity-edge counts match the tabulated results.
