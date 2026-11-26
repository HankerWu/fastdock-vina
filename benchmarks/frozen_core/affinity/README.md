# Binding-affinity ranking with a contact-dependent rotor penalty

On the frozen five-target panel, replacing the empirical Vina torsion divisor with a fixed cost for two-sided, contact-constrained rotors improved the paired correct count in both budgets. The larger external human Factor Xa series gave a three-pair gain when the original top pose supplied the affinity score and a six-pair loss when the corrected score selected its own pose. These results support an exploratory affinity annotation for the original top pose. They provide insufficient evidence to change the default ranking or search.

## Scorer and physical assumptions

Native Vina evaluates U = E(inter + intra) − E(unbound) and reports U/(1 + 0.05846 N). The unbound reference is shared across candidates in a frozen pool. The engine computes N from the ligand's inferred bonds, heavy-neighbor counts and rotatability flags. Running `native_rotors.cpp` on each prepared input showed three of 32 ligands for which N differs from the declared PDBQT TORSDOF (N=2, 7.5 and 3.5 for 126565, 3062316 and 5287969). The largest discrepancy reconstructing the archived native score using this N and the two printed energy terms is 0.001345 kcal/mol over 4,940 poses. The printed components have precision 0.001 kcal/mol. [Original Vina](https://doi.org/10.1002/jcc.21334) · [Scoring implementation](../../../src/lib/vina.cpp) · [Rotor counting](../../../src/lib/conf_independent.cpp)

The primary corrected affinity score is U + RT ln(3) × K. We use 298.15 K as a reference temperature, giving 0.65091 kcal/mol per constrained rotor. This would be the torsional free-energy cost if each independent free rotor had three equiprobable angular wells and only one of these wells remained accessible upon binding. Real wells differ in occupancy, rotors can be correlated and intramolecular hydrogen bonds may restrict the unbound ligand. Counting contacting rotors is an operational proxy. The output is an empirical energy-like ranking score.

For each declared PDBQT BRANCH, split the ligand into the complete two subgraphs obtained by cutting the axis bond. Consider bonds with at least two heavy atoms on each side. A heavy ligand atom qualifies when its nearest receptor heavy atom is at most 4.0 Å away and its perpendicular distance from the infinite axis line is at least 1.0 Å. K counts bonds with a qualifying atom on both sides. A rotation displaces atoms in proportion to the distance from the axis. Both sides must have off-axis contact because a ligand can rotate as a whole about a singly anchored fragment. The fixed cutoffs and three-state model were specified before calculating corrected outcomes. They are implementation choices. No thresholds were swept, no experimental values select poses, no coordinates change and no new docking searches take place. Continuous receptor relaxation, solvent, proton-linked equilibria and the lifetimes or populations of these contacts remain unresolved.

The two-sided contact criterion has a precedent in ChemScore. Its original scoring model fitted multiple terms on measured complexes, and the GOLD implementation describes polarity weighting and empirical coefficients. This test substitutes the stated off-axis geometry and an analytic three-state penalty. Earlier repository experiments already evaluated a uniform RT ln(3) and an angular partition approximation on other targets. The controls here are native, U alone and U + RT ln(3) × all eligible heavy-atom rotors. The comparison with uniform isolates the added value of contact information. [Eldridge et al.](https://doi.org/10.1023/A:1007996124545) · [CCDC GOLD §7.4.5](https://www.ccdc.cam.ac.uk/media/Documentation/GOLD_User_Guide_2020_1.pdf) · [Earlier rotor study](../../affinity_partition/README.md)

## Frozen five-target panel

The previous panel has 47 assay-compound rows, 32 distinct compounds, 27 Bemis–Murcko scaffolds, seven receptors and 70 distinct receptor–compound combinations. Two archived Late refiners (Armijo and geometry), two numerical seeds and two search caps supply 560 original pools with 4,940 poses. Native refers to the unmodified score on these coordinates. For each compound take its median across both seeds within a receptor. Only same-target, same-assay Ki or Kd ratios of at least three enter pair credit. Higher concentration denotes weaker binding, lower score predicts stronger binding, and a prediction tie receives half credit. The primary analysis uses the original Vina rank one. A secondary comparison selects the minimum of the same corrected score among all archived candidates per seed and then aggregates. The total of 187 receptor-specific pair decisions includes repeated chemistry across receptors. A second summary averages each compound's receptor predictions inside the same assay and counts each of 122 assay-specific ligand pairs once. [Frozen inputs and experimental provenance](../../geometry_comparison/affinity/README.md) · [Davis et al.](https://doi.org/10.1038/nbt.1990) · [Cumming et al.](https://doi.org/10.1016/j.bmcl.2012.02.013) · [Stamford et al.](https://doi.org/10.1021/ml3001165)

| Source of frozen poses | Score selection | 12k native / raw / uniform / contact | 60k native / raw / uniform / contact | 12k / 60k contact after averaging receptor predictions |
|---|---|---:|---:|---:|
| Geometry | Fixed first | 87 / 71 / 95 / 99 | 84 / 61 / 86 / 89 | 69 / 68 of 122 |
| Geometry | Select from pool | 87 / 71 / 95 / 98 | 84 / 61 / 86 / 89 | 66 / 70 of 122 |
| Armijo | Fixed first | 93 / 73 / 95 / 103 | 86 / 61 / 88 / 90 | 71 / 71 of 122 |
| Armijo | Select from pool | 93 / 73 / 95 / 102 | 86 / 61 / 88 / 92 | 69 / 70 of 122 |

After averaging the receptor predictions, the native reference is 63/122 at 12k and 59/122 at 60k for geometry, or 66 and 63 for Armijo. Below each cell gives correct directions for native, fixed-first corrected affinity and corrected-score selection from the full geometry pool.

| Target and assay | Receptor | Pairs | 12k native / fixed / selected | 60k native / fixed / selected |
|---|---|---:|---:|---:|
| ALK, Davis Kd | 7BTT | 19 | 16 / 17 / 17 | 15 / 15 / 15 |
| MAPK14, Davis Kd | 1KV2 | 29 | 19 / 20 / 20 | 6 / 6 / 6 |
| MAPK14, Davis Kd | 3GCP | 29 | 12 / 13 / 13 | 13 / 14 / 14 |
| RET, Davis Kd | 7DUA | 19 | 8 / 9 / 9 | 14 / 16 / 16 |
| TNK2, Davis Kd | 7KP6 | 19 | 12 / 11 / 10 | 12 / 11 / 11 |
| BACE1, CHEMBL2015090 Ki | 4DJW | 14 | 2 / 5 / 5 | 2 / 3 / 3 |
| BACE1, CHEMBL2015090 Ki | 4GID | 14 | 4 / 4 / 5 | 11 / 12 / 12 |
| BACE1, CHEMBL2155063 Ki | 4DJW | 22 | 1 / 3 / 2 | 2 / 2 / 2 |
| BACE1, CHEMBL2155063 Ki | 4GID | 22 | 13 / 17 / 17 | 9 / 10 / 10 |

For geometry the fixed-first 12k improvement comprises 15 newly correct and three newly wrong directions. At 60k there are six gains and one loss. Adding contact information to the uniform control gives only four and three additional net credits. The geometric classification leaves every eligible rotor constrained in 4,056 of 4,940 poses and 488 of 560 native top poses. Selecting with the correction changes 83 of 560 model indices. At equal heavy-atom count and formal charge, native versus fixed-first contact becomes 8→9/17 at 12k and 11→11/17 at 60k. Dependence on ligand size and composition limits interpretation. The small, fixed, retrospectively selected assay series cannot establish generalization to unseen chemical families. Sensitivity and per-pair results expose dependencies without treating edges as independent draws.

## External human Factor Xa Ki series

After freezing the scorer on the five-target panel, we applied it unchanged to prior coordinates from two additional receptors belonging to the same human Factor Xa target. These 28 inhibitors from one experimental series have 295 qualifying unique pairs. The previous benchmark's native Late runs, with two seeds and a cap of 50,000, supplied 112 pools and 1,000 models. All 112 top scores match the historical results term for term at their printed precision. Following the original FXa analysis, each pair uses the median of its predicted differences across matched receptor–seed blocks. This four-block median of differences differs from the five-target panel's mean-of-receptor-predictions summary. Residue and protonation assignments remain fixed. [Zhao et al.](https://doi.org/10.1016/S0960-894X(00)00139-6) · [Experimental assay](https://www.ebi.ac.uk/chembl/api/data/activity.json?assay_chembl_id=CHEMBL661804&limit=1000) · [FXa protocol, source identities and chemistry exclusions](../../affinity_partition/README.md)

| FXa panel | Pairs | Native | Raw | Uniform | Contact, fixed first | Contact, select from pool |
|---|---:|---:|---:|---:|---:|---:|
| 2P16 | 295 | 187 | 131 | 190 | 191 | 181 |
| 1FJS | 295 | 231.5 | 157 | 236 | 246 | 242 |
| Both receptors, matched blocks | 295 | 223 | 146 | 218 | 226 | 217 |
| Both, equal heavy count and formal charge | 17 | 9 | 8 | 7 | 8 | 7 |

The three-credit fixed-first gain combines 15 gains and 12 losses against native. Deleting one ligand at a time gives differences ranging from −0.37 to +2.53 percentage points. Corrected-score pose selection changes 35 of 112 winning models and gives 14 gains against 20 losses. Its accuracy difference remains negative for every single-ligand deletion (−3.30 to −0.73 percentage points). The apparent transfer for fixed-first is small, changes sign in the size/charge-matched subset, and permits reversals in ligand-deletion sensitivity. Neither protocol establishes an absolute free energy. Only the native score determines the reference pose for the fixed-first protocol.

## Replay and independent audit

Run from the repository root with Python, NumPy, SciPy, pandas, RDKit and xlrd. The scripts read the archived ligand, receptor and pool coordinates directly. The FXa frozen pool archive has SHA-256 `c7786f8f37562a339298f390574ff8aed15131f5aaabc17c11562e80989a5475`. `fxa_activity.json.xz` contains the full 34-entry original ChEMBL assay export. The existing eligibility rule retained 28 stereochemically assigned compounds.

```sh
A=benchmarks/frozen_core/affinity
python "$A/test_entropy.py"
python "$A/external_fxa.py"
python "$A/sensitivity.py"
python "$A/audit.py"
```

The separate audit reconstructs each split from the RDKit SMILES graph and the Meeko serial correspondence, calculates receptor distances with a dense matrix and calculates off-axis displacement by scalar projection. It checks all 5,940 per-pose rotor counts, all five-target native and corrected pair credits, the 47 source measurements against the frozen Davis spreadsheet and ChEMBL exports, and external FXa pair credits against the separate FXa ChEMBL assay export. Per-pose, selected-pose, paired-gap, sensitivity and per-stratum Spearman results are provided. Every geometry check uses printed PDBQT coordinates at 0.001 Å precision. An implementation evaluated before coordinate rounding can differ at the boundary of the two hard tests. Primary geometry and external FXa fixed-first pair gaps all exceed twice the printed raw-energy unit in magnitude on their respective receptor-specific (geometry) and combined (FXa) comparisons.

The archived `native_rotors.csv` and `fxa_native_rotors.csv` come from a separate C++ extraction using the repository's actual model reader. Regenerate them by linking `native_rotors.cpp` to the native build's non-main objects and its Boost libraries, then passing the unique prepared ligand PDBQT files as arguments. This probe reads ligand files only. It never estimates N by inverting observed score ratios.
