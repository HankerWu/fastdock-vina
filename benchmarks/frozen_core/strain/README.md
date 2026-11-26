# Full intraligand MMFF94s replacement on frozen Vina pools

The tested hypothesis is that replacing the empirical internal energy difference in a saved Vina pose by a unit-weight conformational strain estimate can improve relative ranking. For reported score S, Vina flexibility multiplier f, saved intraligand energy I and common unbound reference U, use

    S_MMFF = S + f [(E_bound - E_free) - (I - U)].

E_bound uses MMFF94s with the original prepared microstate and every heavy atom held fixed; only H coordinates relax. E_free is the best converged vacuum minimum from eight independently embedded ETKDGv3 conformers under the repository's existing physical-pool protocol. E_free cancels when ranking conformations of the same ligand and matters when ranking compounds. The assay analysis also repeats with 32 independently embedded conformers, separate random seed 49043. These approximate energy differences omit solvent, binding entropy and changes in protonation. The compound outputs are empirical scores in model units. Ligand chemistry records and graphs were checked against the original Meeko preparation for every model.

The earlier strain study tested 0.1-weight replacement, screened electrostatics and local references in T4L, hCAII and BRD4. This study retains its 0.1 reference and a control that just subtracts the native intraligand difference. The unit-weight hypothesis and comparison arms were set before these outcomes. Relevant prior investigations include a 150-complex strain analysis and torsion-profile filtering in a large docking library. [Previous in-repository study](../../strain_transfer/README.md) · [Perola and Charifson, 2004](https://doi.org/10.1021/jm030563w) · [Gu et al., 2021](https://doi.org/10.1021/acs.jcim.1c00368)

## Pose recovery

1,214 already-emitted models from 136 pools and 13 human targets were rescored without any search or coordinate moves. Success means the selected pose is ≤2 Å from the mapped experimental heavy-atom coordinates in the receptor frame. Repeated seeds, budgets, ligands and receptor structures make the rows statistically dependent.

| Frozen source | Pools | Native top1 | Interaction-only | 0.1 replacement | Unit replacement | Any emitted |
|---|---:|---:|---:|---:|---:|---:|
| Armijo 12k, two original seeds | 28 | 0 | 0 | 0 | 1 | 1 |
| Armijo 120k, same seeds | 28 | 1 | 1 | 2 | 2 | 4 |
| Armijo 12k, other seeds, cognate | 20 | 1 | 1 | 1 | 0 | 3 |
| Armijo 12k, other seeds, reciprocal | 20 | 7 | 7 | 7 | 3 | 9 |
| Higher matched-budget Native, same other seeds, cognate | 20 | 3 | 3 | 3 | 2 | 10 |
| Higher matched-budget Native, same other seeds, reciprocal | 20 | 10 | 10 | 9 | 6 | 14 |

The original-seed gain in both rows comes from the same BPTF 7RWO ligand/seed. At the high budget its native first pose (7.34 Å) has 10.57 kcal/mol estimated strain; a rank-six pose (1.55 Å) has 2.88 kcal/mol. In DYRK1A 7AKL→3ANR on seed 85231, a 0.78 Å native winner has 6.23 kcal/mol of strain and was displaced by a 7.08 Å pose with 1.94 kcal/mol. BPTF 7F5D on seed 85231 at matched budget likewise went from 0.47 to 6.23 Å. These observations reject default structural reranking by the replacement score.

## Measured same-assay affinity ordering

The experimental panel has four Davis kinase Kd assays and two BACE1 Ki assay series, 32 distinct ligands, five human targets and seven receptors. Each row conditions on one target, assay, endpoint and receptor. For each independent search seed and compound, score either the saved native rank-one pose or the minimum across its emitted models. Take the median across the two seeds, then compare experimentally ≥3-fold-separated pairs with one-half credit for ties. Source measurements, conditions and input structures are documented with the original panel. [Panel provenance](../../geometry_comparison/affinity/README.md)

| Pool and policy | Native | Eight free conformers | Independent 32 | Ligand-only prior, 8 / 32 | Pairs |
|---|---:|---:|---:|---:|---:|
| Earlier 12k Armijo, pool minimum | 93 | 124 | 123 | 76 / 91 | 187 |
| Earlier 60k Armijo, pool minimum | 86 | 112 | 115 | 78 / 91 | 187 |
| Fresh native search, keep native pose | 89 | 132 | 132 | 91 / 98 | 187 |
| Fresh native search, pool minimum | 89 | 133 | 131 | 91 / 98 | 187 |

The ligand-only control adds the normalized strain of the independently generated input conformer to native Vina; it uses no bound conformation. The heavy-atom count baseline scores 68.5/187. The 187 comparisons reuse ligands and, in some strata, the same assay across two receptors. On strictly equal-heavy-atom and equal-formal-charge pairs only 17 comparisons per budget survive. Aggregate native vs unit credit is 8→9 at 12k, 11→8 at 60k, and 8→9 on fresh Native for either pose policy. Chemical series and composition remain confounders.

For BACE1 CHEMBL2155063 on 4DJW, all eight assayed inhibitors use a +1 prepared microstate. Full replacement changes native 2→13/22 (12k), 2→9/22 (60k) and 4→14/22 (fresh Native, pool minimum). The second receptor, 4GID, yields 14→19, 9→14 and 10→19. A post-hoc sensitivity on both +1 BACE1 series reran MMFF with dielectric 4 in the bound and reference states. For CHEMBL2155063 this gives 14/22 and 13/22 on 4DJW and 19/22 and 15/22 on 4GID at the two earlier budgets. For the separate BACE1 assay CHEMBL2015090 on 4GID, 60k worsens from 12/14 to 9/14 using either dielectric. The RET Kd improvement can also be obtained with the 32-conformer ligand-only prior; that target alone does not identify a bound-pose mechanism.

Concrete fresh-seed 4DJW example, keeping native rank-one poses: CHEMBL2151138 (57 nM) scores −8.0210 natively and −5.9213 after replacement. CHEMBL2151136 (9,500 nM) scores −8.9710 and −4.0084. The predicted weak-minus-strong gap moves from −0.9500 to +1.9129. Free references based on 32 rather than eight conformers shifted individual free energies by −11.935 to +0.191 kcal/mol across the panel. Absolute score magnitudes should not be interpreted as measured binding free energies.

A retrospective external transfer uses 28 human Factor Xa inhibitors in a common Ki assay, two structures and two seeds. The pre-existing endpoint is the median of within-receptor/seed score differences. The native baseline was reproduced at 223/295. The unit correction drops to 206/295 keeping native pose and 202/295 selecting by corrected score; for the 17 equal-size/charge pairs it drops from 9/17 to 7/17 under both policies. With pool minimum, 2P16 changes from 187 to 197/295 but 1FJS drops from 231.5 to 205/295. This is a negative transfer result for a universal correction. [FXa assay and estimators](../../affinity_partition/README.md)

## Running and auditing

`tools/score_physical_pool.py --mmff-affinity` emits a separate standalone affinity annotation, with per-model scores, the saved native rank-one score and the minimum over the emitted pool. Existing structural orders use the existing rules. The regression exercises the normal CLI on a real BACE1 pool, compares every model against the frozen archive and verifies both selection conventions and the legacy orders. [Reproducible command and example](../../../docs/FIXED_PHYSICAL_POOL.md)

Use the same Python/RDKit/Meeko environment as the archived studies and an executable supporting fixed-coordinate net-polar scoring. The experiments extract frozen inputs from repository archives. `pose_summary.csv` and related plain CSVs show per-stratum results. Compressed CSVs retain every model and pair. The direct coordinate-independent audit reconstructs the replacement using the raw saved INTER line; all 1,214 pose and 2,488 initial affinity model scores agree within 0.00150 units (rounding of PDBQT terms) and all 416 selected model IDs match. A fresh full pose rerun reproduced all per-mode numerical values.

```sh
D=benchmarks/frozen_core/strain
python "$D/pose.py" "$VINA" work/strain-poses --jobs 4
python "$D/affinity.py" "$VINA" work/strain-affinity --jobs 4
python "$D/references.py" work/strain-affinity --jobs 4
python "$D/charged.py" work/strain-affinity --jobs 4
python "$D/native.py" "$VINA" work/strain-fresh --jobs 4
python "$D/matched_pairs.py" work/strain-affinity work/strain-fresh work/strain-fresh
python "$D/factor_xa.py" "$VINA" work/strain-fxa --jobs 4
python "$D/audit.py" work/strain-poses work/strain-affinity work/strain-audit.json
python "$D/check_cli.py" "$VINA"
```
