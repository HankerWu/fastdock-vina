# Soft amide geometry in search and pose ranking

An ordinary neutral carboxamide usually favors a planar O–C–N–X dihedral. A geometry penalty can down-rank a strained decoy. Applying it to every search step also changes the decoys that reach the endpoint. We varied these two choices independently. Torsional strain has affected the rank of experimentally tested docking decoys in other systems. The test amplitude here is an empirical parameter and is not calibrated to a solution or binding free energy. [Gu et al., 2021](https://doi.org/10.1021/acs.jcim.1c00368)

## Using the geometry term

Generate sites from the prepared Meeko graph and its actual PDBQT torsion tree. Each row contains the PDBQT serials for O, C, N and all one or two N-carbon substituents, followed by a nonnegative barrier. For each active C–N bond, the energy contribution is the barrier times the mean of the arm-specific squared sines. Two equivalent N substituents give one contribution per bond. The two planar wells have equal energy. Rotating between them is allowed. Primary amides, lactams, charged nitrogens, imides, carbamates, ureas, and bonds fixed in the prepared tree are excluded by the site generator. The PDBQT graph alone lacks bond orders and formal charges; generate the site declaration from the same prepared input used for docking.

```sh
python tools/amide_rotamers.py ligand.pdbqt --search-amide ligand.sites --barrier 5
build/vina --receptor receptor.pdbqt --ligand ligand.pdbqt --config box.txt \
    --rank_amide ligand.sites --out ranked.pdbqt
python tools/rank_amide_pool.py ligand.pdbqt native_pool.pdbqt --barrier 5 --geometry all-arms
```

If the prepared nitrogen is pyramidal, its rigid internal geometry can prevent both substituents from being simultaneously coplanar. The mean then has a preparation-dependent floor. For example, the prepared 6WQ9 arm angles are 7.88° and −151.41° while the crystal gives −179.52° and −1.81°. This limits comparisons across prepared states.

The native executable can rank its own refined poses through `--rank_amide`. The PDBQT opening VINA RESULT retains the native affinity; AMIDE RESULT and AMIDE PENALTY record the ranking value and its additive term. The native affinity window set by `--energy_range` determines which poses can be emitted. If it excludes the highest-ranked pose, the output RMSDs are referenced to the first pose actually written. The offline Python command ranks Cartesian coordinates and native scores at PDBQT output precision and writes a CSV. With multiple independent searches or starting states, use the common-reference protocol below before comparing scores.

The separate `--search_amide ligand.sites` option adds analytic Cartesian derivatives during Monte Carlo local optimization. The term is released for explicit-receptor refinement and native affinity calculation. The option can be used with existing receptor or scaffold anchors. A search barrier of five is applied before Vina's conformation-independent normalization; a rank barrier of five is applied to the final normalized docking score. They are separate experimental choices.

## Paired experiments on new seeds

Five human CA II ligands, two receivers and two seeds form 20 cells. Sixteen cells transfer the ligand to a different receptor structure. Each cell has four trajectories covering two previously prepared ammonium geometries with two starts each. All comparisons use the same seeds, box, graph, receiver, independent scaffold anchors, one CPU per trajectory, a 50,000-evaluation search limit and at most nine modes per trajectory. Each pose must have registered whole-ligand heavy-atom RMSD ≤2 Å and the deposited amine handedness to count as a hit. Each query’s own deposited pose is withheld from the scaffold restraint and torsion calculation. The protocol assumes a known pocket and uses aligned crystal coordinates from other ligands as search priors. The seeds, 1618033 and 1414213, were outside the two previous panels; the five compounds and two receivers were already part of development. [Deposited hCA II series](https://doi.org/10.1021/acs.jmedchem.0c00733)

| Search with scaffold anchors | Final all-arms barrier | Top poses / 20 | Cross-docking top / 16 | Pools with a hit / 20 | Same-assay Ki pairs / 7 |
|---|---:|---:|---:|---:|---:|
| Baseline | 0 | 1 | 1 | 12 | 3 |
| Baseline | 5 | 9 | 7 | 12 | 1 |
| Baseline | 15 | 7 | 6 | 12 | 1 |
| Amide-guided, B=5 | 0 | 1 | 0 | 12 | 2 |
| Amide-guided, B=5 | 5 | 1 | 0 | 12 | 1 |
| Amide-guided, B=5 | 15 | 2 | 1 | 12 | 1 |

With rank barrier five and baseline search, winners for 6WQ4/5/7/8/9 counted 1/4, 4/4, 0/4, 4/4 and 0/4. The two receivers contributed 4/10 and 5/10. Guiding only slot zero or only slot one gave 2/20 and 4/20, respectively, for the same total budget. Full guidance raised the number of acceptable poses from 93/691 to 102/687. It also increased the number of incorrect poses with mean squared sine below 0.10 from 84/598 to 510/585. This explains how both candidate classes became more planar while a planarity penalty lost discrimination. For 6WQ7 and 6WQ9, the best RMSDs improved from 4.40 and 4.28 Å to 4.15 and 4.06 Å. The deposited terminal tails of these two ligands have relatively weak omit density, an uncertainty in the distal reference coordinates. [Original structural interpretation](https://doi.org/10.1021/acs.jmedchem.0c00733)

For factor Xa we used 28 compounds from the same Ki assay, two receptor structures and the same two seeds. The threefold cutoff gives 295 unique compound pairs. For each pair the sign is evaluated on its median predicted gap over matched receptor/seed blocks. Ties receive half credit. The baseline, baseline plus rank B=5, and B=5 search plus B=5 rank scored 218/295, 221/295 and 218/295. Three compounds have eligible amides. Of the 66 eligible-to-other edges, the corresponding scores were 63, 66 and 63. The other 229 edges contributed 155 in each arm. No pair among the three eligible compounds reaches the threefold cutoff. The hCA II and FXa results give different directions for the same change in affinity ordering. [FXa assay](https://doi.org/10.1016/S0960-894X(00)00139-6)

## Individual CLI jobs and pooled results

Running the production CLI with rank barrier five on each of the 80 new-seed, baseline-search hCA II jobs gave 15 correct first poses against 13 native, with a candidate ceiling of 21/80. For cross-docking the counts were 11 against 10, with a ceiling of 14/64. The 9/20 result above combines four outputs per cell under a common energy reference. A single CLI search exposes only its own candidates. All 691 poses matched the prior native output by atom records, native affinity and energy components. The independent Python ranker returned the same complete order for 79/80 pools. In the remaining pool, the first two candidates differed by 0.0002634 kcal/mol inside Vina. Rounding coordinates and affinities to PDBQT precision reversed them; both met the structural endpoint. The largest per-model discrepancy in the amide penalty from rounding was 0.00716 kcal/mol. [Per-job comparison](results/guidance/production.csv)

## Generating missing tail poses with a second donor

The first donor for each query is the other ligand with the largest common substructure under the existing mapping and tie-breaking rules. We removed that donor and reapplied the same rules to choose the second, without access to the query's bound pose or affinity. Each donor supplies its own registered, symmetry-aware set of anchor assignments. Retaining both starting amine states, we gave one of the two search slots to each donor, then exchanged their slot assignments. Both mixed arms have the same four trajectories per cell as the first-donor baseline. There is no amide bias during these searches. Using only the second donor in both slots is a separate equal-cost control.

| Seeds | Donors assigned to slots 0 / 1 | Cells containing a joint hit / 20 | First poses with B=5 / 20 | 6WQ9 cells containing a joint hit / 4 |
|---|---|---:|---:|---:|
| 1618033, 1414213 | first / first | 12 | 9 | 0 |
| 1618033, 1414213 | second / first | 13 | 6 | 1 |
| 1618033, 1414213 | first / second | 15 | 1 | 3 |
| 1618033, 1414213 | second / second | 9 | 2 | 3 |
| 314159, 271828 | first / first | 12 | 7 | 0 |
| 314159, 271828 | second / first | 14 | 7 | 2 |
| 314159, 271828 | first / second | 13 | 8 | 2 |
| 314159, 271828 | second / second | 11 | 3 | 3 |

For 6WQ9, the first donor (6WQ7) matches 35 atoms and the second (6WQ8) matches 32. With only the first donor, the best correct-handed RMSDs were 4.277 and 4.515 Å in the two seed panels. With only the second donor they were 1.386 and 1.485 Å, with a joint hit in three of four cross-docking cells in each panel. For those six cells the first correct candidate ranked 4th, 2nd, 10th, 2nd, 3rd and 10th at B=5. The mixed arm recovered a first-ranked 6WQ9 pose in one cell in the replication panel. Neither donor yielded a joint hit for 6WQ7. These results support spending some search budget on a different reference when the goal is to recover alternative candidate geometries. Selection of the first pose and prediction of Ki still require separate validation. [New seed results](results/donor_diversity/new_seeds) · [Replication](results/donor_diversity/replication)

## Reproduction

Extract and verify the [prepared fixture](fixtures.md), build the binary, and run the following commands from the repository root. `BASE` contains the extracted fixture, and `WORK` is a fresh output directory.

```sh
python benchmarks/amide_search/guidance.py dock build/vina "$BASE" "$WORK" --jobs 4
python benchmarks/amide_search/guidance.py score build/vina "$BASE" "$WORK" --jobs 4
python benchmarks/amide_search/verify_native_rank.py build/vina "$BASE" "$WORK" "$WORK/production" --jobs 4
python benchmarks/amide_search/donor_diversity.py build/vina "$BASE" "$WORK/donor-new" --jobs 4
python benchmarks/amide_search/donor_diversity.py build/vina "$BASE" "$WORK/donor-repeat" --seeds 314159 271828 --jobs 4
repo_dir="$PWD"
(cd "$WORK" && sha256sum -c "$repo_dir/benchmarks/amide_search/results/guidance/pools.sha256")
```

The 304 baseline and 104 additional guided searches yielded 3,484 poses. A separate audit recomputed all 265 distinct first-choice poses with the previous binary’s ordinary score-only entry point and an explicit zero unbound energy. The largest discrepancy against the recorded native score at that same reference was 0.000500 kcal/mol. The 25 FXa ligands without a movable eligible amide share their matching baseline outputs. For each receptor, ligand and seed, the baseline pose with the lowest explicit INTER + INTRA supplies the native unbound INTRA reference for every treatment. All candidate coordinates are scored using one reference ligand topology and explicit receptor. The new [mode records](results/guidance/modes.csv), [cell winners](results/guidance/cells.csv) and [summary](results/guidance/summary.csv) include both one-slot crossover arms. The mode file contains one row per native/guided pose before reranking. No trajectory count or seed replication creates additional independent ligands.
