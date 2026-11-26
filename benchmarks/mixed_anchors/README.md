# Native and anchor-guided trajectories at equal budget

The optional `--search_anchor_trajectories all|even|odd` applies an existing point, shell or multi-assignment search prior to all trajectories (the default) or to the selected 1-based trajectory indices. For example,

```sh
vina --receptor rec.pdbqt --ligand lig.pdbqt --config box.txt --search_anchor sites.txt \
     --search_anchor_trajectories even --exhaustiveness 4 --max_evals 50000 \
     --search_counters --seed 314159 --out mixed.pdbqt
```

The two split settings require at least two trajectories. With four trajectories, even uses the anchor on trajectories 2 and 4 and odd uses it on 1 and 3. Odd exhaustiveness gives the odd arm the extra guided trajectory. The same trajectory index receives the same random stream independent of worker assignment. The native and guided trajectories differ only in their receptor search grid. The engine pools their saved minima, scores them against the same explicit native receptor, performs the existing Late selection and refinement, then ranks the surviving poses without the anchor. The option reuses the existing alternate-search-landscape dispatch. Split and all accept the identical anchor file.

## Hypothesis and protocol

If independent native and guided trajectories generate a desired pose with respective probabilities p and q, then n native and g guided trials miss it with probability (1−p)^n (1−q)^g. This explains how reserving native trajectories can limit the damage from an imperfect prior, at a cost when the prior is reliable. The actual engine prunes a pooled set before refinement and removes duplicates. Its returned-pool coverage can therefore depart from this independent-trial formula. We compare a fixed four trajectories with the same per-trajectory soft evaluation limit. Both complementary assignments are scored separately on the same seeds, ligands and receptors. Pooling the two output files would double the search budget. [Arcon et al., experimental use of interaction priors](https://doi.org/10.1093/bioinformatics/btz152)

Each assignment contains 224 human FXa Ki dockings, 40 bovine trypsin Kd dockings and eight cross-dockings of experimental FXa ligand DX-9065a into two different receptor structures. The wrong-end stress test reuses the 28 FXa ligands at two seeds (112 dockings) and deliberately anchors a second basic group in DX-9065a (eight more). Three external scaffold-anchor matrices add 48 BRD4, 72 T4L and 40 human CA II dockings. In total, 552 half-guided calculations per assignment were compared with pure-arm references. For each mixed run all four search-counter lines must match those from the corresponding full-native/full-guided index. Pure-arm controls and prepared inputs come from the prior branch benchmark. The initial native FXa files lacked trajectory logs; 112 full-native searches were replayed for that audit. All cells use CPU 1, four trajectories, 50,000 soft evaluations per trajectory, up to nine models, a 100 kcal/mol energy range and native Late ranking. Equal limits allow small deterministic terminal-step overshoots.

The DX9 and trypsin crystal labels use symmetry-corrected heavy-atom RMSD ≤2 Å. A known Asp189 interaction supplies the ion-pair prior. Of the 28 Zhao FXa affinity compounds, none supplies a whole-ligand experimental pose here. DX9 is a separate ligand with an independent full-pose label. Trypsin has a full-pose label for benzamidine; its other four ligands contribute to the affinity test. On CA II the fixed preparation in this scaffold matrix contains only the original local nitrogen state. Use the accompanying plain-RMSD and centre-geometry results separately. [FXa crystal](https://www.rcsb.org/structure/1FAX) · [FXa experimental series](https://doi.org/10.1016/S0960-894X(00)00139-6) · [Trypsin ITC](https://doi.org/10.1046/j.1432-1327.2001.01991.x)

## Results

**Assay and coordinates.** A unique chemical edge has at least a threefold observed within-assay difference. For each arm, calculate the weak-minus-strong native score difference within the same receptor and seed, then take the median across matched blocks and give exact ties half credit. The constitutive ligands are shared among edges. The ligand bootstrap resamples whole ligand identities and carries the weights to their edges. Four-seed FXa and trypsin data have eight blocks per chemical edge. The wrong-FXa control has four blocks per edge and is compared to the corresponding two-seed native baseline. A second calculation fixes the unbound reference energy for each ligand/seed/receptor at its native top-pose value across all search arms. Sizes and formal charges are the measured-panel ligand preparations from the prior experiment. Neither seed allocation uses the experimental pose or Ki/Kd to select a score. Structural counts below refer to retained output poses after common refinement and pruning.

| Affinity, one vote per chemical pair | Native | All guided | Even guided | Odd guided |
|---|---:|---:|---:|---:|
| FXa Ki, /295 | 221 | 234 | 238 | 229 |
| FXa Ki with common native unbound reference, /295 | 221 | 225 | 232 | 223 |
| FXa Ki, equal heavy atom count and charge, /17 | 8 | 13.5 | 14 | 12 |
| FXa Ki, every receptor/seed block correct, /295 | 35 | 170 | 93 | 119 |
| FXa Ki, measured ratio ≥10, /231 | 178 | 192 | 195 | 187 |
| Wrong FXa head, original two seeds, /295 | 223 | 199 | 207 | 190 |
| Trypsin Kd, /3 | 0 | 0 | 0 | 0 |
| BRD4 Kd, /5 | 1 | 1 | 1 | 0 |
| T4L Kd, /19 | 12 | 12 | 12 | 12 |
| CA II Ki, original local state, /7 | 2 | 3 | 1 | 3 |

The per-receptor FXa Ki credits at 1FJS were 244.5/234/240/237, in the column order above; at 2P16 they were 187/228/231/231. On the reported median sign the even arm beat all-guided by four pair credits and the odd arm lost five. Their whole-ligand bootstrap intervals relative to native span −5.3 to +17.4 and −9.0 to +13.6 percentage points. Every-block wins favor full guidance. Reversing the wrong head recovered correctly bound aryl contacts in some mixed pools but did not reliably improve Ki rank. An independently written Decimal parser reproduced all reported chemical pair counts directly from the first MODEL of each native, guided and mixed file (1,128 file reads per allocation).

| Crystal pose or contact | Native | All guided | Even guided | Odd guided |
|---|---:|---:|---:|---:|
| DX9 in two FXa receptors, top / pool coverage, /8 | 0 / 0 | 5 / 8 | 4 / 8 | 6 / 7 |
| DX9 deliberately reversed head, top / pool, /8 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| Trypsin benzamidine, top / pool, /8 | 4 / 4 | 8 / 8 | 8 / 8 | 8 / 8 |
| FXa aryl amidinium two contacts, top / pool, /224 | 33 / 54 | 223 / 224 | 208 / 224 | 207 / 224 |
| Wrong-head FXa experiment, correct head two contacts, top / pool, /112 | 14 / 26 | 0 / 0 | 7 / 11 | 6 / 17 |
| Trypsin two contacts, top / pool, /40 | 24 / 24 | 40 / 40 | 39 / 39 | 40 / 40 |
| BRD4 scaffold, top / pool, /48 | 1 / 2 | 19 / 24 | 14 / 23 | 12 / 23 |
| T4L scaffold, top / pool, /72 | 55 / 72 | 58 / 71 | 54 / 72 | 56 / 71 |
| CA II scaffold, RMSD alone, top / pool, /40 | 0 / 0 | 11 / 16 | 3 / 15 | 4 / 13 |

Both the top pose and best retained pose passed the stricter CA II centre-geometry plus RMSD criterion in zero cells in all four arms. These 40 cells use the single starting state and are distinct from the larger pool combining two local states described in the previous research. Choosing the parity after seeing these outcomes would overstate generalization. The optional split can help hedge a questionable prior and retain chemistry missed by fully guided trajectories. The experiments give no consistent pose or affinity advantage over full guidance. Search cost across 552 elementary comparisons was 110,464,071 evaluations for even and 110,463,275 for odd, a difference of 796 evaluations (about 0.00072%). Files under results/even and results/odd use the arm label mixed inside their corresponding allocation directory.

## Reproduction

The prior physics evidence archive supplies the input fdv-* directories, including the prepared FXa, trypsin and cross-docking ligands, their template-free native controls, original full-anchor controls and crystal coordinates. Set BASE to its unpacked directory, BIN to the built fastdock-vina-pro CLI, and R to a fresh output directory. Scripts require Python, NumPy, SciPy, RDKit, Meeko and Gemmi as used by the directional-ion-pair benchmark. Each parity has its own output directory.

```sh
for parity in even odd; do
  dest="$R/$parity"
  python benchmarks/mixed_anchors/dock.py "$BASE" "$BIN" "$dest" --parity "$parity" --jobs 4
  python benchmarks/mixed_anchors/scaffolds.py "$BASE" "$BIN" "$dest" --parity "$parity" --jobs 4
  python benchmarks/mixed_anchors/cross.py "$BASE" "$BIN" "$dest" --parity "$parity" --jobs 4
  python benchmarks/mixed_anchors/summarize.py "$BASE" "$dest" "$dest/results"
  python benchmarks/mixed_anchors/verify.py "$BASE" "$dest" "$dest/results"
done
```

The two result directories retain per-cell scores, per-edge differences, per-pose-set RMSDs, numerical summaries and trajectory counters. In the wrong-head stress test both electron donors are selected from the known terminal cyclic amidine rather than the aryl amidine. The same intentionally reversed ligand head is used at both parities. The core regression `benchmarks/check_search_anchor.py BIN BASELINE_BIN` verifies a three-trajectory case with odd exhaustiveness, CPU-count invariance, the legacy all-guided setting, native default equivalence and selected/released pose energies.
