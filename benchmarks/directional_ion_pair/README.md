# Two-contact searches for aryl amidinium binding to a protein carboxylate

An optional search prior gives two chemically equivalent nitrogens on a protonated aryl amidinium access to both oxygens on a named Asp or Glu. The linked distance-shell version restores native-ranked benzamidine poses within 2 Å of their crystal structures in all eight repeats across two trypsin receptors, up from four with ordinary search. Human FXa affinity ranks improve on one receptor and decline on the other. The small trypsin affinity panel fails across all four seeds. These findings support a search option when independent structural information identifies a primary aryl amidinium–carboxylate site. They provide no basis for a general affinity correction.

## Normal use

The ligand must have a declared protonated primary aryl amidinium with a Meeko SMILES-to-PDBQT map. Select the residue and charge states using information independent of the ligand's affinity. Here Asp189 was identified by earlier structures. A prepared single receptor conformation and the prepared ligand suffice to generate a search file. An experimental ligand pose is unnecessary. The ordinary Vina late-selection path releases the search guide before explicit refinement and ranking. See the [engine equations](../../docs/ANCHOR_AND_MICROSTATE_RESEARCH.md).

```sh
python tools/make_ion_pair_anchors.py receptor.pdbqt ligand.pdbqt pair.txt \
  --chain A --residue 189 --geometry shell
build/vina --receptor receptor.pdbqt --ligand ligand.pdbqt --config box.txt \
  --pose_selection late --search_anchor pair.txt --out poses.pdbqt
```

For a carboxylate center C and oxygens O1, O2, define the unit vector u from C to their midpoint. The point control steers either amidinium nitrogen toward either target Oi + 2.90 u Å. The two-point prior steers both nitrogens as one assignment with both permutations allowed. The double shell instead requests two N–O distances of 2.90 Å and preserves angular freedom. Each complete assignment is evaluated before combining the alternatives. Multiple matching primary amidine groups produce additional whole-group alternatives. A single point can localize the correct end without achieving the two distances. Scoring uses the unchanged Vina release path.

## Experiments and structural predictions

The main comparison has four native/point/pair/shell arms, four seeds (271828, 314159, 161803, 141421), two receptor conformations per target, 28 same-assay FXa compounds and five matched trypsin compounds. The initial native FXa outputs came from the branch's earlier holdout. The new seeds were run after the original two-seed analysis. The two-point arm was frozen first. The shell was a follow-up inspired by the crystal audit and was then held fixed on the new seeds. Every arm used four trajectories, one CPU per docking, a per-trajectory soft budget of 50,000 evaluations, late selection and up to nine returned modes. The 1,056 main pose files comprise 9,263 modes. We also ran 112 reversed-head controls and 32 input-repair dockings. Each experimental pair contributes one vote using the median of within-receptor, within-seed gaps (half credit for a tie). Ki and Kd edges are always evaluated within their own target and assay.

The external test of the bisector locates two independent nitrogens in 3PTB, 1BTY and neutron structure 5MO0 (benzamidine), 1FJS (Z34), and 1FAX (DX9). The best assignment to the two ideal points has 2.05–3.24 Å RMS displacement across these structures. The two actual N–O distances range from 2.71 to 3.01 Å in four complexes. In 1FAX the distances are 2.80 and 4.17 Å. This directly contradicts a universal two-short-contact requirement. For docking, 1BTY and 3PTB have complete trypsin binding sites for this preparation. 5MO0 has unmodeled protein side chains and supplies an independent geometry measurement only. [Crystal measurements](results/crystals.csv) · [Neutron experiments](https://doi.org/10.1038/s41467-018-05769-2) · [FXa–DX-9065a](https://www.rcsb.org/structure/1FAX) · [FXa–ZK-807834](https://www.rcsb.org/structure/1FJS)

| Four-seed result | Original search | One point | Two points | Two shells |
|---|---:|---:|---:|---:|
| FXa top with both N–O distances ≤3.5 Å, /224 | 33 | 218 | 224 | 223 |
| FXa pool containing such a pose, /224 | 54 | 222 | 224 | 224 |
| Trypsin top with both N–O distances ≤3.5 Å, /40 | 24 | 40 | 40 | 40 |
| Benzamidine crystal RMSD ≤2 Å in top, /8 | 4 | 8 | 8 | 8 |
| FXa ≥3-fold Ki pairs, /295, initial seeds | 223 | 220 | 233 | 235 |
| FXa ≥3-fold Ki pairs, /295, new seeds | 210 | 249 | 226 | 233.5 |
| FXa ≥3-fold Ki pairs, /295, four seeds | 221 | 230.5 | 229 | 234 |
| FXa, equal heavy atoms and formal charge, /17 | 8 | 14 | 11 | 13.5 |
| Trypsin ≥3-fold Kd pairs, /3 | 0 | 0 | 0 | 0 |

The cognate redocking endpoint in the main panel is benzamidine in the two trypsin receptors. The eight repeats share one ligand. FXa contacts report local geometry around the specified head. The 28 query ligands have no cognate crystal references in this panel. FXa shell ranks are 234/295 against 244.5 native in the 1FJS receptor and 228/295 against 187 native in the apixaban-bound 2P16 receptor. Across both receptors 39 chemical pairs improve, 26 deteriorate and the other 230 keep their credit. A paired ligand-vertex bootstrap (3,000 replicates) gives a 95% sensitivity interval of −7.9 to +16.8 percentage points around the +4.4-point difference. Deleting each one of the 28 ligands gives +1.5 to +6.5 points. Both calculations preserve the sharing of ligands between pair comparisons. The two related proteases share an S1-site motif. Transfer to a distant target family or a charge-neutral ligand is unresolved. [Pairs and structures](results/contrast.csv) · [Mode-level outputs](results/modes.csv)

## Independent FXa cross-docking on the one-sided crystal contact

We tested a structure-defined exception using the chiral ligand DX-9065a in 1FAX, whose measured Asp189 distances are 2.80 and 4.17 Å. Its aryl amidinium and acetimidoylpyrrolidine are protonated and its carboxyl group is deprotonated in the assumed +1 state. We checked the two stereocenters against the wwPDB CCD's canonical CACTVS SMILES, discarded all experimental ligand coordinates, then made a fresh conformer from the graph. Search used the prepared 1FJS and 2P16 receptors and their existing boxes. Aggregation compares the 33-heavy-atom DX9 poses with the 1FAX ligand after superimposing 74 shared receptor Cα positions within 15 Å of the 1FAX Asp189-Cγ. The receptor fits have RMSD 0.561 and 0.583 Å. The eight allowed graph permutations include the resonance-equivalent amidinium and carboxylate atoms. Every arm used the same four seeds and search budget.

| ≤2 Å full-ligand cross-pose recovery | Ordinary | One point | Two points | Two shells |
|---|---:|---:|---:|---:|
| Top on 1FJS, /4 | 0 | 0 | 0 | 1 |
| Top on 2P16, /4 | 0 | 1 | 4 | 4 |
| Any of the returned candidates, /8 combined | 0 | 8 | 8 | 8 |
| Median top RMSD, Å, /8 combined | 8.75 | 2.05 | 2.01 | 1.95 |

Thus the shell gives 5/8 first-pose successes on this distinct scaffold even though it shortens the long experimental contact (median top 3.07 Å in the two cross receptors versus 4.17 Å in 1FAX). On 1FJS three top poses still miss 2 Å despite a ≤2 Å candidate in every run. This separates approximate localization of the molecule from faithful reproduction of the salt bridge and from recognition of the best saved candidate. Repeats still represent one chemical identity. An independent positional graph search and SciPy rotation fit reproduced every top and candidate-minimum RMSD in 32 additional files (285 modes) to 1.1 × 10⁻¹⁴ Å. [Actual DX9 cross-docking results](results/dx9.csv) · [Independent audit](results/verification_cross.txt) · [1FAX coordinates](https://www.rcsb.org/structure/1FAX) · [DX9 chemistry](https://www.rcsb.org/ligand/DX9)

## Separating search quality from a moving energy reference

The native Vina report uses the intramolecular energy of the first native pose retained by each search as its common unbound reference. The first retained conformation can change between search arms. For each ligand, receptor and seed, we recomputed all top-pose scores with the original native arm's unbound U using (E_inter + E_intra − U)/(1 + 0.05846 Nrot). We recovered the native torsion factor from all reported component ratios for that ligand, constrained it to its half-integer rotor count, and verified replay of the arm-specific scores to 0.0015 kcal/mol. At the PDBQT component print precision, reconstructing each arm with its *own* U gives 221, 230, 229 and 235 credits, against printed totals of 221, 230.5, 229 and 234. With the *same* U the four totals are 221, 224, 220 and 225. Shell versus native with common U is 227 versus 223 on the initial seeds and 217 versus 210 on the fresh seeds. A common native U controls a computational discrepancy. The reference still lacks a sampled, solvated ligand ensemble. [Reference reconstruction](results/reference_summary.csv)

## The 173-fold FXa contrast contains two chemical changes

Number the N-naphthylmethyl benzimidazole from substituted N1. The two changes are n-propyl to isopropyl at C2 and migration of the piperidin-4-yloxy attachment between C5 and C6. All four molecules have 36 heavy atoms and formal charge +2 in the docked state. The manifest retains the published identities, neutral graph SMILES and prepared protonation assumptions. [Four complete graphs](results/substitutions.csv) · [Original assay and series](../affinity_partition/README.md) · [Zhao et al.](https://doi.org/10.1016/S0960-894X(00)00139-6)

| C5/C6 attachment | C2 side chain | ChEMBL identifier | Ki (nM) |
|---|---|---|---:|
| C6 | isopropyl | CHEMBL177142 | 0.026 |
| C6 | n-propyl | CHEMBL176813 | 0.23 |
| C5 | isopropyl | CHEMBL176744 | 2.0 |
| C5 | n-propyl | CHEMBL369161 | 4.5 |

Switching C5 to C6 improves measured affinity 19.57-fold with n-propyl and 76.92-fold with isopropyl. Branching the side chain improves it 8.85-fold at C6 and 2.25-fold at C5. Thus the difference of the two branching free-energy gaps is 0.811 kcal/mol (298.15 K nominal for energy comparisons because the FXa assay temperature is undocumented). The shell predicts 0.088 kcal/mol for this interaction between substitutions. Its 0.944 kcal/mol diagonal underestimates the experimental 3.054 by 2.11 kcal/mol. The ordinary result is 0.064 kcal/mol. Holding U fixed gives 0.556 versus 0.065 kcal/mol. On the original seeds, the deliberately incorrect distal acetimidoylpiperidine head also gives the correct diagonal sign and 1.853 kcal/mol, while placing the primary aryl head in a double contact in 0/112 first poses. That control earns only 199/295 overall versus 235 for the declared shell and 223 for the ordinary search. The single largest potency difference is insufficient to identify a microscopic mechanism. [Substitution cycle](results/substitution_cycle.csv) · [Reversed-head control](results/inverted.txt)

One of the other FXa ligands, CHEMBL369480, converged under MMFF with seven distorted C–C/C–N bonds (1.875–2.094 Å) in the original seed-2026 input. Sixteen independent ETKDG initializations gave a converged low-energy conformer at seed 2036 (52.17 vs 693.48 kcal/mol, longest heavy-atom bond 1.546 Å). Rerunning all 32 cells for this ligand yields 220/295 native and 234/295 shell over the full panel. Both versions are retained. [Conformers](results/conformer_selection.csv) · [Original/repaired contrasts](results/contrast.csv)

## Solution and ligand-partner competition

As a diagnostic, let B and I be Gaussian overlap factors for the two matched donor nitrogens relative to the receptor carboxylate or the ligand's own carboxylate, respectively. We used exp(−d²/2) for d in Å, where d is the best two-point assignment RMS. With a fixed maximum stabilization of 2RT, q = exp(2) − 1, and mutually exclusive partners, a dimensionless ratio of contact-state weights is (1 + q I + q B)/(1 + q I). Taking minus RT times its logarithm gives a correction that saturates and declines when the intramolecular contact is already formed. The candidate bonus uses each bound conformer, implicitly assuming that its internal geometry represents a competing free-ligand state. This is a diagnostic approximation. Applied with and without I, the native-pool rank drops from 221 to 216/295 in both cases. On the shell pool the two variants give the same 233/295, versus 234 without the bonus. For the five carboxylate-bearing ligands, 32 separate ETKDG starts each, minimized under MMFF94s at dielectric 4 or 78.4, gave 320 converged geometries and zero internal bidentate contacts under the 3.5 Å criterion. The lowest larger N–O distances were 3.69 and 4.60 Å. These discrete minimizations do not measure aqueous populations. [Per-conformer distances](results/free_ligand.csv) · [Ablation](results/summary.csv)

The bovine trypsin panel comes from one matched ITC series at 25.1 °C, 50 mM Tris and 10 mM calcium chloride, pH 8.0. We inverted the measured association constants to Kd. The authors' tables give CONH2 for one para substituent while the materials text says 4-amidinobenzamidine. The calculation uses the table formula and separately excludes that row. This leaves one ≥3-fold comparison, amino versus methoxy. The shell recovers the same local Asp contact in both ligands, yet predicts −0.019 kcal/mol against +0.928 experimentally (ordinary −0.157). ITC enthalpy favors the amino compound by 11.2 kJ/mol, opposed by 7.4 kJ/mol in TΔS. Ser195 had been proposed as a second contact. Reading the actual top poses gives amino N···Ser195-Oγ median 2.94 Å and methoxy O···Ser195-Oγ 2.85 Å, both ≤3.5 Å in 8/8 shell cells. Counting either short distance cannot separate the pair. [Series and Kd](results/series.csv) · [Ser195 measurements](results/ser195.csv) · [Talhout and Engberts, Table 1](https://doi.org/10.1046/j.1432-1327.2001.01991.x)

An independent polarizable explicit-solvent perturbation study tested benzamidine → 4-aminobenzamidine. The reported complex leg was −1.74 and the aqueous leg −1.46 kcal/mol, giving −0.28 compared with −0.615 kcal/mol from this ITC series. Our shell predicts −0.437 kcal/mol for this 2.82-fold matched edge, below the primary 3-fold cutoff. That simulation found the extra 4-amino group contacts Ser195 and one of the initial Asp189 hydrogen bonds could exchange with water. No methoxy transformation was reported there. Neutron data independently show structured waters around the bound head and mobile, incompletely hydrogen-bonded waters in the unoccupied pocket. A direct two-bond guide can recover a crystallographic motif while overrepresenting it in the solution ensemble. Testing ligand-dependent water exchange requires sampling both aqueous-ligand and solvated-complex legs, including one-contact and water-bridged configurations. [Explicit-solvent thermodynamic cycle and structural analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC2752704/) · [Neutron structures and water occupancies](https://doi.org/10.1038/s41467-018-05769-2)

## Reproduction

Run from the repository root with Python, NumPy, RDKit, Meeko, the normal receptor-preparation dependencies, and the compiled branch binary. The experimental coordinates needed under raw are public PDB structures. The two fast FXa baseline files per ligand and receptor can be regenerated by the repository's native benchmark. Raw directories can instead be restored from the accompanying evidence. The archived results keep the same stems as the commands below.

```sh
d=benchmarks/directional_ion_pair
w=/tmp/ionpair-replay
f1=$w/fxa-1FJS
f2=$w/fxa-2P16
bin=build/vina
mkdir -p "$w/raw" "$w/results"
for code in 1BTY 3PTB 5MO0 1FAX; do
  curl -fsSL "https://files.rcsb.org/download/$code.pdb" -o "$w/raw/$code.pdb"
done
curl -fsSL https://files.rcsb.org/ligands/download/DX9.cif -o "$w/raw/DX9.cif"
for code in 1FJS 2P16; do
  python benchmarks/affinity_partition/fxa_holdout.py prepare "$w/fxa-$code" --receptor "$code"
  python benchmarks/affinity_partition/fxa_holdout.py dock "$w/fxa-$code" --receptor "$code" --binary "$bin"
done
python "$d/experiment.py" "$w/fxa" "$f1" "$f2" "$bin" --arms point pair shell
python "$d/experiment.py" "$w/fxa" "$f1" "$f2" "$bin" --arms native point pair shell --seeds 161803 141421
python "$d/trypsin.py" "$w/trypsin" "$w/raw" "$bin" --seeds 271828 314159 161803 141421
python "$d/reversed_head.py" "$w/inverted" "$f1" "$f2" "$bin"
python "$d/repair.py" "$w/repair" "$f1" "$f2" "$bin"
python "$d/evaluate.py" "$w" "$f1" "$f2" "$w/results" --seeds 271828 314159 161803 141421
python "$d/contrasts.py" "$w" "$w/results" "$w/results"
python "$d/control.py" "$w" "$f1" "$f2" "$w/results"
python "$d/reference.py" "$w" "$f1" "$f2" "$w/results"
python "$d/free_ligand.py" "$w/results/free_ligand.csv"
python "$d/ser195.py" "$w" "$w/results/ser195.csv"
python "$d/verify.py" "$w" "$f1" "$f2" "$w/results"
python "$d/cross_dock.py" "$w" "$f1" "$f2" "$bin" "$w/results/dx9.csv"
python "$d/verify_cross.py" "$w" "$f1" "$f2" "$w/results/dx9.csv"
```

The audit independently parses every PDBQT model by column, infers protonated aryl amidinium nitrogens from bonds and explicit donor hydrogens, reconstructs both N–O permutations, and enumerates element and connectivity preserving atom mappings for the entire benzamidine crystal ligand. It reproduced all 9,263 raw scores and contact distances, 1,536 distinct-pair release rows (including those below threefold), and all benzamidine RMSDs to within 4 × 10⁻¹⁵ Å. Replacing the source PDB carboxylate with the normal prepared receptor PDBQT reproduced all 198 non-native search files byte for byte. [Audit output](results/verification.txt)
