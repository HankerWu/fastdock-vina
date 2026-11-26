# Small-ring conformer selection on eight human drug targets

This benchmark tests whether committing half a fixed search budget to a distinct ring shape improves recovery of an experimental pose. In two paired experiments the strategy decreased top-1 recovery compared with repeating the lowest-energy input. These results support keeping state diversification optional. A fivefold increase in the maximum evaluation count helped more than the tested conformer rules. The latter comparison spends additional compute.

## Receptors and experimental observations

One query molecule per protein was drawn from human drug-target structures outside the earlier hCA II, BRD4, ABL and FXa panels. Each query starts from the wwPDB chemical graph. Deposited coordinates enter only pocket localization, structural evaluation and retrospective mechanism checks. The other-ligand receptor box is centered using its own ligand, never the experimental pose of the query.

| Human target | Query PDB and ligand | Other receptor | Pocket Cα count / fit RMSD (Å) |
|---|---|---|---|
| HSP90α | [2VCI](https://www.rcsb.org/structure/2VCI), 2GJ | [2XAB](https://www.rcsb.org/structure/2XAB), different ligand | 76 / 0.161 |
| DPP4 | [1X70](https://www.rcsb.org/structure/1X70), sitagliptin | [2RGU](https://www.rcsb.org/structure/2RGU), different ligand | 81 / 0.228 |
| p38α | [1KV2](https://www.rcsb.org/structure/1KV2), BIRB 796 | [3GCP](https://www.rcsb.org/structure/3GCP), different ligand | 61 / 1.054 |
| PARP1 | [5DS3](https://www.rcsb.org/structure/5DS3), olaparib | [4ZZZ](https://www.rcsb.org/structure/4ZZZ), different ligand | 61 / 0.796 |
| PDE5 | [1UDT](https://www.rcsb.org/structure/1UDT), sildenafil | [2H42](https://www.rcsb.org/structure/2H42), same ligand in a separate deposit | 77 / 0.624 |
| CDK6 | [5L2I](https://www.rcsb.org/structure/5L2I), palbociclib | [2EUF](https://www.rcsb.org/structure/2EUF) excluded | See below |
| Glucocorticoid receptor | [1M2Z](https://www.rcsb.org/structure/1M2Z), dexamethasone | — | — |
| Progesterone receptor | [1A28](https://www.rcsb.org/structure/1A28), progesterone | — | — |

Rigid receptor preparation used observed target-chain atoms, a 26 Å cubic box centered on the experimental binding site, and retained deposited Zn/Mg ions in PDE5. Waters and ancillary protein chains were omitted. Missing sidechains on observed residues within 8 Å of the deposited ligand were repaired using PDBFixer when Meeko omitted the residue. The 2EUF receiver still lost an entire residue 4.20 Å from the ligand after repair, so it was excluded before docking. Full inputs, reconstruction logs and receptor transforms are in `fixtures`. The eight query molecules provide chemically distinct examples, including two steroid hormones. This design has one ligand per target, so target/family and ligand chemistry are confounded. It contains no within-target affinity benchmark.

## Physical hypothesis and protocol

The native small-ring PDBQT representation fixes ring internal geometry. Rotations about permitted tree branches can orient a ring and its substituents in the binding site, but cannot switch between arbitrary endocyclic puckers. McNutt and colleagues found multi-conformer docking improved recovery over matched sampling from a single input in their protocol. Kubincová and colleagues found limited small-ring pucker recovery in selected alternative-mode structures. The mechanism allows a gain if an additional input samples a productive, previously unreachable shape. Diverting a trajectory from an already productive state also has an opportunity cost. [McNutt et al.](https://doi.org/10.1021/acs.jcim.3c01245) · [Kubincová et al.](https://doi.org/10.1021/acs.jcim.6c00356) · [Cremer–Pople analysis of drug-like rings](https://doi.org/10.1021/acs.jcim.0c01144)

`tools/ring_conformers.py` discards supplied 3D coordinates, embeds 64 conformers using ETKDGv3 with small-ring torsions (seed 20261126), and minimizes each with MMFF94s. The state with the lowest gas-phase energy among converged conformers supplies `base.pdbqt`. The lowest-energy alternative within 12 kcal/mol must change the heavy-atom geometry of at least one nonaromatic five- or six-membered ring including its first heavy-atom neighbors by 0.35 Å or more. Local proper-fit RMSDs are maximized over ring neighborhoods and minimized over chirality-preserving whole-molecule graph symmetries. Including first-shell neighbors captures changes in substituent orientation as well as endocyclic pucker. It also admits geometry changes involving the first shell. A sham start changes full-ligand geometry by at least 0.5 Å but all measured ring neighborhoods by at most 0.15 Å, up to graph symmetry. If no alternative qualifies the slot reuses `base`. Atom types, charges, atom identities and the torsion tree must match in all inputs. The CCD ionization and stereochemistry were held fixed for each ligand. Solution populations and assay-specific protonation were not measured.

Every arm has two independent Late search slots, each running three trajectories on one CPU, refining up to twelve candidates and outputting up to nine. The two slots have identical seeds across arms. Four pairs of seeds per receiver produce 52 paired cells per arm. Reused identical slots count only once in the number of actual binary invocations. Two evaluation caps, 12,000 and 60,000 per trajectory, were tested independently. The arms are (i) two searches from the first converged conformer, (ii) two searches from `base`, (iii) one from `base` and one from the alternative, and (iv) one from `base` and one from the sham. The first-conformer control was added after the primary alternative-vs-repeat experiment to examine the effect of energy selection. The docking budgets match within each cap. The 64-conformer preprocessing has an additional cost over generating the first conformer alone.

The primary pool merge ranks the printed `INTER + INTRA` values. Within one ligand, identical torsional topology gives the same monotone conformation-independent normalization, so this ranks candidates at a common unbound reference. It avoids importing run-specific best-pose offsets from separately printed `VINA RESULT` lines. The standard potential omits internal energy that is constant within a rigid fragment. We therefore report both the uncorrected common-reference ranking and a sensitivity adding the input MMFF energy difference. This full-input, gas-phase difference includes acyclic conformational effects and lacks solvent, receptor relaxation and intrabasin entropy. It should not be interpreted as a binding free energy or a pure ring-strain penalty. The merged set is scored against heavy-atom experimental RMSD in the receptor frame with graph symmetries; no ligand superposition is allowed. The oracle is the lowest RMSD among returned modes and is unavailable to a prospective selection algorithm.

## Paired results

| Search cap | Input strategy | Top-1 ≤2 Å / 52 | Oracle ≤2 Å / 52 | Median top-1 RMSD (Å) |
|---:|---|---:|---:|---:|
| 12,000 | First conformer twice | 15 | 15 | 6.47 |
| 12,000 | Lowest energy twice | 19 | 24 | 3.96 |
| 12,000 | Lowest energy + distinct ring | 17 | 19 | 5.99 |
| 12,000 | Lowest energy + sham | 14 | 17 | 5.90 |
| 60,000 | First conformer twice | 28 | 32 | 1.47 |
| 60,000 | Lowest energy twice | 25 | 31 | 2.99 |
| 60,000 | Lowest energy + distinct ring | 22 | 29 | 3.04 |
| 60,000 | Lowest energy + sham | 19 | 25 | 3.26 |

At 12,000 evaluations the distinct-ring arm gains one top-1 success and loses three relative to repeating `base`. At 60,000 it gains two and loses five. Weighting the eight targets equally yields changes of −3.125 and −4.6875 percentage points. At the larger cap, cognate top-1 recovery changes from 21 to 18 of 32, different-ligand cross-docking changes from 4 to 3 of 16, and the independent sildenafil structure changes from 0 to 1 of 4 (oracle 2 to 4). Each seed repeats an existing physical complex. These observations quantify variability on a selected panel and do not give 52 independent examples of transfer.

Six of eight ligands generated a distinct alternative. Only sildenafil had improved local experimental geometry with that alternative, from 0.531 to 0.386 Å ring-neighborhood error. For sitagliptin it went from 0.174 to 0.523 Å despite a gas-phase difference of only 0.022 kcal/mol, and for the HSP90α query it rose from 0.196 to 0.488 Å. In the p38α query the first conformer cost 15.80 kcal/mol more than `base` and had a larger local ring error (0.374 versus 0.257 Å). At the larger cap it nevertheless recovered the cognate top pose in all four seeds, where `base` recovered zero. A full-ligand high gas-phase energy is insufficient as a rule for discarding an input geometry used to propose docking poses. The productive difference may involve other frozen coordinates or finite search. This pairing does not isolate a causal ring distortion at p38α.

For the 60,000-evaluation alternative arm, native run-specific offsets, the primary common reference and the extra MMFF sensitivity yield 21, 22 and 19 top-1 hits, respectively. The paired deterioration of the candidate-pool ceiling under the primary mechanism also appears at both caps. Re-ranking a fixed pool cannot account for that decrease. The low-energy preprocessing result itself reverses between caps, with +4 then −3 top-1 hits relative to the first-conformer control. These data do not support a default replacement of either the conformer or the scoring protocol.

## Reproduction and audit

The archived inputs were produced with Python 3.12, RDKit 2026.9.2, Meeko 0.8.0, Gemmi 0.7.5, NumPy 2.3.5, SciPy 1.17.0, OpenMM 8.6.1 and PDBFixer 1.12.0. The FastDock binary was built from the `79c793a` source tree (`f04fa353`). A standalone ligand preparation uses, for example, `python tools/ring_conformers.py --smiles 'O=C(Nc1ccccc1)C1CCCCC1' --out ring-inputs`. It writes the four named PDBQT inputs and `states.json`, including energies, conformer IDs, and whether a distinct state was found.

Run from the repository root, after compiling `build/vina` and installing the Python packages above. `prepare.py` can replace the input-fixture extraction if a fresh download and protein preparation are wanted. Each cap has its own work directory. Running the full benchmark without `--evaluate-existing` regenerates every docking pool. Preserve source/target paths when unpacking.

```bash
B=benchmarks/ring_ensemble
(cd "$B" && sha256sum -c fixtures.sha256)
for cap in 12000 60000; do
  work="/tmp/ring-docking-$cap"
  mkdir -p "$work"
  tar -xzf "$B/fixtures/structures-prepared.tar.gz" -C "$work"
  tar -xzf "$B/fixtures/poses-$cap.tar.gz" -C "$work"
  python "$B/run.py" "$work" --evaluate-existing --out "$work/results" --max-evals "$cap"
  python "$B/summarize.py" "$work" "$work/results"
  python "$B/verify.py" "$work" "$work/results"
done
# Regenerate one budget from the prepared inputs
python "$B/run.py" /tmp/ring-docking-60000 --binary build/vina --max-evals 60000 --out /tmp/ring-redocked
```

The independent script fetches or reads the original RCSB ligand-instance SDF, reconstructs models with Meeko, superposes independently parsed receptor coordinates with SciPy and computes uncompensated-frame `RDKit CalcRMS`. It verified all 304 unique runs at each cap, totaling 5,426 models. The maximum disagreement with the primary RMSD was 0.000000969 Å. Per-case scores, all modes and per-file digests appear under `results/` and `results/60000/`. The prepared archive contains source PDB/CCD/SDF, alternate-receptor exclusion details, state metadata, prepared molecules, boxes and receptor-preparation logs. `fixtures/poses-*.tar.gz` contains every original native output. The SHA-256 files distinguish unique runs from the intentional reuse of identical inputs in some experimental arms.
