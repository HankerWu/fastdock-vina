# Human carbonic anhydrase II paired docking study

Bonardi et al. reported five crystal complexes and stopped-flow CO₂ hydration inhibition constants for the same human CA II assay at pH 7.5. The paper gives means of at least three determinations. The authors estimate errors of 5–10%. Correspondence comes from Tables 1 and 2. [Primary paper](https://doi.org/10.1021/acs.jmedchem.0c00733)

| PDB | Compound | CCD | Ki (nM) |
|---|---:|---|---:|
| 6WQ4 | 34 | U7A | 368.5 |
| 6WQ5 | 41 | U7G | 57.0 |
| 6WQ7 | 42 | U6V | 30.4 |
| 6WQ8 | 46 | U77 | 2.4 |
| 6WQ9 | 48 | U7J | 72.5 |

The main panel docks each of the five ligands into 6WQ4 and 6WQ8. Three extra redocks complete cognate controls for 6WQ5, 6WQ7 and 6WQ9. Seeds are 42 and 1729. Each job uses one CPU, exhaustiveness 4, up to 9 modes, an energy range of 100 and a nominal 50,000 evaluations per trajectory. The paired command adds only `--smooth_cutoff`. The initial ligand geometry comes from RDKit ETKDGv3 (seed 9033) and MMFF minimization of the CCD graph.

## Preparation and interpretation

The receptor retains chain A, altloc A or blank, and active-site Zn. Protein hydrogen placement and AutoDock typing use Meeko. His94/His96 use HID and His119 uses HIE, leaving the observed Zn donors free. Water, buffer, additional ligands and crystal contacts are omitted. The ligand graph uses a deprotonated coordinating sulfonamide, ionized carboxylate and protonated free aliphatic amines as applicable. These are assigned microstates. Neither a microstate free energy correction nor induced fit enters the score.

The box center is 5 Å from Zn in the open direction defined by the three coordinating histidines. The first 26 Å cube put portions of aligned reference ligands up to 1.24 Å beyond the box. The primary rerun enlarged the cube to 36 Å for all queries, with all reference atoms at least 3.76 Å inside. Report the 26 Å runs as a box sensitivity analysis. The two main receptor conformations have a common chain-A Cα RMSD of 0.149 Å.

Reference poses use the active-site ligand in chain A at residue 303. For 6WQ8, this distinguishes the Zn-bound instance from the second copy at residue 304. Pose labels use all ligand heavy atoms, mapped from CCD atom names through the prepared ligand’s PDBQT index map. Compare coordinates after aligning common receptor Cα atoms by Kabsch. RDKit CalcRMS permits graph symmetry and conjugated terminal-group equivalences while retaining the receptor frame. Some peripheral tails have weak experimental density.

The rank endpoint counts correctly ordered pairs with a Ki ratio of at least 3 within each receptor and seed. A tied predicted score receives half credit. Four repetitions share the same seven unique pairs. The geometry diagnostic reads the coordinating nitrogen from the original crystallographic LINK records and checks its distance to receptor Zn for every MODEL.

## Reproduce

The recorded environment used Python 3.12.14, RDKit 2026.09.2, Meeko 0.8.0, Gemmi 0.7.5, NumPy 2.3.5 and SciPy 1.17.0. The preparation command downloads or verifies the input digests in [sources.json](sources.json), including the full-text XML used to extract Ki. Preparation and docking create fresh output directories.

```bash
scripts=benchmarks/cutoff_switch/hcaii
work=/tmp/fastdock-hcaii-primary
python "$scripts/prep.py" --work "$work"
python "$scripts/run.py" --work "$work" --binary build/vina --jobs 4
python "$scripts/run.py" --work "$work" --binary build/vina --jobs 4 --extra-cognate
python "$scripts/diagnose.py" --work "$work"
```

For the original small-box comparison use a separate work directory and add `--box-size 26` to preparation. Existing outputs can be re-evaluated with `run.py --evaluate-only` (also `--extra-cognate` for the three extra redocks). [Rows and PDBQT hashes](rows.csv) · [Group summaries](groups.csv) · [Pairwise decisions](pairs.csv) · [Every recovered primary-box pose](all_modes.csv)

## Crystal-initialized local refinement

```bash
python "$scripts/crystal_start.py" --work "$work" --binary build/vina --jobs 3
```

This diagnostic uses the known cognate crystal coordinates. It projects each reference onto the same independently prepared ligand by fitting its rigid placement and PDBQT branch torsions. Internal bond lengths, bond angles and rigid-fragment geometries match the original search input apart from PDBQT coordinate rounding. RDKit supplies the receptor-frame symmetry-aware RMSD evaluation. Meeko topology and atom-index remarks supply an explicit checked bond list for refinement and every fixed-coordinate score.

The projected poses lie 0.27–0.44 Å from the crystal references. All five remain within 2 Å after local minimization under either potential. Native-refined poses span 0.57–0.87 Å and smooth-refined poses span 0.80–1.09 Å. For every ligand at both seeds, the best rescored member of the original native docking pool obtains a lower inter+intra score than the refined crystal-initialized candidate. The latter score is higher by 4.18–6.60 native units or 3.20–5.23 smooth units. The two seed-specific comparisons in each arm share one refined reference pose. These are comparisons within one ligand, in each arm's own objective. Other near-native local minima may exist. [Matched energy comparisons](crystal_start.csv)

Both components favor the selected decoy in all ten comparisons per arm. Under native scoring, the refined near-crystal ligand has 1.54–4.42 higher intermolecular and 1.09–3.08 higher intramolecular energy. Under smooth scoring these intervals are 1.51–3.80 and 1.05–2.35. The intermolecular difference is larger in eight of ten comparisons under either potential. This separates the two computed contributions but gives no independent estimate of their physical accuracy.
