# Galectin Kd ranking and transfer of a conserved sugar core

Peterson et al. measured inhibitors 5–17 against human galectin-1 and galectin-3 in one competitive fluorescence polarization study. The 26 dissociation constants and reported errors in Table 1 are transcribed in [table1.csv](table1.csv). Every central value and exact-equality operator was independently checked against the GFA release. This gives 13 compounds per target, 23 unique galectin-1 pairs and 35 galectin-3 pairs with a difference of at least threefold. Both lectins recognize a related sugar scaffold and the two target rows share all their compounds. [Peterson et al., original Table 1](https://doi.org/10.1021/acs.jmedchem.7b01626) · [GFA dataset, pinned revision](https://github.com/choutkaj/GFA-ligand-set/tree/04dd9b0024ba3a81730731d256b6e7b285977aec)

Each of the five arms uses the same two trajectories, 25,000-evaluation per-trajectory cap, 24 Å cube, native Late refinement, nine requested modes and one CPU. Four matched seeds and two experimental receptor states give 208 searches per arm, 1,040 total, and 9,193 returned poses. All 13 ligands start from SMILES. ETKDGv3 generated eight conformers per ligand. The lowest-energy converged MMFF conformer supplied the initial internal coordinates. No ligand in Table 1 served as the donor. The receptor-bound galactose or talose located the docking region in each crystal.

- Native and balanced randomize initial ligand torsions. Balanced distributes each individual torsional proposal between both sides of the bond with a compensating rigid transformation.
- Prepared allocates its second trajectory to the initial MMFF-derived torsions and random rigid placement. Both trajectories then explore torsions.
- Anchored uses the existing search-stage positional potential on six sugar-ring atoms in both trajectories. Balanced + anchored combines the same restraints with balanced torsional proposals. These two diagnostic controls were added after observing the unanchored core coverage. Every final refinement and reported score uses the native unrestrained objective.

## Affinity

For each compound and receptor state, take the median of the four top-pose scores. The primary summary uses the arithmetic mean of the two state medians as a consensus score. An equilibrium association-constant average would require state populations and calibrated state energies. Exponential pooling with equal weights, each individual state, two separate seed cohorts, and the lowest native score among core-supported returned poses are provided in [summary.csv](results/summary.csv). The latter uses the native best score when the pool lacks any core-supported pose.

| Target and compounds | Unique pairs | Native | Balanced | Prepared | Anchored | Balanced + anchored | More heavy atoms wins |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Gal-1, complete Table 1 | 23 | 17 | 20 | 16 | 20 | 19 | 16.5 |
| Gal-1, triazoles 5–16 | 13 | 7 | 10 | 6 | 10 | 9 | 6.5 |
| Gal-3, complete Table 1 | 35 | 28 | 28 | 25 | 28 | 28 | 34 |
| Gal-3, triazoles 5–16 | 25 | 18 | 18 | 15 | 18 | 18 | 24 |

Compound 17 lacks the pendant aryltriazole and is the weakest binder for both proteins. Excluding it leaves the 12-member substitution series. On Gal-1 the two seed cohorts change the native versus balanced results from 10 to 7 and from 6 to 11 correct out of 13. The four-seed 7 to 10 gain has a compound-bootstrap descriptive 95% interval of −22 to +100 percentage points. On Gal-3 the same comparison is 17 to 12 and 20 to 21 out of 25. The receptor state alone changes its direction, from 14 to 17 in 1KJL and 21 to 13 in 3T1L. Four thousand whole-compound bootstrap draws were attempted per comparison. Draws without any eligible pair were excluded. The small number of compounds and shared endpoints make the uncertainty large. [Paired edges](results/pairs.csv) · [Compound-level bootstrap](results/paired_bootstrap.csv)

Compounds 9 and 13 reverse order between the proteins. Gal-1 prefers 13 over 9 (100 vs 530 µM), while Gal-3 prefers 9 over 13 (8.8 vs 44 µM). Their measured weak-minus-strong gaps at the conventional conversion temperature of 298.15 K are +0.988 and +0.954 kcal/mol. Native gives +0.048 and +0.215, balanced +0.480 and +0.298, and the anchored comparator −0.049 and −0.128. The crystal experiments in the paper concern different, disaccharide analogues 36 and 37. Their multipolar fluoro-amide contacts motivate a hypothesis about the aryl substituents. They do not establish the exact binding mode of monogalactosides 9 and 13.

## What the donor coordinates test

The homology measurement uses the mapped C1–C5 and O5 ring atoms in receptor coordinates with no additional least-squares fit. For three donors this is a galactosyl ring. In 3T1L the donor is a taloside with different C2 stereochemistry. Only the six ring atoms are transferred. Counts below refer to a ring RMSD of at most 2 Å, with 104 receptor–compound–seed conditions per target and arm. This diagnostic has no query-specific full-ligand crystal coordinates.

| Target | Native top / any | Balanced top / any | Prepared top / any | Anchored top / any | Balanced + anchored top / any |
| --- | ---: | ---: | ---: | ---: | ---: |
| Gal-1 | 32 / 38 | 30 / 40 | 29 / 41 | 104 / 104 | 103 / 103 |
| Gal-3 | 4 / 8 | 6 / 11 | 9 / 13 | 89 / 104 | 94 / 104 |

On Gal-3 the existing anchor raises common-core coverage from 8 to 104 pools. Its top scores still correctly rank 18 of 25 threefold triazole pairs. The size control gets 24. Substituent orientation, multipolar recognition, solvent displacement, flexibility and the reference energy remain candidate explanations. None has been identified as the cause by this experiment. [Per-state counts](results/core_counts.csv) · [All returned poses](results/poses.csv)

## Receptor and energetic assumptions

The prepared rigid receptor chain and independent donor are 1GZW B / lactose D, 4Q26 A / LacNAc C, 1KJL A / LacNAc B, and 3T1L A / taloside MQT A. Alternative coordinates use blank or A. The 4Q26 cysteine adducts at residues 16, 88 and 130 are reduced to standard Cys, and PDBFixer completes the unobserved SG of Cys2. The 1GZW buffer-linked cysteines are represented in their reduced state. The rigid models omit solvent, ions and other chains. For Gal-3, His158 uses the neutral Nδ1-H tautomer reported for the binding site. Except for the known-site donor, all test poses and all scores are generated by the branch executable. [1GZW](https://www.rcsb.org/structure/1GZW) · [4Q26](https://www.rcsb.org/structure/4Q26) · [1KJL](https://www.rcsb.org/structure/1KJL) · [3T1L](https://www.rcsb.org/structure/3T1L) · [Histidine assignment](https://doi.org/10.1002/pro.70143)

A fixed-coordinate rescore recomputed all 1,040 top energies from the serialized geometry using an independent CLI path. Maximum absolute discrepancy was 0.0113 kcal/mol. Every one of the 58 target–chemical-pair decisions in each of the five arms matched the original sign. Giving each ligand a common prepared-conformer internal reference changes some rankings. For Gal-3 triazoles the counts become native 14, balanced 13, prepared 12, anchored 11, and balanced + anchored 11 out of 25. This reference tests sensitivity to the estimated intramolecular energy of the pose pool. It omits solution entropy and solvation. A second preparation from the frozen source files reproduced all 74 prepared inputs, boxes, anchor files and the panel table byte for byte. All 1,040 raw PDBQT output files were then reproduced byte for byte with the final binary. [Independent energy checks](results/readback_pairs.csv) · [Reproduction hashes](results/repeated.json)

## Reproduction

Use Python 3.12 with RDKit 2026.9.2, Meeko 0.8.0, PDBFixer 1.12.0, OpenMM 8.6.1, Gemmi 0.7.5, SciPy 1.17.0, and NumPy 2.3.5. From the repository root, point VINA to the branch build with the two opt-in torsion switches.

~~~bash
B=benchmarks/galectin_affinity
W=$(mktemp -d)
export OPENMM_CPU_THREADS=1 OPENBLAS_NUM_THREADS=1
python "$B/prepare_inputs.py" "$W"
python "$B/dock.py" "$W" "$VINA" --jobs 3
python "$B/evaluate.py" "$W" "$W/results"
python "$B/audit.py" "$W" "$VINA" "$W/results" --jobs 3
~~~

For offline replay, verify the checksum of [evidence.tar.xz](evidence.tar.xz) with [evidence.sha256](evidence.sha256), then extract the archive into the work directory. It contains the pinned source rows and crystals, input structures, complete PDBQT pose pools, logs, and per-run hashes. The first four arms used the binary recorded as early_binary_sha256 in [repeated.json](results/repeated.json). Prepared and the independent five-arm replay used its final two-switch build. Per-process wall times from the concurrent run are descriptive.
