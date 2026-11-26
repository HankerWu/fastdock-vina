# Alkylbenzene cross-docking on deposited T4 lysozyme substates

Can enumerating symmetry-related scaffold assignments improve first-pose accuracy when a homologous ligand series binds to three discrete protein conformations? A separate n-butylbenzene complex supplies the transfer geometry and all receptor structures. Each search arm sees the same prepared molecules and locations of the protein atoms for that ligand, seed and state. The six query depositions are external to the original Morton structure panel. Benzene and ethylbenzene occur as molecules in both series. The new panel tests transfer within the same protein and scaffold family.

## Structures and thermodynamics

| Query | Ligand | Complete deposited ligand alternatives and occupancies | Calorimetric Kd (μM) | Binding ΔG (kcal/mol) | Distinct donor mappings |
|---|---|---|---:|---:|---:|
| 4W52 | benzene | single, 0.70 | 175 | −5.19 | 12 |
| 4W53 | toluene | single, 0.80 | 102 | −5.52 | 2 |
| 4W54 | ethylbenzene | A 0.60, B 0.40 | 68 | −5.76 | 2 |
| 4W55 | n-propylbenzene | A 0.60, B 0.40 | 18 | −6.55 | 2 |
| 4W58 | n-pentylbenzene | A 0.70 | unavailable | unavailable | 2 |
| 4W59 | n-hexylbenzene | A 0.70 | unavailable | unavailable | 2 |

Merski et al. determined these structures and reported the same-assay calorimetric Kd series in Table 1. The published benzene ΔG and relative energies give the four binding free energies above. Solubility limited the two longer-chain measurements. The inclusion criterion is Kd(weaker)/Kd(stronger) ≥ 3, verified again with the printed ΔG differences at 302 K. The three qualifying pairs each compare n-propylbenzene with one shorter homolog. Six receptor/seed conditions repeatedly query those same three pairs. The 18 observations per arm share compounds.

[Merski et al., PNAS, Table 1](https://doi.org/10.1073/pnas.1500806112) · [Morton et al., Biochemistry, calorimetry](https://doi.org/10.1021/bi00027a006) · [4W57 experimental donor](https://www.rcsb.org/structure/4W57)

## Preparation and controls

- For all six queries, the donor is 4W57 ligand altloc A, occupancy 0.60. The intermediate (A, 0.60), open (B, 0.30), and closed (C, 0.10) receptor variants use the corresponding deposited altlocs for residues 109–113, alongside the same A/blank protein coordinates elsewhere. These occupancy labels describe the donor crystal under n-butylbenzene binding. The experiments compare three fixed receptor conditions without using the occupancies as unbound-state statistical weights.
- PDBFixer built unresolved side-chain atoms at residues 16, 53, 60, 83, 122, 125, 147, 162 and 163 once in the common background (seed 2026). For the unobserved OXT of terminal Leu164 the program places an oxygen 1.25 Å from terminal carbon in the carboxylate plane, keeping the deposited C, CA and O coordinates. The same background and terminal oxygen are used for all three states. The four peptide bonds joining residues 108/109 and 113/114 in states B/C range from 1.263 to 1.331 Å. Meeko prepares receptor hydrogens and atom types. Protein chains are rigid. Deposited solvent is omitted.
- Each query starts from a SMILES-derived RDKit conformer (embedding seed 2026 and MMFF minimization). Each ligand was assigned using a chemical graph without query crystal coordinates. All searches use the same 20 Å cube centered on the mean Cα coordinates of residues 78, 99, 102, 118, 121, 133 and 153. The nearest experimental heavy atom to a box face is 4.72 Å inside.
- Search settings use seeds 161803 and 141421, four trajectories, one CPU per run, nominal maximum 50,000 evaluations per trajectory, nine requested modes, 100 kcal/mol energy window and late native pose selection. In the three arms, generation uses unconstrained native Vina, the single canonical MCS assignment, or the complete set of graph-supported assignments. The anchored runs all use 4W57 ligand A. The emitted order and affinities use the final native Vina score. Observed search evaluations (hunt + second optimization) ranged from 50,001 to 50,112 per trajectory across arms, as the stopping criterion is checked between search iterations.
- The target frame is fitted using the 72 shared Cα positions between residues 80 and 160 inclusive after excluding 107–115. All query crystal poses are used exclusively for evaluation. Ligand RMSD is the minimum over exact chemical graph symmetries and complete deposited ligand alternatives in the receptor frame. It includes no subsequent ligand superposition. The per-altloc RMSDs and major-altloc sensitivity are also recorded.
- Environment used RDKit 2026.9.2, Meeko 0.8.0, gemmi 0.7.5, PDBFixer 1.12.0, OpenMM 8.6.1, NumPy 2.3.5 and SciPy 1.17.0. The executable used for this run had SHA-256 e4263aa6ec18dad85be5f9095c04f3cf5a6efdc01b34667667f9720bfb7cc011.

## Results

Each entry shows first-ranked successes / successes anywhere in the returned pool at ≤2 Å. A pool success is an oracle for returned geometry and requires access to the experimental structure. For each arm the denominator is 12 cells per state (six queries × two seeds) and 36 overall.

| Donor receptor state | Native | Single mapping | All mappings |
|---|---:|---:|---:|
| A, intermediate | 8 / 8 | 8 / 11 | 8 / 11 |
| B, open | 11 / 12 | 11 / 12 | 10 / 12 |
| C, closed | 8 / 9 | 8 / 9 | 8 / 9 |
| Combined | **27 / 29** | **27 / 32** | **26 / 32** |

Median first-pose RMSD across the 36 cells was 0.793, 1.087 and 1.012 Å for native, single and all mappings. Restricting evaluation to the major deposited ligand coordinates gave the same first-pose success counts. Per-ligand first-pose outcomes were 6/6 in each arm for benzene through propylbenzene. Pentylbenzene scored 1/6, 1/6, 0/6 and hexylbenzene 2/6 in each arm. For open-state pentylbenzene seed 141421, the single assignment returned a 1.945 Å first pose (−7.155 kcal/mol), compared with 2.528 Å (−7.105 kcal/mol) for all assignments. Each arm ranked the three eligible affinity pairs correctly in all six conditions.

The distinct donor improved the availability of near-native geometry for larger ligands in the intermediate receptor state. Enumerating assignments did not further increase that pool coverage or improve first-pose accuracy on the held-out depositions. The affinity subset was already perfectly ordered by native scoring and lacks the long-chain analogues responsible for most structural failures. The results therefore constrain separate explanations for sampling, final pose selection, and intratarget affinity ordering.

## Reproduction

Run from the repository root with the research branch Vina build and the Python dependencies above.

~~~bash
export OPENMM_CPU_THREADS=1 OPENBLAS_NUM_THREADS=1
python benchmarks/alkyl_holdout/prepare.py /tmp/fdv-alkyl
python benchmarks/alkyl_holdout/run.py build/ensemble/vina /tmp/fdv-alkyl --jobs 4
python benchmarks/alkyl_holdout/evaluate.py /tmp/fdv-alkyl /tmp/fdv-alkyl/results
python benchmarks/alkyl_holdout/independent.py /tmp/fdv-alkyl /tmp/fdv-alkyl/results
~~~

The independent verifier re-parses raw crystal rows, uses a separate SciPy rotational fit, reconstructs each all-carbon molecular graph from bond lengths, enumerates graph isomorphisms, and re-reads all native scores and accepted affinity-pair margins. It reproduces every crystal-alternative RMSD for all 868 emitted poses to 7.2 × 10⁻¹⁵ Å and every native score exactly. A fresh preparation and all 108 searches in a second directory reproduced the prepared-input and raw-output file hashes.

[State and combined counts](results/summary.csv) · [Ligand-level results](results/cells.csv) · [All emitted poses](results/poses.csv) · [Matched affinity pairs](results/pairs.csv) · [Independent verification](results/verification.json) · [Source hashes](results/sources.sha256) · [Prepared input hashes](results/prepared.sha256) · [Output hashes](results/runs.sha256)
