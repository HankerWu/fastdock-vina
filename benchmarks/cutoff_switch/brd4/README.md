# Held-out human BRD4 first bromodomain

These four ligands and three common receptors were selected before the mixed-search results. Ciceri et al. measured all four Kd values with the same BRD4(1) isothermal titration calorimetry protocol at 15 °C in 50 mM HEPES, pH 7.4, and 150 mM NaCl, using a single-site fit. 4OGI and 4OGJ are the original authors' co-crystals. 4O76 and 4O77 are independent BRD4(1) experimental structures from Ember et al. [ITC measurements and methods](https://doi.org/10.1038/nchembio.1471) · [Additional crystal structures](https://doi.org/10.1021/cb500072z) · [Third receptor, JQ1-bound BRD4](https://doi.org/10.2210/pdb3mxf/pdb)

| Reference | Compound | Kd, nM | Docking microstate |
|---|---|---:|---|
| 4OGI A, R78 | BI-2536 | 37 ± 3 | 7R, protonated methylpiperidine |
| 4OGJ B, 2TA | TG-101348 | 164 ± 10 | protonated pyrrolidine |
| 4O76 A, 1M3 | TG-101209 | 123 ± 18 | protonated non-aryl methylpiperazine N |
| 4O77 A, 2RE | SB-202190 | 3400 ± 130 | neutral phenol and CCD imidazole tautomer |

Both arms use one common preparation per receptor and ligand. The input charge and protonation assignments are modeled approximations. RDKit prepared ligands from explicit isomeric SMILES, using ETKDG seed 2026 and MMFF optimization. Scoring used rigid protein structures with waters and additives removed. PDBFixer 1.12.0 with OpenMM 8.6.1 and seed 2026 filled missing side-chain atoms on observed residues. No unresolved residues were added. Chain B of 4OGJ has a more complete set of observed side chains than chain A. The remaining protein atoms were prepared by Meeko. See the T4L directory for shared dependency versions.

Three common receptor structures were applied to every ligand, 3MXF A (an external JQ1-bound structure), 4OGI A and 4OGJ B. Additional cognate redocks used 4O76 A and 4O77 A. The primary 28 Å cube is centered on the mean Cα of residues 81, 82, 83, 85, 97, 140 and 146. All heavy-atom references are inside. The narrowest face margin is 0.44 Å for 4OGI ligand in the 4OGJ frame. The wider 34 Å sensitivity used all four ligands on all three common receptors, with all remaining inputs and budgets fixed. Its smallest face margin is 3.44 Å.

The evaluation aligns 123 protein Cα atoms, residues 44–166, followed by heavy-atom graph symmetry minimization in the receptor frame. Receptor fits range from 0 to 0.726 Å RMSD. Seeds 314159 and 271828, one CPU, exhaustiveness 4, max_evals 50000 per trajectory, num_modes 9 and energy_range 100 are shared within each pair. Both arms use Late selection. The mixed arm alternates search surface across trajectories via `--mixed_cutoff_search`; final refinement and reported score use the native surface.

| 28 Å endpoint | Native | Mixed |
|---|---:|---:|
| Cognate top-1 / pool ≤2 Å, 8 cells | 2 / 2 | 1 / 1 |
| Cross-dock top-1 / pool ≤2 Å, 20 cells | 1 / 2 | 1 / 1 |
| Concordant at least threefold Kd, 5 ligand pairs × 3 common receptors × 2 seeds | 10/30 | 12/30 |
| All six unique Kd pairs across the same strata | 15/36 | 17/36 |

The cognate SB-202190 case at seed 271828 shifted from top-1 1.713 Å and best 1.038 Å to top-1 6.786 Å and best 6.397 Å. The alternative method independently checked every primary and sensitivity pose against RCSB ligand-instance SDF and matched model-specific RMSDs within 0.000006 Å. The 34 Å common-template sensitivity had native/mixed top-1 success 0/24 and 0/24, pool coverage 1/24 and 2/24, and threefold-pair agreement 11/30 and 13/30. Each set of 30 rank observations reuses five experimental pairs.

```bash
export OPENMM_CPU_THREADS=1 OPENBLAS_NUM_THREADS=1
p=benchmarks/cutoff_switch/cutoff_benchmark.py
python "$p" prepare /tmp/brd4-run --panel brd4
python "$p" dock /tmp/brd4-run --panel brd4 --binary build/vina
python "$p" evaluate /tmp/brd4-run --panel brd4
python benchmarks/cutoff_switch/verify_cutoff_poses.py /tmp/brd4-run
```

Use a fresh directory and add `--box-size 34 --shared-only` to the prepare, dock and evaluate stages to reproduce the sensitivity. The `cells*.csv` and `pairs*.csv` files retain the individual strata and do not treat repeated seeds and templates as independent molecule pairs. See `alignment.json` for the rigid transformations and face margins and the T4L directory for digest conventions. The independent structure parser also compares isomeric SMILES after applying the declared protonation states.
