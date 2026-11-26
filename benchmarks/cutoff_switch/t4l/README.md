# Matched T4 lysozyme cavity experiments

The primary experiment used the C54T/C97A/L99A construct throughout. The nine complexes are the 181L–188L and 1NHB series of Morton and Matthews. Each ligand was docked to its own receptor and to two common receptors, apo 1L90 and n-butylbenzene-bound 186L. The 186L cognate case appears only once per seed. The PDB headers and chain coordinates provide construct identity. All inputs use chain A, with waters and other additives removed.

| PDB | Ligand | ITC Kd, µM |
|---|---|---:|
| 181L | Benzene | 175 |
| 182L | Benzofuran | 112 |
| 183L | Indene | 193 |
| 184L | Isobutylbenzene | 19 |
| 185L | Indole | 290 |
| 186L | n-Butylbenzene | 14 |
| 187L | p-Xylene | 422 |
| 188L | o-Xylene | 470 |
| 1NHB | Ethylbenzene | 68 |

The affinity experiment is the single ITC series of Morton, Baase and Matthews. The tabulated Kd values use the BindingDB transcription of the Graves et al. collection. Among the 36 unordered pairs, 19 differ by at least threefold. Ligand poses from the companion crystallographic study are the structural endpoint. [Original measurements](https://doi.org/10.1021/bi00027a006) · [Crystal structures](https://doi.org/10.1021/bi00027a007) · [Graves et al.](https://doi.org/10.1021/jm0491187) · [BindingDB records](https://www.bindingdb.org/rwd/jsp/dbsearch/PrimarySearch_ki.jsp?tag=tg&target=Endolysin+%5BL99A%5D&kiunit=nM&icunit=nM&column=ki&submit=Search)

As an independent label check, the nine matched records in the Mobley Lab T4L compilation all identify the same original ITC DOI. Recomputing RT ln(Kd/1 M) at 302.15 K gives all nine reported free energies to within 0.018 kcal/mol, consistent with tabular rounding. The individual values and reported errors are in `affinity_check.csv`. [Machine-readable cross-check](https://github.com/MobleyLab/lysozyme_binding/blob/master/L99A/T4_ligand_information.json)

Meeko prepared each complete protein and the neutral ligands. Each ligand starts from its SMILES in a fresh ETKDG conformer, seed 2026, with MMFF optimization. The 16 Å cubic box is centered on the mean Cα of residues 78, 99, 102, 118, 121, 133 and 153. The C-terminal domain framework, 72 Cα atoms at 80–106 and 116–160, supplies the rigid transform for the experimental reference. The primary evaluator uses the CCD heavy-atom graph and enumerates symmetry-equivalent atom maps without a ligand fit. The independent check uses RCSB ligand-instance SDF, Meeko reconstruction, SciPy rotation and RDKit CalcRMS.

The paired configurations were standard Vina and the all-stage smooth cutoff, seeds 42 and 1729, exhaustiveness 4, 50,000 maximum evaluations per trajectory, one CPU, up to nine modes, 100 kcal/mol energy window. The output files cover 52 pairs of runs. Among 18 cognate cells, top-1 RMSD ≤2 Å fell from 18 to 15, with at least one acceptable pose in all 18 pools in either arm. In 34 cross-docking cells, top-1 remained 25 and pool coverage rose from 31 to 33. For the 19 threefold Kd pairs evaluated separately in each of two conformers and seeds, the totals were 50/76 and 49/76. The 76 observations reuse the same 19 molecule pairs.

From the repository root, use an empty work directory and the same executable for both arms. The tested Python stack includes RDKit 2026.09.2, Meeko 0.8.0, Gemmi 0.7.5, NumPy 2.3.5 and SciPy 1.17.0.

```bash
export OPENBLAS_NUM_THREADS=1
p=benchmarks/cutoff_switch/cutoff_benchmark.py
python "$p" prepare /tmp/t4l-run
python "$p" dock /tmp/t4l-run --binary build/vina
python "$p" evaluate /tmp/t4l-run
python benchmarks/cutoff_switch/verify_cutoff_poses.py /tmp/t4l-run
```

Run `sha256sum -c` with the absolute path to `sources.sha256`, `prepared.sha256` or `poses.sha256` from the generated `raw/`, `prepared/` or `runs/` directory respectively. The per-arm CSV includes each output model's RMSD in ranked order, the original PDBQT top score, actual gradient evaluation counts and an explicit-refinement counter. `pairs.csv` records every within-receptor, within-seed pair. A positive score gap means the lower experimentally measured Kd received the more favorable prediction. Ties receive one half.
