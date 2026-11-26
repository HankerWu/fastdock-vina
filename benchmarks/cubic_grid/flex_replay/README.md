# Selective side-chain flexibility in two diagnosed pockets

A crystal-informed diagnostic identified two side chains with short contacts to a projected ligand. FABP4 Thr60 and RBP4 Gln98 were made flexible in both receivers of each target. The 48 searches used upstream Vina bb64e5e, three fresh seeds, two ligands and two receivers per target, a 24 Å known-site box and three trajectories capped at 24,000 evaluations each. The baseline and rigid/flexible input partitions contain precisely the same atom identities, coordinates, atom types and partial charges. Both arms use default Vina scoring and refinement. Site selection used the experimental pose, making this a retrospective mechanism test.

| Target | Rigid first / top 5 / all | Flexible first / top 5 / all |
|---|---:|---:|
| FABP4, 12 cells | 4 / 6 / 6 | 4 / 5 / 6 |
| RBP4, 12 cells | 1 / 3 / 3 | 1 / 2 / 2 |

The total number of first-pose hits stayed at 5/24. Half the cells were cognate, half cross-ligand. In the RBP4 cross subset the top-five count rose from 0 to 1 of 6. Gln98 and Thr60 actually moved, with first-pose median displacements of 1.99 and 0.45 Å for side-chain heavy atoms excluding Cα. [425 modes](modes.csv) · [Independent atom audit](../official_audit/README.md) · [Upstream flexible docking](https://autodock-vina.readthedocs.io/en/latest/docking_flexible.html)

Extract the human_panel prepared archive into fixtures. The included archive holds rigid/flexible PDBQT, raw pools, docking and preparation logs. Use the human panel Python environment.

```bash
python benchmarks/cubic_grid/flex_replay/run_flex.py fixtures flex_results /path/to/official/vina
```
