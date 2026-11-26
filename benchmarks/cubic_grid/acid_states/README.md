# Carboxylic acids in upstream Vina and AD4

Two FABP4 ligands and RBP4 ligand A1120 were each tested in two receivers. Each of the four arms used the same three new seeds, a 24 Å known-site box, 0.375 Å maps and three trajectories capped at 24,000 evaluations each. Both chemical states share every initial heavy-atom coordinate. Removing the acid proton changes the donor, terminal torsion, partial charges and Vina atom types. Neither chemical state has a measured population in these preparations. Crystal evidence for arginine recognition of acid head groups motivated the comparison. [RBP4](https://doi.org/10.2210/pdb3FMZ/pdb) · [FABP4](https://pmc.ncbi.nlm.nih.gov/articles/PMC4834654/) · [Force-field definitions](https://doi.org/10.1021/acs.jcim.1c00203)

| Official score and chemical state | Top-1 / 18 | Top-5 / 18 | Any returned mode / 18 |
|---|---:|---:|---:|
| Vina, neutral acid | 6 | 7 | 7 |
| Vina, anion | 5 | 7 | 7 |
| AD4, neutral acid | 4 | 7 | 7 |
| AD4, anion | 6 | 7 | 7 |

On RBP4 the anion produced zero successful candidates in six cells with either force field. On FABP4 3FR2, AD4 top-1 improved from 4/6 to 5/6. Both chemical graphs were evaluated using a common RMSD with oxygen exchange allowed, with a separate analysis retaining atom-labeled oxygens. No threshold decisions changed. [637 modes](modes.csv) · [Independent checks](../official_audit/README.md)

Upstream binary bb64e5e and Debian autogrid 4.2.6-10 amd64 (package SHA-256 751e5256ce2d8452abdddff83abb2ce92d43e629ccd9cf84daa0e354f24149b7). Grid recipes, logs, ligand states and output pools are archived here. Regenerate maps from the recipes and the human panel receptors and compare to maps.sha256. Extract human_panel/prepared.tar.xz into fixtures. Use the human panel Python environment.

```bash
python benchmarks/cubic_grid/acid_states/run_states.py fixtures acid_results /path/to/official/vina /path/to/autogrid4
```
