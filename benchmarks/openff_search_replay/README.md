# Source-stratified experimental Ki and two search-potential perturbations

This replays all three affinity partitions in the previous mixed-Vinardo research workflow on pinned OpenFF inputs (60 ligands, 3 seeds, 3 search arms, 540 docking calls). The input and docking protocols of that workflow are retained. Source measurements are grouped by target, endpoint, DOI and source table. Native and mixed use the binary at 4de709a; soft uses the formal branch with the new receptor-only steric proposal. All final scores are native Late energy estimates. Fixed receptor and prepared ligand coordinates, full PDBQT output, per-run counters and logs are in the three target archives. Reported run seconds include map generation and are sensitive to concurrent CPU work.

The frozen OpenFF source provides one starting receptor and one chemical graph per ligand. Article/table identifiers provide experimental provenance but do not guarantee equivalent measurement conditions between tables. Aggregation is within each table. Every ligand in the table must have all arms and seeds to enter the comparison. Summarize.py takes the per-ligand median over seeds, then counts every unique measured pair separated by at least threefold, once. Pairs share ligands. Different source tables are never pooled to inflate the denominator. Error values reported by the sources are preserved in the raw rows.

| Source and table | Ligands | Unique ≥3-fold pairs | Native | Mixed Vinardo | Soft receptor |
|---|---:|---:|---:|---:|---:|
| MCL1, Friberg Table 3 | 25 | 218 | 108 | 103 | 85 |
| TYK2, Liang Table 4 | 10 | 28 | 8 | 16 | 14 |
| TYK2, jm400266t Table 2 | 3 | 1 | 0 | 0 | 0 |
| PTP1B, jm0702478 Table 1 | 10 | 22 | 11 | 17 | 17 |
| PTP1B, jm0702478 Table 2 | 4 | 1 | 0 | 0 | 0 |
| PTP1B, jm0702478 Table 3 | 8 | 25 | 10 | 10 | 11 |

Friberg Table 3 reports MCL1 Ki from fluorescence-polarization competition measured in triplicate. The rows here are an OpenFF subset. The two TYK2 experiments and the three PTP1B tables remain identifiable in the downloadable edges. The full conditions have not been independently reconciled for every compound in these two targets; treat those rows as exploratory. Neither scoring method is established as a general affinity improvement. [MCL1](https://pmc.ncbi.nlm.nih.gov/articles/PMC3646517/) · [TYK2](https://doi.org/10.1016/j.ejmech.2013.03.070) · [PTP1B](https://doi.org/10.1021/jm0702478) · [frozen OpenFF records](https://github.com/openforcefield/protein-ligand-benchmark/tree/fd88824f9114244f95a14b485e6d6c96c1de716d)

## Reproduce

Use Python 3.11 with RDKit 2025.3.5, Meeko 0.6.1, Gemmi 0.7.3, PoseBusters 0.6.3, NumPy 2.4.6, SciPy 1.17.1 and PyYAML 6.0.4. Activate this environment so mk_prepare_receptor.py is on PATH. Provide a checkout of the repository at commit 4de709a6312c6424b0a78ec856f419c17c08b619 with its compiled binary, and a compiled new binary. replay.py imports the originally frozen panel implementation directly from that checkout.

```bash
B=benchmarks/openff_search_replay
W=$(mktemp -d)
python "$B/replay.py" "$W" /path/to/research-checkout /path/to/mixed/vina /path/to/new/vina
mkdir "$W/tables"
for t in mcl1 ptp1b tyk2; do cp "$W/$t/three_arms.csv" "$W/tables/${t}_rows.csv"; done
python "$B/summarize.py" "$W/tables"
diff "$B/results/summary.csv" "$W/tables/summary.csv"
```

For the pose subset use the PoseBusters authors' archive with MD5 f004ac7c4e68317b5348497d2bb6bee6. Use the Python 3.12 preparation environment described in [paired cross-docking](../paired_crossdocking/README.md) for `pose_subset.py prepare DATA OUT`, followed by the pinned Python 3.11 environment for `pose_subset.py dock DATA OUT --prior-source CHECKOUT --mixed-vina MIXED --soft-vina SOFT`. All physical checks use the original conditioning structure. Ground-truth molecules define the known-site boxes. Three runs from the same ligand are repeat measurements. This is a named human-target subset of the previous 15-case list. [Archive](https://zenodo.org/records/8278563) · [PoseBusters](https://doi.org/10.1039/D3SC04185A)

The known-site human subset has 27 completed runs on FABP3 (7WPW), JAK2 JH2 (8EX2) and PPARD (7VWF). Three scoring arms × three seeds on each ligand return zero top-1 or top-5 joint ≤2 Å/PoseBusters successes. Two prespecified preparations are exclusions: 5SB2 has an iodide 6.17 Å from the query and 7NR8 has a sulfate 2.70 Å away. The dry-protein preparation rejects a non-water omitted heavy atom within 8 Å. Raw traces, the excluded groups and their actual distances, prepared structures and individual PoseBusters verdicts are in posebusters-human.tar.xz. The other ten cases in the earlier 15-case list were outside this human subset. The original GitHub workflow failed at installation because Meeko 0.6.1 requires Python <3.12. The interpreter pin has been corrected on its research branch; all scientific measurements reported here were reproduced locally.
