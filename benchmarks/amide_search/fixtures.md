# Reproducing the amide search benchmark

From the repository root, build the FastDock Vina Pro binary with the build prerequisites in the main README. The benchmark scripts require Python 3.10 or later with NumPy, SciPy, RDKit, Meeko and Gemmi. Once these are installed, the four commands below use the checked-in data offline.

```sh
set -eu
repo_dir="$PWD"
bench="$repo_dir/benchmarks/amide_search"
base=$(mktemp -d)
work=$(mktemp -d)
python3 tools/build.py --jobs 4
binary="$repo_dir/build/vina"

(cd "$bench" && echo '0fa23896ffbf6c39602b293ee55267eabc020ec45f84f62f8c40ca2ac48c69d2  prepared-fixtures.tar.gz' | sha256sum -c -)
tar -xzf "$bench/prepared-fixtures.tar.gz" -C "$base"
(cd "$base" && sha256sum -c "$bench/prepared-fixtures.sha256")

python3 "$bench/run.py" "$binary" "$base" "$work/matrix" --jobs 4
python3 "$bench/evaluate.py" "$binary" "$base" "$work/matrix" "$work/results" --jobs 4
python3 "$bench/archive.py" "$base" "$work/results"
python3 "$bench/verify.py" "$base" "$work/matrix" "$work/results" --binary "$binary" --jobs 4
```

The archive contains 215 regular files (1,174,340 compressed bytes). The checksum manifest names each file relative to the extraction directory. The run enumerates 69 receptor–ligand combinations and 548 fixed-budget searches. The independent verifier checks the prepared search trees, aligns the experimental hCA II structures, and reconstructs pose selection and pair counts. Its optional `--binary` check recomputes the energies through ordinary single-pose scoring. Compare the generated CSVs with the committed files in `benchmarks/amide_search/results/`.

## Contents and provenance

- **Human carbonic anhydrase II.** Five deposited structures (6WQ4, 6WQ5, 6WQ7, 6WQ8 and 6WQ9), their experimental ligand instances and CCD definitions, prepared graphs and both declared ammonium geometries. Search uses 6WQ4 and 6WQ8 receptor conformations. Forty archived pose pools contain the 359 historical poses, paired to the corresponding rows of `fdv-window-fine/hcaii.csv`. The measured Ki series is reported in [the original study](https://doi.org/10.1021/acs.jmedchem.0c00733).
- **BRD4.** The original preparation manifest and four prepared ligands permit eligibility to be recomputed. Ligand 4OGI qualifies for the amide comparison on receptors 3MXF, 4OGI and 4OGJ. The measured Kd series is from [Ciceri et al.](https://doi.org/10.1038/nchembio.1471).
- **Human factor Xa.** The 28 selected ligands, prepared receptors 1FJS and 2P16, and receptor-specific search restraints reproduce the original selection in `benchmarks/affinity_partition/fxa_series.csv`. The selected start for CHEMBL369480 is the repaired geometry. The source measurements are from [Zhao et al.](https://doi.org/10.1016/S0960-894X(00)00139-6). FXa ligand inputs were adapted from EMBL-EBI ChEMBL by preparing docking coordinates and restraints. ChEMBL-derived data are distributed under [Creative Commons Attribution-ShareAlike 3.0 Unported](https://creativecommons.org/licenses/by-sa/3.0/). [Source assay CHEMBL661804](https://www.ebi.ac.uk/chembl/explore/assay/CHEMBL661804)
- The crystallographic coordinates and CCD entries are from the wwPDB archive, available under [CC0 1.0](https://www.rcsb.org/pages/usage-policy). Protein alignments, search restraints, input microstates and sampled poses are benchmark preparations.

In a clean extraction, all 24 searches for ligand 6WQ8 in receptor 6WQ4 reproduced the reference pose pools byte for byte. Common-topology rescoring reproduced all 209 mode rows for that case. Replaying the historical pools reproduced both historical CSVs byte for byte.
