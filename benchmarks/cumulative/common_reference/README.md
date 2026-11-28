# Does one energy reference choose a better pose across two searches?

Each cell pools two already paid blind searches of the same prepared query and receptor. Ordinary Vina stores an unbound energy from its own best pose. In the rigid-receptor Vina scoring function, the torsion transform is monotone and constant across poses of this same molecule. Minimizing the published INTER + INTRA therefore removes the run-specific zero. We compared this choice with the lowest published VINA RESULT, with fixed zero-reference batch re-parsing, and with scoring all saved coordinates on the original ligand topology. These are selectors on an existing pool. They do not estimate a solution-state population, binding free energy or the accuracy of one twice-as-long blind search.

For the final engine all **376/376** published modes in **21/21** two-seed comparisons accepted the common-reference rescore. For the initial Pro the corresponding numbers were **371/376** and **17/21**. Official and parent each accepted **369/374** and completed **17/21**. One out-of-domain input makes the complete-pool selector undefined; all source pools, partial rescoring, and failed runs remain in the denominator. All four engines were rescored by the very same fixed-score implementation (4dd4e37). The visible difference in coverage is upstream of the scorer.

| Final Pro cohort | Two-seed pools | Lowest printed score, mean error (Å) | Fixed topology and common zero (Å) | Released union oracle (Å) | Better / tied / worse |
|---|---:|---:|---:|---:|---:|
| Historical cognate | 10 | 6.384 | 6.777 | 4.401 | 0 / 9 / 1 |
| Historical reciprocal | 10 | 3.719 | 3.719 | 2.986 | 0 / 10 / 0 |
| External 4OCZ | 1 | 10.310 | 10.310 | 1.762 | 0 / 1 / 0 |

Only BRDT 7MRH changed winner in the final engine. The error increased from 3.606 to 7.531 Å. On 4OCZ, both the original unrounded-search component sums as printed (−12.051 vs −12.055) and the common scorer on the serialized fixed graph favored the distant seed's top pose. The two dry total energies are separated by only a few thousandths of a score unit. This is an ordering under the specified empirical potential, without a calibrated physical uncertainty. In every fully scorable pool the batch and fixed-topology routes selected the same pose. Individual lower-ranked energies can still differ. The archive check independently reconstructed 316 selections from 1,500 raw models and the actual CLI/API logs. Each version/structure/seed is a repeated condition, not an independent patient, scaffold or family.

Extra work across all four arms was 6.084 CPU s / 11.823 summed wall s for the batch zero-reference replay, 2.708 / 5.958 s for the fixed-graph replay, and 0.393 / 0.931 s for the four external minus-three controls. The two underlying dockings per condition and their CPU/wall costs are listed separately in [pairs.csv](results/pairs.csv). [All mode energies and displacements](results/scores.csv) · [all four arms and coverage](results/summary.csv) · [independent reconstruction](results/archive_audit.json) · [registration](PROTOCOL.md) · [topology control](TOPOLOGY.md).

## Reproduction

Use score-only candidate 4dd4e37 for both methods. Build it with the same tools/build.py and Boost prefix used for the four production builds. Compile the probe with `bash build.sh /path/to/4dd4e37/build/obj "$boost" "$probe"`. The C++ probe uses public Vina APIs to replace coordinates on the prepared ligand, preserving connectivity. PDBQT scoring inputs are deterministically split from the archived blind pools. Every single-model byte hash is in scores.csv.

From the repository root, with an empty work directory and absolute paths to the intermediate CLI and the probe:

```sh
D=benchmarks/cumulative/common_reference
P=benchmarks/geometry_comparison/upstream/four_versions/results
E=benchmarks/cumulative/external
mkdir -p "$work/main/raw" "$work/external/raw" "$work/fx/cognate" "$work/fx/reciprocal" "$work/fx/external"
cp "$P"/{cells,ranks}.csv "$work/main/"
cp "$E"/results/{cells,ranks}.csv "$work/external/"
tar xf "$P/raw.tar.xz" -C "$work/main/raw"
tar xf "$E/results/raw.tar.xz" -C "$work/external/raw"
tar xf benchmarks/geometry_comparison/pose/inputs.tar.xz -C "$work/fx/cognate"
tar xf benchmarks/cubic_grid/human_panel/prepared.tar.xz -C "$work/fx/reciprocal"
tar xf "$E/inputs.tar.xz" -C "$work/fx/external"
python "$D/run.py" "$work/main" "$work/fx" "$work/replay-main" "$intermediate_vina" "$probe" --fixture-subdirectories
python "$D/run.py" "$work/external" "$work/fx/external" "$work/replay-external" "$intermediate_vina" "$probe"
python "$D/collect.py" "$work/combined" "$work/replay-main" "$work/replay-external"
python "$D/verify_archives.py"
```

The pre-publication-fix candidate's 21 two-seed comparisons are archived separately in intermediate/. No selector was added to the default prediction path.
