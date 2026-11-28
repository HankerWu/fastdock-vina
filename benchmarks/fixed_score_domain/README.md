# Fixed-coordinate scores share one domain

A caller who supplies an unbound reference changes the zero of the empirical score. The grid on which the receptor field is defined still limits the valid pose coordinates. Before `4dd4e37`, the single-ligand CLI passed the supplied reference to `Vina::score(double)` without the ligand, map, or box checks used by `Vina::score()`. A batch performed the calculation but discarded its returned energy vector. These defects affect downstream users exporting fixed scores, for example before calculating a ranking. They do not explain a blind search that already produced an in-domain but wrongly ranked pose.

## Registered result

The explicit XS receptor implementation evaluates its field at the componentwise projection of the atom coordinates onto the box. It additionally charges the L1 distance beyond the box, times a slope of 10⁶ score units per Å. With all scored atoms past the same face, moving one more Å outward leaves the projected coordinates, ligand internal geometry, and common reference fixed. The raw ligand–receptor term must increase by exactly one slope per atom. This prediction concerns the unnormalized receptor term; the displayed total also contains Vina's conformation-independent normalization. [Explicit receptor field](../../src/lib/non_cache.cpp) · [Original Vina scoring definition](https://doi.org/10.1002/jcc.21334) · [Registration](PROTOCOL.md)

| Receptor and ligand | Bounded atoms | Measured old raw receptor change, +40 Å to +41 Å | Valid in-box total, native / 0 / −3 reference |
|---|---:|---:|---|
| Mouse ABL 1IEP, imatinib | 37 | +37,000,000 | −12.513 / −12.857 / −10.728 |
| Human BPTF 7F5D, DC-BPi-03 | 21 | +21,000,000 | −6.870 / −7.264 / −4.712 |

The increments replicated at both fixed references. Before repair, the four exterior automatic-reference queries failed and the eight exterior fixed-reference queries returned large positive numbers. All twelve now fail with the grid-domain error. All six interior reference/system combinations reproduced every printed energy component, including after shifting the box center by 1 Å while keeping the pose inside. 1IEP uses the repository's prepared input, 7F5D uses the first published blind pose of the registered Armijo search. Both are historical controls; neither is a new pose test. The registration commit title described both targets as human in error. [1IEP deposition](https://www.rcsb.org/structure/1IEP) · [7F5D deposition](https://www.rcsb.org/structure/7F5D) · [Machine-readable results](results/summary.csv)

The public C++ overloads now validate setup and geometry before calculating scores. Both call the same numerical implementation as the private global-search finalization path, which retains its own selection rules and wall treatment. A batch prints the input path followed by the ten energy lines for every successful ligand. Six two-pose forward batches now return twelve identifiable scores; the previous version emitted zero. In every reverse-order batch the two blocks reproduced the single-input results in reversed order. For an input outside the domain between two valid inputs, the CLI reported the failed identity, scored both valid inputs and returned a nonzero process status. The extra CLI regression controls cover mixed parse failures and explicitly specified NaN or infinity. The direct C++ regression checks missing ligand/maps, exact recovery of the native-reference vector, all three nonfinite references and both exterior overloads. [CLI and C++ checks](run.py) · [API test](check_api.cpp) · [API receipt](results/api.txt)

No structures were sampled or optimized. No temperature, standard state or chemical potential was specified; the three reference choices are Vina score conventions. A single successful rigid explicit score evaluates the receptor–ligand term once, without a gradient. The automatic reference computes the ligand internal term before scoring as well. Box validation involves coordinate comparisons. No pair-table lookup counter was collected. These measurements supply neither a pose error nor a Ki/Kd accuracy estimate.

## Cost and replay

All 92 complete CLI processes used one CPU. For the same 36 registered process configurations the previous version used 1.061 CPU s / 2.144 elapsed s and the candidate 1.046 / 1.985 s. Across all 92, including additional candidate failure and order checks, totals were 2.816 / 5.500 s; the largest per-process RSS was 16,384 KiB. Timings include process startup, PDBQT reads, receptor indexing, pair-table initialization, validation, scoring and printing. Both programs used the same explicit receptor pathway and deferred grid population. Input unpacking and coordinate translation, shared across arms, took 0.098 CPU s / 0.223 elapsed s on replay. No search, gradient, rescoring of alternative modes or further conformer preparation was performed. The second complete run reproduced all 92 merged stdout/stderr logs byte for byte, every exit status and the summary. Individual commands and process costs are in runs.json inside the archive. [Inputs, coordinates and executable hashes](results/manifest.json) · [Original logs](results/raw.tar.xz) · [Replay receipt](results/replay.json)

From the repository root, build the two pinned versions with `tools/build.py`. Pass their paths to the runner and specify a fresh output directory. All structural inputs are in this repository. The C++ test links against the candidate library objects without main.o; use the same compiler and Boost configuration as the executable build.

```sh
python benchmarks/fixed_score_domain/run.py /path/to/f46aa37/vina /path/to/4dd4e37/vina /tmp/fixed-domain-fresh
diff benchmarks/fixed_score_domain/results/summary.csv /tmp/fixed-domain-fresh/summary.csv
```
