# Pair derivatives on human adenosine receptor A2A and renin

One changed search found a 1.790 Å first pose for ZM241385 in the receptor crystallized with T4G. The second seed missed this pose. No candidate recovered a first pose within 2 Å on the twelve held-out, time-calibrated seed cases. These data support keeping the parent algorithm as the production default.

## Structures and controls

The [protocol](PROTOCOL.md) fixes common inputs, the two trajectory caps and the seed used for timing calibration. This panel contains six distinct chemical graphs, four successfully prepared receivers and twelve query/receiver combinations. Four are cognate and eight are cross-docks. Each has two seeds and two caps. The 48 cells per arm are repeated measurements of two protein targets. See [data sources, scaffolds and chemical states](results/sources.json) and [reference-pose geometry](results/reference_geometry.csv).

| Target | Deposited complex, chain, resolution | Role |
|---|---|---|
| Human ADORA2A | [3EML](https://www.rcsb.org/structure/3EML), A, 2.60 Å and [3UZA](https://www.rcsb.org/structure/3UZA), A, 3.27 Å | Two queries and two receivers, including reciprocal cross-docks |
| Human REN | [1RNE](https://www.rcsb.org/structure/1RNE), A, 2.40 Å and [2V0Z](https://www.rcsb.org/structure/2V0Z), C, 2.20 Å | Two queries. Meeko found invalid oxygen valence in rebuilt glutamate residues, so both receivers were excluded before REN docking |
| Human REN | [3G6Z](https://www.rcsb.org/structure/3G6Z), A, 2.00 Å and [3SFC](https://www.rcsb.org/structure/3SFC), B, 2.10 Å | Structural fallbacks, used as two receivers and two additional queries after the first two preparations failed |

Each query starts from one ETKDGv3/MMFF conformer generated from the CCD graph. All six CCD graphs have formal charge zero. These are controlled chemical inputs. They do not establish the population of each protonation state under binding conditions. The crystallographic pocket center is supplied to every arm. All four retained receptor preparations have zero rebuilt nonterminal heavy atoms recorded by PDBFixer. Reproduction from locally cached raw structures took 24.83 s for the six preparations combined, including the two failed receivers. Every ligand and eligible receptor PDBQT reproduced byte-for-byte. The docking times below include per-run grid and table construction.

## End-to-end results

| Arm | Top-1 ≤2 Å | Top-5 ≤2 Å | Returned pool ≤2 Å | Top-1 given pool hit | Returned modes | Joint energy/gradient calls | Local optimizations | Total elapsed, s |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Upstream fork-base Vina | 0/48 | 0/48 | 0/48 | — | 230 | 375,512 | 7,602 | 71.96 |
| Parent Pro | 0/48 | 0/48 | 0/48 | — | 233 | 375,560 | 7,600 | 10.54 |
| Exact piecewise-linear derivative | 1/48 | 1/48 | 1/48 | 1/1 | 223 | 377,845 | 7,110 | 10.78 |
| Monotone cubic | 0/48 | 1/48 | 1/48 | 0/1 | 185 | 375,532 | 7,052 | 11.67 |
| Exact derivative plus shifted receptor cutoff | 1/48 | 1/48 | 1/48 | 1/1 | 225 | 373,142 | 7,106 | 11.28 |

Upstream means the original source at commit 3c65c0b compiled with the same toolchain. The independently downloaded [official v1.2.7 release](https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.7) reproduced every output PDBQT byte-for-byte across these 48 jobs. An upstream build instrumented only for counts also reproduced all 48. Timings in the table come from the uninstrumented upstream build. Pro arms report native counters. Each joint call returns an objective and gradient together. Counts include hunt, second optimization and explicit refinement. The logs give their individual counts as well as the number of explicit scoring inputs. Standalone native score decompositions and interpolation operations are included in elapsed time, outside the joint-call column.

At cap 3,000, all methods scored 0/24 for all three pose metrics. Both favorable cells below occurred at cap 12,000, in cross-docking, on seed 75467. Success denotes agreement with an aligned experimental pose; the query's actual structure in the alternate protein is unobserved. All 16 cognate cells per arm missed 2 Å. The median top-1 error was 8.91 Å for parent, 9.11 Å for linear, 8.83 Å for cubic and 8.84 Å for shifted. Full per-complex, per-target and per-cap results, including the larger regressions, are in [summary.csv](results/summary.csv).

* **ADORA2A, 3EML into 3UZA.** The parent first pose was 7.393 Å from the transferred ZMA pose, and no returned mode was closer. Linear recovered 1.790 Å (12,313 calls, 0.304 s), and shifted recovered 1.813 Å (12,201 calls, 0.277 s). Parent used 12,267 calls and 0.206 s. The transferred label has a 2.73 Å closest receptor contact, no severe clashes and an alignment RMSD of 0.568 Å over 58 pocket Cα atoms. The second seed did not reproduce the hit.
* **REN, 3SFC into 3G6Z.** Cubic returned a 1.233 Å mode behind a 3.603 Å first pose. Parent's best returned mode was 8.079 Å. The cubic and parent jobs took 0.253 and 0.225 s. The transferred label has no severe contact, but pocket Cα RMSD is 1.792 Å and a crystal mate approaches the ligand in the source structure to 4.877 Å. The reverse, 3G6Z into 3SFC, is an especially poor structural label with a minimum heavy-atom distance of 0.694 Å and 37 severe interatomic contacts.
* **Regressions.** On 2V0Z into 3G6Z, seed 98251 at cap 12,000, parent top-1 was 9.284 Å and linear/shifted were 14.423/14.439 Å. On 2V0Z into 3SFC, seed 75467 at cap 3,000, parent was 5.319 Å and cubic was 11.046 Å. Both are failures at the 2 Å threshold. Binary success alone would hide these degradations.
* **Geometry and cost.** The number of first poses with any protein–ligand contact below 0.75 of summed element van der Waals radii was 1/48 for upstream, 2/48 for parent, 2/48 for linear, 3/48 for cubic and 4/48 for shifted. The 48 shifted runs had exactly the same hunt and second-stage counts as linear in every cell. Only the explicit refinement count changed from 15,974 to 11,271 in aggregate. Faster refinement did not secure faster total elapsed time. All favorable primary hits lay outside the subset where both elapsed time and joint calls did not exceed parent.

## Held-out time control

The cap for each query/receiver and treatment was calibrated using only seed 75467. Each calibrated job on seed 98251 was run once or reused if the selected cap was already measured. Cost predictions had substantial error. Three further randomly ordered timing-only replays of all 48 held-out jobs verified byte-identical coordinates for every repetition. No replay contributes a new accuracy observation.

| Arm | Top-1 / top-5 / returned pool hits of 12 | First elapsed / parent | Sum of median replay elapsed / parent | Sum of median replay CPU / parent | Joint calls / parent |
|---|---|---:|---:|---:|---:|
| Parent, cap 12,000 | 0 / 0 / 0 | 1.000 | 1.000 | 1.000 | 1.000 |
| Linear, calibrated | 0 / 0 / 0 | 0.957 | 0.912 | 0.907 | 0.697 |
| Cubic, calibrated | 0 / 0 / 0 | 1.230 | 0.977 | 0.977 | 0.584 |
| Shifted, calibrated | 0 / 0 / 0 | 0.940 | 0.863 | 0.860 | 0.579 |

Only 1, 3 and 3 of 12 cells respectively matched their parent to within ±5% of median replay wall time. All per-pair discrepancies remain in [matched_pairs.csv](results/matched_pairs.csv). These short timings support a cost range and a negative replication, with limited precision for direct speed comparisons. The small panel and low baseline recovery give little power for a general accuracy conclusion. Missing water, the frozen neutral ligand graph, shuffled branches after gradient changes, finite sampling and unobserved cross-docked poses are plausible contributors. Their individual effects are unmeasured here.

## Reproduction

The [input archive](inputs.tar.xz) contains raw PDB and CCD files, prepared queries and receivers, preparation metadata and failed receptor logs. The [run archive](runs.tar.xz) contains raw output coordinates, stdout/stderr, timings and upstream counter replays. The [checksums](SHA256SUMS) cover these files and all result tables. The independent audit reconstructs reference geometry from raw PDB/CCD, checks graph symmetries and chirality, and uses a second sequence alignment and rigid transform implementation. Across 1,427 distinct modes the largest RMSD discrepancy was 2.31 × 10⁻¹⁴ Å, no input stereocenter changed tag, and the largest covalent bond-length change from its prepared value was 0.00136 Å.

The archived analysis used Python 3.12, RDKit 2025.9.3, Meeko 0.7.1, Gemmi 0.7.4, OpenMM 8.6.1, PDBFixer 1.12.0, Biopython 1.88, NumPy, SciPy and pandas. A valid single-core affinity index depends on the host. For the measured primary comparisons all arms ran sequentially on logical CPU 8 of an AMD EPYC 9V74, with one docking thread. The official-release portability check ran separately on logical CPU 5.

From the repository root, set `B` to the directory of the study binaries built by [build_arms.py](../build_arms.py). The following recovers the archived audit without redocking.

```bash
P=benchmarks/pair_gradient/pose
W=/tmp/fastdock-pose-audit
mkdir -p "$W"
(cd "$P" && sha256sum -c SHA256SUMS)
tar -xJf "$P/inputs.tar.xz" -C "$W"
tar -xJf "$P/runs.tar.xz" -C "$W"
cp "$P"/results/* "$W/"
python "$P/audit.py" "$W"
python "$P/sources.py" "$W"
python "$P/summarize.py" "$W"
```

For a fresh execution, use an empty work directory, extract only the inputs and run the commands below on one pinned core. Supply the official release binary to check serialization across distributions. Preparation from the original sources can be replayed with `sources.py --reprepare` once inputs are extracted.

```bash
export OPENMM_CPU_THREADS=1 OPENBLAS_NUM_THREADS=1
python "$P/paired.py" "$W" --upstream "$B/upstream/build/vina" --parent "$B/parent/build/vina" \
  --linear "$B/linear/build/vina" --hermite "$B/hermite/build/vina" --counts-replay "$B/upstream-counts/build/vina"
python "$P/paired.py" "$W" --phase shift --shift "$B/shifted/build/vina"
python "$P/paired.py" "$W" --phase matched --linear "$B/linear/build/vina" \
  --hermite "$B/hermite/build/vina" --shift "$B/shifted/build/vina"
python "$P/paired.py" "$W" --phase release --release /path/to/official-vina-1.2.7
python "$P/time_replay.py" "$W" --parent "$B/parent/build/vina" --linear "$B/linear/build/vina" \
  --hermite "$B/hermite/build/vina" --shift "$B/shifted/build/vina"
python "$P/audit.py" "$W"
python "$P/summarize.py" "$W"
```
