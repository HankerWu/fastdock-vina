# Native and official controls at fixed evaluations and matched CPU cost

Five human targets with complete reciprocal structures from [the primary panel](../human_panel) were selected before reading these test seeds. Each contributes both cross-ligand directions. Seeds 68041 and 95087 give 20 paired cells. The receptors and ligands overlap the primary panel; the two seeds are new. Crystal data define a known pocket and the evaluation frame.

| Search | Evaluation cap per trajectory | CPU seconds, all 20 | Top 1 / top 5 / all within 2 Å | Median top-1 RMSD, Å |
|---|---:|---:|---:|---:|
| Official Vina, 0.375 Å | 24,000 | 19.18 | 10 / 12 / 12 | 1.675 |
| FastDock Late, trilinear, 0.375 Å | 24,000 | 9.32 | 10 / 13 / 13 | 1.679 |
| FastDock Late, tricubic, 0.375 Å | 24,000 | 12.86 | 9 / 11 / 11 | 2.371 |
| FastDock Late, trilinear, 0.25 Å | 24,000 | 15.81 | 8 / 14 / 14 | 3.397 |
| FastDock Late, trilinear, 0.375 Å | 41,000 | 12.78 | 10 / 15 / 15 | 1.680 |

All runs had one CPU and three trajectories. An individual minimization can run past its evaluation cap. Official Vina used its own final selection and refinement. The four FastDock arms requested 15 Late refinement candidates and nine outputs. Both source trees were built with the same compiler, flags and Boost libraries. The upstream comparison applies to these commands and inputs. Official source revision bb64e5e26b520569a86360d69248069b5e665b3f. FastDock source revision [23b24f6](https://github.com/HankerWu/fastdock-vina/commit/23b24f6c53892e3edef8fbfe726cd2b34a080a89).

An independent calibration seed on one direction per target measured aggregate CPU time for native at 24k and 48k and cubic at 24k. Linear interpolation of the measured native cost, rounded to the nearest thousand, gave 41k. The 20 confirmation cells were then run with the frozen cap. CPU includes process startup, grid construction, search, explicit refinement and writing. Pooled elapsed time for 41k native and 24k cubic was 12.901 and 12.903 seconds. The machine, low caps and known-site boxes limit extrapolation. [Calibration](calibration.json) · [Confirmation](confirmation.csv)

## Crystal-informed diagnostic

Each experimental pose was transformed using local protein Cα correspondences. The projected start retains the prepared molecule's fixed geometry, transplants movable dihedrals and fits the ligand to the reference before native local relaxation. The raw energy sums native intermolecular and intramolecular contributions within one molecule and receiver. Each relaxation is one trial using the known answer.

For EPAS1 5TBM into 4GHI at seed 68041 the cubic first pose was 14.00 Å from the reference and all nine modes missed 2 Å. A 0.52 Å projected local pose was favored by 1.90 raw energy units. For FABP4 2NNQ into 3FR2 the 1.08 Å projected geometry was about 7.8 units above the selected poses. For RBP4 3FMZ into 1RBP the 1.42 Å projected geometry still had a 0.57 Å heavy-atom clash with Gln98 NE2. Its energy sat over 66 units above the selected poses. [Projected poses](crystal_starts.csv) · [Gaps and atom identities](cross_diagnosis.csv)

For the wall correction, 37 of 40 native/cubic confirmation pools were byte-identical. All 40 top-1 and best-pool success classifications persisted. A separate replay covers 132 primary pose and 192 FKBP12 runs. [Correction replay](wall_replay.csv) · [Primary replay](../human_panel/wall_replay_summary.json)

## Reproduction

Extract the human panel prepared archive to fixtures. Use separate builds of official Vina, the interpolation commit and the wall correction. For the original timing control use the archived 41,000 cap; a fresh calibration on different hardware may select a different value. The scripts require NumPy and SciPy. The projected replay additionally uses the human panel chemistry dependencies. Offline inspection requires expanding confirmation-pools.tar.xz into the same directory as confirmation.csv.

```bash
B=benchmarks/cubic_grid/confirmation
python "$B/run.py" calibrate fixtures confirmation /path/to/interpolation/vina
python "$B/run.py" confirm fixtures confirmation /path/to/interpolation/vina --official /path/to/official/vina
python "$B/check_wall_outputs.py" fixtures confirmation wall /path/to/corrected/vina
python "$B/projected.py" fixtures confirmation diagnostics /path/to/interpolation/vina
```
