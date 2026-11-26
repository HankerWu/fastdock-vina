# Elapsed-time control for thrombin cross-docking

A second receptor grid adds preprocessing time to the soft-repulsion arm. The previous study held the evaluation cap constant. This additional test uses measured process elapsed time to choose a higher cap for the native arm and compares structural recovery on the same five crystallographic ligands and frozen confirmation seeds. The additional calibration seed is 101921. All commands use the same executable, one worker thread and inherited affinity to Linux logical CPU 8. The three initial warm-ups were discarded. Timed runs were serial. Reported elapsed time includes process launch, receptor/ligand I/O, map creation, Monte Carlo, explicit refinement and writing the outputs. Child user plus system CPU time provides a parallel measurement. These timings apply to this executable, grid and host.

Calibration used a fully crossed three-receptor × three-ligand panel (compounds 5, 1c, 4 with 26, 27 and 29 heavy atoms) in two order-balanced blocks, with three policies for each of the resulting 18 observations. This deliberately spans the size range in the structural test panel. Calibration did not read the crystal labels or any pose results. The second block reverses receiver/ligand order; within triplets the run order rotates and is reversed between blocks. All comparable runs in calibration reuse the same input seed. The hard-75k reference gives the observed marginal cost of additional MC evaluations when each trajectory's ordinary 50k budget changes.

Let H, S and P be the total measured time across the calibration panel with hard 50k, soft 50k and hard 75k, respectively. Choose a common hard cap by rounding 50000 + 25000(S−H)/(P−H) to the nearest thousand, floored at 50000. The wall-time measurements were 20.111, 24.347 and 28.914 seconds, yielding 62000. Child CPU totals were 19.926, 24.132 and 28.629 seconds and independently yielded 62000. Relative to the hard arm, the 50k soft arm took 21.1% more total elapsed time and 22.4% more time at the median paired observation. The cost may differ for other box sizes, atom-type sets, chemistry, compiler and hardware.

Validation uses the structurally labeled compounds 5, 1a, 1b, 3a and 4 with three receiving structures and seeds 88091 and 302123. All 30 cases ran three arms (hard 50k, soft 50k and calibrated hard 62k). Run order was rotated and reversed across consecutive test cases; each case's arms used the same input seed and receptor. Every model and search counter is archived. The independent crystal RMSD evaluator reads PDBQT coordinates, reconstructs atom correspondence from the CCD and rotates the crystal upon a local protein alignment.

| Validation outcome | Native 50k | Native 62k | Soft 50k |
|---|---:|---:|---:|
| All 30 cases, end-to-end seconds | 34.91 | 41.76 | 40.86 |
| 26 distinct-deposit cases, end-to-end seconds | 30.13 | 36.48 | 35.65 |
| 26 distinct-deposit cases, child CPU seconds | 29.87 | 36.08 | 35.33 |
| 26 distinct-deposit cases, first pose ≤2 Å | 9 | 11 | 15 |
| 26 distinct-deposit cases, best pose ≤2 Å | 15 | 16 | 20 |
| All 30 cases, first pose ≤2 Å | 11 | 13 | 17 |
| All 30 cases, best pose ≤2 Å | 17 | 18 | 23 |

On the 26 cross-docking cases, the calibrated hard baseline consumed 2.34% more aggregate elapsed time than the soft arm. The median per-case hard-62k/soft-50k duration ratio was 1.020 (interquartile range 0.987–1.070). Against the higher-budget hard arm, soft search gains eight hits and loses four. Against the ordinary hard arm, the corresponding count is eight gains and two losses. Of the six original net additional hits, two are recovered by the additional hard-search budget. This panel has five distinct molecules on a single target and scaffold. It does not measure the time-matched performance of inter-ligand affinity ranking. The 30 baseline and 30 soft files at 50k were identical by SHA-256 to their corresponding previous confirmation outputs. Absolute elapsed times depend on the host; costs are measured here with a fresh CLI process and map preparation for every search.

## Reproduction

From the parent benchmark's work directory holding `prepared/`, `raw/` and `series.csv`, use the full path to this directory as `B` and the new Vina binary as `VINA`. The sibling `independent_audit.py` must be available in the bench tree above this directory as in the checkout. Time all cases afresh or inspect the archived timed outputs. The first command recalibrates the native budget on the current machine.

```bash
python "$B/equal_time.py" "$VINA" "$PWD" timed_again --core 8
python "$B/summarize_time.py" "$PWD" timed_again
```

Choose a logical CPU available on the host or omit `--core` to pin to the highest available CPU. For an offline audit with the frozen budget, unzip `timed_evidence.zip` in the work directory and run `python "$B/summarize_time.py" "$PWD" timed`. Supply `--previous-confirm confirmation` if the original soft-search confirmation archive has been extracted to `confirmation/`; this additionally verifies both 50k outputs in each of the 30 cases against the previous raw files by SHA-256. Numerical tables and individual paired timings are in `expected/`.
