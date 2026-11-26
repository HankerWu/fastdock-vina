# Fixed-budget trajectory breadth and depth

The study used three new seeds 76231, 76243 and 76249 on the ten prepared
cognates and four reciprocal directions. Both arms used the analytic metric
in Late refinement, nine output modes, one CPU thread and 144,000 nominal
global evaluations. Twelve trajectories capped at 12,000 evaluations represent
breadth. Three capped at 48,000 represent depth. An MC iteration can complete
past the cap. The final explicit refinement and wall time were measured
separately from global evaluations.

| Set | Breadth top1 / pool | Depth top1 / pool | Actual global evaluations breadth / depth |
|---|---:|---:|---:|
| 10 cognates × 3 seeds | 3 / 8 of 30 | 4 / 7 of 30 | 4,333,903 / 4,322,997 |
| 4 reciprocal directions × 3 seeds | 0 / 0 of 12 | 0 / 0 of 12 | 1,734,537 / 1,729,468 |

Neither allocation dominates these endpoints. Output pools are the nine
emitted candidates. `python allocation/run.py PREPARED VINA OUT --jobs 2`
repeats the experiment from the parent directory. The input bundle and case
eligibility are in ../pose. Read the CPU-calibrated experiment in ../upstream
for a separate comparison against the official executable.
