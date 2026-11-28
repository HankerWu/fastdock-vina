# Where did the near-reference HUP positions survive?

The two ordinary searches at seed 1103503 use the same ligand and budget. Their published first poses differ by 5.329 Å. A passive single-CPU replay with no crystal input reproduced the full production PDBQT and each search-counter line for both receivers. The stage records locate a difference before either final ranking or explicit-receptor refinement. [Observation plan](PROTOCOL.md) · [Native commands, release outputs and costs](../blind/results/raw) · [Observer patch](observer.patch)

| Stage | 4EY6 records ≤2 Å / all | Best RMSD (Å) | 4EY7 records ≤2 Å / all | Best RMSD (Å) |
|---|---:|---:|---:|---:|
| Randomized trajectory start | 0/3 | 10.447 | 0/3 | 10.409 |
| Actual proposal after mutation | 731/2,431 | 0.207 | 1/2,828 | 1.990 |
| Hunt minimization returns | 768/2,431 | 0.289 | 0/2,828 | 2.121 |
| Accepted hunt returns | 620/1,955 | 0.349 | 0/2,160 | 2.131 |
| Second minimization returns | 229/780 | 0.347 | 0/204 | 2.960 |
| Entered trajectory container | 16/120 | 0.369 | 0/84 | 2.960 |
| Subsequently evicted from that container | 14/94 | 0.369 | 0/57 | 2.960 |
| Still present at trajectory end | 2/26 | 0.414 | 0/27 | 2.981 |
| After trajectory merge | 1/9 | 0.414 | 0/9 | 2.981 |
| Explicit refinement inputs | 1/9 | 0.414 | 0/9 | 2.981 |
| Explicit refinement returns | 1/9 | 0.414 | 0/9 | 2.981 |
| Published, after coordinate rounding | 1/9 | 0.414 | 0/9 | 2.981 |

Counts refer to visits, including strongly correlated steps within a trajectory. The best across a stage can change identity; an individual trajectory's transformations can be read from the per-step and lineage tables. These data do not count BFGS line-search trials inside a minimization, and cannot rule out transient visits there. [Each attempt](results/steps.csv.gz) · [Per-trajectory distributions of near visits and minima](results/stages.csv) · [Every container candidate's lifetime](results/trajectory_lifetimes.csv)

On 4EY7 the sole sub-2 Å proposal is trajectory 1 step 163. It starts at 1.989812 Å and leaves hunt at 2.287681 Å. Monte Carlo accepts the return, but the native save gate does not run the second minimization on it. The best hunt return over all steps is 2.121356 Å (trajectory 1, step 183), and fails acceptance. The best published mode is 2.981410 Å versus 5.743813 Å for the published first mode. Here the search generated one borderline geometry, the first minimization moved it farther from the crystal, and final ranking adds a separate 2.762402 Å regret. The one-position observation is conditional on this empirical search, preparation and seed.

On 4EY6 the identity published first comes from trajectory 3 step 447. Mutation puts the ligand at 1.569063 Å and hunt returns 0.414348 Å. Second minimization and explicit refinement preserve this coordinate to the recorded precision, and final output rounds it to 0.414354 Å. Trajectory 3 step 421 had already entered the container at 0.369112 Å, only to be evicted at step 425. Across whole runs 45 initially >2 Å proposals crossed into ≤2 Å during hunt, and eight moved the other way. Total near counts across 4EY6 are consequentially different from improvement counts. [Origins of merged and published poses](results/origins.csv) · [Published first mode back to its proposal](results/winner_path.csv)

The binary stream reuses the repository's basin observer to write all 18 ordered heavy atom serials and full binary64 coordinates. Membership is observed after each native insertion attempt, using the exact vector of stored heavy coordinates plus energy to maintain identities. The final trajectory, merged and pre-explicit stages are joined using exact heavy coordinates and objective value; each join was unique. A published pose must have a unique explicit-refined precursor within 0.000501 Å per component. The initial and proposal energy payloads are placeholders; the observer does not score them. Released geometry is measured in the frozen protein frame with the one verified stereochemical mapping. The observer does not record any Boltzmann population.

4EY6 uses 60,973 hunt, 11,101 second-minimization and 135 explicit-refinement joint energy/gradient evaluations. 4EY7 uses 69,377, 2,688 and 151. The first two traced native invocations consumed 0.357 and 0.369 CPU seconds, respectively. Compression consumed an additional 0.697 and 0.598 CPU seconds. Offline linkage and RMSD measurement used 0.402 CPU seconds. A second execution in a clean directory reproduced the compressed coordinate/membership streams and evaluation tables byte for byte. All streams use lossless XZ; the .00 suffix is a segment ordinal, and the manifest records the SHA-256 of the reassembled uncompressed stream. [Executables, product hashes and replay costs](results/runs.json) · [4EY6 log](results/4EY6.log) · [4EY7 log](results/4EY7.log)

From the repository root, set BUILD to the native production objects directory's parent, BOOST to the Boost prefix and WORK to an empty temporary directory. The run step only opens the two noncognate receptors, graph-only ligand, prior CLI arguments, counters and released PDBQTs for parity. The separate evaluation step consumes the archived diagnostic registration and crystal mapping.

```bash
D=benchmarks/ache_receptor_state/trace
bash "$D/build.sh" "$BUILD" "$BOOST" "$WORK/tracer"
python "$D/replay.py" "$WORK/tracer/vina" "$WORK/run" "$WORK/results"
OPENBLAS_NUM_THREADS=1 python "$D/evaluate.py" "$WORK/results" \
    benchmarks/ache_receptor_state/diagnostic/inputs.json "$WORK/results"
for f in "$D"/results/{*.xz.*,*.csv,*.csv.gz}; do cmp "$f" "$WORK/results/${f##*/}"; done
```
