# Unseen androgen receptor transfer and explicit-refinement stage audit

Five engines recovered the same geometric hits on this frozen wild-type human androgen receptor panel. Each obtained 1/8 Top-1, 1/8 Top-5 and 1/8 oracle at 3,000 and again at 12,000 nominal evaluations. Both hits are the **same reciprocal direction and same seed**. No cognate run passed 2 Å. This steroid series supports no accuracy gain over upstream or the branch parent.

## Inputs and scope

R18 (metribolone) was prepared from human AR 1E3G chain A, TES (testosterone) from 2AM9 chain A. Receptor sites use ligand centroids, and the independent RDKit embed uses each chemical graph without crystal coordinates. Two closely related rigid steroid cores occupy two similar receptor structures. Seventy-five pocket Cα atoms align at 0.299 Å in either cross direction. The original PDB and CCD files and checksums are included with both dock-ready fixtures. 1E3G had no omitted nonwater heavy atoms or repaired heavy atoms. In 2AM9 the closest omitted nonwater atom was 11.77 Å from the ligand and the closest repaired atom was 14.61 Å. Identical one-thread preparation from the saved raw files took 3.54 s and 2.92 s respectively, reproducing both prepared structures and metadata byte-for-byte. [Frozen protocol](PROTOCOL.md) · [1E3G deposition](https://www.rcsb.org/structure/1E3G) · [2AM9 deposition](https://www.rcsb.org/structure/2AM9)

We used the unmodified upstream 1.2.7 binary, branch parent, exact linear slope, monotone Hermite, and shifted explicit-receptor potential. Each cell runs with one CPU, one trajectory, at most nine modes, a 24 Å cube, energy range 100, and seeds 75467 and 98251. The broad energy window retains poor poses for diagnosing failures. Every arm used the same prepared input bytes. Upstream timing uses its uninstrumented executable. All 16 passive upstream-count replays and all 32 linear/shift stage replays returned byte-identical pose files against their respective production runs. Time spent verifying results is outside the predictive measurements.

| Search cap | Engine | Top-1 / Top-5 / oracle | Joint objective and gradient calls | Explicit calls / joint evaluations | Sum of wall seconds, eight runs |
|---|---|---|---:|---|---:|
| 3,000 | Upstream | 1 / 1 / 1 | 25,148 | 17 / 785 | 8.343 |
| 3,000 | Parent | 1 / 1 / 1 | 25,150 | 17 / 787 | 1.431 |
| 3,000 | Linear | 1 / 1 / 1 | 25,315 | 19 / 952 | 1.375 |
| 3,000 | Hermite | 1 / 1 / 1 | 25,578 | 21 / 1,215 | 1.687 |
| 3,000 | Shift | 1 / 1 / 1 | 25,305 | 21 / 942 | 1.444 |
| 12,000 | Upstream | 1 / 1 / 1 | 97,239 | 21 / 889 | 8.316 |
| 12,000 | Parent | 1 / 1 / 1 | 97,258 | 21 / 908 | 1.977 |
| 12,000 | Linear | 1 / 1 / 1 | 97,343 | 21 / 993 | 1.491 |
| 12,000 | Hermite | 1 / 1 / 1 | 97,324 | 21 / 974 | 1.470 |
| 12,000 | Shift | 1 / 1 / 1 | 97,493 | 25 / 1,143 | 1.444 |

The count sums are hunt + second optimization + explicit-refinement joint objective/gradient requests. Separate score-only operations, map construction and startup are included in wall time. The parent already contains substantial implementation speedups against upstream. Wall results include shared-host scheduling noise and show no consistent additional speed advantage for shift. [Per-run actual counts, times and hit criteria](cells.csv) · [Per-mode scores, receptor distances and clashes](modes.csv)

| Query → receiver | 3k hits, two seeds | 12k hits, two seeds | Diagnosis shared by all five engines |
|---|---:|---:|---|
| 1E3G → 1E3G | 0 | 0 | Top poses stay near 2.76–2.77 Å. |
| 1E3G → 2AM9 | 0 | 0 | One seed stays near 2.72–2.76 Å, the other gives a distant pose. |
| 2AM9 → 2AM9 | 0 | 0 | At 12k, the first poses lie near 6.29 Å. |
| 2AM9 → 1E3G | 1 | 1 | Only seed 98251 recovers the aligned cross reference. |

The aligned TES crystal reference has one contact at 1.23 Å against the alternate 1E3G receptor. Its ≤2 Å recovery measures agreement with the aligned proxy; the cognate crystal structure of TES is 2AM9. Both cognate crystal labels have zero contacts under the stated test, as does R18 aligned to 2AM9. Across 187 emitted poses, the permissive output window admits 116 with at least one receptor contact below 0.75 × the sum of the tabulated radii. Each engine has three clashing top poses out of eight at 3k and two at 12k. Independent stereo-sensitive RDKit correspondences, gemmi coordinates and SciPy rotations verified every mode's RMSD, with largest discrepancy 5.33e−15 Å. No predicted stereocenter inverted; maximum heavy bond-length change from its prepared ligand was 0.00130 Å. Optimal rigid fits of prepared R18 and TES to their cognate crystal ligands give 0.427 and 0.068 Å respectively. Their single prepared conformers therefore permit a near-native heavy-atom geometry. [Per-pose independent geometry measurements](geometry.csv) · [Reference packing and audit](verification.json)

## Where the linear and shifted runs diverged

The passive patch serializes full-precision Cartesian heavy-atom coordinates and the stored energy at three locations. It records the ordered grid-search pool, the actual argument to each explicit optimization, and the rescored poses after explicit optimization but before final sorting. Each record also includes a serial-identified PDBQT geometry for computing receptor-frame RMSD. Coordinate-only and joint energy/coordinate SHA-256 digests are separate. Every grid-search counter and every coordinate/energy digest agrees across the 38 common search candidates and all 38 inputs to explicit refinement. After refinement, all 38 Cartesian digests differ. Both arms preserve the same two hits and leave all 36 other inputs outside 2 Å. There are two emitted winners conditional on the two successful pools. The other 14 pools already lack a near-native candidate when the grid search returns. This identifies a limitation upstream of explicit refinement in these runs; the trace cannot distinguish an unvisited basin from a visited basin discarded inside the grid search. [Full stage ledger](stages.csv) · [Passive trace patch](stage_trace.patch)

## Reproduction

Use the same Python chemistry environment and five binary builds as the parent pair-gradient experiment. The stage patch applies on top of each of the two production source trees for linear and shift, respectively; rebuild with the same options used for its production version. The upstream count build applies the parent experiment's `upstream_counts.patch` to upstream 1.2.7. For example, after building the executables in separate directories:

```sh
D=benchmarks/pair_gradient/holdout
W=$(mktemp -d)
tar -xJf "$D/inputs.tar.xz" -C "$W"
OPENBLAS_NUM_THREADS=1 python "$D/run_holdout.py" "$W" \
  --upstream "$UPSTREAM" --parent "$PARENT" --linear "$LINEAR" --hermite "$HERMITE" --shift "$SHIFT" \
  --counts-replay "$UPSTREAM_COUNT" --trace-linear "$TRACED_LINEAR" --trace-shift "$TRACED_SHIFT"
python "$D/verify_holdout.py" "$W"
```

For offline pose/geometry and stage verification from the original raw outputs, extract `runs.tar.xz` beside the fixture, copy `modes.csv`, `cells.csv` and `stages.csv` there, then run `verify_holdout.py`. To reconstruct dock-ready inputs from deposited structures and CCD graphs, run `OPENMM_CPU_THREADS=1 OPENBLAS_NUM_THREADS=1 python "$D/prepare_holdout.py" "$W"`. Each archive contains the original inputs or production outputs and diagnostic sidecar logs, including the round-trippable full-precision stage coordinates.
