# Galectin Kd ranking after elapsed-time calibration

At almost identical total predictive runtime per target, the exact radial-table derivative improves the 13-compound Gal-3 consensus from **25/35 to 31/35** direction credits. It leaves Gal-1 at **9/23**. One Gal-3 ligand accounts for the sign of the aggregate change: omit compound 11 and the result becomes **24/26 for parent, 22/26 for linear**. A descending heavy-atom-count rule alone scores 34/35 for Gal-3. The result supports a timing-controlled, series-specific observation and leaves the evidence for general improvement insufficient.

The same 13 ligands, chemical states, Kd table and four receptor fixtures were used in the original nominal-12k experiment. That study averaged seeds 75467 and 98251; its Gal-3 result was 23/35 and 30/35. This follow-up keeps seed 75467 for timing calibration and scores seed 98251 alone. On that seed at the original fixed cap the matching controls are 25/35 and 30/35. Both seeds and both chemical series were already used in the earlier study. [Frozen protocol](PROTOCOL.md) · [Source study and assay table](../../galectin_affinity/README.md) · [Original four-structure comparison](../galectins/summary.csv)

## Accuracy by receptor

Each arrow goes from parent at 12k to linear at its calibrated cap. Each prediction uses the first emitted pose. Receptor consensus is the unweighted mean of the two scores per compound. All pair counts use a threefold or larger Kd ratio. Half-credit handles score ties. Spearman and Kendall use all 13 compounds irrespective of pair gap.

| Target | Receptor | Pair credits | Spearman | Kendall tau-b | Equal-heavy pair credits |
|---|---|---|---|---|---|
| Gal-1 | 1GZW | 12 → 8.5 / 23 | −.143 → −.061 | −.090 → −.052 | 2 → 1 / 7 |
| Gal-1 | 4Q26 | 14 → 14 / 23 | .143 → −.030 | .090 → −.039 | 5 → 2 / 7 |
| Gal-1 | Consensus | 9 → 9 / 23 | −.209 → −.162 | −.142 → −.168 | 1 → 1 / 7 |
| Gal-3 | 1KJL | 25 → 28 / 35 | .256 → .459 | .194 → .323 | 0 → 2 / 2 |
| Gal-3 | 3T1L | 18 → 25 / 35 | .165 → .404 | .142 → .374 | 2 → 2 / 2 |
| Gal-3 | Consensus | 25 → 31 / 35 | .377 → .680 | .323 → .477 | 2 → 2 / 2 |

The Gal-3 triazole consensus changes 15→21/25, while the corresponding Gal-1 subset changes 3→1/13. Eight previously missed Gal-3 comparisons involving compound 11 are recovered, with a net loss of two on the remaining edges. Across all 13 single-ligand deletions the Gal-3 change ranges from −7.7 to +25 percentage points. There are only two threefold-separated Gal-3 pairs with equal heavy-atom count. The pre-existing heavy-atom rule gives 16.5/23 for Gal-1 as well. Pair edges share compounds and the two receptor runs share experimental measurements. [All pair decisions](pairs.csv) · [Subsets and correlations](summary.csv) · [Compound and size controls](influence.csv)

## What the clock and evaluation counters show

The complete validation ran serially with CPU affinity `{7}`. Caps preserve one CPU, one trajectory, an identical search box and nine output slots. Every elapsed measurement includes starting the docking subprocess, grid work, search, explicit refinement and writing the result. The per-compound/receiver caps were frozen from seed-75467 timings only. Of 52 cells, 28 kept 12,000 and 24 were reduced to 7,950–11,993. The generated cap-ledger SHA-256 was `0d8618b6aa851d2973d7105fe8be08d49e08360e28949e312f4a31e884a88eb3`.

| Target | Parent seconds | Linear seconds | Parent joint evaluations | Linear joint evaluations | Within ±5% time per cell | Range of linear / parent time |
|---|---:|---:|---:|---:|---:|---|
| Gal-1 | 3.415 | 3.352 | 320,562 | 302,400 | 17/26 | 0.796–1.083 |
| Gal-3 | 3.564 | 3.529 | 321,565 | 305,798 | 20/26 | 0.888–1.103 |
| All four receivers | 6.978 | 6.881 | 642,127 | 608,198 | 37/52 | 0.796–1.103 |

Both aggregate target-level elapsed ratios meet the prespecified 5% criterion. Individual runs show wider variation and 24 of 52 linear cells were slower than parent. The 208 calibration runs consumed another **25.539 s** of measured subprocess time, plus 13.860 s for the 104 validation runs. The entire comparative procedure took 39.51 s wall clock. Repeating this calibration at deployment would be an additional cost. These timings do not support a claim of intrinsic throughput improvement. [Per-cell costs](paired_costs.csv) · [Receiver and aggregate costs](costs.csv) · [Calibration measurements](calibration.csv) · [Frozen caps](budgets.csv)

## Reproduction

Use the parent/linear executables from the pair-gradient build instructions. Rebuild timing data on an available isolated CPU. Wall times and the derived calibrated caps depend on the host. For exact offline replay of the published numbers, extract the committed raw logs and poses and use the committed cap/timing records instead. The audit reads integer-nanosecond timings, the source Kd table in µM, and first scores/counters from raw files. All 184 runs made at 12k in this trial reproduce their historical score and all six search/explicit counters.

```sh
D=benchmarks/pair_gradient/galectin_time
W=$(mktemp -d)
tar -xJf benchmarks/galectin_affinity/evidence.tar.xz -C "$W" prepared
(cd "$W/prepared" && sha256sum -c "$OLDPWD/$D/inputs.sha256")
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python "$D/run.py" "$W/prepared" "$W/fresh" \
  --cpu 7 --parent "$PARENT" --linear "$LINEAR"
python "$D/evaluate.py" "$W/fresh"

mkdir "$W/offline"
tar -xJf "$D/raw.tar.xz" -C "$W/offline"
cp "$D"/{calibration.csv,validation.csv,budgets.csv,run.json} "$W/offline/"
python "$D/evaluate.py" "$W/offline"
for f in summary.csv pairs.csv influence.csv paired_costs.csv costs.csv; do
  cmp "$D/$f" "$W/offline/$f"
done
```
