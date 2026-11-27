# Gal3 rank transfer at lower actual cost

At a fresh seed and the full 12,000-evaluation cap, one-slot retention gained 3/403 pair credits and a recurrence-based choice gained 4/403 against the parent. At the smaller receiver-specific caps, both gains disappeared. All three lower-cap methods scored 183/403, compared with 201/403 for the full-cap parent. A same-cap parent ablation attributes this 18-credit shortfall to shortening search on this seed. The molecular evidence for the full-cap gain concentrates on compound 59 in one receptor.

## Design and computation

The earlier seed 614039 gave 221→225/403 for energy-gated retention, with one ligand responsible for the four flips and a 3.97% higher measured wall time. Before obtaining any fresh Gal3 outputs, we froze one validation seed (411119), the exact same 47 ligands and both receptors, a 12,000-evaluation parent, and two candidate caps. From the earlier receiver-specific CLI times, take 12,000 × parent time / candidate time, clip at 12,000 and round down to the next 500 evaluations. The raw estimates of 11,372.49 and 11,700.99 gave 11,000 for 1KJL and 11,500 for 3T1L. Neither receiver saturated the 12,000 ceiling. The linear approximation includes a conservative rounding margin and treats fixed per-process work only approximately. It imports no affinity labels. [Frozen controller and primary design](PROTOCOL.md) · [Previous series and inputs](../affinity/README.md)

We ran the full-cap parent, short-cap parent, short-cap energy gate and short-cap recurrence gate in alternating order within randomized receiver–compound cells. Each CLI used CPU 7, one trajectory, seed 411119 and max nine modes with the common input box, ligand and receptor. This block contains all direct timing controls. After observing its results, a separately documented follow-up ran both policies at 12,000 with the identical seed and 94 inputs. This completes the comparison of three policies at two caps. The same full-cap parent result is used in both comparisons. Wall times for the two added arms are descriptive due to their later timing block. [Full-cap extension recorded before extra runs](FULL_CAP.md)

## Same-assay ranking

The labels are the frozen Sörme et al. human Gal3 competitive-FP Kd estimates at the common 8 µM inhibitor concentration. The primary metric gives half credit to score ties over the same 403 unique, ≥3-fold-separated pairs. Predictions average the two first native Vina scores per ligand. Analyses use exact decimal Kd ratios; pKd = 6 − log10(Kd in µM). The table uses per-compound Spearman ρ and Kendall τ-b between pKd and the negative native score. The 35 identically sized and charged pairs overlap the 403. The 63 tenfold-separated pairs provide a second affinity-gap check. [Article, Fig. 4 and supplement p. S6](https://doi.org/10.1021/ja043475p) · [Per-compound Kd, pKd and native scores](predictions.csv) · [All metrics](metrics.csv)

| Arm | 1KJL /403 | 3T1L /403 | Mean /403 | Matched /35 | ≥10-fold /63 | Mean ρ | Mean τ-b |
|---|---:|---:|---:|---:|---:|---:|---:|
| Parent, 12,000 | 192 | 204.5 | 201 | 17 | 29 | −0.0091 | 0.0019 |
| Energy gate, 12,000 (follow-up) | 196 | 204.5 | 204 | 17 | 29 | 0.0061 | 0.0148 |
| Recurrence gate, 12,000 (follow-up) | 196 | 204.5 | 205 | 17 | 29 | 0.0094 | 0.0167 |
| Parent, lower caps | 187 | 188.5 | 183 | 16 | 25 | −0.0934 | −0.0557 |
| Energy gate, lower caps | 187 | 188.5 | 183 | 16 | 25 | −0.0928 | −0.0538 |
| Recurrence gate, lower caps | 187 | 188.5 | 183 | 16 | 25 | −0.0928 | −0.0538 |

Under the full cap, Kd 3.4 µM compound 59 changes from −4.150 to −4.903 kcal/mol on 1KJL in both new arms. This corrects four consensus orderings against 8, 16, 27 and 28. The energy gate additionally changes compound 8 on 3T1L (−4.357 to −4.383), losing its comparison with 21. Removing any one ligand gives a range of −1 to +4 consensus credits for the energy gate and 0 to +4 for recurrence. Removing 59 leaves −1 and 0 respectively. The chemically correlated comparisons do not give four independent successes. The earlier seed's affected ligand was 32. The tenfold and equal-heavy/charge subsets receive no added full-cap credit on the fresh seed. Both full-cap consensus correlations remain close to zero.

Under the lower caps, both policies change the short-parent first coordinates only for compound 33 on 1KJL (−3.796 to −3.516) and change zero decisions on any receptor or the mean at the threefold threshold. Relative to the full-cap parent, the short energy-gate consensus loses 21 decisions and gains 3. Omitting any one compound leaves a net shortfall of 10 to 19 credits. Removing 49 or 58 gives the smallest shortfall. Both the small positive policy effect at full cap and its disappearance at lower caps occur at the same seed. Across two seeds in this one existing chemical series, there is no evidence of a broadly transferable improvement in relative affinity at matched or reduced actual cost. [Per-pair flips](flips.csv) · [All leave-one-compound-out contrasts](leave_one_out.csv)

## Actual cost and controller overhead

| Arm | Block | Wall s | Joint evaluations | Explicit inputs | Explicit optimizer calls | Replacements |
|---|---:|---:|---:|---:|---:|---:|
| Parent, 12,000 | 1 | 16.514 | 1,165,612 | 518 | 530 | 0 |
| Parent, lower caps | 1 | 15.466 | 1,095,490 | 510 | 522 | 0 |
| Energy gate, lower caps | 1 | 15.644 | 1,095,098 | 510 | 521 | 44 |
| Recurrence gate, lower caps | 1 | 15.864 | 1,095,218 | 510 | 521 | 50 |
| Energy gate, 12,000 | 2 | 15.549 | 1,164,473 | 518 | 525 | 48 |
| Recurrence gate, 12,000 | 2 | 15.207 | 1,164,463 | 518 | 525 | 55 |

Within the originally blocked comparison, the lower-cap energy gate took 5.27% less total wall time and 6.05% fewer joint evaluations than the full-cap parent. The recurrence gate took 3.94% less wall time and 6.04% fewer evaluations. Each lower-cap method used less wall time and fewer evaluations on **both receiver totals**. Cell wall time was at or below the full-cap parent in 57/94 energy-gate and 59/94 recurrence-gate runs, and all 94 in each arm used fewer evaluations. All 564 jobs reached the nominal cap in their hunt-plus-second stage. Because the cap is checked between MC steps, this stage exceeded the nominal budget by up to 201 evaluations for the 12,000 arm and 268 for the short arms. Explicit refinement follows the cap. Per-cell hunt and second BFGS counts and evaluations, and the input count to explicit handling, agree across the three policies at each common cap. The retained slots sometimes lead to different explicit optimizer retries, scalar rescoring and final output counts. CLI wall times include map setup, geometric checks, file I/O and these other operations. [Per-receiver costs and saturation](cost.csv) · [Per-cell counters, pose hashes and times](cells.csv)

The controller read and processed the pre-existing calibration table in 0.397 ms. Constructing that table from scratch would require 32.470 s of earlier parent-plus-gated docking across these receivers; this cost is separate from the 564 validation runs. Reusing the earlier Gal3 molecular inputs required no new target or ligand preparation. Preparing the original shared ligands cost 105.50 s summed over processes in the input study. Statistical evaluation and passive reruns for verification are outside predictive timing. The labelling source reports duplicate single-concentration inhibition measurements, no ligand-specific errors and no FP temperature for these Kd estimates. Thus pKd and pair ranks can be checked without adding a temperature-dependent experimental free-energy calibration.

## Reproduce

Use the same parent, final gated and recurrence research binaries used for the accepted-retention core benchmark. Run from repository root in a scientific Python environment with NumPy and SciPy. The script extracts only the prepared Gal3 fixtures from their existing archive. Use a reserved CPU number available on the host. The historical production used CPU 7; wall ratios on a different host can vary. Both raw timing blocks and their executable hashes, pose files, logs, actual controller outputs and CPU IDs are archived. [Production archive](runs.tar.xz)

```bash
D=benchmarks/accepted_retention/affinity_time
W=$(mktemp -d)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python "$D/transfer.py" "$W" --cpu 7 \
  --parent "$PARENT" --retention "$GATED" --recurrence "$RECURRENCE"
python "$D/transfer.py" "$W" --cpu 7 --full-cap \
  --parent "$PARENT" --retention "$GATED" --recurrence "$RECURRENCE"
python "$D/analyze.py" "$W" "$W/tables"
```

For offline verification, extract `runs.tar.xz` into a fresh directory and run only `analyze.py` on it. This checks the 564 raw PDBQT scores and counting logs before writing the six result tables. The cap controller hashes the earlier calibration file. Source ligand identities and experimental values are joined using the committed 47-graph manifest.
